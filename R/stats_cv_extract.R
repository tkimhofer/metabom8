
#' @title OPLS Component Estimation via Cross-Validation
#' @description
#' Performs orthogonal projections to latent structures (OPLS) modeling on
#' training folds of a cross-validation (CV) scheme.
#' The function fits one orthogonal and one predictive component per fold,
#' tracking scores, predictions, and residuals.
#' For the first orthogonal component (\code{nc = 1}), the full model structure
#' is initialized.
#' For later components (\code{nc > 1}), modeling proceeds on residual matrices
#' carried over from previous iterations.
#'
#' @param X Numeric matrix. Input data matrix (samples x features); only
#' required for \code{nc = 1}.
#' @param Ycs_fold List containing the response values (`Y`) for each
#' cross-validation split.
#' @param cv.set List of integer vectors. Each element contains training indices
#' for one CV round.
#'   These indices are adjusted internally by -1 due to 0-based indexing in the
#'   underlying C++ (Rcpp) routines.
#' @param nc Integer. Component number currently being estimated (\code{nc = 1}
#' for first orthogonal component).
#' @param mod.cv List. Model state to be updated with results from each CV iteration.
#'
#' @return A list of the same length as \code{cv.set}. Each list element contains:
#' \itemize{
#'   \item \code{t_xo}: Orthogonal component scores (matrix, samples x components).
#'   \item \code{t_xp}: Predictive component scores (matrix, samples x 1).
#'   \item \code{y_pred_train}: Predicted responses for training samples.
#'   \item \code{y_pred_test}: Predicted responses for held-out test samples.
#'   \item \code{x_res}: Residual matrix after filtering out orthogonal structure.
#' }
#'
#' @note Training indices in \code{cv.set} are internally adjusted by -1 before
#' being passed to Rcpp routines.
#' This is required because C++ uses zero-based indexing.
#'
#' @keywords internal
.oplsComponentCv <- function(X, Ycs_fold, cv.set, nc, mod.cv, acc) {

  n <- nrow(Ycs_fold[[1]]);

  if (nc == 1) {
    q <- ncol(Ycs_fold[[1]]);
    p <- ncol(X)

    mod.cv <- vector("list", length(cv.set))
    for(i in seq_along(cv.set)) {
      mod.cv[[i]] <- list(
        t_xo = matrix(0, n, 1),
        t_xp = matrix(0, n, 1),
        y_pred_train = matrix(0, n, q),
        y_pred_test = matrix(0, n, q),
        x_res = matrix(0, n, p)
      )
    }

  } else {
    for(i in seq_along(mod.cv)) {
      mod.cv[[i]]$t_xo <- cbind(mod.cv[[i]]$t_xo, 0)
      mod.cv[[i]]$t_xp <- cbind(mod.cv[[i]]$t_xp, 0)
    }
  }

  for (k in seq_along(cv.set)) {
    idc <- cv.set[[k]]
    not_idc <- setdiff(seq_len(n), idc)

    Xcs <- if (nc == 1) .scaleMatRcpp(X, idc - 1, center = TRUE, scale_type = 1)[[1]] else mod.cv[[k]]$x_res
    Ycs <- Ycs_fold[[k]]

    opls_filt  <- .nipOplsRcpp(X = Xcs[idc, , drop = FALSE], Y = Ycs[idc, , drop = FALSE])
    opls_filt$w_o <- .ensure_matrix(opls_filt$w_o)
    opls_filt$p_o <- .ensure_matrix(opls_filt$p_o)

    pred_comp  <- .nipPlsCompRcpp(opls_filt$X_res, Ycs[idc, , drop = FALSE], it_max = 800, eps = 1e-8)

    pred_comp$w_p <- .ensure_matrix(pred_comp$w_x)
    pred_comp$p_y <- .ensure_matrix(pred_comp$p_y)

    pred_comp$b   <- as.numeric(pred_comp$b)

    Xte <- as.matrix(Xcs[not_idc, , drop = FALSE])
    storage.mode(Xte) <- "double"

    opls_pred <- .oplsPredRcpp(opls_mod = opls_filt, pred_mod = pred_comp, Xnew = Xte)

    mod.cv[[k]]$t_xo[not_idc, nc]      <- opls_pred$t_xo_new
    mod.cv[[k]]$t_xp[not_idc, nc]      <- opls_pred$t_pred
    mod.cv[[k]]$y_pred_train[idc, ]    <- pred_comp$y_pred
    mod.cv[[k]]$y_pred_test[not_idc, ] <- opls_pred$y_pred
    mod.cv[[k]]$x_res[idc, ]           <- opls_filt$X_res
    mod.cv[[k]]$x_res[not_idc, ]       <- opls_pred$Xres

    acc$sum_test[not_idc, ] <- acc$sum_test[not_idc, ] + opls_pred$y_pred
    acc$n_test[not_idc]     <- acc$n_test[not_idc] + 1

    acc$sum_train[idc, ] <- acc$sum_train[idc, ] + pred_comp$y_pred
    acc$n_train[idc]     <- acc$n_train[idc] + 1
  }

  return(list(mod.cv = mod.cv, acc = acc))
}


.ensure_matrix <- function(x, ncol = 1L) {
  if (is.null(dim(x)))
    matrix(x, ncol = ncol)
  else
    x
}

#' @title Extract Cross-Validated OPLS Features
#'
#' @description
#' Extracts a specific data structure (e.g., scores or predictions) from cross-validation output returned by `.oplsComponentCv()`.
#' Depending on the cross-validation type and model structure, the output is aggregated across folds by computing the mean,
#' standard deviation (SD), and coverage.
#'
#' Here, **Coverage** refers to the proportion of cross-validation folds in which a given sample's value was observed
#' (i.e., not missing), serving as an estimate of how frequently that value was evaluated across folds.
#'
#' @param cv_obj Named list. Output of `.oplsComponentCv()` for a specific OPLS component.
#' @param feat Character. Feature name to extract from each CV fold object (e.g., `"t_xp"` or `"y_pred_test"`).
#'
#' @return A numeric matrix or array of the requested feature:
#' \itemize{
#'   \item For Monte Carlo CV and single-column Y: matrix with rows `"mean"`, `"sd"`, and `"coverage"` x samples.
#'   \item For Monte Carlo CV and multi-column Y: 3D array with shape \code{[3, samples, outcomes]}.
#'   \item For k-fold CV: matrix or array containing values for each sample (no aggregation, just fold values).
#' }
#'
#' @importFrom abind abind
#' @importFrom stats sd
#' @keywords internal
.extMeanCvFeat <- function(cv_obj, feat = "t_xp") {

  if (!feat %in% names(cv_obj[[1]]))
    stop("Feature name not in feature list.")

  inter <- lapply(cv_obj, `[[`, feat)

  f_mat <- abind::abind(inter, along = 3)

  apply(f_mat, c(1,2), function(x) {
    nz <- !is.na(x)
    if (any(nz)) sum(x[nz]) / sum(nz) else NA_real_
  })
}

