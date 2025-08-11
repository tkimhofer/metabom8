#' @title Check Y for missing values (PLS context)
#' @description Checks input Y for NA/NaN/Inf values, determines analysis type
#' (Regression or Discriminant Analysis), and returns cleaned Y matrix,
#' associated levels, and type.
#' @param Y A vector, matrix, or data.frame. Target variable.
#' @return A list: cleaned Y matrix, levels (if applicable), and type ('R' or 'DA')
#' @keywords internal
.checkYclassNas <- function(Y) {
  if (!inherits(Y, c("matrix", "data.frame"))) {
    Y <- matrix(Y, nrow = length(Y), ncol = 1)
  } else {
    Y <- as.matrix(Y)
  }

  if (anyNA(Y) || any(is.nan(Y)) || any(is.infinite(Y))) {
    stop("Input Y contains NA, NaN, or Inf values.")
  }

  if (is.numeric(Y)) {
    type <- "R"
    Y_out <- .prepareY(Y)
    Y <- Y_out[[1]]
    levs <- data.frame()
  } else {
    type <- "DA"
    Y_out <- .prepareY(Y)
    Y <- Y_out[[1]]

    levs <- unique(apply(Y, 2, function(x) length(unique(x))))
    if (length(levs) == 1 && levs[1] == 1) {
      stop("Input Y has only a single level.")
    }
    message("Performing discriminant analysis.")
  }

  if (ncol(Y) > 1) {
    type <- paste0(type, "-mY")
  }

  return(c(Y_out, type))
}


#' @title Check X matrix (PLS context)
#' @description Validates that input \code{X} is a numeric matrix and does not contain missing, NaN, or infinite values.
#' @param X Input data (univariate or multivariate), formatted as a matrix.
#' @return \code{NULL} if all checks pass; otherwise throws an informative error.
#' @keywords internal
#' @rdname checkXclassNas
.checkXclassNas <- function(X) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop(sprintf("Input X must be a numeric matrix. Use matrix() or as.matrix() to convert."))
  }

  if (any(is.na(X) | is.nan(X) | is.infinite(X))) {
    stop("Input X contains missing (NA), NaN, or infinite values.")
  }

  return(invisible(NULL))
}

#' @title Check Dimensions of X and Y Matrices (PLS Context)
#' @description
#' Checks that the input matrices \code{X} and \code{Y} have compatible dimensions
#' for Partial Least Squares (PLS) analysis.
#' Specifically, the number of observations (rows) must match, and
#' the number of variables (columns) in \code{X} must be greater than in \code{Y}.
#'
#' @param X Numeric matrix. Predictor data matrix (observations x variables).
#' @param Y Numeric matrix. Response data matrix (observations x variables).
#'
#' @return
#' Invisibly returns \code{NULL} if checks pass.
#' Throws an error if dimensions are incompatible.
#' @keywords internal
.checkDimXY <- function(X, Y) {
  if (nrow(X) != nrow(Y)) {
    stop("Dimensions of input X and Y do not match.")
  }
  if (ncol(X) <= ncol(Y)) {
    stop("Number of variables (columns) in X should be higher than in Y.")
  }
  invisible(NULL)
}

#' @title Evaluate Model Fit Progression Based on Generalization Performance
#'
#' @description
#' Determines whether an additional component improves model generalization based on cross-validated performance indices:
#' \eqn{Q^2} for regression (R) and AUROC for discriminant analysis (DA).
#' AUROC is computed on cross-validation test folds to ensure it reflects out-of-sample classification performance,
#' analogous to \eqn{Q^2} in regression.
#'
#' @param type Character. Either 'R' for regression or 'DA' for discriminant analysis. May be prefixed with '-mY' to indicate multi-column Y.
#' @param q2s Numeric vector of \eqn{Q^2} values for fitted components (used for regression).
#' @param auroc_cv Numeric vector of AUROC values computed on test folds during cross-validation (used for classification).
#' @param pc_max Integer. Maximum number of components allowed for fitting.
#'
#' @return Logical. Returns \code{FALSE} if another component should be fitted;
#' \code{TRUE} if the current component starts to overfit or provides minimal gain.
#'
#' @keywords internal
#'
.evalFit <- function(type, q2s, auroc_cv, pc_max) {
  type <- strsplit(type, '-')[[1]][1]

  if (type == "R") {
    cv_est <- q2s
  } else {
    cv_est <- auroc_cv
  }

  nc <- length(cv_est)
  if (nc == pc_max)
    return(TRUE)

  switch(type,
         R = {
           if (any(is.na(cv_est))) {
             stop("Something went wrong: Q2 is NA.")
           }
           if (nc > 1 && (diff(cv_est[(nc - 1):nc]) < 0.05 || cv_est[nc] > 0.98)) return(TRUE)
         },
         DA = {
           if (any(is.na(cv_est))) {
             stop("Something went wrong: AUROC is NA.")
           }
           if (nc > 1 && (diff(cv_est[(nc - 1):nc]) < 0.05 || cv_est[nc] > 0.97)) return(TRUE)
         })

  return(FALSE)
}

#' @title k-fold cross-validation index generator
#'
#' @description
#' Creates a list of training set indices for k-fold cross-validation (CV),
#' without considering class balance in Y. Each element of the list represents
#' the row indices used for training in one CV fold.
#'
#' @param k Integer. Number of CV folds. If \code{k} is not valid or too high relative to number of rows in Y,
#' it defaults to leave-one-out CV (LOO-CV).
#' @param Y Matrix. Outcome matrix (observations x variables).
#'
#' @return
#' A list of length \code{k}, where each element contains the training indices (integers)
#' for one CV fold.
#'
#' @keywords internal
.kFold <- function(k, Y){
  n <- nrow(Y)

  if (!is.numeric(k) || is.na(k) || is.infinite(k) || k < 2 || k > n) {
    message("Invalid k; defaulting to leave-one-out CV (k = n).")
    k <- n
  }

  if (k != n && n / k < 2) {
    message("Fold size is small; consider using fewer folds.")
  }

  sets <- sample(rep_len(seq_len(k), nrow(Y)))
  sets_list <- lapply(seq_len(k), function(i) {
    which(sets != i)
  })

  return(sets_list)
}


#' @title Stratified k-fold cross-validation index generator
#'
#' @description
#' Generates a list of training set indices for stratified k-fold cross-validation (CV).
#' Stratification is performed based on the first column of \code{Y}. For regression tasks,
#' \code{Y} is binned by quantiles to emulate class balance. Ensures class proportions are
#' preserved across folds when possible.
#'
#' @param k Integer. Number of folds.
#' @param stratified List with three elements:
#' \itemize{
#'   \item \code{type}: Character. Either \code{"R"} for regression or \code{"DA"} for discriminant analysis.
#'   \item \code{Y}: Matrix. Outcome matrix with a single column.
#'   \item \code{probs}: Numeric vector. Probabilities used for stratification of numeric \code{Y} (regression only).
#' }
#'
#' @return
#' A list of length \code{k}, each element containing the training set row indices (integers) for one CV fold.
#' Returns \code{NULL} and a warning if stratification is not feasible due to class imbalance.
#'
#' @keywords internal
.kFoldStratified <- function(k, stratified){

  if (!is.list(stratified) || length(stratified) != 3) {
    stop("'stratified' must be a list of length 3: type, Y, probs.")
  }

  Y <- stratified[[2]]
  if (anyNA(Y)) {
    stop("Stratified CV not possible: Y contains NA values.")
  }


  if (anyNA(Y)) {
    stop("Stratified CV not possible: Y contains NA values.")
  }

  if (grepl('R', stratified[[1]])) {
    Yori <- Y
    breaks <- unique(quantile(Yori[, 1], stratified[[3]]))
    if (length(breaks) < 3) {
      warning("Insufficient variability in Y for stratification - use unstratified k-fold or MC.")
      return(NULL)
    }
    Y <- cbind(cut(Yori[, 1], breaks = breaks, include.lowest = TRUE))
  } else {
    Y <- Y[, 1]
  }

  ct <- table(Y)
  if (min(ct) / max(ct) < 0.1) {
    message("Skewed class distribution detected. Consider using non-stratified k-fold CV.")
    return(NULL)
  }

  levs <- names(ct)
  if (min(ct) <= k || min(ct) / k < 2 || is.na(k) || is.infinite(k) || k < 2) {
    message(sprintf("Reducing k to %d due to small group size (min n = %d).", min(ct), min(ct)))
    k <- min(ct)
  }

  nobs <- floor(min(ct) / k)
  if (nobs == 0) stop("Too few observations per class for stratified CV.")

  set_lev <- lapply(levs, function(x, ks = k, y = Y, nob = nobs) {
    idx <- sample(which(y == x))
    k1 <- rep(seq_len(ks), length.out = length(idx))
    names(k1) <- idx
    return(k1)
  })

  idc <- unlist(set_lev)
  sets_list <- lapply(seq_len(k), function(i) {
    idx <- which(idc != i)
    as.numeric(names(idc[idx]))
  })

  return(sets_list)
}


#' @title Generate Monte Carlo Cross-Validation (MCCV) Training Indices
#' @description Generates a list of training set indices for Monte Carlo Cross-Validation.
#' @param k Integer. Number of MCCV iterations (i.e., training sets to generate).
#' @param Y A matrix (n x p), where rows are observations and columns are outcomes. Only the number of rows is used here.
#' @param split Numeric. Fraction of the data to include in each training set. Should be between 0 and 1.
#' @return A list of length \code{k}, each containing a vector of training set indices.
#' @details For each fold, a random sample (with replacement) of size \code{floor(nrow(Y) * split)} is drawn.
#' @keywords internal
.mc <- function(k, Y, split) {

  k <- ceiling(k)
  # # Check input: k must be a positive integer
  # if (!is.integer(k)) {
  #   k <- ceiling(k)
  #   warning(sprintf("The 'k' parameter should be an integer. Rounding up to k = %d.", k))
  # }

  if (k <= 0 || k > 1e6 || is.na(k) || is.infinite(k)) {
    stop("Invalid value for 'k'. Must be a positive integer and <= 1e6.")
  }

  # Check split parameter
  if (split <= 0 || split >= 1 || is.na(split) || is.infinite(split)) {
    stop("Invalid value for 'split'. Must be between 0 and 1 (exclusive).")
  }

  if (split > 0.9 || split < 0.3) {
    message("Unusual value for 'split'. Typical values are between 0.5 and 0.9.")
  }

  n <- nrow(Y)
  n_train <- floor(n * split)
  if (n_train < 1) {
    stop("Split results in zero training samples. Increase split or use a larger dataset.")
  }

  sets_list <- lapply(seq_len(k), function(i) {
    sample(seq_len(n), size = n_train, replace = TRUE)
  })

  return(sets_list)
}


#' @title Class-balanced Monte Carlo Cross-Validation (MCCV)
#' @description Generates class/group-balanced training set indices for each MCCV round.
#' @param k Integer. Number of MCCV training sets to generate.
#' @param split Numeric. Fraction of observations per class to include in each training set. Must be between 0 and 1.
#' @param stratified List of three elements:
#'   \itemize{
#'     \item{1: } Character. Outcome type: 'R' (regression) or 'DA' (discriminant analysis), optionally with suffix '-mY' for multi-column Y.
#'     \item{2: } Matrix. Response matrix \code{Y}, with \code{ncol(Y) == 1}.
#'     \item{3: } Numeric vector of probabilities used to stratify numeric Y if regression.
#'   }
#' @return A list of length \code{k}, each element containing a vector of row indices for the training set.
#' @details Each round samples a class-balanced subset with replacement. Useful for high-variance, imbalanced-class modeling.
#' @keywords internal
.mcBalanced <- function(k, split, stratified) {

  k <- ceiling(k)
  # # Validate k
  # if (!is.integer(k)) {
  #   k <- ceiling(k)
  #   warning(sprintf("The k-fold parameter should be an integer. Rounding to k = %d", k))
  # }
  if (k <= 2 || k > 1e6 || is.na(k) || is.infinite(k)) {
    stop("Invalid 'k' value. Must be between 3 and 1e6.")
  }

  # Validate split
  if (split <= 0 || split >= 1 || is.na(split) || is.infinite(split)) {
    stop("Split must be between 0 and 1 (exclusive).")
  }
  if (split < 0.3 || split > 0.9) {
    message("Unusual value for split. Typical values are between 0.5 and 0.9.")
  }

  # Prepare outcome Y
  type <- stratified[[1]]
  Ymat <- stratified[[2]]
  probs <- stratified[[3]]

  if (grepl("R", type)) {
    breaks <- quantile(Ymat[, 1], probs = probs)
    if (length(unique(breaks)) <= 2) {
      stop("Y has insufficient variability to create stratified bins - switch to unstratified CV.")
    }
    Y <- cut(Ymat[, 1], breaks = breaks, include.lowest = TRUE)
  } else {
    Y <- Ymat[, 1]
  }

  tabY <- table(Y)

  # Check class imbalance
  if (min(tabY) / max(tabY) < 0.1) {
    stop("Skewed class frequencies. Results will be biased. Consider using non-stratified CV.")
  }

  # Determine number of training samples per class
  n_obs <- floor(min(tabY) * split)
  if (n_obs == 0) {
    stop("Insufficient observations in minority class. Consider increasing split or aggregating classes.")
  }

  # Generate training sets
  sets_list <- lapply(seq_len(k), function(i) {
    idx_k <- unlist(lapply(names(tabY), function(class) {
      sample(which(Y == class), size = n_obs, replace = TRUE)
    }))
    return(idx_k)
  })

  return(sets_list)
}



#' @title Generate Cross-Validation (CV) Training Set Indices
#' @description Generates a list of row indices representing training sets for various CV strategies including stratified and Monte Carlo methods. Handles both regression (R) and discriminant analysis (DA), and multi-column Y.
#'
#' @param Y matrix or data.frame. Outcome matrix. If multi-column, only the first column will be used.
#' @param type character. Type of analysis: "R" for regression or "DA" for discriminant analysis. Use suffix "-mY" to indicate multi-column Y.
#' @param method character. One of: `"k-fold"`, `"k-fold_stratified"`, `"MC"`, `"MC_balanced"`.
#' @param k integer. Number of CV iterations/folds.
#' @param split numeric. For Monte Carlo methods, the fraction of samples to use in the training set (0 < split < 1).
#'
#' @return A list of integer vectors. Each vector contains row indices representing a training set.
#'
#' @details
#' - In `"k-fold_stratified"` and `"MC_balanced"` modes, stratification is applied to preserve class proportions or outcome distribution.
#' - If `type` is regression (`"R"`), quantile-based binning is applied to stratify continuous Y values.
#' - If `type` is discriminant analysis (`"DA"`), class labels are used directly.
#' @keywords internal
.cvSetsMethod <- function(Y, type, method = "k-fold_stratified", k = 7, split = 2/3) {

  valid_methods <- c("k-fold", "k-fold_stratified", "MC", "MC_balanced")
  if (!method %in% valid_methods) {
    stop("Check method argument. Valid methods are: ", paste(valid_methods, collapse = ", "))
  }

  Y <- as.matrix(Y)
  if (ncol(Y) > 1) {
    Y <- Y[, 1, drop = FALSE]
  }

  valid_methods <- c("k-fold", "k-fold_stratified", "MC", "MC_balanced")
  if (!method %in% valid_methods) {
    stop("Check method argument. Valid methods are: ", paste(valid_methods, collapse = ", "))
  }

  if (k %% 1 == 0) {
    k <- as.integer(k)
  }

  if (grepl("MC", method) && (split <= 0 || split >= 1)) {
    stop("Parameter 'split' must be between 0 and 1 (exclusive).")
  }

  k <- as.integer(k)
  if (k < 2) stop("Number of folds or iterations (k) must be at least 2.")


  sets_list <- switch(method,
                      `k-fold_stratified` = .kFoldStratified(k, stratified = list(type, Y, probs = c(0, 0.33, 0.66, 1))),
                      `k-fold` = .kFold(k, Y),
                      MC = .mc(k, Y, split),
                      MC_balanced = .mcBalanced(k, split, stratified = list(type, Y, probs = c(0, 0.33, 0.66, 1)))
  )

  if (is.null(sets_list)) {
    stop("Cross-validation set generation failed. Check input arguments.")
  }

  return(sets_list)
}

#' @title Prepare Response Vector for OPLS/OPLS-DA
#' @description
#' Converts a response vector \code{Y} into a numeric matrix suitable for regression or classification models.
#' - Numeric input is returned as a column matrix.
#' - Factor or character input:
#'   - If binary (2 levels): returns a single numeric column (values 1 and 2).
#'   - If multiclass (>2 levels): returns a dummy matrix with 1 for presence and -1 for absence.
#' Also returns a mapping between original labels and numeric values (only for categorical input).
#'
#' @param Y A numeric, factor, or character vector.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{[[1]]}}{Numeric matrix of processed response.}
#'   \item{\code{[[2]]}}{Data frame mapping class labels to encoding (empty for numeric input).}
#' }
#' @keywords internal
.prepareY <- function(Y) {
  if (!is.numeric(Y)) {
    Y_factor <- factor(Y)
    Y_levels <- levels(Y_factor)
    n_levels <- length(Y_levels)

    if (n_levels == 2) {
      # Binary classification: encode as 1/2
      Y_num <- as.numeric(Y_factor)
      mapping <- data.frame(Level = Y_levels,
                            Code = sort(unique(Y_num)),
                            stringsAsFactors = FALSE)
      return(list(matrix(Y_num, ncol = 1), mapping))

    } else {
      # Multiclass classification: dummy matrix with 1/-1
      Y_mat <- vapply(
        Y_levels,
        function(lvl) ifelse(Y_factor == lvl, 1, -1),
        numeric(length(Y_factor))
      )
      colnames(Y_mat) <- Y_levels

      mapping <- data.frame(Level = Y_levels,
                            Column = seq_along(Y_levels),
                            stringsAsFactors = FALSE)

      return(list(Y_mat, mapping))
    }

  } else {
    # Check if multi-column numeric matrix
    if (is.matrix(Y) && ncol(Y) > 1) {
      stop("Only single-column numeric input is supported for Y at this time.")
    }
    # Numeric (regression)
    return(list(matrix(Y, ncol = 1), data.frame()))
  }
}

#' @title OPLS Component Estimation via Cross-Validation
#' @description
#' Performs orthogonal projections to latent structures (OPLS) modeling on training folds of a cross-validation (CV) scheme.
#' The function fits one orthogonal and one predictive component per fold, tracking scores, predictions, and residuals.
#' For the first orthogonal component (\code{nc = 1}), the full model structure is initialized.
#' For later components (\code{nc > 1}), modeling proceeds on residual matrices carried over from previous iterations.
#'
#' @param X Numeric matrix. Input data matrix (samples x features); only required for \code{nc = 1}.
#' @param Y Numeric matrix. Response variable (can be dummy-coded for classification).
#' @param cv.set List of integer vectors. Each element contains training indices for one CV round.
#'   These indices are adjusted internally by -1 due to 0-based indexing in the underlying C++ (Rcpp) routines.
#' @param nc Integer. Component number currently being estimated (\code{nc = 1} for first orthogonal component).
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
#' @note Training indices in \code{cv.set} are internally adjusted by -1 before being passed to Rcpp routines.
#' This is required because C++ uses zero-based indexing.
#'
#' @keywords internal
.oplsComponentCv <- function(X, Y, cv.set, nc, mod.cv) {

  if (nc == 1) {
    mod.cv <- lapply(seq_along(cv.set), function(i) {
      list(
        t_xo = matrix(NA_real_, nrow = nrow(Y), ncol = 1),  # orth scores
        t_xp = matrix(NA_real_, nrow = nrow(Y), ncol = 1),  # predictive scores
        y_pred_train = matrix(NA_real_, nrow = nrow(Y), ncol = ncol(Y)),
        y_pred_test = matrix(NA_real_, nrow = nrow(Y), ncol = ncol(Y)),
        x_res = matrix(NA_real_, nrow = nrow(Y), ncol = ncol(X))  # residuals
      )
    })
  } else {
    mod.cv <- lapply(mod.cv, function(fold) {
      fold$t_xo <- cbind(fold$t_xo, NA_real_)
      fold$t_xp <- cbind(fold$t_xp, NA_real_) # single predictive component
      fold
    })
  }

  mod.cv <- lapply(seq_along(cv.set), function(k) {
    idc <- cv.set[[k]]

    Xcs <- if (nc == 1) .scaleMatRcpp(X, idc - 1, center = TRUE, scale_type = 1)[[1]] else mod.cv[[k]]$x_res
    Ycs <- .scaleMatRcpp(Y, idc - 1, center = TRUE, scale_type = 1)[[1]]

    opls_filt  <- .nipOplsRcpp(X = Xcs[idc, , drop = FALSE], Y = Ycs[idc, , drop = FALSE])
    pred_comp  <- .nipPlsCompRcpp(opls_filt$X_res, Ycs[idc, , drop = FALSE], it_max = 800, eps = 1e-8)

    Xte <- Xcs[-idc, , drop = FALSE]
    opls_pred <- .oplsPredRcpp(opls_mod = opls_filt, pred_mod = pred_comp, Xnew = Xte)

    # Xtr <- Xcs[idc, , drop = FALSE]
    # opls_pred_tr <- .oplsPredRcpp(opls_mod = opls_filt, pred_mod = pred_comp, Xnew = Xtr)

    mod.cv[[k]]$t_xo[-idc, nc]      <- opls_pred$t_xo_new
    mod.cv[[k]]$t_xp[-idc, nc]      <- opls_pred$t_pred
    mod.cv[[k]]$y_pred_train[idc, ] <- pred_comp$y_pred
    mod.cv[[k]]$y_pred_test[-idc, ] <- opls_pred$y_pred
    mod.cv[[k]]$x_res[idc, ]        <- opls_filt$X_res
    mod.cv[[k]]$x_res[-idc, ]       <- opls_pred$Xres

    mod.cv[[k]]
  })

  return(mod.cv)
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
#' @param cv_type Character. Cross-validation type, one of `"k-fold"`, `"k-fold_stratified"`, `"MC"`, or `"MC_balanced"`.
#' @param model_type Character. Indicates the output structure, typically `"mY"` for multi-column response matrix.
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
.extMeanCvFeat <- function(cv_obj, feat = "t_xp", cv_type, model_type) {

  ### feat t_xo -> can be more than single orthogonal component -> keep structure intact

  if (!feat %in% names(cv_obj[[1]])) stop("Feature name not in feature list.")
  oc_ind <- (feat == 't_xo') && (ncol(cv_obj[[1]][[feat]]) > 1)

  inter <- lapply(cv_obj, `[[`, feat)

  if (grepl("MC", cv_type) || (grepl("k-fold", cv_type) && grepl("train", feat))) {
    f_mat <- abind::abind(inter, along = 3)
    fout <- apply(f_mat, c(1, 2), function(x) {
      valid <- which(!is.na(x))
      mean_val <- sum(x, na.rm = TRUE) / length(valid)
      return(mean_val)
    })
    # if (grepl("mY", model_type)) {
    #   f_mat <- abind::abind(inter, along = 3)
    #   fout <- apply(f_mat, c(1, 2), function(x) {
    #     valid <- which(!is.na(x))
    #     mean_val <- sum(x, na.rm = TRUE) / length(valid)
    #     return(mean_val)
    #   })
    # } else {
    #   ### this is not bivariate Y, but it has multiple components for t_xo -> third dimension
    #   # ned to make sure this is always column vector, even in only single column
    #   f_mat <- do.call(cbind, inter)
    #   fout <- apply(f_mat, 1, function(x) {
    #     idx <- which(!is.na(x))
    #     mean = mean(x[idx], na.rm = TRUE)
    #     return(mean)
    #   })
    # }
  } else {
    if (grepl("mY", model_type) || !grepl("mY", model_type)) {
      # no mean required
      f_mat <- abind::abind(inter, along = 3)
      fout <- apply(f_mat, c(1, 2), function(x) {
        idx <- which(!is.na(x))
        return(x[idx])
      })
    }
  }

  return(fout)
}

.extMeanCvFeat_alt <- function(cv_obj, feat = "t_xp", cv_type, model_type) {
  if (!feat %in% names(cv_obj[[1]])) stop("Feature name not in feature list.")
  inter <- lapply(cv_obj, `[[`, feat)
  if (grepl("MC", cv_type) || (grepl("k-fold", cv_type) && grepl("train", feat))) {
    if (grepl("mY", model_type)) {
      # cbind(inter)
      f_mat <- abind::abind(inter, along = 3)
      fout <- apply(f_mat, c(1, 2), function(x) {
        valid <- which(!is.na(x))
        mean_val <- sum(x, na.rm = TRUE) / length(valid)
        return(mean_val)
      })
    } else {
      f_mat <- do.call(cbind, inter)
      fout <- apply(f_mat, 1, function(x) {
        idx <- which(!is.na(x))
        mean <- mean(x[idx], na.rm = TRUE)
        return(as.vector(mean))
      })
    }
  } else {
    if (grepl("mY", model_type) || !grepl("mY", model_type)) {
      f_mat <- abind::abind(inter, along = 3)
      fout <- apply(f_mat, c(1, 2), function(x) {
        idx <- which(!is.na(x))
        return(as.vector(x[idx]))
      })
    }
  }

  return(fout)
}



#' @title \eqn{R^2} and \eqn{Q^2} Calculation for OPLS Models
#'
#' @description
#' Computes the coefficient of determination (\eqn{R^2} or \eqn{Q^2}) for OPLS regression models using the formula:
#' \deqn{R^2 = 1 - \frac{PRESS}{TSS}}, where PRESS is the prediction error sum of squares,
#' and TSS is the total sum of squares. If `ytss` is not provided, it is computed directly from `Y`.
#'
#' This function supports matrix inputs (e.g., for multi-class outcomes) and averages across columns
#' to return a single summary \eqn{R^2}/\eqn{Q^2} value.
#'
#' @param Y Numeric matrix. True response values. For classification, this should be a dummy matrix.
#' @param Yhat Numeric matrix. Predicted response values from the model.
#' @param ytss Numeric (optional). Total sum of squares of `Y`. If not provided, it will be computed internally.
#'
#' @return A single numeric value representing \eqn{R^2} or \eqn{Q^2}.
#'
#' @keywords internal
.r2 <- function(Y, Yhat, ytss = NULL) {
  if (!is.numeric(Y)) stop("Y must be numeric")
  if (!is.numeric(Yhat)) stop("Yhat must be numeric")

  if (is.null(dim(Y))) {
    if (length(Y) != length(Yhat)) stop("Lengths of Y and Yhat must match")
    ncol_Y <- 1
  } else {
    if (!all(dim(Y) == dim(Yhat))) stop("Dimensions of Y and Yhat must match")
    ncol_Y <- ncol(Y)
  }

  if (!is.null(ytss) && !is.numeric(ytss)) stop("ytss must be numeric or NULL")

  press <- sum((Y - Yhat)^2, na.rm = TRUE) / ncol_Y

  if (is.null(ytss)) {
    ytss <- sum((Y - mean(Y))^2, na.rm = TRUE) / ncol_Y
  }

  return(1 - (press / ytss))
}

#' @title Summary of OPLS Model Components
#'
#' @description
#' Generates a tabular and graphical summary of orthogonal and predictive components from an OPLS model.
#' Includes explained variance (\eqn{R^2}), predictive accuracy (\eqn{Q^2}), and AUROC metrics (for classification models).
#'
#' For discriminant models (DA), AUROC values are included. For regression models (R), \eqn{R^2} and \eqn{Q^2} are reported.
#'
#' @param type Character. Model type: either "DA" (discriminant analysis) or "R" (regression). May include prefix "-mY" for multi-response models.
#' @param r2x_comp Numeric vector. Proportion of X variance explained by each orthogonal component.
#' @param r2_comp Numeric vector. Proportion of Y variance explained by the predictive component.
#' @param q2_comp Numeric vector. Cross-validated predictive accuracy per component (\eqn{Q^2}).
#' @param aucs_tr Numeric vector. AUROC values for training sets (only for "DA").
#' @param aucs_te Numeric vector. AUROC values for test sets (only for "DA").
#' @param cv List. Cross-validation configuration as used in the `opls()` function.
#'
#' @return A list with:
#' \itemize{
#'   \item A \code{data.frame} summarizing the metrics per component.
#'   \item A \code{ggplot2} object visualizing the component summary.
#' }
#'
#' @importFrom scales breaks_pretty
#' @importFrom ggplot2 ggplot geom_bar aes_string scale_fill_manual scale_alpha labs theme_bw theme
#' @importFrom ggplot2 element_line element_blank element_text element_rect
#' @importFrom reshape2 melt
#' @keywords internal
.orthModelCompSummary <- function(type, r2x_comp, r2_comp, q2_comp, aucs_tr, aucs_te, cv) {
  type <- strsplit(type, '-')[[1]][1]

  # Construct summary table
  model_summary <- switch(type,
                          'DA' = data.frame(
                            PC_pred = 1,
                            PC_orth = seq_along(aucs_tr),
                            R2X = round(r2x_comp, 2),
                            AUROC = round(aucs_tr, 2),
                            AUROC_CV = round(aucs_te, 2)
                          ),
                          'R' = data.frame(
                            PC_pred = 1,
                            PC_orth = seq_along(q2_comp),
                            R2X = round(r2x_comp, 2),
                            R2Y = round(r2_comp, 2),
                            Q2 = round(q2_comp, 2)
                          )
  )

  model_summary$PC_orth <- seq_len(nrow(model_summary))

  mm <- reshape2::melt(model_summary, id.vars = c("PC_orth", "PC_pred"))
  mm$PC <- paste0(mm$PC_pred, '+', mm$PC_orth)
  mm$alph <- 1
  mm$alph[mm$PC_orth == max(mm$PC_orth)] <- 0.7

  # Avoid very negative bars pulling down y-axis
  idx <- which(mm$value < (-0.015) & mm$variable == 'Q2')
  mm$value[idx] <- -0.01

  mm[is.na(mm$value),"value"] <- 0 # r2x not calc for overfitted component

  if(grepl('k-fold', cv$method)){
    cv_text <- paste0('\nCross validation: ', cv$method, ' (k=', cv$k, ')')
  } else{
    cv_text <- paste0('\nCross validation: ', cv$method, ' (k=', cv$k, ', ratio test/training set=', round(cv$split, 2), ')')
  }

  g <- ggplot2::ggplot(mm, ggplot2::aes(x = !!sym("PC"), y = !!sym("value"), fill = !!sym("variable"))) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", colour = NA, ggplot2::aes(alpha = !!sym("alph"))) +
    ggplot2::scale_fill_manual(
      values = c(R2X = "#1b9e77", R2Y = "#d95f02", Q2 = "#7570b3", AUROC = "#d95f02", AUROC_CV = "#7570b3"),
      labels = c(
        R2X = expression(R^2 * X),
        R2Y = expression(R^2 * Y),
        Q2 = expression(Q^2),
        AUROC = expression(AUROC),
        AUROC_CV = expression(AUROC[cv])
      ),
      name = ""
    ) +
    ggplot2::scale_alpha(guide = "none", limits = c(0, 1)) +
    ggplot2::labs(
      x = "Predictive + Orthogonal Component(s)",
      y = "",
      title = paste("O-PLS-", type, "  - Component Summary", sep = ""),
      caption = cv_text
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.text = element_text(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "black", linewidth = 0.15),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.55),
      axis.line.y = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_rect(colour = "white"),
      text = element_text(family = "Helvetica")
    )

  # Set y-axis limits
  if (any(mm$value < 0, na.rm = TRUE)) {
    g <- g + ggplot2::scale_y_continuous(limits = c(-0.015, 1), breaks = breaks_pretty(), expand = c(0, 0))
  } else {
    g <- g + ggplot2::scale_y_continuous(limits = c(0, 1), breaks = breaks_pretty(), expand = c(0, 0))
  }

  return(list(model_summary, g))
}

#' @title Summary of PLS Model Components
#'
#' @description
#' Generates a tabular and graphical summary of PLS model components.
#' Includes explained variance in X Space (\eqn{R^2X}), predictive accuracy (\eqn{Q^2Y}), and AUROC metrics, depending on model type.
#'
#' For discriminant models (DA), AUROC values are included. For regression models (R), \eqn{R^2Y} and \eqn{Q^2} are reported.
#'
#' @param type Character. Model type: either "DA" (discriminant analysis) or "R" (regression). May include prefix "-mY" for multi-response models.
#' @param r2x_comp Numeric vector. Proportion of X variance explained by each component.
#' @param r2_comp Numeric vector. Proportion of Y variance explained by each component.
#' @param q2_comp Numeric vector. Cross-validated predictive accuracy per component (\eqn{Q^2}).
#' @param aucs_tr Numeric vector. in-sample prediction accuracy for categorical Y ("DA").
#' @param aucs_te Numeric vector. out-of-sample prediction accuracy for categorical Y ("DA").
#' @param cv List. Cross-validation parameters.
#'
#' @return A list with:
#' \itemize{
#'   \item A \code{data.frame} summarizing the metrics per component.
#'   \item A \code{ggplot2} object visualizing the component summary.
#' }
#'
#' @importFrom scales breaks_pretty
#' @importFrom ggplot2 ggplot geom_bar aes_string scale_fill_manual scale_alpha labs theme_bw theme
#' @importFrom ggplot2 element_line element_blank element_text element_rect
#' @importFrom reshape2 melt
#' @keywords internal
.plsModelCompSummary <- function(type, r2x_comp, r2_comp, q2_comp, aucs_tr, aucs_te, cv) {
  type <- strsplit(type, '-')[[1]][1]

  # Construct summary table
  model_summary <- switch(type,
                          'DA' = data.frame(
                            PC_pred = seq_along(r2x_comp),
                            PC_orth = 0,
                            R2X = round(r2x_comp, 2),
                            AUROC = round(aucs_tr, 2),
                            AUROC_CV = round(aucs_te, 2)
                          ),
                          'R' = data.frame(
                            PC_pred = 0,
                            PC_orth = seq_along(q2_comp),
                            R2X = round(r2x_comp, 2),
                            R2Y = round(r2_comp, 2),
                            Q2 = round(q2_comp, 2)
                          )
  )

  # Reshape for plotting
  mm <- reshape2::melt(model_summary, id.vars = c("PC_orth", "PC_pred"))
  mm$PC <- paste0('PC', mm$PC_pred)
  mm$alph <- 1
  mm$alph[mm$PC_orth == max(mm$PC_orth)] <- 0.7

  # Avoid very negative bars pulling down y-axis
  idx <- which(mm$value < (-0.015) & mm$variable == 'Q2')
  mm$value[idx] <- -0.01

  if(grepl('k-fold', cv$method)){
    cv_text <- paste0('\nCross validation: ', cv$method, ' (k=', cv$k, ')')
  } else{
    cv_text <- paste0('\nCross validation: ', cv$method, ' (k=', cv$k, ', ratio test/training set=', round(cv$split, 2), ')')
  }

  g <- ggplot2::ggplot(mm, ggplot2::aes(x = !!sym("PC"), y = !!sym("value"), fill = !!sym("variable"))) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", colour = NA, ggplot2::aes(alpha = !!sym("alph"))) +
    ggplot2::scale_fill_manual(
      values = c(R2X = "lightgreen", R2Y = "lightblue", Q2 = "red", AUROC = "black", AUROC_CV = "red"),
      labels = c(
        R2X = expression(R^2 * X),
        R2Y = expression(R^2 * Y),
        Q2 = expression(Q^2),
        AUROC = expression(AUROC),
        AUROC_CV = expression(AUROC[cv])
      ),
      name = ""
    ) +
    ggplot2::scale_alpha(guide = "none", limits = c(0, 1)) +
    ggplot2::labs(
      x = "Component(s)",
      y = "",
      title = paste("PLS-", type, "  - Component Summary", sep = ""),
      caption =  cv_text
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.text = element_text(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "black", linewidth = 0.15),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.55),
      axis.line.y = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_rect(colour = "white"),
      text = element_text(family = "Helvetica")
    )

  # Set y-axis limits
  if (any(mm$value < 0, na.rm = TRUE)) {
    g <- g + ggplot2::scale_y_continuous(limits = c(-0.015, 1), breaks = breaks_pretty(), expand = c(0, 0))
  } else {
    g <- g + ggplot2::scale_y_continuous(limits = c(0, 1), breaks = breaks_pretty(), expand = c(0, 0))
  }

  return(list(model_summary, g))
}



#' @title Check Consistency of Spectral Data and ppm Vector
#'
#' @description
#' Validates that the NMR data matrix and the chemical shift reference vector (`ppm`) are consistent.
#' Specifically, it checks that:
#' \itemize{
#'   \item `ppm` contains no NA or infinite values.
#'   \item The number of columns in `X` matches the length of `ppm`.
#' }
#'
#' This function is typically used as a preliminary data validation step before downstream spectral analysis.
#'
#' @param X Numeric matrix. NMR spectral matrix with samples in rows and variables in columns (e.g., intensities).
#' @param ppm Numeric vector. Chemical shift reference values corresponding to columns of `X`.
#'
#' @return Logical. Returns \code{TRUE} if inputs are valid, \code{FALSE} otherwise.
#'
#' @keywords internal
.check_X_ppm <- function(X, ppm) {
  if (any(is.na(ppm) | is.infinite(ppm))) return(FALSE)

  if (ncol(X) != length(ppm)) return(FALSE)

  return(TRUE)
}


#' @title Ensure Input is a Row Matrix
#'
#' @description
#' Ensures that the input object `X` is in row-matrix form. If `X` is a vector (i.e., no `ncol` attribute),
#' it is converted to a matrix with 1 row and \code{length(X)} columns.
#' This is helpful for standardizing inputs before applying matrix operations.
#'
#' @param X Numeric vector, matrix, or array. Typically NMR data where rows represent spectra.
#'
#' @return A numeric matrix with one row if input was a vector; otherwise, returns `X` unchanged.
#'
#' @keywords internal
.dimX <- function(X) {
  if (is.null(ncol(X))) return(t(X))
  return(X)
}

#' @title Select Indices for a Chemical Shift Region
#'
#' @description
#' Returns the indices of the chemical shift vector (\code{ppm}) that fall within the specified range.
#'
#' @param range Numeric vector of length 2. Specifies the chemical shift region (in ppm) of interest.
#' @param ppm Numeric vector. The full chemical shift axis (in ppm).
#'
#' @return Integer vector of indices corresponding to \code{ppm} values within the given range.
#'
#' @export
#'
#' @examples
#' data(covid_raw)
#' X <- covid_raw$X
#' ppm <- covid_raw$ppm
#' idx_tsp <- get_idx(c(-0.1, 0.1), ppm)
#' ppm[range(idx_tsp)]
#' plot(ppm[idx_tsp], X[1, idx_tsp], type = 'l')
#'
#' @family NMR
get_idx <- function(range = c(1, 5), ppm) {
  range <- sort(range, decreasing = TRUE)
  which(ppm <= range[1] & ppm >= range[2])
}

#' @rdname get_idx
#' @export
get.idx <- function(range = c(1, 5), ppm) {
  warning("`get.idx` is deprecated and will be removed in future versions. Use `get_idx` instead.", call. = FALSE)
  get_idx(range, ppm)
}


#' @title Min-Max Scaling to Arbitrary Range
#' @description
#' Rescales a numeric vector to a specified range using min-max scaling.
#' This is a generalized form of min-max normalization allowing any output range.
#'
#' @param x Numeric vector. Input values to be scaled.
#' @param ra Numeric vector of length 2. Desired output range (e.g., \code{c(5, 10)}).
#'
#' @return A numeric vector of the same length as \code{x}, scaled to the range \code{ra}.
#'
#' @details
#' The scaled values are computed as:
#' \deqn{x_{scaled} = \frac{x - \min(x)}{\max(x) - \min(x)} \cdot (r_{max} - r_{min}) + r_{min}}
#'
#' @examples
#' x <- rnorm(20)
#' plot(x, type = 'l'); abline(h = range(x), lty = 2)
#' points(scRange(x, ra = c(5, 10)), type = 'l', col = 'red'); abline(h = c(5, 10), col = 'red', lty = 2)
#'
#' @seealso [minmax()]
#' @export
scRange <- function(x, ra) {
  ra <- sort(ra)
  (ra[2] - ra[1]) * (x - min(x)) / (max(x) - min(x)) + ra[1]
}

#' @title Min-Max Scaling to \eqn{[0,1]}
#' @description
#' Scales a numeric vector to the range [0, 1] using min-max normalization.
#' This is a special case of \code{\link{scRange}}.
#'
#' @param x Numeric vector. Input values to be scaled.
#'
#' @return A numeric vector of the same length as \code{x}, scaled to the range \code{[0, 1]}.
#'
#' @details
#' The scaled values are computed as:
#' \deqn{x_{scaled} = \frac{x - \min(x)}{\max(x) - \min(x)}}
#'
#' Equivalent to \code{scRange(x, ra = c(0, 1))}.
#'
#' @examples
#' x <- rnorm(20)
#' plot(x, type = 'l'); abline(h = range(x), lty = 2)
#' points(minmax(x), type = 'l', col = 'blue'); abline(h = c(0, 1), col = 'blue', lty = 2)
#'
#' @seealso [scRange()] for flexible output ranges.
#' @family NMR ++
#' @export
minmax <- function(x) {
  stopifnot(is.numeric(x))
  r <- range(x, na.rm = na.rm)
  d <- r[2] - r[1]
  if (is.na(d) || d == 0) {
    return(ifelse(is.na(x), NA_real_, 0))
  }
  (x - r[1]) / d
}

#' @title Full Width at Half Maximum (FWHM) Estimation
#'
#' @description
#' Calculates the full width at half maximum (FWHM, or line width) of a peak
#' within a specified chemical shift range in each NMR spectrum.
#'
#' @param X Numeric matrix. NMR spectral data with rows as spectra and columns as chemical shift (ppm) variables.
#' @param ppm Numeric vector. Chemical shift axis corresponding to columns of \code{X}.
#' @param shift Numeric vector of length 2. Chemical shift range that includes the singlet peak (e.g., \code{c(-0.1, 0.1)} for TSP).
#' @param sf Numeric scalar. Spectrometer frequency in MHz (e.g., 600 for 600 MHz).
#'
#' @return Numeric vector. FWHM in Hz, one value per row in \code{X}.
#'
#' @details
#' For each spectrum, the function extracts the region defined by \code{shift},
#' interpolates the signal, identifies the full width at half maximum (FWHM),
#' and converts the width from ppm to Hz using \code{sf}.
#'
#' This is commonly used as part of technical QC for NMR data, especially for evaluating the TSP reference peak,
#' where a sharp peak implies good resolution and field homogeneity.
#'
#' @note
#' The chemical shift axis direction is auto-detected. If a spectrum does not have a clear half-maximum width,
#' \code{NA} is returned for that spectrum.
#'
#' @examples
#' # Simulated data example
#' ppm <- seq(-0.2, 0.2, length.out = 1000)
#' spec <- dnorm(ppm, mean = 0, sd = 0.02)
#' X <- matrix(rep(spec, each = 10), nrow = 10, byrow = TRUE)
#' sf <- 600  # MHz
#' fwhm_vals <- lw(X, ppm, shift = c(-0.1, 0.1), sf = sf)
#' hist(fwhm_vals, main = "Line Width (Hz)", xlab = "FWHM [Hz]")
#'
#' @seealso \code{\link{get_idx}}
#' @family NMR
#' @importFrom stats approxfun
#' @export
lw <- function(X, ppm, shift = c(-0.1, 0.1), sf) {
  idx <- get_idx(shift, ppm)
  asign <- sign(ppm[2] - ppm[1])  # direction of ppm axis

  fwhm <- apply(X[, idx, drop = FALSE], 1, function(x) {
    pp <- ppm[idx]
    if (min(x, na.rm = TRUE) < 0) {
      x <- x + abs(min(x, na.rm = TRUE))
      baseline <- 0
    } else {
      baseline <- min(x, na.rm = TRUE)
    }

    height <- max(x, na.rm = TRUE) - baseline
    hw <- baseline + 0.5 * height

    f <- approxfun(pp, x)
    x_new <- seq(pp[1], pp[length(pp)], by = 1e-05 * asign)
    y_new <- f(x_new)

    above_hw <- which(y_new > hw)
    if (length(above_hw) < 2) {
      return(NA_real_)
    }

    fwhm_val <- diff(range(x_new[above_hw]))
    return(fwhm_val)
  })

  return(fwhm * sf)
}



#' @title Estimate Noise Level in 1D NMR Spectra
#'
#' @description
#' Estimates the noise level for each spectrum in a 1D NMR dataset by analyzing a signal-free region.
#' Useful for quality control (QC) before or after spectral processing such as normalization or baseline correction.
#'
#' @param X Numeric matrix. NMR data with spectra represented in rows (samples x chemical shift bins).
#' @param ppm Numeric vector. Chemical shift positions in ppm. Must match the number of columns in \code{X}.
#' @param where Numeric vector of length 2. Specifies the ppm range to use for noise estimation.
#'   This range should be free of peaks (i.e., no expected metabolite signal). Default is \code{c(14.6, 14.7)}.
#'
#' @details
#' Noise is estimated after baseline correction using asymmetric least squares smoothing
#' (via \code{\link[ptw]{asysm}}). The output is based on the 95th percentile of the signal in the specified
#' region and assumes a minimum of 50 points in that region. If fewer points are available, the function stops.
#'
#' @return Numeric vector of noise levels, one value per spectrum (length equals \code{nrow(X)}).
#'
#' @note Choose the \code{where} region carefully based on acquisition setup.
#' Signal-free regions are typically at the high or low ppm ends (e.g., >10 ppm or <0 ppm).
#'
#' @importFrom ptw asysm
#' @export
#'
#' @examples
#' data(covid_raw)
#' X <- covid_raw$X
#' ppm <- covid_raw$ppm
#' noise_vals <- noise.est(X, ppm, where = c(8.8, 8.9))
#' hist(noise_vals, main = "Estimated Noise", xlab = "Noise Level")
#'
#' @family NMR
noise.est <- function(X, ppm, where = c(14.6, 14.7)) {
  idx <- get_idx(where, ppm)

  if (length(idx) < 50) {
    stop("Insufficient number of points in selected region for noise estimation (minimum: 50).")
  }

  noise <- apply(X[, idx, drop = FALSE], 1, function(x) {
    corrected <- x - asysm(x, lambda = 10000)
    quantile(corrected, p = 0.95, na.rm = TRUE)
  })

  return(noise)
}



#' @title Prepare Data Frame for Projection or Score Plotting
#'
#' @description
#' Helper function to extract and structure projection (`p`) or score (`t`) components from a PCA or OPLS object for downstream visualization (e.g., scatter plots).
#'
#' @param obj An object of class `PCA_metabom8` or `OPLS_metabom8`.
#' @param pc Character vector indicating which components to extract, e.g., \code{c("1", "o1")} for predictive and orthogonal components.
#' @param an Named list of annotation vectors to append to the output data frame. Each vector should be length 1 or \code{nrow(X)}.
#' @param type Character. One of `"p"` (projection), `"t"` (scores), or `"t_cv"` (cross-validated scores).
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{df}: A data frame containing extracted components and annotations.
#'   \item \code{an_le}: Number of annotation columns.
#' }
#'
#' @details
#' Used internally for constructing input to plotting functions. The \code{pc} vector can contain `"o"`-prefixed entries to access orthogonal components in OPLS.
#'
#' @keywords internal
.viz_df_helper <- function(obj, pc, an, type = "p") {
  an <- .check_an_viz(an, obj)
  idx_orth <- grepl("o", pc)
  pc1 <- if (any(idx_orth)) as.numeric(gsub("o", "", pc)) else pc

  if (any(is.na(pc1) | is.infinite(pc1))) {
    stop("Invalid 'pc' argument: contains NA or Inf.")
  }

  if (any(idx_orth) && any(pc1[idx_orth] > nrow(obj@p_orth))) {
    stop("Orthogonal component index exceeds number of components.")
  }

  com <- list()
  cls <- class(obj)[1]

  extract_matrix <- function(source, comp_index) {
    if (length(comp_index) == 1) return(source[, comp_index])
    lapply(seq_along(comp_index), function(i) source[, comp_index[i]])
  }

  if (type == "p") {
    for (i in seq_along(pc)) {
      com[[i]] <- if (cls == "PCA_metabom8") {
        obj@p[, pc[i]]
      } else if (cls == "PLS_metabom8") {
        obj@p_pred[ pc[i],]
      } else {
        if (grepl("o", pc[i])) obj@p_orth[pc1[i], ] else obj@p_pred[1, ]
      }
    }
  } else if (type == "t") {
    for (i in seq_along(pc)) {
      com[[i]] <- if (cls == "PCA_metabom8") {
        obj@t[, pc[i]]
      } else if (cls == "PLS_metabom8") {
        obj@t_pred[, pc[i]]
      } else {
        if (grepl("o", pc[i])) obj@t_orth[, pc1[i]] else obj@t_pred[, 1]
      }
    }
  } else if (type == "t_cv") {
    if (cls == "PCA_metabom8") stop("t_cv not defined for PCA objects.")
    if (cls == "PLS_metabom8") stop("t_cv not implemented for PLS objects.")
    for (i in seq_along(pc)) {
      com[[i]] <- if (grepl("o", pc[i])) obj@t_orth_cv[, pc1[i]] else obj@t_pred_cv[, 1]
    }

  }

  melted <- as.data.frame(com)
  colnames(melted) <- paste0("pc_", pc)

  # Process annotations
  an_full <- lapply(an, function(x, le_m = nrow(melted)) {
    len_x <- length(x)
    if (len_x == 1) return(rep(x, le_m))
    if (len_x == le_m) return(x)
    stop("Mismatch between annotation length and data frame rows.")
  })

  an_df <- as.data.frame(an_full)
  colnames(an_df) <- names(an)
  melted <- cbind(melted, an_df)

  return(list(df = melted, an_le = ncol(an_df)))
}


#' @title Check Validity of Selected Principal or Orthogonal Component
#' @description
#' Internal consistency check for principal (or orthogonal) component selection used in plotting or extraction functions.
#'
#' @param pc Character or numeric. Selected component, e.g., \code{1} or \code{"o1"}.
#' @param mod An object of class \code{PCA_metabom8} or \code{OPLS_metabom8}.
#' @param le Integer. Expected maximum length of \code{pc}; usually 1.
#' @param type Character. Type of component to check: \code{"p"} for projection or \code{"t"} for scores. Currently only used for logic branching.
#'
#' @return Nothing. Function is used for its side effect (i.e., error checking).
#'
#' @keywords internal
.check_pc_model <- function(pc, mod, le = 1, type = 'p') {

  if (is.na(pc) || is.infinite(pc) || length(pc) > le) {
    stop("Check 'pc' argument: It must be finite and of expected length.")
  }

  mod_class <- class(mod)[1]

  if (mod_class == "PCA_metabom8") {
    if (!is.numeric(pc)) stop("PC value must be numeric for PCA.")
    if (max(pc) > ncol(mod@p)) stop("PC value exceeds number of PCA components.")
  }

  if (mod_class == "OPLS_metabom8") {
    idx_orth <- grepl("o", pc)
    if (any(idx_orth)) {
      pc1 <- as.numeric(gsub("o", "", pc))
      if (any(is.na(pc1)) || any(is.infinite(pc1))) {
        stop("Invalid orthogonal component identifier in 'pc'.")
      }
      if (any(pc1 > nrow(mod@p_orth))) {
        stop("Requested orthogonal component exceeds available components.")
      }
    } else {
      if (!is.numeric(as.numeric(pc))) {
        stop("Predictive component index should be numeric or coercible to numeric.")
      }
      if (length(mod@p_pred) == 0) {
        stop("Predictive component appears to be missing in model object.")
      }
    }
  }
}



#' @title Backscaling for NMR Loadings
#' @description
#' Computes backscaled loadings for visualization in NMR spectral space. The loadings are rescaled to original variable units using standard deviations from preprocessed input data. This is often used to interpret model loadings in terms of chemical shift.
#'
#' @param mod An object of class \code{PCA_metabom8} or \code{OPLS_metabom8}.
#' @param pc Character or numeric. Principal or orthogonal component identifier (e.g., \code{1}, \code{"o1"}).
#' @param idx Integer vector. Indices defining the region of interest in the \code{ppm} scale.
#' @param ppm Numeric vector. Chemical shift axis in parts per million (ppm).
#'
#' @return A \code{data.frame} containing:
#' \describe{
#'   \item{p_bs}{Backscaled loadings for the selected component.}
#'   \item{p_abs}{Normalized absolute loadings (min-max scaled).}
#'   \item{ppm}{Chemical shift values (ppm).}
#' }
#'
#' @details
#' Backscaling is computed as: \code{p_bs = p * sd(X)}, where \code{p} is the loading vector and \code{sd(X)} the original variable standard deviations.
#'
#' @keywords internal
#' @examples
#' # Example usage:
#' # df <- .load_backscaled_nmr(model, "1", get_idx(c(0.5, 4.5), ppm), ppm)
.load_backscaled_nmr <- function(mod, pc, idx, ppm) {
  p_mod <- .viz_df_helper(mod, pc, an = NA, type = 'p')

  p_bs <- p_mod$df[, 1] * mod@X_sd         # backscaling using variable standard deviations
  p_abs <- minmax(abs(p_mod$df[, 1]))      # normalized (absolute) loadings for visualization

  df <- data.frame(p_bs = p_bs, p_abs = p_abs, ppm = ppm)
  df <- df[idx, ]

  return(df)
}


#' @title Correlation and Covariance of Scores with NMR Data
#' @description
#' Computes the correlation and covariance between model component scores and NMR spectral variables.
#' This is typically used to reconstruct statistical relevance maps (e.g., statistical total correlation spectroscopy, STOCSY-like plots).
#'
#' @param mod An object of class \code{PCA_metabom8} or \code{OPLS_metabom8}.
#' @param pc Character or numeric. Component identifier (e.g., \code{1}, \code{"o1"}).
#' @param X Numeric matrix. NMR spectra matrix (samples x variables).
#' @param idx Integer vector. Indices specifying the region of interest in the chemical shift vector.
#' @param ppm Numeric vector. Chemical shift values (ppm).
#'
#' @return A \code{data.frame} with the following columns:
#' \describe{
#'   \item{cor}{Absolute correlation between model scores and spectral intensities.}
#'   \item{cov}{Covariance between model scores and spectral intensities.}
#'   \item{ppm}{Chemical shift values (ppm).}
#' }
#'
#' @details
#' Correlations are computed as Pearson coefficients between the selected component's scores and each spectral variable across samples.
#' Covariances are similarly computed. This provides insight into variables driving component separation.
#'
#' @examples
#' # df <- .load_stat_reconstr_nmr(model, "1", X, get_idx(c(0.5, 4.5), ppm), ppm)
#'
#' @keywords internal
.load_stat_reconstr_nmr <- function(mod, pc, X, idx, ppm) {
  t_mod <- .viz_df_helper(mod, pc, an = NA, type = 't')
  cc <- cor(t_mod$df[, 1], X)[1, ]
  cv <- cov(t_mod$df[, 1], X)[1, ]

  df <- data.frame(cor = abs(cc), cov = cv, ppm = ppm)
  df <- df[idx, ]

  return(df)
}


#' @title Validate and Format Annotation Data for Plotting
#' @description
#' Helper function to validate and format annotation metadata used for plotting PCA or OPLS scores.
#' Ensures that the input is a named list of up to three elements, and assigns default names if necessary.
#'
#' @param an A list of annotation variables (e.g., group labels or sample metadata).
#'   If missing, defaults to the original response \code{Y} from the model object (for OPLS) or a placeholder for PCA.
#' @param obj A \code{PCA_metabom8} or \code{OPLS_metabom8} object containing the model and optionally the annotation data.
#'
#' @return A cleaned and named list suitable for use with plotting functions (e.g., \code{.viz_df_helper}).
#'
#' @details
#' If annotation names are missing or duplicated, default names are assigned. Only non-\code{NULL} entries are retained.
#' A message is printed if more than three annotation variables are supplied.
#'
#' @examples
#' # .check_an_viz(list(Group = c("A", "B")), model)
#'
#' @keywords internal
.check_an_viz <- function(an, obj) {
  if (missing(an)) {
    if (inherits(obj, "OPLS_metabom8")) {
      an <- list(Y = obj@Y$ori)
    } else {
      an <- list("All samples")
    }
  }

  if (length(an) > 3) {
    message("Annotation list `an` should contain a maximum of three elements.")
  }

  # remove NULLs
  an <- an[!vapply(an, is.null, logical(1))]

  # assign default names if none exist
  if (is.null(names(an))) {
    names(an) <- paste0("Var", seq_along(an))
  }

  # fix duplicated names
  if (any(duplicated(names(an)))) {
    dup_names <- duplicated(names(an))
    names(an)[dup_names] <- paste0(names(an)[dup_names], seq_len(sum(dup_names)))
  }

  # fix missing names
  missing_names <- is.na(names(an)) | names(an) == ""
  if (any(missing_names)) {
    names(an)[missing_names] <- paste0("Var", seq_len(sum(missing_names)))
  }

  return(an)
}


#' @title OPLS Y-permutation Modeling
#' @description
#' Performs OPLS modeling on permuted Y to estimate cross-validated model performance under the null hypothesis.
#' Extracts predictive performance (R2, Q2, AUROC) from cross-validation for either regression or classification.
#'
#' @param Xs Numeric matrix. Input data matrix (samples x features).
#' @param Y Numeric matrix. Response variable (numeric or dummy-coded).
#' @param cv List. Cross-validation parameters, including \code{method} and \code{cv_sets}.
#' @param type Character. Model type: \code{"DA"}, \code{"R"}, possibly with \code{"-mY"} suffix for multi-column Y.
#' @param nc_o Integer. Number of orthogonal components to be fitted.
#'
#' @return List with elements:
#' \itemize{
#'   \item \code{r2_comp}: Numeric, R2 statistic for test data.
#'   \item \code{q2_comp}: Numeric, Q2 statistic from cross-validation.
#'   \item \code{aucs_tr}: Numeric, AUROC on training data (for classification only).
#'   \item \code{aucs_te}: Numeric, AUROC on test data (for classification only).
#' }
#'
#' @details
#' When \code{nc_o > 1}, the function fits additional orthogonal components sequentially, updating the model object.
#' Classification performance is computed using AUROC (via \code{pROC}). Regression performance uses R2 and Q2.
#'
#' @importFrom pROC roc multiclass.roc
#' @keywords internal
.permYmod <- function(Xs, Y, cv, type, nc_o) {

  if(is.vector(Y) && is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
  }

  for (i in seq_len(nc_o)) {
    if (i == 1) {
      tt <- .oplsComponentCv(Xs, Y = Y, cv$cv_sets, nc_o[i], mod.cv = NULL)
    } else {
      tt <- .oplsComponentCv(X = NA, Y = Y, cv$cv_sets, nc_o[i], mod.cv = tt)
    }
  }

  # Extract CV-based predictions
  preds_test <- .extMeanCvFeat(cv_obj = tt, feat = "y_pred_test", cv_type = cv$method, model_type = type)
  preds_train <- .extMeanCvFeat(tt, feat = "y_pred_train", cv_type = cv$method, model_type = type)

  r2_comp <- q2_comp <- aucs_tr <- aucs_te <- numeric(1)
  nc <- 1
  tssy <- .tssRcpp(Y) / ncol(Y)

  cv_method <- strsplit(cv$method, "_")[[1]][1]




  # Performance evaluation
  if (cv_method == 'MC') {
    if (grepl('DA', type)) {
      if (grepl("mY", type)) {
        Y_vec <- factor(apply(Y, 1, function(x, nam=colnames(Y)){nam[which(x==max(x))]}))
        colnames(preds_test) <- colnames(preds_train) <- colnames(Y)
        auc_te <- multiclass.roc(response = Y_vec, predictor = preds_test, quiet = TRUE)$auc
        auc_tr <- multiclass.roc(response = Y_vec, predictor = preds_train, quiet = TRUE)$auc
      } else {
        Y_vec <- as.factor(apply(Y, 1, function(x) {x>0}))
        auc_te <- roc(response = Y_vec, predictor = as.vector(preds_test), quiet = TRUE)$auc
        auc_tr <- roc(response = Y_vec, predictor = as.vector(preds_train), quiet = TRUE)$auc
      }
      list(q2 = NA, r2 = NA, aucs_te = auc_te, aucs_tr = auc_tr)
    } else {
      list(q2 = .r2(Y, preds_test, NULL), r2 = .r2(Y, preds_train, NULL), aucs_te = NA, aucs_tr = NA)
    }
  } else {
    if (grepl('DA', type)) {
      if (grepl("mY", type)) {
        Y_vec <- factor(apply(Y, 1, function(x, nam=colnames(Y)){nam[which(x==max(x))]}))
        colnames(preds_test) <- colnames(preds_train) <- class_memb
        auc_te <- multiclass.roc(response = Y_vec, predictor = preds_test, quiet = TRUE)$auc
        auc_tr <- multiclass.roc(response = Y_vec, predictor = preds_train, quiet = TRUE)$auc
      } else {
        Y_vec <- as.factor(apply(Y, 1, function(x) {x>0}))
        auc_te <- roc(response = Y_vec, predictor = preds_test, quiet = TRUE)$auc
        auc_tr <- roc(response = Y_vec, predictor = as.vector(preds_train), quiet = TRUE)$auc
      }
      list(q2 = NA, r2 = NA, aucs_te = auc_te, aucs_tr = auc_tr)
    } else {
      list(q2 = .r2(Y, preds_test, NULL), r2 = .r2(Y, as.vector(preds_train), NULL), aucs_te = NA, aucs_tr = NA)
    }
  }

  # if (cv_method == "MC") {
  #   if (grepl("DA", type)) {
  #     if (grepl("mY", type)) {
  #       pred_mean <- preds_test[1, , ]
  #       colnames(pred_mean) <- colnames(Y)
  #       mod <- multiclass.roc(response = factor(Y), predictor = pred_mean)
  #       aucs_te[nc] <- mod$auc
  #       pred_tr_mean <- preds_train[1, , ]
  #       colnames(pred_tr_mean) <- colnames(Y)
  #       mod <- multiclass.roc(response = factor(Y), predictor = pred_tr_mean)
  #       aucs_tr[nc] <- mod$auc
  #     } else {
  #       mod <- roc(response = Y, predictor = preds_test[1, ], quiet = TRUE)
  #       aucs_te[nc] <- mod$auc
  #       mod <- roc(response = Y, predictor = preds_train[1, ], quiet = TRUE)
  #       aucs_tr[nc] <- mod$auc
  #     }
  #   } else {
  #     if (grepl("mY", type)) {
  #       r2_comp[nc] <- .r2(Y, preds_test[1, , ], NULL)
  #       q2_comp[nc] <- .r2(Y, preds_test[1, , ], tssy)
  #     } else {
  #       r2_comp[nc] <- .r2(Y, preds_test[1, ], NULL)
  #       q2_comp[nc] <- .r2(Y, preds_test[1, ], tssy)
  #     }
  #   }
  # } else if (cv_method == "k-fold") {
  #   if (grepl("DA", type)) {
  #     if (grepl("mY", type)) {
  #       colnames(preds_test) <- colnames(Y)
  #       mod <- multiclass.roc(response = Y, predictor = apply(preds_test, 2, as.numeric))
  #       aucs_te[nc] <- mod$auc
  #       preds_te <- preds_train[1, , ]
  #       colnames(preds_te) <- colnames(Y)
  #       mod <- multiclass.roc(response = Y, predictor = preds_te, quiet = TRUE)
  #       aucs_tr[nc] <- mod$auc
  #     } else {
  #       mod <- roc(response = as.vector(Y), predictor = as.vector(preds_test), quiet = TRUE)
  #       aucs_te[nc] <- mod$auc
  #       mod <- roc(response = as.vector(Y), predictor = preds_train[1, ], quiet = TRUE)
  #       aucs_tr[nc] <- mod$auc
  #     }
  #   } else {
  #     if (grepl("mY", type)) {
  #       r2_comp[nc] <- .r2(Y, preds_test[1, , ], NULL)
  #       q2_comp[nc] <- .r2(Y, preds_test[4, , ], tssy)
  #     } else {
  #       r2_comp[nc] <- .r2(Y, t(t(preds_train[1,])), NULL)
  #       q2_comp[nc] <- .r2(Y, as.array(preds_test), tssy)
  #
  #       # r2_comp[nc] <- .r2(Y, preds_test[1, ], NULL)
  #       # q2_comp[nc] <- .r2(Y, preds_test[4, ], tssy)
  #     }
  #   }
  # }

  # return(list(r2_comp = r2_comp, q2_comp = q2_comp, aucs_tr = auc_tr, aucs_te = auc_te))
}


#' @title Check if input is numeric-like using tryCatch
#' @description
#' Tests whether a vector \code{x} can be safely coerced to numeric without
#' producing \code{NA}s (excluding original \code{NA}s) or warnings.
#'
#' This function attempts to coerce the input to numeric and returns \code{TRUE}
#' if all non-missing elements convert without coercion failure; otherwise \code{FALSE}.
#' It catches warnings and errors during coercion and treats those as non-numeric.
#' @param x A vector to test for numeric coercion.
#' @return Logical \code{TRUE} if \code{x} is numeric-like, \code{FALSE} otherwise.
#' @keywords internal
.is_numeric_trycatch <- function(x) {
  tryCatch({
    num_x <- as.numeric(x)
    # Check if coercion produced NA for non-NA inputs
    if (any(is.na(num_x) & !is.na(x))) {
      FALSE
    } else {
      TRUE
    }
  }, warning = function(w) {
    # On warning -> not numeric
    FALSE
  }, error = function(e) {
    # On error -> not numeric
    FALSE
  })
}

