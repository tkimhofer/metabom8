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
    stop("Input Y contains NA, NaN, or Inf values.", call. = FALSE)
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
      stop("Input Y has only a single level.", call. = FALSE)
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
    stop(sprintf("Input X must be a numeric matrix. Use matrix() or as.matrix() to convert."), call. = FALSE)
  }

  if (any(is.na(X) | is.nan(X) | is.infinite(X))) {
    stop("Input X contains missing (NA), NaN, or infinite values.", call. = FALSE)
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
    stop("Dimensions of input X and Y do not match.", call. = FALSE)
  }
  if (ncol(X) <= ncol(Y)) {
    stop("Number of variables (columns) in X should be higher than in Y.", call. = FALSE)
  }
  invisible(NULL)
}

#' Evaluate model fit progression based on cross-validated performance
#'
#' Determines whether an additional component should be fitted based on
#' cross-validated generalization performance. For regression models
#' (`type = "R"`), the decision is based on the \eqn{Q^2} statistic.
#' For discriminant analysis models (`type = "DA"`), the decision is
#' based on the cross-validated AUROC (area under the ROC curve).
#'
#' AUROC is computed on cross-validation test folds to ensure that it
#' reflects out-of-sample classification performance, analogous to the
#' use of \eqn{Q^2} in regression.
#'
#' @param type Character. Model type: `"R"` for regression or `"DA"` for
#'   discriminant analysis. The value may optionally include a suffix
#'   such as `"-mY"` to indicate multi-response models.
#' @param q2s Numeric vector of \eqn{Q^2} values for the fitted
#'   components (used for regression models).
#' @param cv_auc Numeric vector of cross-validated AUROC values for
#'   the fitted components (used for classification models).
#' @param pc_max Integer. Maximum number of components allowed.
#' @param min_delta Numeric. Minimum improvement in the cross-validated
#'   performance metric required to justify fitting an additional
#'   component.
#' @param sat_q2 Numeric. Saturation threshold for \eqn{Q^2} in
#'   regression models.
#' @param sat_auc Numeric. Saturation threshold for cross-validated
#'   AUROC in discriminant analysis models.
#' @param min_q2 Numeric. Minimum acceptable \eqn{Q^2} value for the
#'   first component.
#' @param min_auc Numeric. Minimum acceptable cross-validated AUROC for
#'   the first component.
#'
#' @return A list containing:
#' \describe{
#'   \item{stop}{Logical indicating whether model fitting should stop.}
#'   \item{reason}{Character string describing the stopping condition
#'     (e.g., `"pc_max"`, `"cv_decreased"`, `"cv_improvement_negligible"`).}
#'   \item{metric}{The current cross-validated performance value.}
#'   \item{delta}{Improvement relative to the previous component.}
#' }
#'
#' @keywords internal
.evalFit <- function(type, q2s, cv_auc,
                     pc_max = 5,
                     min_delta = 0.05,
                     sat_q2  = 0.98,
                     sat_auc = 0.97,
                     min_q2    = 0.05,
                     min_auc   = 0.60
) {

  type <- strsplit(type, '-')[[1]][1]

  cv_est <- if (type == "R") q2s else cv_auc

  if (any(is.na(cv_est))) {
    stop("Something went wrong: CV metric contains NA.", call. = FALSE)
  }

  nc <- length(cv_est)

  if (type == "R" && nc ==1 && cv_est[nc] <= min_q2) {
    return(list(
      stop   = TRUE,
      reason = sprintf("q2_below_%s", min_q2),
      metric = cv_est[nc],
      delta  = if (nc > 1) cv_est[nc] - cv_est[nc - 1] else NA_real_
    ))
  }

  if (type == "DA" && nc ==1 && cv_est[nc] <= min_auc) {
    return(list(
      stop   = TRUE,
      reason = sprintf("auc_below_%s", min_auc),
      metric = cv_est[nc],
      delta  = if (nc > 1) cv_est[nc] - cv_est[nc - 1] else NA_real_
    ))
  }

  if (nc < 2) {
    return(list(
      stop   = FALSE,
      reason = NULL,
      metric = cv_est[nc],
      delta  = NA_real_
    ))
  }

  if (nc == pc_max) {
    return(list(
      stop   = TRUE,
      reason = "pc_max",
      metric = cv_est[nc],
      delta  = if (nc > 1) cv_est[nc] - cv_est[nc - 1] else NA_real_
    ))
  }

  delta <- cv_est[nc] - cv_est[nc - 1]

  thresh <- if (type == "R") sat_q2 else sat_auc

  if (delta < 0) {
    return(list(
      stop   = TRUE,
      reason = "cv_decreased",
      metric = cv_est[nc],
      delta  = delta
    ))
  }

  if (delta < min_delta) {
    return(list(
      stop   = TRUE,
      reason = "cv_improvement_negligible",
      metric = cv_est[nc],
      delta  = delta
    ))
  }

  if (cv_est[nc] > thresh) {
    return(list(
      stop   = TRUE,
      reason = "saturation",
      metric = cv_est[nc],
      delta  = delta
    ))
  }

  return(list(
    stop   = FALSE,
    reason = NULL,
    metric = cv_est[nc],
    delta  = delta
  ))
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
      Y_num <- as.numeric(Y_factor)
      mapping <- data.frame(Level = Y_levels,
                            Code = sort(unique(Y_num)),
                            stringsAsFactors = FALSE)
      return(list(matrix(Y_num, ncol = 1), mapping))

    } else {
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
    if (is.matrix(Y) && ncol(Y) > 1) {
      stop("Only single-column numeric input is supported for Y at this time.", call. = FALSE)
    }
    return(list(matrix(Y, ncol = 1), data.frame()))
  }
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
    if (any(is.na(num_x) & !is.na(x))) {
      FALSE
    } else {
      TRUE
    }
  }, warning = function(w) {
    FALSE
  }, error = function(e) {
    FALSE
  })
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
    stop("Check 'pc' argument: It must be finite and of expected length.", call. = FALSE)
  }

  mod_class <- class(mod)[1]

  if (mod_class == "PCA_metabom8") {
    if (!is.numeric(pc)) stop("PC value must be numeric for PCA.", call. = FALSE)
    if (max(pc) > ncol(mod@p)) stop("PC value exceeds number of PCA components.", call. = FALSE)
  }

  if (mod_class == "OPLS_metabom8") {
    idx_orth <- grepl("o", pc)
    if (any(idx_orth)) {
      pc1 <- as.numeric(gsub("o", "", pc))
      if (any(is.na(pc1)) || any(is.infinite(pc1))) {
        stop("Invalid orthogonal component identifier in 'pc'.", call. = FALSE)
      }
      if (any(pc1 > nrow(mod@p_orth))) {
        stop("Requested orthogonal component exceeds available components.", call. = FALSE)
      }
    } else {
      if (!is.numeric(as.numeric(pc))) {
        stop("Predictive component index should be numeric or coercible to numeric.", call. = FALSE)
      }
      if (length(mod@p_pred) == 0) {
        stop("Predictive component appears to be missing in model object.", call. = FALSE)
      }
    }
  }
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


#' @title Check Consistency of Spectral Data and ppm Vector
#'
#' @description
#' Validates that the NMR data matrix and the chemical shift reference vector
#' (\code{ppm}) are consistent.
#' Specifically, it checks that
#' \itemize{
#'   \item \code{ppm} contains no NA or infinite values.
#'   \item The number of columns in \code{X} matches the length of \code{ppm}.
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

  if (is.null(ppm) || !is.numeric(ppm)) return(FALSE)

  if (anyNA(ppm) || any(!is.finite(ppm))) return(FALSE)

  if (is.null(dim(X)) || ncol(X) != length(ppm)) return(FALSE)

  if (anyDuplicated(ppm)) return(FALSE)
  TRUE
}


#' @title Ensure Input is a Row Matrix
#'
#' @description
#' Ensures that the input object `X` is in row-matrix form. If `X` is a numeric
#' vector  (i.e., has no dimensions), it is converted to a matrix with 1 row
#' and \code{length(X)} columns. This is helpful for standardizing inputs before
#' applying matrix operations.
#'
#' @param X Numeric vector, matrix, or array. Typically NMR data where rows represent spectra.
#'
#' @return A numeric matrix with one row if input was a vector; otherwise, returns `X` unchanged.
#'
#' @keywords internal
.dimX <- function(X) {

  if (!is.numeric(X))
    stop("X must be numeric.", call. = FALSE)

  if (is.atomic(X) && is.null(dim(X))) {
    return(matrix(X, nrow = 1L))
  }

  if (is.matrix(X) || is.data.frame(X)) {
    return(as.matrix(X))
  }

  stop("X must be a numeric vector, matrix, or data.frame.", call. = FALSE)
}
