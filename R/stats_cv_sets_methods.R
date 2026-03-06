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
    stop("Invalid k value.")
  }

  if (k != n && n / k < 2) {
    message("Fold size is small; consider using fewer folds.")
  }

  sets <- sample(rep_len(seq_len(k), n))
  sets_list <- lapply(seq_len(k), function(i) {
    which(sets != i)
  })

  sets_list
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
#' @importFrom stats quantile
#' @keywords internal
.kFoldStratified <- function(k, stratified){

  if (!is.list(stratified) || length(stratified) != 3) {
    stop("'stratified' must be a list of length 3: type, Y, probs.")
  }

  Y <- stratified[[2]]
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
    warning("Skewed class distribution detected. Consider using non-stratified k-fold CV.")
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

  sets_list
}



#' @title Generate Monte Carlo Cross-Validation (MCCV) Training Indices
#' @description Generates a list of training set indices for Monte Carlo Cross-Validation.
#' @param k Integer. Number of MCCV iterations (i.e., training sets to generate).
#' @param Y A matrix (n x p), where rows are observations and columns are outcomes. Only the number of rows is used here.
#' @param split Numeric. Fraction of the data to include in each training set. Should be between 0 and 1.
#' @return A list of length \code{k}, each containing a vector of training set indices.
#' @details For each fold, a random sample (without replacement) of size \code{floor(nrow(Y) * split)} is drawn.
#' @keywords internal
.mc <- function(k, Y, split) {

  k <- ceiling(k)

  if (k <= 0 || k > 1e6 || is.na(k) || is.infinite(k)) {
    stop("Invalid value for 'k'. Must be a positive integer and <= 1e6.")
  }

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
    sample(seq_len(n), size = n_train, replace = FALSE)
  })

  sets_list
}


#' @title Class-balanced Monte Carlo Cross-Validation (MCCV) splits
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
#' @details Each round samples a class-balanced subset without replacement. Useful for high-variance, imbalanced-class modeling.
#'
#' @importFrom stats quantile
#' @keywords internal
.mcBalanced <- function(k, split, stratified) {

  k <- ceiling(k)

  if (k <= 1 || k > 1e6 || is.na(k) || is.infinite(k)) {
    stop("Invalid 'k' value. Must be between 2 and 1e6.")
  }

  if (split <= 0 || split >= 1 || is.na(split) || is.infinite(split)) {
    stop("Split must be between 0 and 1 (exclusive).")
  }
  if (split < 0.3 || split > 0.9) {
    message("Unusual value for split. Typical values are between 0.5 and 0.9.")
  }

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

  if (min(tabY) / max(tabY) < 0.1) {
    stop("Skewed class frequencies. Results will be biased. Consider using non-stratified CV.")
  }

  n_obs <- floor(min(tabY) * split)
  if (n_obs == 0) {
    stop("Insufficient observations in minority class. Consider increasing split or aggregating classes.")
  }

  sets_list <- lapply(seq_len(k), function(i) {
    idx_k <- unlist(lapply(names(tabY), function(class) {
      sample(which(Y == class), size = n_obs, replace = FALSE)
    }))
    return(idx_k)
  })

  sets_list
}

#' @title Class-balanced resampling (with replacement)
#'
#' @description
#' Generates k training index sets. Total train size per set is floor(N * split),
#' where N is the total number of samples. Each stratum contributes (approximately)
#' equally, sampling WITH replacement within each stratum.
#'
#' @param k Integer. Number of resampling sets.
#' @param split Numeric in (0,1). Fraction of total N used as training size.
#' @param stratified List of length 3:
#'   1) type: "DA" for classification or "R" for regression
#'   2) Y: matrix with one column (n x 1)
#'   3) probs: numeric vector of quantile probs for regression binning (e.g. c(0, .33, .66, 1))
#' @param remainder Character. How to handle n_train %% G remainder:
#'   - "drop": keep perfect balance, total size becomes G * floor(n_train/G)
#'   - "distribute": distribute leftover 1-by-1 to random strata (still near-balanced)
#'
#' @return List of length k. Each element is an integer vector of training indices.
#' @keywords internal
.balanced_bootstrap_mc <- function(k, split, stratified, remainder = c("distribute", "drop")) {

  remainder <- match.arg(remainder)

  k <- ceiling(k)
  if (!is.numeric(k) || is.na(k) || is.infinite(k) || k <= 0 || k > 1e6) {
    stop("Invalid 'k'. Must be a positive integer <= 1e6.")
  }

  if (!is.numeric(split) || is.na(split) || is.infinite(split) || split <= 0 || split >= 1) {
    stop("Invalid 'split'. Must be between 0 and 1 (exclusive).")
  }

  if (!is.list(stratified) || length(stratified) != 3) {
    stop("'stratified' must be a list of length 3: type, Y, probs.")
  }

  type  <- stratified[[1]]
  Ymat  <- as.matrix(stratified[[2]])
  probs <- stratified[[3]]

  if (ncol(Ymat) != 1) stop("Y must be a single-column matrix.")
  if (anyNA(Ymat)) stop("Y contains NA values; cannot stratify.")

  if (grepl("R", type)) {
    if (!is.numeric(probs) || length(probs) < 2) stop("For regression, 'probs' must be a numeric vector of quantile probs.")
    breaks <- unique(stats::quantile(Ymat[, 1], probs = probs, na.rm = TRUE))
    if (length(breaks) <= 2) stop("Insufficient variability in Y to create stratified bins.")
    Ygrp <- cut(Ymat[, 1], breaks = breaks, include.lowest = TRUE)
  } else {
    Ygrp <- Ymat[, 1]
  }

  N <- length(Ygrp)
  n_train <- floor(N * split)
  if (n_train < 1) stop("Split results in zero training samples.")

  tab <- table(Ygrp)
  groups <- names(tab)
  G <- length(groups)
  if (G < 2) stop("Need at least 2 strata/classes to balance.")

  n_per <- floor(n_train / G)
  if (n_per < 1) stop("Training size too small to sample at least 1 per group.")

  extra <- n_train - (n_per * G)  # remainder

  sets_list <- lapply(seq_len(k), function(i) {
    alloc <- rep(n_per, G)
    names(alloc) <- groups

    if (extra > 0 && remainder == "distribute") {
      add_to <- sample(groups, size = extra, replace = TRUE)
      alloc[add_to] <- alloc[add_to] + 1
    }

    unlist(lapply(groups, function(g) {
      idx_g <- which(Ygrp == g)
      sample(idx_g, size = alloc[[g]], replace = TRUE)
    }), use.names = FALSE)
  })

  sets_list
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
#' - If `Y` has multiple columns, only the first column is used to create resampling sets.

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

  sets_list
}
