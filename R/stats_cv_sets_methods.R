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
    stop("Invalid k value.", call. = FALSE)
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
    stop("'stratified' must be a list of length 3: type, Y, probs.", call. = FALSE)
  }

  Y <- stratified[[2]]
  if (anyNA(Y)) {
    stop("Stratified CV not possible: Y contains NA values.", call. = FALSE)
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
  if (nobs == 0) stop("Too few observations per class for stratified CV.", call. = FALSE)

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
    stop("Invalid value for 'k'. Must be a positive integer and <= 1e6.", call. = FALSE)
  }

  if (split <= 0 || split >= 1 || is.na(split) || is.infinite(split)) {
    stop("Invalid value for 'split'. Must be between 0 and 1 (exclusive).", call. = FALSE)
  }

  if (split > 0.9 || split < 0.3) {
    message("Unusual value for 'split'. Typical values are between 0.5 and 0.9.")
  }

  n <- nrow(Y)
  n_train <- floor(n * split)
  if (n_train < 1) {
    stop("Split results in zero training samples. Increase split or use a larger dataset.", call. = FALSE)
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
    stop("Invalid 'k' value. Must be between 2 and 1e6.", call. = FALSE)
  }

  if (split <= 0 || split >= 1 || is.na(split) || is.infinite(split)) {
    stop("Split must be between 0 and 1 (exclusive).", call. = FALSE)
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
      stop("Y has insufficient variability to create stratified bins - switch to unstratified CV.", call. = FALSE)
    }
    Y <- cut(Ymat[, 1], breaks = breaks, include.lowest = TRUE)
  } else {
    Y <- Ymat[, 1]
  }

  tabY <- table(Y)

  if (min(tabY) / max(tabY) < 0.1) {
    stop("Skewed class frequencies. Results will be biased. Consider using non-stratified CV.", call. = FALSE)
  }

  n_obs <- floor(min(tabY) * split)
  if (n_obs == 0) {
    stop("Insufficient observations in minority class. Consider increasing split or aggregating classes.", call. = FALSE)
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
    stop("Invalid 'k'. Must be a positive integer <= 1e6.", call. = FALSE)
  }

  if (!is.numeric(split) || is.na(split) || is.infinite(split) || split <= 0 || split >= 1) {
    stop("Invalid 'split'. Must be between 0 and 1 (exclusive).", call. = FALSE)
  }

  if (!is.list(stratified) || length(stratified) != 3) {
    stop("'stratified' must be a list of length 3: type, Y, probs.", call. = FALSE)
  }

  type  <- stratified[[1]]
  Ymat  <- as.matrix(stratified[[2]])
  probs <- stratified[[3]]

  if (ncol(Ymat) != 1) stop("Y must be a single-column matrix.", call. = FALSE)
  if (anyNA(Ymat)) stop("Y contains NA values; cannot stratify.", call. = FALSE)

  if (grepl("R", type)) {
    if (!is.numeric(probs) || length(probs) < 2) stop("For regression, 'probs' must be a numeric vector of quantile probs.", call. = FALSE)
    breaks <- unique(stats::quantile(Ymat[, 1], probs = probs, na.rm = TRUE))
    if (length(breaks) <= 2) stop("Insufficient variability in Y to create stratified bins.", call. = FALSE)
    Ygrp <- cut(Ymat[, 1], breaks = breaks, include.lowest = TRUE)
  } else {
    Ygrp <- Ymat[, 1]
  }

  N <- length(Ygrp)
  n_train <- floor(N * split)
  if (n_train < 1) stop("Split results in zero training samples.", call. = FALSE)

  tab <- table(Ygrp)
  groups <- names(tab)
  G <- length(groups)
  if (G < 2) stop("Need at least 2 strata/classes to balance.", call. = FALSE)

  n_per <- floor(n_train / G)
  if (n_per < 1) stop("Training size too small to sample at least 1 per group.", call. = FALSE)

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
#' @param probs Numeric vector. Probabilities used for stratification of numeric \code{Y} (regression only).
#' @return A list of integer vectors. Each vector contains row indices representing a training set.
#'
#' @details
#' - In `"k-fold_stratified"` and `"MC_balanced"` modes, stratification is applied to preserve class proportions or outcome distribution.
#' - If `type` is regression (`"R"`), quantile-based binning is applied to stratify continuous Y values.
#' - If `type` is discriminant analysis (`"DA"`), class labels are used directly.
#' - If `Y` has multiple columns, only the first column is used to create resampling sets.

#' @keywords internal
.cvSetsMethod <- function(Y, type, method = "k-fold_stratified", k = 7, split = 2/3, probs=NULL) {

  valid_methods <- c("k-fold", "k-fold_stratified", "MC", "MC_balanced", "Boot_balanced")
  if (!method %in% valid_methods) {
    stop("Check method argument. Valid methods are: ", paste(valid_methods, collapse = ", "), call. = FALSE)
  }

  Y <- as.matrix(Y)
  if (ncol(Y) > 1) {
    Y <- Y[, 1, drop = FALSE]
  }

  if (k %% 1 == 0) {
    k <- as.integer(k)
  }

  if (grepl("MC", method) && (split <= 0 || split >= 1)) {
    stop("Parameter 'split' must be between 0 and 1 (exclusive).", call. = FALSE)
  }

  if (k < 2) stop("Number of folds or iterations (k) must be at least 2.", call. = FALSE)


  sets_list <- switch(method,
                      `k-fold_stratified` = .kFoldStratified(k, stratified = list(type, Y, probs=probs)),
                      `k-fold` = .kFold(k, Y),
                      MC = .mc(k, Y, split),
                      MC_balanced = .mcBalanced(k, split, stratified = list(type, Y, probs = probs)),
                      Boot_balanced = .balanced_bootstrap_mc(k, split, stratified = list(type, Y, probs = probs), remainder = "distribute")
  )

  if (is.null(sets_list)) {
    stop("Cross-validation set generation failed. Check input arguments.", call. = FALSE)
  }

  sets_list
}


#' K-fold cross-validation strategy
#' @param k Integer number of folds.
#' @details
#' Partitions the data into \code{k} folds. Each fold is used once as a test set,
#' with the remaining folds used for training.
#' No stratification is applied; folds are created by random partitioning.
#' @return A named \code{list} with elements:
#' \describe{
#'   \item{train}{List of integer vectors containing training set indices
#'   for each resampling iteration.}
#'   \item{strategy}{Character string indicating the resampling strategy.}
#'   \item{n}{Integer. Number of samples in the dataset.}
#'   \item{seed}{Integer. Random seed used to generate the resampling splits,
#'   ensuring reproducibility.}
#' }
#' @family resampling strategies
#' @examples
#' n <- 100
#' thr <- 1.5
#' Y <- c(rnorm(80, thr - 3, 0.3), rnorm(20, thr + 3, 0.3))  # unbalanced outcome
#' mean(Y > thr)
#' cv_k  <- kfold(k = 10)
#' k_inst <- metabom8:::.arg_check_cv(cv_pars=cv_k, model_type='R', n=n, Y_prepped=cbind(Y))
#' sapply(k_inst$train, function(i) length(i))
#' @export
kfold <- function(k){

  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 2L)
    stop("k must be an integer >= 2.", call. = FALSE)

  list(method = 'kfold', k = k)
}


#' Y-stratified k-fold cross-validation strategy
#'
#' @param k Integer. Number of folds.
#' @param type Character. Either \code{"DA"} (classification)
#'   or \code{"R"} (regression).
#' @param probs Numeric vector of quantile probabilities used
#'   to stratify continuous \code{Y} when \code{type = "R"}.
#'
#' @details
#' For classification (\code{type = "DA"}), folds are generated such that
#' class proportions are approximately preserved in each fold.
#'
#' For regression (\code{type = "R"}), \code{Y} is discretised into bins
#' defined by \code{probs}, and folds are generated to approximately
#' preserve the bin distribution.
#'
#' @return A named \code{list} with elements:
#' \describe{
#'   \item{train}{List of integer vectors containing training set indices
#'   for each resampling iteration.}
#'   \item{strategy}{Character string indicating the resampling strategy.}
#'   \item{n}{Integer. Number of samples in the dataset.}
#'   \item{seed}{Integer. Random seed used to generate the resampling splits,
#'   ensuring reproducibility.}
#' }
#' @family resampling strategies
#' @examples
#' set.seed(1)
#' n <- 100
#' thr <- 1.5
#' Y <- c(rnorm(80, thr - 3, 0.3), rnorm(20, thr + 3, 0.3))  # unbalanced outcome
#' mean(Y > thr)
#'
#' q80 <- quantile(Y, 0.8)  # defines the rare "high" stratum (top 20%)
#'
#' cv_k  <- kfold(k = 10)
#' cv_sk  <- stratified_kfold(k = 10, type = "R", probs = c(0, 0.8, 1))
#'
#' k_inst <-  metabom8:::.arg_check_cv(cv_pars=cv_k, model_type='R', n=n, Y_prepped=cbind(Y))
#' sk_inst <-  metabom8:::.arg_check_cv(cv_pars=cv_sk, model_type='R', n=n, Y_prepped=cbind(Y))
#'
#' round(sapply(k_inst$train,  function(i) mean(Y[i] > q80)), 2)  # reflects imbalance
#' round(sapply(sk_inst$train, function(i) mean(Y[i] > q80)), 2)  # balanced across strata
#' @export
stratified_kfold <- function(k,
                             type = c("DA", "R"),
                             probs = NULL
                             ) {

  type <- match.arg(type)

  if (type == "R") {
    if (is.null(probs)) {
      stop("'probs' must be provided when type = 'R'", call. = FALSE)
    }
  } else {
    if (!is.null(probs)) {
      warning("'probs' ignored when type = 'DA'", call. = FALSE)
    }
  }

  k <- as.integer(k)
  if (is.na(k) || k < 2L)
    stop("k must be an integer >= 2.", call. = FALSE)

  list(method = 'StratifiedKFold', k = k, type  = type, probs = probs)

}

#' Monte-Carlo cross-validation strategy
#'
#' @param k Integer. Number of repeated random splits.
#' @param split Numeric. Fraction of samples assigned to the
#'   training set (e.g. \code{2/3}).
#'
#' @details
#' Monte-Carlo cross-validation generates \code{k} random train/test splits
#' without replacement. No stratification is applied; samples are drawn
#' uniformly at random.
#'
#' @return A named \code{list} with elements:
#' \describe{
#'   \item{train}{List of integer vectors containing training set indices
#'   for each resampling iteration.}
#'   \item{strategy}{Character string indicating the resampling strategy.}
#'   \item{n}{Integer. Number of samples in the dataset.}
#'   \item{seed}{Integer. Random seed used to generate the resampling splits,
#'   ensuring reproducibility.}
#' }
#' @family resampling strategies
#' @examples
#' n <- 100
#' # bivariate outcome
#' thr <- 1.5
#' Y <- c(rnorm(80, thr-3, 0.3), rnorm(20, thr+3, 0.3))  # unbalanced low/high outcome
#' mean(Y>thr)
#'
#' cv_mc <- mc(k = 10, split = 2/3)
#' mc_inst <- metabom8:::.arg_check_cv(cv_pars=cv_mc, model_type='R', n=n, Y_prepped=cbind(Y))
#' sapply(mc_inst$train, function(i) length(i))
#' @export
mc <- function(k, split) {

  k <- as.integer(k)
  if (is.na(k) || k < 1L)
    stop("k must be an integer >= 1.", call. = FALSE)

  if (!is.numeric(split) || length(split) != 1L)
    stop("split must be a single numeric value.", call. = FALSE)

  if (split <= 0 || split >= 1)
    stop("split must be strictly between 0 and 1.", call. = FALSE)

  list(method = "MonteCarlo", k = k, split = split)
}

#' Balanced Monte-Carlo resampling strategy
#'
#' @param k Integer. Number of repeated random splits.
#' @param split Numeric. Fraction of samples assigned to the training set
#'   (e.g. \code{2/3}).
#' @param type Character. Either \code{"DA"} (classification) or \code{"R"} (regression).
#' @param probs Numeric vector of quantile probabilities used to stratify
#'   continuous \code{Y} when \code{type = "R"}.
#'
#' @details
#' Generates \code{k} Monte-Carlo resampling splits by randomly partitioning
#' the data into training and test sets without replacement.
#'
#' Balancing ensures equal representation of strata in the training data:
#' \describe{
#'   \item{\code{type = "DA"}}{
#'     Class labels define the strata, and sampling is balanced across classes.
#'   }
#'   \item{\code{type = "R"}}{
#'     The response is discretised into bins using quantiles defined by
#'     \code{probs}, and each bin contributes equally to the training set.
#'   }
#' }
#'
#' This strategy can improve robustness of model evaluation in settings with
#' limited samples size and imbalanced or unevenly distributed outcome variables.
#'
#' @return A named \code{list} with elements:
#' \describe{
#'   \item{train}{List of integer vectors containing training set indices
#'   for each resampling iteration.}
#'   \item{strategy}{Character string indicating the resampling strategy.}
#'   \item{n}{Integer. Number of samples in the dataset.}
#'   \item{seed}{Integer. Random seed used to generate the resampling splits,
#'   ensuring reproducibility.}
#' }
#' @family resampling strategies
#' @examples
#' n <- 100
#' # bivariate outcome
#' thr <- 1.5
#' Y <- c(rnorm(80, thr-3, 0.3), rnorm(20, thr+3, 0.3))  # unbalanced low/high outcome
#' mean(Y>thr)
#'
#' cv_k <- kfold(k = 10)
#' cv_mc <- balanced_mc(k = 10, split = 2/3, type = "R", probs = c(0, 0.8, 1))
#'
#' k_inst <- metabom8:::.arg_check_cv(cv_pars=cv_k, model_type='R', n=n, Y_prepped=cbind(Y))
#' mc_inst <- metabom8:::.arg_check_cv(cv_pars=cv_mc, model_type='R', n=n, Y_prepped=cbind(Y))
#'
#' # balanced splits: proportion above global median stays ~0.5
#' q80 <- quantile(Y, 0.8)
#' round(sapply(k_inst$train, function(i) mean(Y[i] > q80)), 2) # resembles original Y distr.
#' round(sapply(mc_inst$train, function(i) mean(Y[i] > q80)), 2) # balanced strata (low/high)
#' @export
balanced_mc <- function(k,
                        split,
                        type = c("DA", "R"),
                        probs = NULL
                        ) {

  type <- match.arg(type)

  if (type == "R") {
    if (is.null(probs)) {
      stop("'probs' must be provided when type = 'R'", call. = FALSE)
    }
  } else {
    if (!is.null(probs)) {
      warning("'probs' ignored when type = 'DA'", call. = FALSE)
    }
  }

  k <- as.integer(k)
  if (is.na(k) || k < 1L)
    stop("k must be an integer >= 1.", call. = FALSE)

  if (!is.numeric(split) || length(split) != 1L)
    stop("split must be a single numeric value.", call. = FALSE)

  if (split <= 0 || split >= 1)
    stop("split must be strictly between 0 and 1.", call. = FALSE)

  list(method = "BalancedMonteCarlo", k = k, split = split,
       type  = type,  probs = probs)

}


#' Balanced bootstrap resampling strategy
#'
#' @param k Integer. Number of bootstrap resamples.
#' @param split Numeric. Fraction of samples drawn for the training set
#'   (e.g. \code{2/3}). Sampling is performed with replacement.
#' @param type Character. Either \code{"DA"} (classification) or \code{"R"} (regression).
#' @param probs Numeric vector of quantile probabilities used to stratify continuous
#'   \code{Y} when \code{type = "R"}.
#'
#' @details
#' Generates \code{k} bootstrap samples (training sets sampled with replacement).
#' The remaining samples (the out-of-bag set) can be used as a test set.
#'
#' Balancing ensures equal representation of strata in the training data:
#' \describe{
#'   \item{\code{type = "DA"}}{
#'     Class labels define the strata, and sampling is balanced across classes.
#'   }
#'   \item{\code{type = "R"}}{
#'     The response is discretised into bins using quantiles defined by
#'     \code{probs}, and each bin contributes equally to the training set.
#'   }
#' }
#'
#' @return A named \code{list} with elements:
#' \describe{
#'   \item{train}{List of integer vectors containing training set indices
#'   for each resampling iteration.}
#'   \item{strategy}{Character string indicating the resampling strategy.}
#'   \item{n}{Integer. Number of samples in the dataset.}
#'   \item{seed}{Integer. Random seed used to generate the resampling splits,
#'   ensuring reproducibility.}
#' }
#' @family resampling strategies
#' @examples
#' n <- 100
#' # bivariate outcome
#' thr <- 1.5
#' Y <- c(rnorm(80, thr-3, 0.3), rnorm(20, thr+3, 0.3))  # unbalanced low/high outcome
#' mean(Y>thr)
#'
#' cv_k <- kfold(k = 10)
#' cv_boot <- balanced_boot(k = 10, split = 2/3, type = "R", probs = c(0, 0.8, 1))
#'
#' k_inst <- metabom8:::.arg_check_cv(cv_pars=cv_k, model_type='R', n=n, Y_prepped=cbind(Y))
#' b_inst <- metabom8:::.arg_check_cv(cv_pars=cv_boot, model_type='R', n=n, Y_prepped=cbind(Y))
#'
#' # balanced splits: proportion above global median stays ~0.5
#' q80 <- quantile(Y, 0.8)
#' round(sapply(k_inst$train, function(i) mean(Y[i] > q80)), 2) # resembles original Y distr.
#' round(sapply(b_inst$train, function(i) mean(Y[i] > q80)), 2) # balanced strata (low/high)
#' @export
balanced_boot <- function(k,
                          split,
                          type = c("DA", "R"),
                          probs = NULL
                          ) {
  type <- match.arg(type)

  if (type == "R") {
    if (is.null(probs)) {
      stop("'probs' must be provided when type = 'R'", call. = FALSE)
    }
  } else {
    if (!is.null(probs)) {
      warning("'probs' ignored when type = 'DA'", call. = FALSE)
    }
  }

  k <- as.integer(k)
  if (is.na(k) || k < 1L)
    stop("k must be an integer >= 1.", call. = FALSE)

  if (!is.numeric(split) || length(split) != 1L)
    stop("split must be a single numeric value.", call. = FALSE)

  if (split <= 0 || split >= 1)
    stop("split must be strictly between 0 and 1.", call. = FALSE)

  list(method = "BalancedBoot", k = k, split = split, type=type, probs=probs)
}

