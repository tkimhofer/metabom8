#####################
### PREPROCESSING ###
#####################

setClass("m8_preprocess",
         slots = c(
           center = "logical",
           scale  = "character",  # "none", "uv", "pareto"

           X_mean = "numeric",
           X_sd   = "numeric"
         )
)


setClass("PreprocessStrategy",
         contains = "VIRTUAL")


setClass("ScalingStrategy",
         contains = "PreprocessStrategy",
         slots = list(
           center = "logical",
           scale  = "character"   # "none","uv","pareto"
         ))

### constructors

#' No Scaling
#' Leaves variables unscaled. Optional centering.
#' @inheritParams UVScaling
#' @return An object of class \code{m8_preprocess}.
#' @family scaling objects & strategies
#' @examples
#' just_centering <- NoScaling(center=TRUE)
#' X <- matrix(c(10,10, 0,0, 0, 10), ncol=3)
#' X_centered <- prep(just_centering, X)
#' str(X_centered)
#' X_centered$X
#' @export
NoScaling <- function(center = FALSE)
  methods::new("m8_preprocess", center=center, scale="None")

#' Unit Variance Scaling
#' Centers variables and scales each feature to unit variance.
#' @param center Logical. If TRUE, variables are mean-centered before scaling.
#' @return An object of class \code{m8_preprocess}.
#' @details
#' UV scaling divides each variable by its standard deviation.
#' This is the default scaling in many multivariate methods such as PCA and PLS.
#' @family scaling objects & strategies
#' @examples
#' autoscale <- UVScaling(center=TRUE)
#' X <- matrix(c(10,10, 0,0, 0, 10), ncol=3)
#' X_scaled <- prep(autoscale, X)
#' str(X_scaled)
#' X_scaled$X
#' @export
UVScaling <- function(center = TRUE)
  methods::new("m8_preprocess", center=center, scale="UV")


#' Pareto Scaling
#' Scales variables by the square root of their standard deviation.
#' @inheritParams UVScaling
#' @return An object of class \code{m8_preprocess}.
#' @details
#' Pareto scaling reduces the influence of large-variance variables
#' while preserving data structure better than full unit variance scaling.
#' @family scaling objects & strategies
#' @examples
#' paritalUV <-ParetoScaling(center=TRUE)
#' X <- matrix(c(10,10, 0,0, 0, 10, 0, 1000), ncol=4)
#' X_scaled <- prep(paritalUV, X)
#' str(X_scaled)
#' X_scaled$X
#' @export
ParetoScaling <- function(center = FALSE)
  methods::new("m8_preprocess", center=center, scale="Pareto")


###################################################
### MODEL VALIDATION - RESAMPLING INSTANTIATION ###
###################################################

setClass("ResamplingStrategy",
         contains = "VIRTUAL")

setClass("ResamplingInstance",
         slots=list(
           train     = "list",
           strategy  = "ResamplingStrategy",   # generator
           n         = "numeric",              # length(Y)
           seed      = "numeric"               # RNG state
         ))


setClass("KFold",
         contains = "ResamplingStrategy",
         slots = list(
           k = "numeric"
         ))

setClass("StratifiedKFold",
         contains = "ResamplingStrategy",
         slots = list(
           k     = "numeric",
           type  = "character",  # "DA" or "R"
           probs = "numeric"     # only used for regression
         ))

setClass("MonteCarlo", # no replacement within each set ("plain holdout")
         contains = "ResamplingStrategy",
         slots = list(
           k     = "numeric",
           split = "numeric"
         ))

setClass("BalancedMonteCarlo", # this is without resampling to keep distirubitions equal
         contains = "ResamplingStrategy",
         slots = list(
           k     = "numeric",
           split = "numeric",
           type  = "character",  # "DA" or "R"
           probs = "numeric"
         ))


setClass("BalancedBootstrap",  # this is wit resampling to keep distirubitions equal
         contains = "ResamplingStrategy",
         slots = list(
           k     = "numeric",
           split = "numeric",
           type  = "character",  # "DA" or "R"
           probs = "numeric"
         ))


### constructors

#' K-fold cross-validation strategy
#' @param k Integer number of folds.
#' @return A \code{ResamplingStrategy} object.
#' @family resampling strategies
#' @examples
#' n <- 100
#' thr <- 1.5
#' Y <- c(rnorm(80, thr - 3, 0.3), rnorm(20, thr + 3, 0.3))  # unbalanced outcome
#' mean(Y > thr)
#' cv_k  <- kfold(k = 10)
#' k_inst <- instantiate(cv_k,  matrix(Y, ncol = 1))
#' sapply(train_idc(k_inst), function(i) length(i))
#' @export
kfold <- function(k){

  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 2L)
    stop("k must be an integer >= 2.")

  methods::new("KFold", k = k)
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
#' @return An object of class \code{StratifiedKFold}.
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
#' k_inst <- instantiate(cv_k,  matrix(Y, ncol = 1))
#' sk_inst <- instantiate(cv_sk, matrix(Y, ncol = 1))
#'
#' round(sapply(train_idc(k_inst),  function(i) mean(Y[i] > q80)), 2)  # reflects imbalance
#' round(sapply(train_idc(sk_inst), function(i) mean(Y[i] > q80)), 2)  # balanced across strata
#' @export
stratified_kfold <- function(k,
                             type = "DA",
                             probs = c(0, 0.33, 0.66, 1)) {

  k <- as.integer(k)
  if (is.na(k) || k < 2L)
    stop("k must be an integer >= 2.")

  type <- match.arg(type, c("DA", "R"))

  if (type == "R") {
    if (!is.numeric(probs))
      stop("probs must be numeric.")

    if (any(probs < 0 | probs > 1))
      stop("probs must be between 0 and 1.")

    if (!isTRUE(all.equal(sort(probs), probs)))
      stop("probs must be sorted in increasing order.")

    if (probs[1] != 0 || tail(probs, 1) != 1)
      stop("probs must start at 0 and end at 1.")
  }

  methods::new("StratifiedKFold",
      k     = k,
      type  = type,
      probs = probs)
}

#' Monte-Carlo cross-validation strategy
#'
#' @param k Integer. Number of repeated random splits.
#' @param split Numeric. Fraction of samples assigned to the
#'   training set (e.g. \code{2/3}).
#'
#' @details
#' Monte-Carlo cross-validation (also called repeated holdout)
#' generates \code{k} random train/test splits without replacement.
#' Each split uses the proportion defined by \code{split}
#' for the training set.
#'
#' @return An object of class \code{MonteCarlo}.
#' @family resampling strategies
#' @examples
#' n <- 100
#' # bivariate outcome
#' thr <- 1.5
#' Y <- c(rnorm(80, thr-3, 0.3), rnorm(20, thr+3, 0.3))  # unbalanced low/high outcome
#' mean(Y>thr)
#'
#' cv_mc <- mc(k = 10, split = 2/3)
#'
#' mc_inst <- instantiate(cv_mc, matrix(Y, ncol = 1))
#' sapply(train_idc(mc_inst), function(i) length(i))
#' @export
mc <- function(k, split) {

  k <- as.integer(k)
  if (is.na(k) || k < 1L)
    stop("k must be an integer >= 1.")

  if (!is.numeric(split) || length(split) != 1L)
    stop("split must be a single numeric value.")

  if (split <= 0 || split >= 1)
    stop("split must be strictly between 0 and 1.")

  methods::new("MonteCarlo",
      k     = k,
      split = split)
}


#' Balanced Monte-Carlo cross-validation strategy
#'
#' @param k Integer. Number of repeated random splits.
#' @param split Numeric. Fraction of samples assigned to the training set
#'   (e.g. \code{2/3}).
#' @param type Character. Either \code{"DA"} (classification) or \code{"R"} (regression).
#' @param probs Numeric vector of quantile probabilities used to stratify continuous
#'   \code{Y} when \code{type = "R"}.
#'
#' @details
#' Generates \code{k} repeated random train/test splits while attempting to keep the
#' class distribution (DA) or the distribution of \code{Y} bins (R) approximately
#' balanced between training and test sets.
#'
#' For regression (\code{type = "R"}), \code{Y} is discretised into bins defined by
#' \code{probs} (quantiles), and balancing is performed on these bins.
#'
#' @return An object of class \code{BalancedMonteCarlo}.
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
#' k_inst <- instantiate(cv_k, matrix(Y, ncol = 1))
#' mc_inst <- instantiate(cv_mc, matrix(Y, ncol = 1))
#'
#' # balanced splits: proportion above global median stays ~0.5
#' q80 <- quantile(Y, 0.8)
#' round(sapply(train_idc(k_inst), function(i) mean(Y[i] > q80)), 2) # resembles original Y distr.
#' round(sapply(train_idc(mc_inst), function(i) mean(Y[i] > q80)), 2) # balanced strata (low/high)
#' @export
balanced_mc <- function(k,
                        split,
                        type = "DA",
                        probs = c(0, 0.33, 0.66, 1)) {

  k <- as.integer(k)
  if (is.na(k) || k < 1L)
    stop("k must be an integer >= 1.")

  if (!is.numeric(split) || length(split) != 1L)
    stop("split must be a single numeric value.")

  if (split <= 0 || split >= 1)
    stop("split must be strictly between 0 and 1.")

  type <- match.arg(type, c("DA", "R"))

  if (type == "R") {
    if (!is.numeric(probs))
      stop("probs must be numeric.")
    if (any(probs < 0 | probs > 1))
      stop("probs must be between 0 and 1.")
    if (!isTRUE(all.equal(sort(probs), probs)))
      stop("probs must be sorted in increasing order.")
    if (probs[1] != 0 || tail(probs, 1) != 1)
      stop("probs must start at 0 and end at 1.")
  }

  methods::new("BalancedMonteCarlo",
      k     = k,
      split = split,
      type  = type,
      probs = probs)
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
#' Generates \code{k} bootstrap resamples (training sets sampled with replacement).
#' The remaining samples (often the out-of-bag set) can be used as a test set.
#'
#' Balancing attempts to preserve class proportions (DA) or the distribution of
#' binned \code{Y} (R) across resamples.
#'
#' For regression (\code{type = "R"}), \code{Y} is discretised into bins defined by
#' \code{probs} (quantiles), and balancing is performed on these bins.
#'
#' @return An object of class \code{BalancedBootstrap}.
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
#' k_inst <- instantiate(cv_k, matrix(Y, ncol = 1))
#' b_inst <- instantiate(cv_boot, matrix(Y, ncol = 1))
#'
#' # balanced splits: proportion above global median stays ~0.5
#' q80 <- quantile(Y, 0.8)
#' round(sapply(train_idc(k_inst), function(i) mean(Y[i] > q80)), 2) # resembles original Y distr.
#' round(sapply(train_idc(b_inst), function(i) mean(Y[i] > q80)), 2) # balanced strata (low/high)
#' @export
balanced_boot <- function(k,
                          split,
                          type = "DA",
                          probs = c(0, 0.33, 0.66, 1)) {

  k <- as.integer(k)
  if (is.na(k) || k < 1L)
    stop("k must be an integer >= 1.")

  if (!is.numeric(split) || length(split) != 1L)
    stop("split must be a single numeric value.")

  if (split <= 0 || split >= 1)
    stop("split must be strictly between 0 and 1.")

  type <- match.arg(type, c("DA", "R"))

  if (type == "R") {
    if (!is.numeric(probs))
      stop("probs must be numeric.")
    if (any(probs < 0 | probs > 1))
      stop("probs must be between 0 and 1.")
    if (!isTRUE(all.equal(sort(probs), probs)))
      stop("probs must be sorted in increasing order.")
    if (probs[1] != 0 || tail(probs, 1) != 1)
      stop("probs must start at 0 and end at 1.")
  }

  methods::new("BalancedBootstrap",
      k     = k,
      split = split,
      type  = type,
      probs = probs)
}


setClassUnion("ResamplingInstanceOrNULL", c("ResamplingInstance", "NULL"))



#############
### MODEL ###
#############

#' @title m8_model class
#' Model object returned by \code{pca()}, \code{pls()}, and \code{opls()}.
#' @slot engine Character. Model engine ("pca", "pls", "opls").
#' @slot ctrl List. Engine-specific control and performance information.
#' @slot fit List. Fitted data (engine specific).
#' @slot cv Resampling instance (may be NULL if not used).
#' @slot prep Scaling and centering information
#' @slot provenance Preprocessing attributes of spectral matrix X
#' @slot session R-session information
#' @slot call Function call
#' @slot dims List with \code{n} and \code{p}.
#' @return
#' An object of class `m8_model`.
#' @name m8_model-class
#' @examples
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- UVScaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' class(model)
#' show(model)
#' @aliases m8_model
#' @exportClass m8_model
setClass("m8_model",
         slots = c(
           engine = "character",
           ctrl = "list",
           fit = "list",
           cv = "ResamplingInstanceOrNULL",
           prep = "m8_preprocess",
           provenance = "list",
           session = "list",
           call = "call",
           dims = "list"
         )
)


setClass("m8_modelSummary",
         slots = list(
           perf = "data.frame",
           engine = "character",
           y_type =  "character"
         ))

