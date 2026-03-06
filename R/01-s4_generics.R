
#####################
### PREPROCESSING ###
#####################

#' Applies a preprocessing strategy to a numeric matrix.
#' @param object A preprocessing strategy object.
#' @param X Numeric matrix.
#' @return A list containing the processed matrix and parameters.
#' @examples
#' uv <- UVScaling(center=TRUE)
#' X <- matrix(c(10,10, 0,0, 0, 10, 0, 1000), ncol=4)
#' X_scaled <- prep(uv, X)
#' str(X_scaled)
#' @export
setGeneric("prep", function(object, X) standardGeneric("prep"))


###################################################
### MODEL VALIDATION - RESAMPLING INSTANTIATION ###
###################################################

setGeneric(".rngPreflight",
           function(object) standardGeneric(".rngPreflight"))


#' Instantiate a resampling strategy
#'
#' Generates a concrete resampling instance from a resampling strategy.
#' The resulting object contains the training set indices for each
#' resampling split and records the random seed used during generation.
#'
#' @param object A resampling strategy object (e.g. \code{MonteCarlo},
#'   \code{BalancedMonteCarlo}, \code{KFold}, or related classes).
#' @param Y Outcome data used to derive resampling splits. Typically a
#'   numeric vector or matrix with observations in rows. Some strategies
#'   (e.g. balanced or stratified variants) use \code{Y} to preserve
#'   class proportions or outcome distribution.
#'
#' @return An object of class \code{ResamplingInstance}. This object
#'   stores the training indices for each resampling iteration, the
#'   originating strategy, the number of observations, and the RNG seed
#'   used to generate the splits.
#'
#' @details
#' Instantiation converts an abstract resampling strategy into concrete
#' training index sets. The exact behaviour depends on the resampling strategy.
#'
#' The RNG seed is stored in the resulting instance to ensure that
#' resampling can be reproduced if needed.
#'
#' @family resampling strategies
#' @examples
#' n <- 100
#' thr <- 1.5
#' Y <- c(rnorm(80, thr - 3, 0.3), rnorm(20, thr + 3, 0.3))  # unbalanced outcome
#' cv_k  <- kfold(k = 10)
#' k_inst <- instantiate(cv_k,  matrix(Y, ncol = 1))
#' sapply(train_idc(k_inst), function(i) length(i))
#'
#' @export
setGeneric("instantiate",
           function(object,Y)
             standardGeneric("instantiate"))

#' Extract training indices from a resampling instance
#' @details Returned is a list of training set indices stored in a
#' \code{ResamplingInstance} object.
#' @param object A \code{ResamplingInstance}.
#' @return A list where each element contains the row indices
#'   corresponding to the training set of a resampling split.
#' @examples
#' n <- 100
#' thr <- 1.5
#' Y <- c(rnorm(80, thr - 3, 0.3), rnorm(20, thr + 3, 0.3))  # unbalanced outcome
#' cv_k  <- kfold(k = 10)
#' k_inst <- instantiate(cv_k,  matrix(Y, ncol = 1))
#' idc <- train_idc(k_inst)
#' head(idc)
#' @export
setGeneric("train_idc", function(object) standardGeneric("train_idc"))


######################
### MODEL FEATURES ###
######################

#' PLS/OPLS model scores
#' @param cv Logical indicating whether cross-validated scores should be returned
#' @param orth Logical indicating whether orthogonal scores should be returned
#'   (only applicable for OPLS models).
#' @param object An object of class \code{m8_model}.
#' @param ... Additional arguments (currently ignored).
#' @return Numeric vector or matrix containing scores.
#' @examples
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- UVScaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' show(model)
#' scores(model, orth=FALSE)
#' scores(model, orth=TRUE)
#' scores(model, cv=TRUE)
#' @export
setGeneric("scores",
           function(object, ...)
             standardGeneric("scores"))

#' Variable Importance in Projection (VIP)
#' @param object An object of class \code{m8_model}.
#' @return Numeric vector or matrix containing the variable importance in
#' projection (VIP) scores.
#' @examples
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- UVScaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' show(model)
#' vip_scores <- vip(model)
#' dim(vip_scores)
#' @export
setGeneric("vip", function(object)
  standardGeneric("vip"))

#' Extract model weights
#' @rdname weights
#' @param object An object of class \code{m8_model}.
#' @param ... Additional arguments (currently ignored).
#' @return Numeric vector or matrix containing model weights.
#' @examples
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- UVScaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' show(model)
#' W <- weights(model)
#' Wo <- weights(model, orth = TRUE)
#' dim(W)
#' dim(Wo) == dim(W)
#' @export
setGeneric("weights", function(object, ...)
  standardGeneric("weights"))

#' Extract fitted Y values
#' @param object An object of class \code{m8_model}.
#' @param ... Additional arguments (currently ignored).
#' @return Numeric vector or matrix containing the fitted response values.
#' @examples
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- UVScaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' show(model)
#' Y_hat_dummy <- fitted(model)
#' @export
setGeneric("fitted", function(object, ...)
  standardGeneric("fitted"))

#' Compute X residual matrix
#' Returns the residual matrix (E) of an OPLS model.
#' @param object An object.
#' @return Numeric matrix containing the X residuals.
#' @examples
#' data(covid)
#' cv <- balanced_mc(k = 5, split = 2/3)
#' scaling <- UVScaling(center = TRUE)
#' model <- opls(X = covid$X, Y = covid$an$type, scaling, cv)
#' X_res <- xres(model)
#' dim(X_res) == dim(covid$X)
#' @export
setGeneric("xres", function(object) standardGeneric("xres"))
