#' @title OPLS model object from metabom8
#' @description An S4 class to represent an OPLS model constructed with metabom8.
#' @slot type Character string indicating model type (e.g., "OPLS")
#' @slot t_pred Predictive component scores (matrix)
#' @slot p_pred Predictive loadings (matrix)
#' @slot w_pred Predictive weights (matrix)
#' @slot betas_pred Regression coefficients (numeric)
#' @slot Qpc Predicted Y values per component (matrix)
#' @slot t_pred_cv Cross-validated predictive scores (matrix)
#' @slot t_orth_cv Cross-validated orthogonal scores (matrix)
#' @slot t_orth Orthogonal scores (matrix)
#' @slot p_orth Orthogonal loadings (matrix)
#' @slot w_orth Orthogonal weights (matrix)
#' @slot nPC Number of predictive components (numeric)
#' @slot summary Model performance metrics (data.frame)
#' @slot X_orth Orthogonal part of X (matrix)
#' @slot Y_res Residuals of Y (matrix)
#' @slot X_res Residuals of X (matrix)
#' @slot X_mean Mean of each feature in X (numeric)
#' @slot X_sd Standard deviation of each feature in X (numeric)
#' @slot Y_mean Mean of Y (numeric)
#' @slot Y_sd Standard deviation of Y (numeric)
#' @slot Parameters List of parameters used during training
#' @slot X Original input data (matrix)
#' @slot X_scaled Scaled version of X (matrix)
#' @slot Y Target variable (list)
#' @author \email{tkimhofer@@gmail.com}
#' @family NMR & MS
#' @export
setClass("OPLS_metabom8", representation(
  type = "character",
  t_pred = "matrix",
  p_pred = "matrix",
  w_pred = "matrix",
  betas_pred = "numeric",
  Qpc = "matrix",
  t_pred_cv = "matrix",
  t_orth_cv = "matrix",
  t_orth = "matrix",
  p_orth = "matrix",
  w_orth = "matrix",
  nPC = "numeric",
  summary = "data.frame",
  X_orth = "matrix",
  Y_res = "matrix",
  X_res = "matrix",
  X_mean = "numeric",
  X_sd = "numeric",
  Y_mean = "numeric",
  Y_sd = "numeric",
  Parameters = "list",
  X = "matrix",
  X_scaled = "matrix",
  Y = "list"
))

# Define generics and methods for OPLS_metabom8 slots


setGeneric("summary", function(object) standardGeneric("summary"))
setMethod("summary", "OPLS_metabom8", function(object) {
  object@summary
})

#
# # Generic and method for 'type'
# setGeneric("type", function(object) standardGeneric("type"))
# setMethod("type", "OPLS_metabom8", function(object) {
#   object@type
# })
#
# # Generic and method for 't_pred'
# setGeneric("t_pred", function(object) standardGeneric("t_pred"))
# setMethod("t_pred", "OPLS_metabom8", function(object) {
#   object@t_pred
# })
#
# # Generic and method for 'p_pred'
# setGeneric("p_pred", function(object) standardGeneric("p_pred"))
# setMethod("p_pred", "OPLS_metabom8", function(object) {
#   object@p_pred
# })

#
# # Helper function to define generic and method for a slot
# makeAccessor <- function(slotName, genericName = slotName) {
#   setGeneric(genericName, eval(substitute(
#     function(object) standardGeneric(genericName),
#     list(genericName = as.name(genericName))
#   )))
#   setMethod(genericName, "OPLS_metabom8", function(object) slot(object, slotName))
# }
#
# # Accessors for all slots
# makeAccessor("type")
# makeAccessor("t_pred")
# makeAccessor("p_pred")
# makeAccessor("w_pred")
# makeAccessor("betas_pred")
# makeAccessor("Qpc")
# makeAccessor("t_pred_cv")
# makeAccessor("t_orth_cv")
# makeAccessor("t_orth")
# makeAccessor("p_orth")
# makeAccessor("w_orth")
# makeAccessor("nPC")
# makeAccessor("summary")
# makeAccessor("X_orth")
# makeAccessor("Y_res")
# makeAccessor("X_res")
# makeAccessor("X_mean")
# makeAccessor("X_sd")
# makeAccessor("Y_mean")
# makeAccessor("Y_sd")
# makeAccessor("Parameters")
# makeAccessor("X")
# makeAccessor("X_scaled")
# makeAccessor("Y")


#' @title PLS model object from metabom8
#' @description An S4 class to represent a PLS model constructed with metabom8.
#' @slot type Character string indicating model type (e.g., "PLS")
#' @slot t_pred Predictive component scores (matrix)
#' @slot p_pred Predictive loadings (matrix)
#' @slot w_pred Predictive weights (matrix)
#' @slot betas_pred Regression coefficients (numeric)
#' @slot Qpc Predicted Y values per component (matrix)
#' @slot t_pred_cv Cross-validated predictive scores (matrix)
#' @slot nPC Number of predictive components (numeric)
#' @slot summary Model performance metrics (data.frame)
#' @slot Y_res Residuals of Y (matrix)
#' @slot X_res Residuals of X (matrix)
#' @slot X_mean Mean of each feature in X (numeric)
#' @slot X_sd Standard deviation of each feature in X (numeric)
#' @slot Y_mean Mean of Y (numeric)
#' @slot Y_sd Standard deviation of Y (numeric)
#' @slot Parameters List of parameters used during training
#' @slot X Original input data (matrix)
#' @slot X_scaled Scaled version of X (matrix)
#' @slot Y Target variable (list)
#' @author \email{tkimhofer@@gmail.com}
#' @family NMR & MS
#' @export
setClass("PLS_metabom8", representation(
  type = "character",
  t_pred = "matrix",
  p_pred = "matrix",
  w_pred = "matrix",
  betas_pred = "numeric",
  Qpc = "matrix",
  t_pred_cv = "matrix",
  nPC = "numeric",
  summary = "data.frame",
  Y_res = "matrix",
  X_res = "matrix",
  X_mean = "numeric",
  X_sd = "numeric",
  Y_mean = "numeric",
  Y_sd = "numeric",
  Parameters = "list",
  X = "matrix",
  X_scaled = "matrix",
  Y = "list"
))

#' @title STOCSY model object from metabom8
#' @description An S4 class to represent a STOCSY model.
#' @slot version Character string indicating version of the model
#' @slot X Spectral data matrix
#' @slot ppm Chemical shift values (vector)
#' @slot driver Numeric vector indicating driver variable
#' @slot r Correlation coefficients (vector)
#' @slot cov Covariance values (vector)
#' @author \email{tkimhofer@@gmail.com}
#' @family NMR
#' @export
setClass("stocsy1d_metabom8", representation(
  version = "character",
  X = "matrix",
  ppm = "vector",
  driver = "numeric",
  r = "vector",
  cov = "vector"
))

#' @title PCA model object from metabom8
#' @description An S4 class to represent a PCA model constructed using metabom8.
#' @slot type Character string indicating model type (e.g., "PCA")
#' @slot t Scores matrix (PCs)
#' @slot p Loadings matrix
#' @slot nPC Number of principal components
#' @slot X_mean Mean of each feature in X (numeric)
#' @slot X_sd Standard deviation of each feature in X (numeric)
#' @slot Parameters List of parameters used
#' @slot X Original input matrix
#' @author \email{tkimhofer@@gmail.com}
#' @family NMR & MS
#' @export
setClass("PCA_metabom8", representation(
  type = "character",
  t = "matrix",
  p = "matrix",
  nPC = "numeric",
  X_mean = "numeric",
  X_sd = "numeric",
  Parameters = "list",
  X = "matrix"
))


# Define generic and methods for each accessor

setGeneric("type", function(object) standardGeneric("type"))
setMethod("type", "PCA_metabom8", function(object) {
  object@type
})

setGeneric("scores", function(object) standardGeneric("scores"))
setMethod("scores", "PCA_metabom8", function(object) {
  object@t
})
#
# setGeneric("loadings", function(object) standardGeneric("loadings"))
# setMethod("loadings", "PCA_metabom8", function(object) {
#   object@p
# })
#
# setGeneric("numPC", function(object) standardGeneric("numPC"))
# setMethod("numPC", "PCA_metabom8", function(object) {
#   object@nPC
# })
#
# setGeneric("meanX", function(object) standardGeneric("meanX"))
# setMethod("meanX", "PCA_metabom8", function(object) {
#   object@X_mean
# })
#
# setGeneric("sdX", function(object) standardGeneric("sdX"))
# setMethod("sdX", "PCA_metabom8", function(object) {
#   object@X_sd
# })
#
# setGeneric("parameters", function(object) standardGeneric("parameters"))
# setMethod("parameters", "PCA_metabom8", function(object) {
#   object@Parameters
# })
#
# setGeneric("originalX", function(object) standardGeneric("originalX"))
# setMethod("originalX", "PCA_metabom8", function(object) {
#   object@X
# })
