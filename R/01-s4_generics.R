
# --- scores -----------------------------------------------------------

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
#' scaling <- uv_scaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' show(model)
#' scores(model, orth=FALSE)
#' scores(model, orth=TRUE)
#' scores(model, cv=TRUE)
#' @export
setGeneric("scores",
           function(object, ...)
             standardGeneric("scores"))


# --- vips -----------------------------------------------------------

#' Variable Importance in Projection (VIP)
#' @param object An object of class \code{m8_model}.
#' @return Numeric vector or matrix containing the variable importance in
#' projection (VIP) scores.
#' @examples
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- uv_scaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' show(model)
#' vip_scores <- vip(model)
#' dim(vip_scores)
#' @export
setGeneric("vip", function(object)
  standardGeneric("vip"))


# --- weights -----------------------------------------------------------

#' Extract model weights
#' @rdname weights
#' @param object An object of class \code{m8_model}.
#' @param ... Additional arguments (currently ignored).
#' @return Numeric vector or matrix containing model weights.
#' @examples
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- uv_scaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' show(model)
#' W <- weights(model)
#' Wo <- weights(model, orth = TRUE)
#' dim(W)
#' dim(Wo) == dim(W)
#' @export
setGeneric("weights", function(object, ...)
  standardGeneric("weights"))


# --- yfitted -----------------------------------------------------------

#' Extract fitted Y values
#' @param object An object of class \code{m8_model}.
#' @param ... Additional arguments (currently ignored).
#' @return Numeric vector or matrix containing the fitted response values.
#' @examples
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- uv_scaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' show(model)
#' Y_hat_dummy <- fitted(model)
#' @export
setGeneric("fitted", function(object, ...)
  standardGeneric("fitted"))


# --- xres -----------------------------------------------------------

#' Compute X residual matrix
#' Returns the residual matrix (E) of an OPLS model.
#' @param object An object.
#' @return Numeric matrix containing the X residuals.
#' @examples
#' data(covid)
#' cv <- balanced_mc(k = 5, split = 2/3)
#' scaling <- uv_scaling(center = TRUE)
#' model <- opls(X = covid$X, Y = covid$an$type, scaling, cv)
#' X_res <- xres(model)
#' dim(X_res) == dim(covid$X)
#' @export
setGeneric("xres", function(object) standardGeneric("xres"))
