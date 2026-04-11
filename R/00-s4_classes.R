
# --- Model -----------------------------------------------------------
#' @import methods
NULL

setClassUnion("listOrNULL", c("list", "NULL"))

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
#' scaling <- uv_scaling(center=TRUE)
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
           cv = "listOrNULL",
           prep = "list",# "m8_preprocess",
           provenance = "list",
           session = "list",
           call = "call",
           dims = "list"
         )
)

