#' @title An S4 class to represent an OPLS model constructed with MetaboMate
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @section


# define slots for OPLS_Torben object
setClass("OplsMate", representation(
  type = "character",
  t_pred = "matrix",

  p_pred = "matrix",
  w_pred = "matrix",

  betas_pred = "numeric",
  Qpc = "matrix",

  t_cv = "matrix",
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
  Y = "list")
  )




#' @title Plotting function OplsMAte object
#' @param x OplsMAte object
#' @param y string, describing plot type: scores, scores_cv
#' @return graphics output
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @examples
#' data(iris)
#' smod=opls(X=iris[,seq_len(4)], Y=iris$Species=='setosa')
#' plot(smod, 'scores_cv')
#' @section

setMethod("plot", "OplsMate",
          function(x, y)
          {
          if(missing(y)) y='scores'
          if(y=='scores')plot(x@t_pred, x@t_orth, col=factor(x@Y$ori))
          if(y=='scores_cv') plot(x@t_cv, x@t_orth_cv, col=factor(x@Y$or))
          }
)