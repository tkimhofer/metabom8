#' @title An S4 class to represent an OPLS model constructed with MetaboMate
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @section


# define slots for OPLS_Torben object
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
  )
  )


#' @title An S4 class to represent STOCSY model constructed with metabom8
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @section
setClass("stocsy1d_metabom8", representation(
  version = "character",
  X = "matrix",
  ppm = "vector",
  driver ='numeric',
  r = "vector",
  cov="vector"
)
)

# @title Plotting function OplsMAte object
# @param x OplsMAte object
# @param y string, describing plot type: scores, scores_cv
# @return graphics output
# @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#@examples
# \dontrun{
# data(iris)
# smod=opls(X=iris[,seq_len(4)], Y=iris$Species=='setosa')
# plot(smod, 'scores_cv')
# }
# @section
#
# setMethod("plot", "OPLS_metabom8",
#           function(x, y)
#           {
#           if(missing(y)) y='scores'
#           if(y=='scores')plot(x@t_pred, x@t_orth, col=factor(x@Y$ori))
#           if(y=='scores_cv') plot(x@t_cv, x@t_orth_cv, col=factor(x@Y$or))
#           }
# )


#' @title An S4 class to represent an OPLS model constructed with MetaboMate
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @section


# define slots for PCA object
setClass("PCA_metabom8", representation(
  type = "character",

  t = "matrix",
  p = "matrix",

  nPC = "numeric",

  X_mean = "numeric",
  X_sd = "numeric",

  Parameters = "list",

  X = "matrix"
  )
)


# @title Plotting function PCA
# @param x PCA_metabom8 object
# @param y string, describing plot type: scores, scores_cv
# @return ggplot2 output
# @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
# @examples
# data(iris)
# mod=pca(X=iris[,seq_len(4)], Y=iris$Species=='setosa')
# \dontrun{
# plot(mod)
# }
# @section

# setMethod("plot", "PCA_metabom8",
#
#           function(x, an, pc, title='', qc=NA, legend='in', cv.scores=T, ...){
#             plotscores(model=x, pc, an, title, qc, legend, cv.scores, ...)
#           }
# )

