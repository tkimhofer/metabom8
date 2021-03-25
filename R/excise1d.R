#' @title Excise 1H  chemical shift reagions for 1D small molecule reserach
#' @export
#' @param X num matrix:  NMR data, spectra in rows
#' @param ppm ppm, num array - chemical shift positions, length matches to columns in X
#' @details
#' Function removes the following chemical shift reagions from 1H NMR spectra acquired for small molecule research:
#' \itemize{
#'   \item min(ppm) to 0.25 ppm  (no signal upfield)
#'   \item 4.5 - 5.2 ppm (residual water)
#'   \item 5.5 - 6.0 ppm (urea)
#'   \item 9.7 to max(ppm)  (no signal downfield)
#' }
#' @return
#' The function exports the following two objects into the currently active R environment (no variable assignments needed):
#' \itemize{
#'   \item Xc, num matrix: column reduced
#'   \item ppc, num array - chemical shift positions, length is equal to ncol(X_c)
#'
#' }
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @family NMR
#' @section

excise1d<-function(X, ppm){
  if(nrow(X)==0){X=t(X)}
  if(ncol(X)!=length(ppm)){ stop('ppm does not match to X')}
  if(any(is.na(ppm))){stop('ppm contains NA values')}

  idx_rm<- c( get.idx(c(min(ppm), 0.25), ppm), get.idx(c(4.5, 5.2), ppm), get.idx(c(5.5,6), ppm), get.idx(c(9.7, max(ppm)), ppm))

  Xc=X[,-idx_rm]
  ppc=ppm[-idx_rm]

  assign("Xc", Xc, envir = .GlobalEnv)
  assign("ppc", ppc, envir = .GlobalEnv)
}
