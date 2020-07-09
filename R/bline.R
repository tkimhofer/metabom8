#' Baseline correction for 1D NMR spectra
#' @export
#' @aliases bline
#' @param X NMR data matrix or dataframe with rows representing spectra.
#' @description Baseline correction for NMR spectra. This function firstly estimates a baseline trend based on a asymmetric least squares. Secondly the estimated baseline is subtracted from the NMR spectrum. Baseline estimation uses the the \code{asysm} function as implemented in the \code{ptw} package.
#' @return Baseline corrected NMR matrix.
#' @importFrom ptw asysm
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
bline <- function(X) {
  X.bl <- t(apply(X, 1, function(x) {
    x - asysm(x)
  }))
  return(X.bl)
}
