#' @title Baseline correction for 1D NMR spectra
#' @export
#' @aliases blcor
#' @param X num matrix or data.frame, NMR data with rows representing spectra
#' @description Non-linear baseline correction for NMR spectra. This function estimates smooth and non-linear baseline trends for each spectrum based on a asymmetric least squares method. The estimated baseline is subtracted from each spectrum and the remainder is returned.
#' @return Baseline corrected NMR matrix an with same dimension as input matrix.
#' @importFrom ptw asysm
#' @seealso \code{\link[ptw]{asysm}}
#' @family NMR
#' @examples
#' data(covid_raw)
#' X1_bc=bline(X[1,])
#' spec(X[1,], ppm, shift=c(3,4))
#' spec(X1_bc[1,], ppm, shift=c(3,4), add=TRUE)
#' @author \email{torben.kimhofer@@murdoch.edu.au}
bline <- function(X) {
    X <- .dimX(X)
    X.bl <- t(apply(X, 1, function(x) {
        x - asysm(x)
    }))
    return(X.bl)
}
