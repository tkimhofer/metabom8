#' @title Baseline correction for 1D NMR spectra
#' @export
#' @aliases blcor
#' @param X num matrix or data.frame, NMR data with rows representing spectra
#' @param lambda, num smoothing paramter
#' @param iter_max num, maximum number of iterations
#' @description Non-linear baseline correction for NMR spectra based on asymetric least squares. This function estimates smooth and non-linear baseline trends for each spectrum. The baseline is then subtracted from a spectrum and the remainder is returned. For more infor on smoothing parameter *lambda*, check out \code{\link[ptw]{asysm}}.
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
bline <- function(X, lambda=1e7, iter_max=30) {
    if(any(is.na(X))){message('X contains missing values, replacing with zeros.')}
    X <- .dimX(X)
    X.bl <- t(apply(X, 1, function(x) {
        x - asysm(x, iter_max, lambda)
    }))
    return(X.bl)
}
