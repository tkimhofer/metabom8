#' @title Excise Standard Chemical Shift Regions from NMR Spectra
#'
#' @description
#' Removes standard chemical shift regions from 1D \eqn{^1}H NMR spectra typically excluded from metabolomics analysis.
#'
#' @param X Numeric matrix. NMR spectra with rows representing samples and columns representing chemical shift variables.
#' @param ppm Numeric vector. Chemical shift positions (in ppm), must match number of columns in \code{X}.
#'
#' @details
#' The following regions are removed:
#' \itemize{
#'   \item Upfield noise: \code{min(ppm)} to 0.25 ppm
#'   \item Residual water: 4.5 to 5.2 ppm
#'   \item Urea region: 5.5 to 6.0 ppm
#'   \item Downfield noise: 9.7 ppm to \code{max(ppm)}
#' }
#'
#' @return A list with:
#' \itemize{
#'   \item \code{Xc}: Numeric matrix with excised chemical shift regions removed.
#'   \item \code{ppc}: Numeric vector of ppm values corresponding to \code{Xc}.
#' }
#'
#' @examples
#' set.seed(1)
#' ppm <- seq(0, 10, length.out = 1000)
#' X <- matrix(rnorm(100 * length(ppm)), nrow = 100)
#' result <- excise1d(X, ppm)
#' dim(result$Xc)
#' length(result$ppc)
#'
#' @author Torben Kimhofer
#' @export
excise1d <- function(X, ppm) {
  if (is.vector(X)) X <- t(X)
  if (ncol(X) != length(ppm)) stop("ppm length does not match number of columns in X.")
  if (anyNA(ppm)) stop("ppm contains NA values.")

  idx_rm <- c(
    get.idx(c(min(ppm), 0.25), ppm),
    get.idx(c(4.5, 5.2), ppm),
    get.idx(c(5.5, 6.0), ppm),
    get.idx(c(9.7, max(ppm)), ppm)
  )

  keep_idx <- setdiff(seq_along(ppm), idx_rm)
  Xc <- X[, keep_idx, drop = FALSE]
  ppc <- ppm[keep_idx]

  return(list(Xc = Xc, ppc = ppc))
}
