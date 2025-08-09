#' @title Baseline Correction for 1D NMR Spectra
#'
#' @description
#' Non-linear baseline correction for NMR spectra based on asymmetric least squares. This function estimates a smooth, non-linear baseline trend for each spectrum. The estimated baseline is then subtracted, returning the corrected spectrum. See \code{\link[ptw]{asysm}} for more information on the smoothing parameter \code{lambda}.
#'
#' @param X Numeric matrix or data frame. NMR data with spectra in rows.
#' @param lambda Numeric. Smoothing parameter passed to \code{\link[ptw]{asysm}}. Larger values result in smoother baselines.
#' @param iter_max Integer. Maximum number of iterations for the baseline estimation algorithm.
#'
#' @return A numeric matrix of the same dimensions as the input, containing the baseline-corrected spectra.
#'
#' @details
#' Missing values in \code{X} are replaced with zeros before baseline correction.
#'
#' @seealso \code{\link[ptw]{asysm}}, \code{\link{spec}}, \code{\link{matspec}}
#'
#' @examples
#' data(covid_raw)
#' X <- covid_raw$X
#' ppm <- covid_raw$ppm
#' X_bc <- bcor(X[1, ])
#'
#' spec(X[1, ], ppm, shift = c(3, 4), interactive = FALSE)
#' spec(X_bc[1, ], ppm, shift = c(3, 4), col='red', add = TRUE, interactive = FALSE)
#'
#' @importFrom ptw asysm
#' @export
bline <- function(X, lambda = 1e7, iter_max = 30) {
  if (any(is.na(X))) {
    message("X contains missing values, replacing with zeros.")
    X[is.na(X)] <- 0
  }

  if (as.character(match.call()[[1]]) == "bline") {
    warning("`bline` will be removed in future versions, please use `bcor` instead.", call. = FALSE)
  }

  X <- .dimX(X)
  X.bl <- t(apply(X, 1, function(x) {
    x - asysm(x, maxit = iter_max, lambda = lambda)
  }))
  return(X.bl)
}

#' @rdname bline
#' @export
bcor <- bline
