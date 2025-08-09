#' @title Normalize Spectra Using ERETIC Signal
#'
#' @description
#' This function normalizes 1D NMR spectra using the ERETIC reference signal. If `integr = TRUE`, it returns the ERETIC integral only. The reference region is inferred automatically based on the observed ERETIC peak (around 12 ppm for urine or 15 ppm for plasma).
#'
#' @param X Numeric matrix. NMR data with spectra in rows and chemical shift positions as column names.
#' @param integr Logical. If \code{TRUE}, the function returns the ERETIC integral for each spectrum. If \code{FALSE}, returns the normalized spectra.
#'
#' @return
#' \itemize{
#'   \item If \code{integr = TRUE}: numeric vector of ERETIC integrals per spectrum.
#'   \item If \code{integr = FALSE}: numeric matrix of spectra normalized by ERETIC integral.
#' }
#'
#' @details
#' The function identifies the dominant ERETIC peak within the range 10–16 ppm and then uses a narrow window around 12 ppm (urine) or 15 ppm (plasma) to compute the normalization factor.
#'
#' @family NMR
#' @export
#' @examples
#' set.seed(123)
#' n <- 1000
#' ppm <- seq(0, 14, length.out = n)
#' gauss <- function(x, c, h, w=0.05) h * exp(-((x - c)^2) / (2 * w^2))
#' heights <- c(100, 80, 60, 40, 20)
#' spectra <- sapply(heights, function(h) rnorm(n, 0, 0.01) + gauss(ppm, 12, h))
#' spectra <- t(spectra)  # 5 spectra x 1000 points
#' colnames(spectra) = ppm
#' matspec(spectra, ppm, shift=c(11, 13))
#'
#' X_norm = normErectic(spectra)
#' matspec(X_norm, ppm, shift=c(11, 13))
normErectic <- function(X, integr = FALSE) {
  ppm <- as.numeric(colnames(X))
  if (is.null(ppm) || any(is.na(ppm)) || length(ppm) != ncol(X)) stop("ppm values must be present in column names of X.")

  idx_scan <- get_idx(c(10, 16), ppm)
  eretic_peak <- ppm[which.max(X[1, idx_scan]) + idx_scan[1] - 1]

  # Determine ERETIC target window based on closest known peak
  if (abs(eretic_peak - 12) < abs(eretic_peak - 15)) {
    idx_eretic <- get_idx(c(11.9, 12.1), ppm)  # Urine
  } else {
    idx_eretic <- get_idx(c(14.9, 15.1), ppm)  # Plasma
  }

  if (integr) {
    integrals <- rowSums(X[, idx_eretic, drop = FALSE])
    return(integrals)
  } else {
    Xn <- t(apply(X, 1, function(x, idx) x / sum(x[idx]), idx = idx_eretic))
    colnames(Xn) <- colnames(X)
    rownames(Xn) <- rownames(X)
    return(Xn)
  }
}
