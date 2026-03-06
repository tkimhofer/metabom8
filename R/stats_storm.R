#' @title Subset Optimisation by Reference Matching (STORM)
#'
#' @description
#' Selects an optimal subset of spectra that best match a specified target signal
#' region, improving downstream correlation-based structural analysis such as STOCSY.
#'
#' @param X Numeric matrix (or data.frame) of NMR spectra with samples in rows
#'   and spectral variables in columns.
#' @param ppm Numeric vector of chemical shift values corresponding to the columns of \code{X}.
#' @param b Integer. Half-window size expressed as number of spectral variables
#'   (data points). The effective window width therefore depends on the ppm spacing.
#' @param q Numeric. P-value threshold for including spectral variables when
#'   updating the reference region.
#' @param idx.refSpec Integer. Row index of \code{X} defining the initial reference spectrum.
#' @param shift Numeric vector of length 2 giving the ppm range (minimum, maximum)
#'   defining the initial target signal region.
#'
#' @details
#' STORM iteratively refines a subset of spectra exhibiting consistent signal
#' position and multiplicity within a specified ppm region.
#'
#' Starting from an initial reference spectrum and ppm window:
#' \enumerate{
#'   \item Spectra positively correlated with the current reference signal are retained.
#'   \item A driver peak (maximum intensity within the reference window) is identified.
#'   \item Correlation and covariance are evaluated within a local window of size \code{b}
#'         around the driver peak.
#'   \item The reference region is updated using variables that satisfy the p-value
#'         threshold (\code{q}) and show positive correlation.
#' }
#'
#' The procedure continues until the selected subset stabilises. The resulting
#' row indices define spectra that most consistently represent the structural
#' pattern of the target signal.
#'
#' STORM does not perform metabolite identification directly. Instead, it refines
#' the dataset to enhance structural coherence prior to correlation-based
#' interpretation methods such as STOCSY.
#'
#' @return Integer vector of row indices of \code{X} defining the selected spectral subset.
#'
#' @references
#' Posma, J. M., et al. (2012). Subset optimisation by reference matching (STORM):
#' an optimised statistical approach for recovery of metabolic biomarker structural
#' information from \eqn{^1H} NMR spectra of biofluids.
#' \emph{Analytical Chemistry}, 84(24), 10694–10701.
#'
#' @family structural_annotation
#' @importFrom stats cor cov pt
#' @export
#'
#' @examples
#' ## Simulated example with three Gaussian signals
#' set.seed(123)
#'
#' n <- 100          # number of spectra
#' S <- 1000         # number of spectral variables
#' ppm <- seq(10, 0, length.out = S)
#'
#' gauss <- function(x, centre, width, height) {
#'   height * exp(-((x - centre)^2) / (2 * width^2))
#' }
#'
#' ## Generate base signals
#' sig1 <- gauss(ppm, 7,   0.05, 10)
#' sig2 <- gauss(ppm, 3.5, 0.10,  8)
#' sig3 <- gauss(ppm, 1,   0.07,  5)
#'
#' spectra <- matrix(0, n, S)
#'
#' for (i in seq_len(n)) {
#'   spectra[i, ] <- sig1 + sig2 + rnorm(S, 0, 0.1)
#'   if (i <= 25) {
#'     spectra[i, ] <- spectra[i, ] + sig3
#'   }
#' }
#'
#' ## Apply STORM to refine spectra containing the 1 ppm signal
#' idx <- storm(
#'   X = spectra,
#'   ppm = ppm,
#'   b = 30,
#'   q = 0.05,
#'   idx.refSpec = 1,
#'   shift = c(0.75, 1.25)
#' )
#'
#' length(idx)   # number of spectra retained
storm <- function(X, ppm, b=30, q=0.05, idx.refSpec, shift){

  X_sc <- .scaleMatRcpp(X, seq.int(0, nrow(X) - 1), center=TRUE, scale_type=0)
  X <- X_sc[[1]]
  ref.idx <- get_idx(range=shift, ppm)
  ref <- X[idx.refSpec, ref.idx]

  subset <- 0
  subset1 <- seq_len(nrow(X))
  i <- 1

  subset  <- integer(0)
  subset1 <- seq_len(nrow(X))

  while (!setequal(subset, subset1)) {

    if (length(subset1) <= 1L) {
      message("STORM reduced to a single spectrum; stopping loop.")
      break
    }

    subset <- subset1

    Xr <- X[subset, ref.idx, drop = FALSE]
    r  <- stats::cor(t(Xr), ref)
    df1 <- length(ref) - 2L
    a  <- abs(r) * sqrt(df1 / pmax(1 - r^2, .Machine$double.eps))
    pval <- 2 * stats::pt(-a, df = df1)

    keep <- which(r > 0)
    subset1 <- subset[keep]
    pval <- pval[keep]

    ord <- order(pval)
    subset1 <- subset1[ord]

    idx_driver <- which.max(ref)
    driver <- ref.idx[idx_driver]

    ppm_win <- seq.int(max(1L, driver - (b + 1L)), min(ncol(X), driver + (b + 1L)))

    r  <- stats::cor(X[subset1, ppm_win, drop = FALSE], X[subset1, driver])
    co <- stats::cov(X[subset1, ppm_win, drop = FALSE], X[subset1, driver])
    co <- as.numeric(co)

    df2 <- length(subset1) - 2L
    a <- abs(r) * sqrt(df2 / pmax(1 - r^2, .Machine$double.eps))
    pval <- 2 * stats::pt(-a, df = df2)

    iid <- which(pval < q & r > 0)
    if (length(iid) == 0L) {
      message("No variables passed p/q criteria; stopping.")
      break
    }

    idx_max <- iid[which.max(co[iid])]

    ref_center <- ppm_win[idx_max]
    ref.idx <- seq.int(max(1L, ref_center - (b + 1L)), min(ncol(X), ref_center + (b + 1L)))

    ref <- stats::cov(X[subset1, ref.idx, drop = FALSE], X[subset1, driver])[, 1]
  }

  return(subset1)
}
