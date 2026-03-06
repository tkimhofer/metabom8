#' @title Full Width at Half Maximum (FWHM) Estimation
#'
#' @description
#' Estimates the full width at half maximum (FWHM; line width) of a singlet-like peak
#' within a specified chemical-shift range for each spectrum.
#'
#' @param X Numeric matrix (spectra in rows) \emph{or} a named list as returned by
#'   \code{\link{read1d}}/\code{read1d_proc} containing \code{X}, \code{ppm}, and \code{meta}.
#' @param ppm Numeric vector of chemical shift values (ppm) corresponding to columns of \code{X}.
#'   If \code{NULL}, \code{ppm} is inferred in the following order:
#'   \enumerate{
#'     \item \code{attr(X, "m8_axis")$ppm} (if present),
#'     \item numeric \code{colnames(X)} (if present).
#'   }
#' @param shift Numeric vector of length 2. Chemical shift range containing the peak
#'   (e.g., \code{c(-0.1, 0.1)} for TSP).
#' @param sf Spectrometer frequency in MHz. Either a single numeric (recycled across spectra)
#'   or a numeric vector of length \code{nrow(X)} (one value per spectrum).
#'   If \code{X} is a list input and \code{sf} is missing, \code{sf} is taken from
#'   \code{X$meta$a_SFO1} when available. If \code{sf} is missing for a matrix input,
#'   the function attempts to use \code{attr(X, "m8_meta")$a_SFO1} when present.
#'
#' @details
#' For each spectrum, the function:
#' \enumerate{
#'   \item extracts the region defined by \code{shift},
#'   \item finds the peak apex within that region,
#'   \item computes the half-height level relative to the local baseline (minimum in the window),
#'   \item estimates the left and right half-height crossing points by linear interpolation,
#'   \item converts the width from ppm to Hz using \code{sf} (MHz), i.e. \code{Hz = ppm * sf}.
#' }
#'
#' If no valid half-height crossings can be found (e.g., very low SNR or truncated peak),
#' \code{NA} is returned for that spectrum.
#'
#' @return Numeric vector of FWHM values in Hz (length \code{nrow(X)}).
#'
#' @note
#' The ppm axis may be increasing or decreasing; FWHM is computed as an absolute width
#' and is therefore independent of axis direction.
#'
#' @examples
#' # Simulated NMR peaks with different linewidths
#' ppm <- seq(-0.2, 0.2, length.out = 1000)
#'
#' # generate peaks with increasing width
#' sds <- seq(0.01, 0.03, length.out = 10)
#'
#' X <- t(sapply(sds, function(s)
#'   dnorm(ppm, mean = 0, sd = s)
#' ))
#'
#' sf <- 600  # spectrometer frequency in MHz
#'
#' fwhm_vals <- lw(X, ppm = ppm, shift = c(-0.1, 0.1), sf = sf)
#'
#' plot(sds, fwhm_vals,
#'   xlab = "Gaussian sd",
#'   ylab = "Estimated FWHM (Hz)",
#'   pch = 16)
#' @export
lw <- function(X, ppm = NULL, shift = c(-0.1, 0.1), sf) {

  if (is.list(X) && !is.null(X$X)) {
    dat <- X
    X <- dat$X
    if (is.null(ppm)) ppm <- dat$ppm
    if (missing(sf) && !is.null(dat$meta) && "a_SFO1" %in% names(dat$meta))
      sf <- dat$meta$a_SFO1
  }

  X <- .dimX(X)

  if (is.null(ppm)) {
    ax <- attr(X, "m8_axis", exact = TRUE)
    if (is.list(ax) && !is.null(ax$ppm)) {
      ppm <- ax$ppm
    } else if (!is.null(colnames(X))) {
      ppm <- as.numeric(colnames(X))
      if (anyNA(ppm)) ppm <- NULL
    }
  }
  if (is.null(ppm)) stop("`ppm` not provided and not found in attr(X,'m8_axis') or numeric colnames(X).")
  ppm <- as.numeric(ppm)
  if (!.check_X_ppm(X, ppm)) stop("X and ppm mismatch or ppm contains NA/Inf.")

  if (missing(sf)) {
    meta <- attr(X, "m8_meta", exact = TRUE)
    if (is.data.frame(meta) && "a_SFO1" %in% names(meta)) {
      sf <- meta$a_SFO1
    } else {
      stop("`sf` missing. Provide sf (MHz) or supply meta$a_SFO1 via list input or attr(X,'m8_meta').")
    }
  }

  sf <- as.numeric(sf)
  if (any(!is.finite(sf))) stop("`sf` contains non-finite values.")
  if (length(sf) == 1L) sf <- rep(sf, nrow(X))
  if (length(sf) != nrow(X)) stop("`sf` must be length 1 or length nrow(X).")

  idx <- get_idx(shift, ppm)
  if (length(idx) < 3L) stop("Selected `shift` region too small for FWHM estimation.")

  vapply(seq_len(nrow(X)), function(i) {
    x <- X[i, idx]
    x[is.na(x)] <- 0
    pp <- ppm[idx]

    baseline <- min(x, na.rm = TRUE)
    apex_i <- which.max(x)
    apex <- x[apex_i]
    height <- apex - baseline
    if (!is.finite(height) || height <= 0) return(NA_real_)

    hw <- baseline + 0.5 * height

    left_le <- which(x[seq_len(apex_i)] <= hw)
    if (length(left_le) == 0L) return(NA_real_)
    iL2 <- max(left_le)
    if (iL2 >= apex_i) return(NA_real_)
    iL1 <- iL2 + 1L

    right_le <- which(x[apex_i:length(x)] <= hw)
    if (length(right_le) == 0L) return(NA_real_)
    iR1 <- min(right_le) + apex_i - 1L
    if (iR1 <= apex_i) return(NA_real_)
    iR2 <- iR1 - 1L

    xL <- pp[iL2] + (pp[iL1] - pp[iL2]) * (hw - x[iL2]) / (x[iL1] - x[iL2])
    xR <- pp[iR2] + (pp[iR1] - pp[iR2]) * (hw - x[iR2]) / (x[iR1] - x[iR2])

    abs(xR - xL) * sf[i]
  }, numeric(1))
}



#' @title Estimate Noise Standard Deviation in 1D NMR Spectra
#'
#' @description
#' Estimates the noise standard deviation (\eqn{\sigma}) for each spectrum from a signal-free ppm region.
#' The default estimator is robust (MAD-based) and suitable for signal-to-noise calculations.
#'
#' @param X Numeric matrix (spectra in rows), numeric vector (single spectrum), or a named list
#'   as returned by \code{\link{read1d}}/\code{read1d_proc} containing \code{X}, \code{ppm}, and \code{meta}.
#' @param ppm Numeric vector of chemical shift values (ppm) corresponding to columns of \code{X}.
#'   If \code{NULL}, \code{ppm} is inferred in the following order:
#'   \enumerate{
#'     \item \code{X$ppm} if \code{X} is a list input,
#'     \item \code{attr(X, "m8_axis")$ppm} (if present),
#'     \item numeric \code{colnames(X)} (if present).
#'   }
#' @param where Numeric vector of length 2. ppm range used for noise estimation.
#'   Should be free of metabolite signals (default \code{c(14.6, 14.7)}).
#' @param method Character. Noise estimator: \code{"mad"} (default), \code{"sd"}, or \code{"p95"} (legacy amplitude).
#' @param baseline_correct Logical. If \code{TRUE}, subtract a smooth baseline in the noise window
#'   using asymmetric least squares (\code{\link[ptw]{asysm}}). Default \code{FALSE}.
#' @param lambda Numeric. Smoothing parameter for \code{\link[ptw]{asysm}} when \code{baseline_correct = TRUE}.
#' @param min_points Integer. Minimum number of points required in the noise window.
#'   May be a single value (recycled across spectra) or a numeric vector of
#'   length \code{nrow(X)}. If \code{X} is a list input as returned by
#'   \code{\link{read1d}}/\code{read1d_proc} and \code{ns} is \code{NULL},
#'   the function attempts to extract the number of scans from
#'   \code{X$meta$a_NS} when available.
#' @param ns Number of scans / transients (required if normalise_scans=TRUE)
#' @param normalise_scans Logical. If \code{TRUE}, noise estimates are
#'   multiplied by \eqn{\sqrt{NS}} to account for the theoretical scaling
#'   of noise with the number of scans (\eqn{\sigma \propto 1/\sqrt{NS}}).
#'   This is useful when comparing noise levels across spectra acquired
#'   with different numbers of scans. Default is \code{FALSE}.
#' @details
#' In NMR spectroscopy, noise scales predictably with the number of scans (\eqn{NS}).
#' For otherwise identical acquisition settings:
#'
#' \itemize{
#'   \item Signal increases approximately proportional to \eqn{NS}.
#'   \item Noise increases approximately proportional to \eqn{\sqrt{NS}}.
#'   \item Consequently, signal-to-noise ratio (SNR) increases proportional to \eqn{\sqrt{NS}}.
#' }
#'
#' Equivalently, the noise standard deviation scales as:
#'
#' \deqn{\sigma(NS) \propto \frac{1}{\sqrt{NS}},}
#'
#' assuming a fixed underlying signal scale and comparable acquisition conditions.
#'
#' To compare noise levels across datasets acquired with different numbers of scans,
#' a scan-normalised noise estimate may be used:
#'
#' \deqn{\sigma_{\mathrm{norm}} = \sigma \cdot \sqrt{NS}.}
#'
#' Under stable receiver gain and processing conditions, this normalised noise
#' should be approximately constant across runs.
#'
#' @return Numeric vector of noise estimates (length \code{nrow(X)}).
#'
#' @examples
#' data(hiit_raw)
#' X <- hiit_raw$X
#' ppm <- hiit_raw$ppm
#' sigma <- noise_sd(X, ppm, where = c(10,11))
#' plot(hiit_raw$meta$a_NS, sigma,
#'   xlab = "Number of scans (NS)",
#'   ylab = expression(sigma~"(noise estimate)"),
#'    pch = 16)
#'  lines(lowess(hiit_raw$meta$a_NS, sigma), col = "red", lwd = 2)
#' @importFrom stats sd quantile mad
#' @importFrom ptw asysm
#' @export
noise_sd <- function(X, ppm = NULL, where = c(14.6, 14.7),
                     method = c("mad","sd","p95"),
                     baseline_correct = FALSE, lambda = 1e4, min_points = 50L,
                     ns = NULL, normalise_scans = FALSE) {

  if (is.list(X) && !is.null(X$X)) {
    dat <- X
    X <- dat$X
    if (is.null(ppm)) ppm <- dat$ppm
    if (is.null(ns) && !is.null(dat$meta) && "a_NS" %in% names(dat$meta))
      ns <- dat$meta$a_NS
  }

  X <- .dimX(X)
  if (is.null(ppm)) {
    ax <- attr(X, "m8_axis", exact = TRUE)
    if (is.list(ax) && !is.null(ax$ppm)) {
      ppm <- ax$ppm
    } else if (!is.null(colnames(X))) {
      ppm <- as.numeric(colnames(X))
      if (anyNA(ppm)) ppm <- NULL
    }
  }
  if (is.null(ppm)) stop("`ppm` not provided and not found in list input, attr(X,'m8_axis'), or numeric colnames(X).")
  ppm <- as.numeric(ppm)
  if (!.check_X_ppm(X, ppm)) stop("X and ppm mismatch or ppm contains NA/Inf.")

  method <- match.arg(method)

  idx <- get_idx(where, ppm)
  if (length(idx) < as.integer(min_points)) {
    stop(sprintf("Insufficient number of points in selected region (need >= %d).", as.integer(min_points)))
  }

  sigma <- apply(X[, idx, drop = FALSE], 1, function(x) {
    x[is.na(x)] <- 0
    if (baseline_correct) x <- x - ptw::asysm(x, lambda = lambda)
    switch(method,
           mad = stats::mad(x, constant = 1.4826, na.rm = TRUE),
           sd  = stats::sd(x, na.rm = TRUE),
           p95 = as.numeric(stats::quantile(x, probs = 0.95, na.rm = TRUE, names = FALSE))
    )
  })

  if (isTRUE(normalise_scans)) {
    if (is.null(ns)) stop("normalise_scans=TRUE requires `ns` or meta$a_NS.")
    ns <- as.numeric(ns)
    if (length(ns) == 1L) ns <- rep(ns, length(sigma))
    if (length(ns) != length(sigma)) stop("`ns` must be length 1 or nrow(X).")
    if (any(!is.finite(ns)) || any(ns <= 0)) stop("`ns` must be positive finite values.")
    sigma <- sigma / sqrt(ns)
  }

  sigma
}

#' @noRd
noise.est <- function(...) {
  warning(msg = "noise.est() is deprecated; use noise_sd() instead.", call.=FALSE)
  noise_sd(...)
}
