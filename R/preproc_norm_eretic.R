#' @title Normalise Spectra Using ERETIC Signal
#'
#' @description
#' Normalises 1D NMR spectra using an ERETIC reference signal. By default
#' the ERETIC position is discovered in a search window and the spectra are
#' normalised by the integral in a narrow window around the detected position.
#' Power users can supply \code{pos} to enforce a fixed ERETIC position.
#'
#' @param X Numeric matrix or vector. NMR spectra with spectra in rows.
#'   If \code{ppm} is not provided, it is inferred from \code{colnames(X)}.
#' @param integr Logical. If \code{TRUE}, returns the ERETIC integral per spectrum.
#'   If \code{FALSE}, returns the normalised spectra.
#' @param ppm Numeric vector of chemical shift positions. If \code{NULL}, inferred from \code{colnames(X)}.
#' @param pos Numeric scalar or \code{NULL}. If \code{NULL} (default), ERETIC position is discovered
#'   within \code{search}. If provided, the ERETIC window is centered on \code{pos}.
#' @param search Numeric vector of length 2. Ppm window used to discover ERETIC when \code{pos = NULL}.
#' @param width Numeric scalar. Full integration window width in ppm (default 0.2 gives \code{pos ± 0.1}).
#' @param warn_absent Logical. If \code{TRUE}, warn when ERETIC integral is zero/invalid for individual spectra.
#' @param noise_win Numeric vector of length 2. Ppm window for estimating noise intensity.
#' @param snr_factor Numeric scaler. Minimum Signal-to-Noise Ratio required for detecting the ERETIC peak.
#' @details
#' Binning/normalisation history is recorded in \code{attr(X, "m8_prep")} if present.
#' Ppm axis values are stored in \code{attr(X, "m8_axis")$ppm}.
#'
#' @return
#' If \code{integr = TRUE}, numeric vector of ERETIC integrals.
#' If \code{integr = FALSE}, numeric matrix of normalised spectra.
#'
#' @examples
#' set.seed(123)
#'
#' n <- 1000
#' ppm <- seq(0, 14, length.out = n)
#'
#' gauss <- function(x, c, h, w = 0.05) {
#'   h * exp(-((x - c)^2) / (2 * w^2))
#' }
#'
#' heights <- c(100, 80, 60, 40, 20)
#'
#' spectra <- sapply(heights, function(h)
#'   rnorm(n, 0, 0.01) + gauss(ppm, 12, h)
#' )
#'
#' spectra <- t(spectra)   # 5 spectra × 1000 variables
#'
#' plot_spec(spectra, ppm, shift = c(11, 13))
#'
#' X_norm <- norm_eretic(spectra, ppm=ppm)
#' plot_spec(X_norm, ppm, shift=c(11, 13))
#' @export
norm_eretic <- function(X, integr = FALSE,
                        ppm = NULL,
                        pos = NULL,
                        search = c(10, 16),
                        width = 0.2,
                        warn_absent = TRUE,
                        noise_win = c(9.8, 10.2),
                        snr_factor = 10) {

  .return_list <- is.list(X) && !is.null(X$X)
  meta <- NULL

  if (.return_list) {
    u <- .m8_unpack_dat(X, ppm = ppm)
    X   <- u$X
    ppm <- u$ppm
    meta <- u$meta
  }

  X0 <- X
  X  <- .dimX(X)

  if (is.null(ppm) && !is.null(colnames(X))) {
    ppm <- as.numeric(colnames(X))
    if (anyNA(ppm)) stop("ppm could not be inferred from colnames(X): non-numeric column names.", call. = FALSE)
  }
  if (!.check_X_ppm(X, ppm))
    stop("Non-matching dimensions X matrix and ppm vector or invalid ppm values.", call. = FALSE)
  ppm <- as.numeric(ppm)

  if (!is.null(pos)) {
    if (!is.numeric(pos) || length(pos) != 1L || !is.finite(pos))
      stop("`pos` must be a single finite numeric value or NULL.", call. = FALSE)
  }
  if (!is.numeric(search) || length(search) != 2L)
    stop("`search` must be a numeric vector of length 2.", call. = FALSE)
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0)
    stop("`width` must be a single positive numeric value.", call. = FALSE)
  if (!is.numeric(noise_win) || length(noise_win) != 2L)
    stop("`noise_win` must be a numeric vector of length 2.", call. = FALSE)
  if (!is.numeric(snr_factor) || length(snr_factor) != 1L || !is.finite(snr_factor) || snr_factor <= 0)
    stop("`snr_factor` must be a single positive numeric value.", call. = FALSE)

  discovered <- FALSE
  if (is.null(pos)) {
    idx_scan <- get_idx(search, ppm)
    if (length(idx_scan) < 3L)
      stop("Search window not present in ppm axis (too few points).", call. = FALSE)
    mean_spec <- colMeans(X, na.rm = TRUE)
    pos <- ppm[idx_scan][which.max(mean_spec[idx_scan])]
    discovered <- TRUE
  }

  half <- width / 2
  win <- c(pos - half, pos + half)
  idx_eretic <- get_idx(win, ppm)
  if (length(idx_eretic) < 2L)
    stop("ERETIC integration window not present in ppm axis.", call. = FALSE)

  idx_noise <- get_idx(noise_win, ppm)
  if (length(idx_noise) < 2L)
    stop("Noise window not present in ppm axis (adjust `noise_win`).", call. = FALSE)

  integrals <- rowSums(X[, idx_eretic, drop = FALSE], na.rm = TRUE)

  eretic_peak <- apply(X[, idx_eretic, drop = FALSE], 1, max, na.rm = TRUE)
  noise_level <- apply(abs(X[, idx_noise, drop = FALSE]), 1, max, na.rm = TRUE)

  present <- is.finite(eretic_peak) & is.finite(noise_level) & (eretic_peak > snr_factor * noise_level)

  if (warn_absent && any(!present)) {
    show <- paste(head(which(!present), 10L), collapse = ",")
    more <- if (sum(!present) > 10L) sprintf(" (+%d more)", sum(!present) - 10L) else ""
    warning(sprintf(
      "ERETIC likely absent (peak <= noise*%g) for %d spectra: %s%s",
      snr_factor, sum(!present), show, more
    ), call. = FALSE)
  }

  if (integr) {
    return(integrals)
  }

  denom <- integrals
  bad <- !is.finite(denom) | denom == 0 | !present
  denom[bad] <- 1
  Xn <- X / denom

  rownames(Xn) <- rownames(X)
  colnames(Xn) <- colnames(X)

  Xn <- .m8_copy_attrs(X0, Xn)

  ax <- attr(Xn, "m8_axis", exact = TRUE)
  if (!is.list(ax)) ax <- list()
  ax$ppm <- ppm
  attr(Xn, "m8_axis") <- ax

  Xn <- .m8_stamp(
    Xn,
    step = "normErectic",
    params = list(
      pos = pos,
      discovered = discovered,
      search = if (discovered) search else NULL,
      width = width,
      window = win,
      noise_win = noise_win,
      snr_factor = snr_factor,
      mode = "divide_by_integral"
    ),
    notes = "Normalised by ERETIC integral; spectra failing peak>noise*SNR left unchanged."
  )

  if (.return_list) {
    list(X = Xn, ppm = ppm, meta = meta)
  } else {
    Xn
  }
}
