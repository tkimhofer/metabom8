#' @title Find Local Extrema in NMR Spectra (Peak Picking)
#'
#' @description
#' Identifies local maxima, minima, or both from smoothed NMR spectra using Savitzky–Golay filtering.
#'
#' @param X Numeric matrix. NMR data with spectra in rows and chemical shifts in columns.
#' @param ppm Numeric vector. Chemical shift values corresponding to columns in \code{X}.
#' @param fil_p Integer. Polynomial order of the Savitzky–Golay filter.
#' @param fil_n Integer. Filter length (must be odd) of the Savitzky–Golay filter.
#' @param type Character. Type of extrema to return: \code{"max"}, \code{"min"}, or \code{"both"}.
#'
#' @details
#' The spectra are smoothed using a Savitzky–Golay filter to reduce noise. Extrema are then detected by identifying sign changes in the first derivative of the smoothed signal.
#'
#' @return
#' A list of data frames, one per spectrum. Each data frame contains:
#' \itemize{
#'   \item \code{idc}: Index of the detected peak.
#'   \item \code{ppm}: Chemical shift at the peak.
#'   \item \code{Int}: Intensity at the peak.
#'   \item \code{Etype}: Extrema type: \code{1} for minima, \code{-1} for maxima.
#' }
#'
#' @importFrom signal sgolayfilt
#' @seealso \code{\link{ppick2}}
#' @examples
#' data(covid)
#' X <- covid$X
#' ppm <- covid$ppm
#'
#' peaklist <- ppick(X, ppm)
#' plot_spec(X[1, ], ppm, shift = c(1.2, 1.4), backend='base')
#' points(peaklist[[1]]$ppm, peaklist[[1]]$Int, col = 'cyan')
#'
#' @importFrom utils head tail
#' @export
#' @export
ppick <- function(X, ppm, fil_p = 3, fil_n = 5, type = "max") {

  X <- .dimX(X)

  if (is.null(ppm)) stop("`ppm` must be provided.", call. = FALSE)
  ppm <- as.numeric(ppm)
  if (anyNA(ppm)) stop("`ppm` must be numeric and must not contain NA.", call. = FALSE)

  if (!.check_X_ppm(X, ppm))
    stop("Dimensions of X and ppm must match and ppm must not contain NA/Inf.", call. = FALSE)

  type <- match.arg(type, c("max", "min", "both"))

  if (!is.numeric(fil_p) || length(fil_p) != 1L || fil_p < 1)
    stop("`fil_p` must be a positive integer.", call. = FALSE)
  if (!is.numeric(fil_n) || length(fil_n) != 1L || fil_n < 3)
    stop("`fil_n` must be an integer >= 3.", call. = FALSE)
  fil_p <- as.integer(fil_p)
  fil_n <- as.integer(fil_n)

  if (fil_n %% 2L == 0L) stop("`fil_n` must be odd.", call. = FALSE)
  if (fil_n <= fil_p) stop("`fil_n` must be greater than `fil_p`.", call. = FALSE)

  lapply(seq_len(nrow(X)), function(i) {
    x <- X[i, ]
    x[is.na(x)] <- 0

    x_smooth <- signal::sgolayfilt(x, p = fil_p, n = fil_n)

    dsgn <- sign(diff(x_smooth))

    if (any(dsgn == 0)) {
      nz <- which(dsgn != 0)
      if (length(nz) > 0) {
        for (k in seq_along(dsgn)) {
          if (dsgn[k] == 0) dsgn[k] <- if (k == 1) dsgn[nz[1]] else dsgn[k - 1]
        }
      }
    }

    extrema_idx <- which(utils::head(dsgn, -1) != utils::tail(dsgn, -1)) + 1L
    if (length(extrema_idx) == 0L) return(NULL)

    extrema_type <- -dsgn[extrema_idx - 1L]

    res <- data.frame(
      idc = extrema_idx,
      ppm = ppm[extrema_idx],
      Int = x[extrema_idx],
      Etype = extrema_type
    )

    switch(type,
           "max"  = res[res$Etype < 0, , drop = FALSE],
           "min"  = res[res$Etype > 0, , drop = FALSE],
           "both" = res
    )
  })
}


#' @title Peak picking using Savitzky–Golay derivatives
#'
#' @description
#' Finds local extrema in 1D spectra using Savitzky–Golay first and second derivatives.
#' Candidate peaks are identified at zero-crossings of the first derivative and classified
#' by the sign of the second derivative. Optional filters control peak height, prominence,
#' SNR, curvature and minimum separation.
#'
#' @param X Numeric matrix or vector. Spectra in rows.
#' @param ppm Numeric vector or NULL. If NULL, inferred from colnames(X).
#' @param type Character. "max", "min", or "both".
#' @param fil_p Integer. Polynomial order for SG filter.
#' @param fil_n Integer. Window length for SG filter (odd).
#' @param noise_win Numeric length-2 vector or NULL. Ppm window used to estimate noise per spectrum.
#'   If NULL, a robust noise estimate is computed from the full spectrum (MAD of first differences).
#' @param min_snr Numeric. Minimum SNR (peak height / noise). Set NULL to disable.
#' @param min_height Numeric. Minimum absolute peak height (in original intensity units). NULL disables.
#' @param min_prominence Numeric. Minimum local prominence. NULL disables.
#' @param prom_half_window_ppm Numeric. Half-window (ppm) for prominence estimation around each peak.
#' @param min_distance_ppm Numeric. Minimum separation between peaks (ppm). NULL disables.
#' @param min_curvature Numeric. Minimum absolute curvature at peak (|d2|). NULL disables.
#' @param keep_cols Character. Extra columns to keep (default keeps all computed).
#'
#' @return List of data.frames (one per spectrum). Each contains:
#' idc, ppm, Int, Etype (+1 max, -1 min), height, snr, curvature, prominence.
#'
#' @importFrom signal sgolayfilt
#' @examples
#' data(covid)
#' X <- covid$X
#' ppm <- covid$ppm
#'
#' peaklist <- ppick2(X[1,], ppm, min_snr=50)
#'
#' plot_spec(X[1, ], ppm, shift = c(3, 4.5), backend='base')
#' points(peaklist[[1]]$ppm, peaklist[[1]]$Int, col = "cyan")
#' head(peaklist[[1]])
#' @export
ppick2 <- function(X,
                   ppm = NULL,
                   type = c("max", "min", "both"),
                   fil_p = 3L,
                   fil_n = 11L,
                   noise_win = NULL,
                   min_snr = 10,
                   min_height = NULL,
                   min_prominence = NULL,
                   prom_half_window_ppm = 0.02,
                   min_distance_ppm = 0.005,
                   min_curvature = NULL,
                   keep_cols = c("height", "snr", "curvature", "prominence")) {


  u <- .m8_unpack_dat(X, ppm = ppm)
  X <- .dimX(u$X)
  ppm <- u$ppm

  X <- .dimX(X)

  if (is.null(ppm)) {
    if (is.null(colnames(X))) stop("ppm is NULL and colnames(X) are missing.", call. = FALSE)
    ppm <-as.numeric(colnames(X))
    if (anyNA(ppm)) stop("ppm could not be inferred from colnames(X): non-numeric column names.", call. = FALSE)
  }
  ppm <- as.numeric(ppm)
  if (!.check_X_ppm(X, ppm)) stop("X and ppm dimensions mismatch or ppm contains NA/Inf.", call. = FALSE)

  type <- match.arg(type)

  fil_p <- as.integer(fil_p)
  fil_n <- as.integer(fil_n)
  if (fil_n < 5L) stop("fil_n should be >= 5 for stable derivative estimation.", call. = FALSE)
  if (fil_n %% 2L == 0L) stop("fil_n must be odd.", call. = FALSE)
  if (fil_n <= fil_p) stop("fil_n must be greater than fil_p.", call. = FALSE)

  if (!is.null(noise_win)) {
    if (!is.numeric(noise_win) || length(noise_win) != 2L) stop("noise_win must be numeric length 2 or NULL.", call. = FALSE)
    idx_noise <- get_idx(noise_win, ppm)
    if (length(idx_noise) < 5L) stop("noise_win too small or not present in ppm axis.", call. = FALSE)
  } else {
    idx_noise <- NULL
  }

  .noise_mad_diff <- function(x) {
    d <- diff(x)
    s <- stats::mad(d, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s == 0)
      s <- sqrt(.Machine$double.eps)
    s
  }

  .prominence_local <- function(x, ppm, i0, half_ppm) {

    win <- c(ppm[i0] - half_ppm, ppm[i0] + half_ppm)
    idx <- get_idx(win, ppm)
    idx <- idx[idx != i0]
    if (length(idx) < 2L) return(NA_real_)

    left <- idx[idx < i0]
    right <- idx[idx > i0]
    if (length(left) < 1L || length(right) < 1L) return(NA_real_)
    base_left <- min(x[left], na.rm = TRUE)
    base_right <- min(x[right], na.rm = TRUE)
    x[i0] - max(base_left, base_right)
  }


  .enforce_min_distance <- function(tab, min_dist_ppm) {

    if (is.null(min_dist_ppm) || nrow(tab) <= 1L) return(tab)
    tab <- tab[order(-tab$height), , drop = FALSE]
    keep <- rep(TRUE, nrow(tab))
    for (i in seq_len(nrow(tab))) {
      if (!keep[i]) next
      too_close <- which(
        is.finite(tab$ppm) &
          abs(tab$ppm - tab$ppm[i]) < min_dist_ppm
      )
      too_close <- setdiff(too_close, i)
      keep[too_close] <- FALSE
    }
    tab2 <- tab[keep, , drop = FALSE]
    tab2[order(tab2$idc), , drop = FALSE]
  }

  ref <- if (nrow(X) > 1) apply(X, 2, median) else X[1,]
  d  <- diff(ref)
  sigma <- .noise_mad_diff(d)
  mask <- abs(ref) > (3 * sigma)
  k <- 9L
  mask <- stats::filter(mask,
                        rep(1, k),
                        sides = 2) > 0

  lapply(seq_len(nrow(X)), function(r) {
    x <- X[r, ]
    x[!mask] <- 0

    d1 <- signal::sgolayfilt(x, p = fil_p, n = fil_n, m = 1)
    d2 <- signal::sgolayfilt(x, p = fil_p, n = fil_n, m = 2)

    s1 <- sign(d1)
    if (any(s1 == 0)) {
      nz <- which(s1 != 0)
      if (length(nz) > 0) {
        for (k in seq_along(s1)) {
          if (s1[k] == 0) s1[k] <- if (k == 1) s1[nz[1]] else s1[k - 1]
        }
      }
    }

    zc <- which(diff(s1) != 0) + 1L

    if (length(zc) == 0L) return(NULL)

    etype <- ifelse(d2[zc] < 0, 1L, -1L)

    if (type == "max") {
      keep0 <- etype == 1L
    } else if (type == "min") {
      keep0 <- etype == -1L
    } else {
      keep0 <- rep(TRUE, length(zc))
    }
    zc <- zc[keep0]
    etype <- etype[keep0]
    if (length(zc) == 0L) return(NULL)

    height <- x[zc]
    curvature <- d2[zc]

    noise <- if (!is.null(idx_noise)) {
      apply(abs(x[idx_noise, drop = TRUE]), 1, max, na.rm = TRUE)
      max(abs(x[idx_noise]), na.rm = TRUE)
    } else {
      .noise_mad_diff(x)
    }
    snr <- if (is.finite(noise) && noise > 0) abs(height) / noise else NA_real_

    prominence <- vapply(zc, function(i0) .prominence_local(x, ppm, i0, prom_half_window_ppm), numeric(1))

    res <- data.frame(
      idc = zc,
      ppm = ppm[zc],
      Int = x[zc],
      Etype = etype,
      height = height,
      snr = snr,
      curvature = curvature,
      prominence = prominence
    )

    # filters
    if (!is.null(min_height)) {
      res <- res[res$height >= min_height, , drop = FALSE]
    }
    if (!is.null(min_snr)) {
      res <- res[is.finite(res$snr) & res$snr >= min_snr, , drop = FALSE]
    }
    if (!is.null(min_curvature)) {
      res <- res[abs(res$curvature) >= min_curvature, , drop = FALSE]
    }
    if (!is.null(min_prominence)) {
      res <- res[is.finite(res$prominence) & res$prominence >= min_prominence, , drop = FALSE]
    }
    if (nrow(res) == 0L) return(NULL)

    res <- .enforce_min_distance(res, min_distance_ppm)
    if (nrow(res) == 0L) return(NULL)

    if (!is.null(keep_cols)) {
      keep_cols <- intersect(keep_cols, names(res))
      base_cols <- c("idc", "ppm", "Int", "Etype")
      res <- res[, unique(c(base_cols, keep_cols)), drop = FALSE]
    }

    res
  })
}



