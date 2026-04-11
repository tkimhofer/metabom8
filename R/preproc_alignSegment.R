#' @title Align NMR Spectra via Cross-Correlation
#'
#' @description
#' Aligns rows of an NMR spectral matrix to a reference spectrum by maximizing
#' cross-correlation. Optionally uses the row-wise median spectrum as reference.
#' Useful for minor spectral misalignments.
#'
#' @param seg Numeric matrix. Each row is a 1D NMR spectrum segment.
#' @param idx_ref Integer. Row index to use as reference spectrum. Ignored if
#'   \code{med = TRUE}.
#' @param clim Numeric. Minimum cross-correlation threshold. Segments with lower
#'   similarity are not shifted.
#' @param med Logical. If \code{TRUE}, the row-wise median spectrum is used as
#'   the reference.
#' @param norm Logical. If \code{TRUE}, rows are z-scaled prior to alignment.
#'
#' @return A numeric matrix of aligned spectra with the same dimensions as
#'   \code{seg}. If \code{med = TRUE}, the output excludes the median row.
#'
#' @details
#' Each spectrum is aligned to the reference by computing the cross-correlation
#' and applying a lag-based linear interpolated shift. The shift is applied only
#' if the cross-correlation exceeds \code{clim}. This function is intended for
#' fine-tuning alignment across short ppm windows (i.e., segments).
#'
#' @importFrom stats ccf
#' @keywords internal
.alignSegment <- function(seg, idx_ref = 1, clim = 0.7,
                         med = TRUE, norm = FALSE,
                         lag.max = 20) {

  if (!is.matrix(seg) || nrow(seg) < 2)
    stop("Input must be a matrix with at least two rows.", call. = FALSE)

  seg <- .dimX(seg)

  ref <- if (med) apply(seg, 2, median) else seg[idx_ref, ]

  if (norm) {
    seg <- t(apply(seg, 1, scale))
    ref <- as.numeric(scale(ref))
  }

  out <- t(vapply(seq_len(nrow(seg)), function(i) {

    cc <- ccf(ref, seg[i, ],
              type = "correlation",
              lag.max = lag.max,
              plot = FALSE)

    acf <- as.vector(cc$acf)
    lags <- as.vector(cc$lag)
    imax <- which.max(acf)

    cval <- acf[imax]
    k <- lags[imax]

    if (is.finite(cval) && cval >= clim)
      .shift_interp(seg[i, ], k)
    else
      seg[i, ]

  }, FUN.VALUE = numeric(ncol(seg))))

  out
}

#' @title Linear Interpolated Shift (Internal)
#'
#' @description
#' Applies a linear shift to a numeric vector using interpolation.
#' Values outside the original domain are extended
#' using the boundary values (no zero padding).
#'
#' @param x Numeric vector.
#' @param lag Numeric scalar. Shift in index units. Positive values shift the
#'   signal to the right, negative values shift to the left.
#'
#' @return Numeric vector of the same length as \code{x}.
#'
#' @keywords internal
.shift_interp <- function(x, lag, ppm = NULL) {

  n <- length(x)

  if (is.null(ppm)) {
    idx <- seq_len(n)
    f <- stats::approxfun(idx, x, rule = 2)
    f(idx - lag)
  } else {
    f <- stats::approxfun(ppm, x, rule = 2)
    f(ppm - lag)
  }
}


#' @title Align NMR Spectra in a Selected Shift Region
#'
#' @description
#' Aligns spectra within a specified ppm window using cross-correlation.
#' Only the selected region is aligned; the remainder of the spectra
#' is unchanged.
#'
#' @param dat Named list with elements:
#'   \describe{
#'     \item{X}{Numeric matrix (spectra in rows)}
#'     \item{ppm}{Numeric vector of chemical shift values}
#'     \item{meta}{Optional metadata}
#'   }
#' @param shift Numeric vector of length 2 specifying ppm region to align.
#' @param idx_ref Integer. Row index to use as reference spectrum.
#' @param med Logical. Use row-wise median spectrum as reference?
#' @param clim Numeric. Minimum correlation threshold.
#' @param norm Logical. Z-scale before alignment.
#' @param lag.max Integer. Maximum lag allowed.
#'
#' @return Updated \code{dat} list with aligned spectra in selected region.
#'
#' @family preprocessing
#' @export
#' @examples
#' data(hiit_raw)
#' plot_spec(hiit_raw, shift=c(1.3,1.4))
#' hiit_aligned <- align_segment(hiit_raw, c(1.3, 1.35))
#' plot_spec(hiit_aligned, shift=c(1.3,1.4)) # aligned segment
#' plot_spec(hiit_aligned, shift=c(3, 3.1)) # segment not aligned
align_segment <- function(dat,
                         shift,
                         idx_ref = 1,
                         med = TRUE,
                         clim = 0.7,
                         norm = FALSE,
                         lag.max = 20) {

  if (!is.list(dat) || !all(c("X", "ppm") %in% names(dat)))
    stop("dat must be a list with elements 'X' and 'ppm'.", call. = FALSE)

  X   <- dat$X
  ppm <- dat$ppm

  if (!is.matrix(X))
    stop("dat$X must be a matrix.", call. = FALSE)

  if (ncol(X) != length(ppm))
    stop("ppm length must match ncol(X).", call. = FALSE)

  if (ppm[1] < ppm[length(ppm)]) {
    ord <- order(ppm, decreasing = TRUE)
    ppm <- ppm[ord]
    X   <- X[, ord, drop = FALSE]
  }

  idx <- get_idx(shift, ppm)

  if (length(idx) < 3)
    stop("Selected shift window too small.", call. = FALSE)

  # extract segment
  seg <- X[, idx, drop = FALSE]

  aligned_seg <- .alignSegment(
    seg      = seg,
    idx_ref  = idx_ref,
    clim     = clim,
    med      = med,
    norm     = norm,
    lag.max  = lag.max
  )

  X[, idx] <- aligned_seg

  dat$X   <- X
  dat$ppm <- ppm

  dat
}


#' @title Dynamic Interval Segmentation (RSPA-style)
#'
#' @description
#' Defines peak-system intervals from the median spectrum for
#' recursive segment-wise alignment (RSPA).
#'
#' @param X Numeric matrix (spectra in rows).
#' @param ppm Numeric ppm vector.
#' @param half_win_ppm Numeric. Half window size around detected peaks.
#'
#' @return List of index vectors defining alignment intervals.
#'
#' @family NMR
#' @keywords internal
.dynamicIntervalsMedian <- function(X, ppm,
                                    half_win_ppm = 0.02,
                                    min_snr = 8) {

  ref <- matrix(apply(X, 2, median), nrow = 1)

  pk <- ppick2(ref,
               ppm = ppm,
               type = "max",
               min_snr = min_snr,
               min_distance_ppm = 0.005)[[1]]

  if (is.null(pk) || nrow(pk) < 1L)
    return(list(seq_along(ppm)))

  intervals <- lapply(pk$idc, function(i0) {
    get_idx(c(ppm[i0] - half_win_ppm,
              ppm[i0] + half_win_ppm), ppm)
  })

  intervals <- intervals[order(vapply(intervals, `[`, 1L, FUN.VALUE = 1L))]

  merged <- list(intervals[[1]])

  for (i in 2:length(intervals)) {
    cur <- intervals[[i]]
    last <- merged[[length(merged)]]

    if (min(cur) <= max(last)) {
      merged[[length(merged)]] <- sort(unique(c(last, cur)))
    } else {
      merged[[length(merged) + 1L]] <- cur
    }
  }

  merged
}


#' @title Cohort-Guided Interval Alignment for 1D NMR Spectra
#'
#' @description
#' Aligns 1D NMR spectra using signal-adaptive ppm intervals derived from the
#' cohort median spectrum. Intervals are constructed around median peak systems
#' and aligned locally via cross-correlation.
#' The function accepts either a metabom8-style `dat` list or plain matrix/vector with `ppm`.
#'
#' @param x A numeric matrix/vector of spectra or a named list with elements
#'   `X`, `ppm`, and optionally `meta`.
#' @param ppm Numeric vector of chemical shift values (ppm). Ignored if `x` is
#'   a `dat` list.
#' @param half_win_ppm Numeric scalar controlling alignment interval width.
#' @param lag.max Integer. Maximum allowed lag for cross-correlation alignment.
#'
#' @details
#' Alignment intervals are derived from peaks in the cohort median spectrum.
#'
#' The `half_win_ppm` parameter controls interval granularity:
#'
#' \itemize{
#'   \item Smaller values (e.g. 0.005–0.007 ppm) generate more, narrower
#'   intervals (RSPA-like behaviour).
#'   \item Larger values (e.g. 0.008–0.015 ppm) merge nearby multiplets into
#'   broader regions (icoshift-like behaviour).
#' }
#'
#' Typical 600–800 MHz \eqn{^1}H NMR data:
#' \itemize{
#'   \item 0.006–0.008 ppm: stable multiplet-level alignment
#'   \item <0.005 ppm: may split J-coupled systems
#'   \item >0.015 ppm: may merge unrelated resonances
#' }
#'
#' Alignment parameters and interval definitions are recorded in
#' \code{attr(X, "m8_prep")}.
#'
#' @return
#' If input is `dat`, returns updated `dat`. Otherwise returns aligned matrix.
#'
#' @seealso \code{\link{align_segment}}
#' @family preprocessing
#' @examples
#' data(hiit_raw)
#' plot_spec(hiit_raw, shift=c(1.3,1.4))
#' hiit_aligns <- align_spectra(hiit_raw)
#' plot_spec(hiit_aligns, shift=c(1.3,1.4)) # aligned segment
#' plot_spec(hiit_aligns, shift=c(3, 3.1)) # aligned segment
#' @export
align_spectra <- function(x,
                                ppm = NULL,
                                half_win_ppm = 0.007,
                                lag.max = 200) {

  inp <- x
  u   <- .m8_unpack_dat(x, ppm = ppm)

  X   <- .dimX(u$X)
  ppm <- u$ppm
  X0  <- X

  if (!.check_X_ppm(X, ppm))
    stop("Dimensions of X and ppm must match.", call. = FALSE)

  ints <- .dynamicIntervalsMedian(X, ppm,
                                  half_win_ppm = half_win_ppm)

  bounds <- lapply(ints, function(i)
    c(ppm[i[1]], ppm[i[length(i)]]))

  for (int in ints) {

    seg <- X[, int, drop = FALSE]
    seg <- .alignSegment(seg,
                        med = TRUE,
                        lag.max = lag.max)

    X[, int] <- seg
  }

  colnames(X) <- formatC(ppm, format = "fg", digits = 10)
  rownames(X) <- rownames(X0)

  X <- .m8_copy_attrs(X0, X)
  attr(X, "m8_axis") <- list(ppm = ppm)

  X <- .m8_stamp(
    X,
    step = "align spectra",
    params = list(
      half_win_ppm = half_win_ppm,
      n_intervals  = length(ints),
      lag.max      = lag.max,
      intervals_ppm = bounds
    ),
    notes = "Median-derived chemical shift intervals aligned independently using cross-correlation."
  )

  if (is.list(inp) && !is.null(inp$X)) {
    inp$X <- X
    inp$ppm <- ppm
    return(inp)
  }

  X
}
