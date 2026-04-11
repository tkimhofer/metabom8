#' @title Spectral data binning
#' @description
#' Equidistant binning of spectra by summarising intensities within ppm bins.
#'
#' @param X Numeric matrix or data frame with spectra in rows, or a named list as
#'   returned by \code{\link{read1d}}/\code{read1d_proc} containing \code{X} and \code{ppm}.
#' @param ppm Numeric vector of chemical shift positions (length must match \code{ncol(X)}).
#'   If \code{NULL}, \code{ppm} is inferred in the following order:
#'   \enumerate{
#'     \item \code{X$ppm} if \code{X} is a list input,
#'     \item \code{attr(X, "m8_axis")$ppm} (if present),
#'     \item numeric \code{colnames(X)} (if present).
#'   }
#' @param width Numeric. Bin size in ppm, or \code{NULL} if \code{npoints} is specified.
#' @param npoints Integer. Desired number of bins per spectrum, or \code{NULL} if \code{width} is specified.
#'   If both are provided, \code{npoints} is used.
#' @param fun Function. Summary function applied to intensities within each bin.
#'   Must return a single numeric value (e.g. \code{sum}, \code{mean}, \code{max}).
#'
#' @details
#' If present, preprocessing provenance is appended to \code{attr(X, "m8_prep")} using
#' \code{.m8_stamp()}. The ppm axis is updated in \code{attr(X, "m8_axis")$ppm} and
#' column names are set to the bin centres.
#'
#' When \code{width} is specified, spectra are interpolated onto a regular ppm grid and
#' then aggregated within bins (\code{interp = TRUE} in provenance). When \code{npoints}
#' is specified, aggregation is performed by index bins on the original grid
#' (\code{interp = FALSE} in provenance).
#'
#' @return Numeric matrix with spectra in rows and binned ppm variables in columns.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(2 * 100), nrow = 2)
#' ppm <- round(seq(10, 0.5, length.out = 100), 3)
#' colnames(X) <- ppm
#'
#' Xb <- binning(X, ppm, width = 0.5)
#' Xb_mean <- binning(X, ppm, width = 0.5, fun = mean)
#' dim(Xb)
#'
#' @family preprocessing
#' @export
binning <- function(X, ppm = NULL, width = NULL, npoints = NULL, fun = sum) {

  .return_list <- is.list(X) && !is.null(X$X)
  meta <- NULL

  if (.return_list) {
    u <- .m8_unpack_dat(X, ppm = ppm)
    X   <- u$X
    ppm <- u$ppm
    meta <- u$meta
  }

  X <- .dimX(X)

  if (!is.function(fun))
    stop("`fun` must be a function.", call. = FALSE)
  .test <- fun(c(1, 2))
  if (!is.numeric(.test) || length(.test) != 1L)
    stop("`fun` must return a single numeric value.", call. = FALSE)

  X0 <- X

  if (is.null(ppm)) {
    ax <- attr(X, "m8_axis", exact = TRUE)
    if (is.list(ax) && !is.null(ax$ppm)) {
      ppm <- ax$ppm
    } else if (!is.null(colnames(X))) {
      ppm <- as.numeric(colnames(X))
      if (anyNA(ppm)) ppm <- NULL
    }
  }
  if (is.null(ppm))
    stop("`ppm` not provided and not found in list input, attr(X,'m8_axis'), or numeric colnames(X).", call. = FALSE)

  ppm <- as.numeric(ppm)
  if (!.check_X_ppm(X, ppm))
    stop("Non-matching dimensions X matrix and ppm vector or missing values in ppm.", call. = FALSE)

  if (!is.null(width) && !is.null(npoints)) width <- NULL
  if (is.null(width) && is.null(npoints))
    stop("Define bin width in ppm or desired number of bins.", call. = FALSE)

  ord <- order(ppm, decreasing = TRUE)
  ppm <- ppm[ord]
  X <- as.matrix(X)[, ord, drop = FALSE]

  ppm_in <- ppm
  axis_in <- list(n = length(ppm_in), range = range(ppm_in, finite = TRUE))

  .restore_names <- function(out, X_in) {
    rownames(out) <- rownames(X_in)
    out
  }

  if (!is.null(width)) {

    dppm <- abs(diff(ppm))
    res <- stats::median(dppm, na.rm = TRUE)

    if (!is.finite(res) || res <= 0) stop("Cannot determine ppm grid spacing.", call. = FALSE)
    if (width <= res) stop("Bin width is smaller or equal to the ppm grid spacing.", call. = FALSE)

    step <- max(1L, as.integer(round(width / res)))
    new_res <- width / step

    ppm_new <- seq(from = max(ppm), to = min(ppm), by = -new_res)

    nbins <- floor(length(ppm_new) / step)
    if (nbins < 1L) stop("Binning width too large for the ppm range.", call. = FALSE)

    ybin <- rep(seq_len(nbins), each = step)[seq_len(nbins * step)]
    ppm_new <- ppm_new[seq_along(ybin)]

    Xb <- t(apply(X, 1, function(x) {
      sInter <- stats::approxfun(ppm, x, rule = 2)
      s <- sInter(ppm_new)
      vapply(seq_len(nbins), function(i) fun(s[ybin == i]), numeric(1))
    }))

    ppm_bin <- vapply(seq_len(nbins), function(i) mean(ppm_new[ybin == i]), numeric(1))
    colnames(Xb) <- formatC(ppm_bin, format = "fg", digits = 10)
    Xb <- .restore_names(Xb, X)

    axis_out <- list(n = length(ppm_bin), range = range(ppm_bin, finite = TRUE))

    Xb <- .m8_copy_attrs(X0, Xb)
    attr(Xb, "m8_axis") <- list(ppm = ppm_bin)
    Xb <- .m8_stamp(
      Xb,
      step = "binning",
      params = list(
        width = width,
        npoints = NULL,
        fun = deparse(substitute(fun)),
        interp = TRUE,
        grid_res_est = res,
        step_points = step,
        axis_in = axis_in,
        axis_out = axis_out
      ),
      notes = "Binning via interpolation to a regular grid followed by aggregation within bins."
    )

    if (.return_list) {
      return(list(X = Xb, ppm = ppm_bin, meta = meta))
    } else {
      return(Xb)
    }
  }

  if (!is.null(npoints)) {

    npoints <- as.integer(npoints)
    if (npoints >= length(ppm)) stop("npoints cannot be larger or equal than length(ppm).", call. = FALSE)
    if (npoints < 2L) stop("npoints must be >= 2.", call. = FALSE)

    breaks <- floor(seq(0, length(ppm), length.out = npoints + 1L))
    breaks[1] <- 0L
    nbins <- npoints

    Xb <- t(apply(X, 1, function(s) {
      vapply(seq_len(nbins), function(i) {
        idx <- (breaks[i] + 1L):breaks[i + 1L]
        fun(s[idx])
      }, numeric(1))
    }))

    ppm_bin <- vapply(seq_len(nbins), function(i) {
      idx <- (breaks[i] + 1L):breaks[i + 1L]
      mean(ppm[idx])
    }, numeric(1))

    colnames(Xb) <- formatC(ppm_bin, format = "fg", digits = 10)
    Xb <- .restore_names(Xb, X)

    axis_out <- list(n = length(ppm_bin), range = range(ppm_bin, finite = TRUE))

    Xb <- .m8_copy_attrs(X0, Xb)
    attr(Xb, "m8_axis") <- list(ppm = ppm_bin)
    Xb <- .m8_stamp(
      Xb,
      step = "binning",
      params = list(
        width = NULL,
        npoints = npoints,
        fun = deparse(substitute(fun)),
        interp = FALSE,
        axis_in = axis_in,
        axis_out = axis_out
      ),
      notes = "Binning by index aggregation on the original grid."
    )

    if (.return_list) {
      return(list(X = Xb, ppm = ppm_bin, meta = meta))
    } else {
      return(Xb)
    }
  }

  NULL
}
