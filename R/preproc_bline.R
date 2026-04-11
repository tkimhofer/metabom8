#' Baseline Correction for Spectral Data
#'
#' Applies baseline correction to each spectrum (row) of a spectral matrix.
#' Multiple correction algorithms are available and selected via the
#' `method` argument.
#'
#' @param X Numeric matrix or metabom8 `dat` object containing spectra
#'   in rows and spectral variables in columns.
#' @param method Character specifying the baseline correction algorithm.
#'   One of:
#'   \describe{
#'     \item{`"asls"`}{Asymmetric least squares baseline estimation.}
#'     \item{`"linear"`}{Linear baseline estimated from edge regions.}
#'   }
#' @param ... Additional parameters passed to the selected method.
#'
#' **Arguments for `method = "asls"`**
#' \describe{
#'   \item{lambda}{Numeric smoothing parameter controlling baseline stiffness.
#'   Larger values produce smoother baselines. Default `1e7`.}
#'   \item{iter_max}{Maximum number of iterations used in the baseline
#'   estimation procedure. Default `30`.}
#' }
#'
#' **Arguments for `method = "linear"`**
#' \describe{
#'   \item{ppm}{Numeric vector describing the spectral axis (e.g. chemical shift).
#'   Must have length equal to `ncol(X)`.}
#'   \item{edge_frac}{Fraction of points at each spectrum edge used to
#'   estimate the baseline. Default `0.1`.}
#' }
#'
#' @return Numeric matrix containing baseline-corrected spectra with the
#'   same dimensions as `X`. If `X` is a metabom8 `dat` object, the corrected
#'   matrix replaces the `X` component while preserving associated metadata.
#'
#' @details
#' Baseline correction is performed independently for each spectrum.
#'
#' The `"asls"` method estimates a smooth baseline using asymmetric least
#' squares smoothing, which is well suited for spectra with positive peaks
#' such as NMR metabolomics data.
#'
#' The `"linear"` method estimates a straight baseline from the edges of
#' each spectrum and subtracts the fitted trend.
#'
#' The returned object includes a `.stamp` attribute recording the
#' baseline correction method and parameters used.
#'
#' @references
#' Eilers PHC, Boelens HFM (2005).
#' Baseline correction with asymmetric least squares smoothing.
#'
#' @family preprocessing
#'
#' @examples
#' data(hiit_raw)
#'
#' plot_spec(hiit_raw$X[1,], hiit_raw$ppm, shift=c(3.1,4), backend='base')
#' hiit_proc <-
#'   hiit_raw |>
#'   excise() |>
#'   correct_baseline()
#'
#' plot_spec(hiit_proc$X[1,], hiit_proc$ppm, shift=c(3.1,4), backend='base', add=TRUE, col='red')
#'
#' @export
correct_baseline <- function(X, method = c("asls", "linear"), ...) {

  method <- match.arg(method)

  inp <- X
  u   <- .m8_unpack_dat(X)

  X0   <- u$X
  Xmat <- .dimX(X0)
  ppm  <- u$ppm
  meta <- u$meta

  args <- list(...)

  Xbc <- switch(
    method,
    linear = .linear_baseline(Xmat, ppm, ...),
    asls   = .baseline_asls(Xmat, ...)
  )

  Xbc <- .m8_copy_attrs(X0, Xbc)

  if (!is.null(ppm))
    attr(Xbc, "m8_axis") <- list(ppm = ppm)

  Xbc <- .m8_stamp(
    Xbc,
    step   = "baseline",
    params = list(method = method, args = args),
    notes  = "Baseline correction applied."
  )

  if (is.list(inp) && !is.null(inp$X)) {
    inp$X    <- Xbc
    inp$ppm  <- ppm
    inp$meta <- meta
    return(inp)
  }

  Xbc
}

#'NL Baseline Correction for 1D NMR Spectra
#' The implementation relies on \code{\link[ptw]{asysm}} from the \pkg{ptw} package.
#' @noRd
.baseline_asls <- function(X, lambda = 1e7, iter_max = 60) {
  if (any(is.na(X))) {
    message("X contains missing values, replacing with zeros.")
    X[is.na(X)] <- 0
  }

  if (!is.numeric(lambda) || lambda <= 0)
    stop("`lambda` must be positive numeric.", call. = FALSE)

  X0 <- X
  X <- .dimX(X)
  X.bl <- t(vapply(seq_len(nrow(X)), function(i) {
    x <- X[i, ]
    x - ptw::asysm(x, maxit = iter_max, lambda = lambda)
  }, numeric(ncol(X))))

  X.bl <- .m8_copy_attrs(X0, X.bl)
  X.bl <- .m8_stamp(
    X.bl,
    step = "baseline_asysm",
    params = list(lambda = lambda, iter_max = iter_max),
    notes = "Missing values replaced by 0 prior to baseline correction (if present)."
  )

  X.bl
}

#' @noRd
.linear_baseline <- function(X, ppm, edge_frac = 0.1) {

  n <- length(ppm)
  stopifnot(n == ncol(X))
  stopifnot(is.numeric(edge_frac), edge_frac > 0, edge_frac < 0.5)

  k <- floor(n * edge_frac)
  idx <- c(seq_len(k), seq(n - k + 1, n))

  A     <- cbind(1, ppm[idx])
  Afull <- cbind(1, ppm)

  beta <- qr.solve(A, t(X[, idx, drop = FALSE]))
  baseline <- t(Afull %*% beta)
  X - baseline
}

#' @rdname correct_baseline
#' @export
bline <- function(X, ...) {
  warning("`bline()` will be removed in a future version; please use `correct_baseline(method = \"asls\")` instead.",
          call. = FALSE)
  correct_baseline(X, method = "asls", ...)
}

