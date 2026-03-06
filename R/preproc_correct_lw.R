#' Linewidth correction by scaling spectra to a reference linewidth
#'
#' Applies a multiplicative correction to each spectrum to compensate for
#' linewidth-induced peak height variation, using a reference peak region
#' (e.g. TSP). The correction assumes an empirical power-law relationship
#' between peak height and linewidth:
#'
#' \deqn{I \propto \mathrm{FWHM}^{-\beta}}
#'
#' Spectra are scaled such that all samples are adjusted to a common reference
#' linewidth \code{lw_ref}:
#'
#' \deqn{X_{\mathrm{adj}} = X \cdot \left(\frac{\mathrm{FWHM}}{lw_{\mathrm{ref}}}\right)^{\beta}}
#'
#' where FWHM is estimated per spectrum within the specified \code{shift}
#' region (typically containing an internal reference signal).
#'
#' If \code{estimate_beta = TRUE}, the exponent \code{beta} is estimated from a
#' log-log regression between peak height (max within \code{shift}) and FWHM
#' across samples. Otherwise, the user-supplied \code{beta} is used.
#'
#' The correction is only applied if it reduces the coefficient of variation (CV)
#' of the reference peak height across samples when \code{only_if_improves = TRUE}.
#'
#' Input may be either:
#' \itemize{
#'   \item a numeric matrix \code{X} (samples x variables), with ppm axis supplied
#'         via \code{attr(X, "m8_axis")$ppm} or numeric \code{colnames(X)}, or
#'   \item a named list with elements \code{X}, \code{ppm}, and optionally \code{meta}.
#' }
#'
#' Attributes \code{"m8_prep"}, \code{"m8_axis"}, and \code{"m8_meta"} are
#' preserved and the applied correction is appended to the preprocessing log
#' via \code{.m8_stamp()}.
#'
#' @param x A numeric matrix (samples x variables) or a named list containing
#'   \code{X}, \code{ppm}, and optionally \code{meta}.
#' @param lw_ref Numeric reference linewidth (FWHM) to which spectra are
#'   normalized. If \code{NULL}, the median linewidth across samples is used.
#' @param beta Numeric exponent governing the linewidth-to-height relationship.
#'   Ignored if \code{estimate_beta = TRUE} and sufficient samples are available.
#' @param shift Numeric vector of length 2 specifying the ppm region used for
#'   linewidth estimation and peak height measurement (e.g. \code{c(-0.1, 0.1)}
#'   for TSP).
#' @param ppm Optional numeric ppm axis. If \code{NULL}, inferred from
#'   \code{attr(X, "m8_axis")$ppm} or numeric \code{colnames(X)}.
#' @param sf Optional spectrometer frequency (MHz). If \code{NULL}, obtained from
#'   \code{meta$a_SFO1} in list input or \code{attr(X, "m8_meta")$a_SFO1}.
#' @param estimate_beta Logical; if \code{TRUE}, estimate \code{beta} from a
#'   log-log regression of reference peak height versus FWHM.
#' @param beta_min_n Minimum number of valid samples required to estimate
#'   \code{beta}.
#' @param only_if_improves Logical; if \code{TRUE}, apply correction only if the
#'   CV of the reference peak height decreases after correction.
#' @param notes Optional character string appended to the preprocessing log.
#'
#' @return An object of the same type as input:
#'   \itemize{
#'     \item If input is a matrix, returns a corrected matrix with attributes
#'           preserved and updated.
#'     \item If input is a list, returns the same list with corrected \code{$X}.
#'   }
#'
#' @details
#' This correction removes systematic peak height variation due to spectral
#' broadening (e.g. shimming differences) under the assumption that the internal
#' reference signal has constant concentration across samples. The exponent
#' \code{beta} reflects the empirical sensitivity of peak height to linewidth
#' in the current acquisition and processing regime.
#'
#' This adjustment improves comparability of signal intensities across spectra
#' by reducing linewidth-induced amplitude bias, thereby enhancing the
#' detectability of small effects in downstream multivariate analyses
#' (e.g. PLS-type models).
#'
#' @seealso \code{\link{lw}}
#'
#' @examples
#' data(covid_raw)
#' # matrix input with ppm in m8_axis
#' Xcor <- correct_lw(covid_raw, shift = c(-0.1, 0.1))
#'
#' @family preprocessing
#' @export

correct_lw <- function(x,
                       lw_ref = 0.8,
                       beta = 0.6,
                       shift = c(-0.1, 0.1),     # region used for FWHM and intensity (e.g. TSP)
                       ppm = NULL,
                       sf = NULL,
                       estimate_beta = TRUE,
                       beta_min_n = 6,
                       only_if_improves = TRUE,
                       notes = NULL) {

  cv <- function(v) stats::sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE)

  is_list_input <- is.list(x) && !is.null(x$X)
  if (is_list_input) {
    dat <- x
    X <- dat$X
    if (is.null(ppm)) ppm <- dat$ppm
    meta <- dat$meta
  } else {
    X <- x
    meta <- NULL
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
  if (is.null(ppm)) stop("correct_lw(): `ppm` not provided and not found in attr(X,'m8_axis') or numeric colnames(X).")
  ppm <- as.numeric(ppm)
  if (!.check_X_ppm(X, ppm)) stop("correct_lw(): X and ppm mismatch or ppm contains NA/Inf.")

  if (is.null(sf)) {
    if (is_list_input && !is.null(meta) && "a_SFO1" %in% names(meta)) {
      sf <- meta$a_SFO1
    } else {
      m <- attr(X, "m8_meta", exact = TRUE)
      if (is.data.frame(m) && "a_SFO1" %in% names(m)) sf <- m$a_SFO1
    }
  }

  fwhm <- lw(X, ppm = ppm, shift = shift, sf = sf)

  idx <- get_idx(shift, ppm)
  if (length(idx) < 3L) stop("correct_lw(): selected `shift` region too small.")

  I_ref <- apply(X[, idx, drop = FALSE], 1, max, na.rm = TRUE)

  beta_used <- beta
  beta_source <- "user"

  if (isTRUE(estimate_beta)) {
    ok <- is.finite(I_ref) & I_ref > 0 & is.finite(fwhm) & fwhm > 0
    if (sum(ok) >= beta_min_n) {
      fit <- stats::lm(log(I_ref[ok]) ~ log(fwhm[ok]))
      bhat <- -stats::coef(fit)[2]
      if (is.finite(bhat) && bhat > 0) {
        beta_used <- as.numeric(bhat)
        beta_source <- "estimated"
      }
    }
  }

  lw_adj <- (fwhm / lw_ref)^beta_used
  X_new <- sweep(X, 1, lw_adj, `*`)

  I_ref_new <- apply(X_new[, idx, drop = FALSE], 1, max, na.rm = TRUE)

  improved <- is.finite(cv(I_ref)) && is.finite(cv(I_ref_new)) && (cv(I_ref_new) < cv(I_ref))

  if (only_if_improves && !improved) {
    return(x)
  }

  X_new <- .m8_copy_attrs(from = X, to = X_new, keys = c("m8_prep", "m8_axis", "m8_meta"))

  step_notes <- if (!is.null(notes)) notes else {
    "Scaled spectra by (FWHM/lw_ref)^beta using FWHM estimated in the reference shift window."
  }

  X_new <- .m8_stamp(
    X_new,
    step = "correct_lw",
    params = list(
      shift = shift,
      lw_ref = lw_ref,
      beta = beta_used,
      beta_source = beta_source,
      cv_before = cv(I_ref),
      cv_after = cv(I_ref_new)
    ),
    notes = step_notes
  )

  if (is_list_input) {
    out <- x
    out$X <- X_new
    if (is.null(out$ppm)) out$ppm <- ppm
    return(out)
  }

  X_new
}
