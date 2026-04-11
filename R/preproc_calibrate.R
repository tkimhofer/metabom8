#' @title Chemical Shift Calibration
#'
#' @description
#' Aligns 1D \eqn{^1}H NMR spectra to a reference signal on the existing ppm grid.
#' Supports singlet (e.g. TSP or custom range) and predefined doublet references
#' (glucose, alanine).
#'
#' @param X Numeric matrix/vector of spectra or a metabom8 data list.
#' @param ppm Numeric chemical shift vector. If \code{X} is a metabom8 list,
#'   this is taken from \code{X$ppm}.
#' @param type Character (\code{"tsp"}, \code{"glucose"}, \code{"alanine"}),
#' or list.
#'
#' @details
#' In addition to predefined references, custom calibration targets can be supplied
#' as a list with elements:
#' \itemize{
#'   \item \code{mode}: \code{"singlet"} or \code{"doublet"}
#'   \item \code{window}: numeric vector of length 2 defining the ppm search region
#'   \item \code{centre}: target ppm position (optional; defaults to mean(window))
#'   \item \code{j}: numeric vector of length 2 specifying expected J-coupling
#'         (required for doublet calibration)
#' }
#' Custom doublet calibration requires the expected J-coupling range (`j`)
#' in ppm to distinguish the two peaks of the multiplet.
#'
#' For example:
#' \preformatted{
#' calibrate(
#'   X, ppm,
#'   list(
#'     mode   = "doublet",
#'     window = c(1.2, 1.35),
#'     j      = c(0.007, 0.009)
#'   )
#' )
#' }
#' @return Calibrated spectra in the same structure as input.
#'
#' @examples
#' data("covid_raw")
#' X=covid_raw$X
#' ppm=covid_raw$ppm
#' X_tsp <- calibrate(X, ppm, type = "tsp")
#' X_glu <- calibrate(X, ppm, type = "glucose")
#' X_custom <- calibrate(X, ppm, type = c(1.9, 2.1))
#' @family preprocessing
#' @export
calibrate <- function(X, ppm = NULL, type = "tsp") {

  inp <- X
  u   <- .m8_unpack_dat(X, ppm = ppm)

  X0   <- u$X
  Xmat <- .dimX(X0)
  ppm  <- u$ppm
  meta <- u$meta

  if (is.null(ppm))
    stop("ppm must be supplied or inferable.", call. = FALSE)

  ppm <- as.numeric(ppm)

  target <- .resolve_calib_target(type)

  Xc <- switch(
    target$mode,
    singlet = .calibrate_singlet(
      Xmat,
      ppm,
      target$window,
      target$centre
    ),
    doublet = .calibrate_doublet(
      Xmat,
      ppm,
      window   = target$window,
      centre   = target$centre,
      j_const  = target$j
    )
  )
  Xc <- .m8_copy_attrs(X0, Xc)
  attr(Xc, "m8_axis") <- list(ppm = ppm)

  Xc <- .m8_stamp(
    Xc,
    step   = "calibrate",
    params = list(target = target),
    notes  = "Chemical shift calibration applied."
  )

  if (is.list(inp) && !is.null(inp$X)) {
    inp$X    <- Xc
    inp$ppm  <- ppm
    inp$meta <- meta
    return(inp)
  }

  Xc
}

#' Resolve calibration specification into a numeric target definition.
#' @keywords internal
#' @noRd
.resolve_calib_target <- function(type) {

  # tsp
  if (identical(type, "tsp")) {
    return(list(
      mode   = "singlet",
      window = c(-0.2, 0.2),
      centre = 0
    ))
  }

  # glucose
  if (identical(type, "glucose")) {
    return(list(
      mode   = "doublet",
      window = c(5.15, 5.30),
      centre = 5.233,
      j      = c(0.006, 0.007)
    ))
  }

  # alanine
  if (identical(type, "alanine")) {
    return(list(
      mode   = "doublet",
      window = c(1.40, 1.56),
      centre = 1.48,
      j      = c(0.0115, 0.0135)
    ))
  }

  # custom singlet
  if (is.numeric(type) && length(type) == 2L) {
    return(list(
      mode   = "singlet",
      window = type,
      centre = mean(type)
    ))
  }

  # custom full specification
  if (is.list(type) ) {

    if (!all(c("mode","window") %in% names(type)))
      stop("Custom calibration list must include 'mode' and 'window'.", call. = FALSE)

    if (!type$mode %in% c("singlet","doublet"))
      stop("`mode` must be 'singlet' or 'doublet'.", call. = FALSE)

    if (type$mode == "doublet" && is.null(type$j))
      stop("Custom doublet requires 'j'.", call. = FALSE)

    if (is.null(type$centre))
      type$centre <- mean(type$window)


    out <- list(
      mode   = type$mode,
      window = type$window,
      centre = type$centre
    )

    if (type$mode == "doublet")
      out$j <- type$j

    return(out)
  }

  stop("Unsupported calibration target.", call. = FALSE)
}

#' Align spectra to the maximum peak within a ppm interval.
#' @keywords internal
#' @noRd
.calibrate_singlet <- function(X, ppm, range, centre = mean(range)) {

  idx <- get_idx(range, ppm)
  if (length(idx) == 0L){
    warning("No ppm values found in calibration window; returning X unchanged.", call.=FALSE)
    return(X)
    }

  target <- centre
  zeroAbs <- idx[which.min(abs(ppm[idx] - target))]

  Xc <- matrix(0, nrow(X), ncol(X))

  for (i in seq_len(nrow(X))) {

    peakAbs <- idx[which.max(X[i, idx])]
    shift <- zeroAbs - peakAbs

    if (shift == 0L) {
      Xc[i, ] <- X[i, ]
    } else if (shift > 0L) {
      Xc[i, seq.int(shift + 1L, ncol(X))] <- X[i, seq_len(ncol(X) - shift)]
    } else {
      sh <- abs(shift)
      Xc[i, seq_len(ncol(X) - sh)] <- X[i, seq.int(sh + 1L, ncol(X))]
    }
  }

  Xc
}


#' Calibrate spectra to a known doublet using J-coupling constraints.
#' @importFrom signal sgolayfilt
#' @importFrom utils combn
#' @importFrom ptw asysm
#' @keywords internal
#' @noRd
.calibrate_doublet <- function(X,
                               ppm,
                               window,
                               centre,
                               j_const,
                               sg_length = 13) {

  idx <- get_idx(window, ppm)

  if (length(idx) < 3L) {
    warning("Calibration region contains too few points; returning X unchanged.", call. = FALSE)
    return(X)
  }

  picks_list <- lapply(seq_len(nrow(X)), function(i) {
    x <- X[i, idx]
    smoothed <- signal::sgolayfilt(x - ptw::asysm(x, lambda = 100, maxit = 1000),
                                   p = 3, n = sg_length)
    picks <- ppick(smoothed, ppm[idx], type = "max")
    if (is.null(picks) || length(picks) == 0L || is.null(picks[[1]])) return(NULL)
    picks[[1]]
  })

  if (all(vapply(picks_list, is.null, logical(1)))) {
    warning("Could not calibrate: no signals detected in calibration window; returning X unchanged.",
            call. = FALSE)
    return(X)
  }

  dbl <- lapply(seq_along(picks_list), function(i) {
    peaks <- picks_list[[i]]
    if (is.null(peaks) || nrow(peaks) < 2L) return(NULL)

    ptab <- peaks[peaks$Etype > 0, , drop = FALSE]
    if (nrow(ptab) < 2L) return(NULL)

    combs <- utils::combn(nrow(ptab), 2)
    ppm_diffs <- apply(combs, 2, function(id) abs(diff(ptab$ppm[id])))
    valid <- which(ppm_diffs > j_const[1] & ppm_diffs < j_const[2])
    if (length(valid) == 0L) return(NULL)

    if (length(valid) == 1L) {
      out <- ptab[combs[, valid], , drop = FALSE]
    } else {
      avg_int <- vapply(valid, function(v) {
        id <- combs[, v]
        mean(ptab$Int[id])
      }, numeric(1))
      out <- ptab[combs[, valid[which.max(avg_int)]], , drop = FALSE]
    }

    out$mm <- mean(out$ppm)
    out
  })

  Xc <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  for (i in seq_len(nrow(X))) {
    if (is.null(dbl[[i]]) || nrow(dbl[[i]]) < 2L) {
      Xc[i, ] <- X[i, ]
      next
    }

    delta <- (dbl[[i]]$mm - centre)

    shift_fun <- stats::approxfun(x = ppm - delta, y = X[i, ],
                                  rule = 2)  # rule=2 avoids NA at edges
    xi <- shift_fun(ppm)

    if (anyNA(xi)) xi[is.na(xi)] <- 0
    Xc[i, ] <- xi
  }

  rownames(Xc) <- rownames(X)
  colnames(Xc) <- colnames(X)

  Xc
}

