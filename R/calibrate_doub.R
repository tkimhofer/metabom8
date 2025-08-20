#' @title Calibrate spectra to a doublet signal
#' @description
#' Calibrate 1D NMR spectra to a doublet peak (e.g. glucose or alanine) by aligning to a known coupling pattern.
#'
#' @param X Numeric matrix. NMR data with spectra in rows.
#' @param ppm Numeric vector. Chemical shift positions corresponding to columns of X.
#' @param type Character or numeric. Use \code{"glu"} (glucose) or \code{"ala"} (alanine) for predefined doublet regions, or a numeric vector with custom ppm range.
#' @param j_const Numeric vector. Expected J-coupling range in ppm (default for glucose).
#' @param sg_length Integer. Smoothing window length for Savitzky-Golay filter.
#'
#' @return Calibrated numeric matrix of the same dimensions as X.
#'
#' @keywords internal
#'
#' @importFrom signal sgolayfilt
#' @importFrom utils combn
#' @importFrom ptw asysm
.calibrate_doub <- function(X, ppm, type = c("glu", "ala"), j_const = c(0.006, 0.007), sg_length = 13) {

  if (type[1] == "glu") {
    idx <- get_idx(c(5.15, 5.3), ppm)
    cent_loc <- 5.233
    j_const <- c(0.006, 0.007)
  } else if (type[1] == "ala") {
    idx <- get_idx(c(1.4, 1.56), ppm)
    cent_loc <- 1.48
    j_const <- c(0.0115, 0.0135)
  } else if (is.numeric(type[1])) {
    idx <- get_idx(type, ppm)
    cent_loc <- mean(type)
  } else {
    stop("Unknown calibration type. Use 'glu', 'ala', or numeric ppm range.")
  }

  test <- apply(X[, idx, drop = FALSE], 1, function(x, pp = ppm[idx]) {
    smoothed <- sgolayfilt(x - asysm(x, lambda = 100, maxit=1000), p = 3, n = sg_length)

    picks <- ppick(smoothed, pp, type = "max")

    if (is.null(picks) || length(picks) == 0 || is.null(picks[[1]])) {
      return(NULL)
    } else {
      return(picks[[1]])
    }
  })

  if (all(is.null(test))) {
    warning('Could not calibrate - no signals detected');
    return(NULL)
    }

  s <- lapply(seq_along(test), function(i) {
    peaks <- test[[i]]

    if (is.null(peaks) || nrow(peaks) < 2){
      warning('Could not calibrate spectrum')
      return(NULL)
    }
    ptab <- peaks[peaks$Etype > 0, , drop = FALSE]

    if (nrow(ptab) < 2) return(NULL)

    combs <- combn(nrow(ptab), 2)
    ppm_diffs <- apply(combs, 2, function(id) abs(diff(ptab$ppm[id])))
    valid <- which(ppm_diffs > j_const[1] & ppm_diffs < j_const[2])

    if (length(valid) == 0) return(NULL)

    if (length(valid) == 1) {
      out <- ptab[combs[, valid], ]
    } else {
      crit <- apply(combs[, valid], 2, function(id) {
        avg <- mean(ptab$Int[id])
        ratio <- max(ptab$Int) / min(ptab$Int[id])
        c(avg = avg, ratio = ratio)
      })
      idx_best <- which.max(crit["avg", ])
      out <- ptab[combs[, valid[idx_best]], ]
    }

    out$diff <- diff(out$ppm)
    out$mm <- mean(out$ppm)
    return(out)
  })

  Xc <- vapply(seq_len(nrow(X)), function(i) {
    if (is.null(s[[i]]) || nrow(s[[i]]) < 2) {
      warning(sprintf("Could not calibrate spectrum %d", i))
      return(X[i, ])
    } else {
      shift_fun <- approxfun(x = ppm - (max(s[[i]]$ppm) - cent_loc), y = X[i, ])
      return(shift_fun(ppm))
    }
  }, FUN.VALUE = X[1, ])

  return(t(Xc))
}
