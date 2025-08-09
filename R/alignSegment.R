#' @title Align NMR Spectra via Cross-Correlation
#'
#' @description
#' Aligns rows of an NMR spectral matrix to a reference spectrum by maximizing cross-correlation.
#' Optionally uses the row-wise median spectrum as reference. Useful for minor spectral misalignments.
#'
#' @param seg Numeric matrix. Each row is a 1D NMR spectrum segment.
#' @param idx_ref Integer. Row index to use as reference spectrum. Ignored if \code{med = TRUE}.
#' @param clim Numeric. Minimum cross-correlation threshold. Segments with lower similarity are not shifted.
#' @param med Logical. If \code{TRUE}, the row-wise median spectrum is used as the reference.
#' @param norm Logical. If \code{TRUE}, rows are scaled to the \code{[0, 1]} range using min-max scaling before alignment.
#'
#' @return A numeric matrix of aligned spectra with the same dimensions as \code{seg}.
#' If \code{med = TRUE}, the output excludes the median row.
#'
#' @details
#' Each spectrum is aligned to the reference by computing the cross-correlation and applying a lag-based circular shift.
#' The shift is applied only if the cross-correlation exceeds \code{clim}.
#' This function is intended for fine-tuning alignment across short ppm windows (i.e., segments).
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(1000), nrow = 10)
#' x[5, ] <- x[1, ] + 0.1  # Simulate minor shift
#' aligned <- alignSegment(x, clim = 0.6)
#'
#' @seealso \code{\link{minmax}}
#' @family NMR
#' @importFrom stats ccf
#' @export
alignSegment <- function(seg, idx_ref = 1, clim = 0.7, med = TRUE, norm = FALSE) {
  if (!is.matrix(seg) || nrow(seg) < 2) stop("Input must be a matrix with at least two rows.")

  if (med) {
    seg <- rbind(apply(seg, 2, median), seg)
    idx_ref <- 1
  }

  if (norm) {
    seg <- t(apply(seg, 1, minmax))
  }

  out <- vapply(seq_len(nrow(seg)), function(i) {
    cc_mod <- ccf(seg[i, ], seg[idx_ref, ], type = "correlation", plot = FALSE)
    acf <- as.vector(cc_mod$acf)
    xlag <- as.vector(cc_mod$lag)[which.max(acf)]
    cval <- max(acf)

    if (cval > clim) {
      if (xlag > 0) {
        return(c(seg[i, -seq_len(xlag)], rep(0, xlag)))
      } else if (xlag < 0) {
        le <- ncol(seg)
        return(c(rep(0, -xlag), seg[i, seq_len(le + xlag)]))
      }
    }
    seg[i, ]
  }, FUN.VALUE = numeric(ncol(seg)))

  out <- t(out)
  if (med) out <- out[-1, , drop = FALSE]
  return(out)
}
