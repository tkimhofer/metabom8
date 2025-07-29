#' Align NMR Spectra to a Reference Spectrum via Cross-Correlation
#'
#' @description
#' Aligns a matrix of NMR spectral segments by x-shifting each row to maximize cross-correlation with a reference spectrum.
#'
#' @param seg A matrix or data frame with spectra as rows.
#' @param idx_ref Integer. Row index of the reference spectrum (ignored if \code{med = TRUE}).
#' @param clim Numeric. Cross-correlation similarity threshold. Segments below this threshold are not aligned.
#' @param med Logical. Use the median spectrum as reference? If \code{TRUE}, \code{idx_ref} is ignored.
#' @param norm Logical. If \code{TRUE}, min-max normalization is applied to each row before alignment.
#'
#' @return A matrix of aligned spectra (rows).
#'
#' @details
#' Spectra are aligned by computing the cross-correlation with a reference and applying a circular shift based on the lag that maximizes similarity.
#' If \code{med = TRUE}, the reference is the row-wise median spectrum. Segments with a cross-correlation below \code{clim} are left unshifted.
#'
#' @importFrom stats ccf
#'
#' @author \email{t.kimhofer@@gmail.com}
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(1000), nrow = 10)
#' x[5, ] <- x[1, ] + 0.1  # create similar spectrum
#' aligned <- alignSegment(x)
#'
#' @export
alignSegment <- function(seg, idx_ref = 1, clim = 0.7, med = TRUE, norm = FALSE) {
  if (is.null(nrow(seg)) || nrow(seg) < 2) stop("Check input dimensions.")

  if (med[1]) {
    seg <- rbind(apply(seg, 2, median), seg)
    idx_ref <- 1
  }

  if (norm) {
    seg <- t(apply(seg, 1, minmax))
  }

  out <- vapply(seq(nrow(seg)), function(i, le = ncol(seg)) {
    cc_mod <- ccf(seg[i, ], seg[idx_ref, ], type = "correlation", plot = FALSE)
    acf <- as.vector(cc_mod$acf)
    idx <- which.max(acf)
    cval <- acf[idx]
    xlag <- as.vector(cc_mod$lag)[idx]

    if (cval > clim) {
      if (xlag > 0) {
        dd <- c(seg[i, -seq_len(abs(xlag))], rep(0, abs(xlag)))
      } else if (xlag < 0) {
        dd <- c(rep(0, abs(xlag)), seg[i, -seq((le - abs(xlag) + 1), le)])
      } else {
        dd <- seg[i, ]
      }
      return(dd)
    } else {
      seg[i, ]
    }
  }, FUN.VALUE = seg[1, ])

  out <- t(out)
  if (med) out <- out[-1, ]

  return(out)
}
