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
#' @seealso \code{\link[signal]{sgolayfilt}}
#' @examples
#' data(covid)
#' X <- covid$X
#' ppm <- covid$ppm
#'
#' peaks <- ppick(X, ppm)
#' spec(X[1, ], ppm, shift = c(4.2, 5.3))
#' points(peaks[[1]]$ppm, peaks[[1]]$Int, col = factor(peaks[[1]]$Etype))
#' 
#' @importFrom utils head tail
#' @export
ppick <- function(X, ppm, fil_p = 3, fil_n = 5, type = "max") {
  X <- .dimX(X)

  if (!.check_X_ppm(X, ppm)) {
    stop("Dimensions of X and ppm must match and contain only numeric values.")
  }

  if (!(type %in% c("max", "min", "both"))) {
    stop("Parameter 'type' must be one of: 'max', 'min', 'both'.")
  }

  apply(X, 1, function(x, pp = ppm) {
    x[is.na(x)] <- 0
    x_smooth <- sgolayfilt(x, p = fil_p, n = fil_n)
    x_diff <- sign(diff(x_smooth))

    extrema_idx <- which(head(x_diff, -1) != tail(x_diff, -1)) + 1
    if (length(extrema_idx) == 0) return(NULL)

    # Correct indexing here
    extrema_type <- x_diff[extrema_idx - 1]

    res <- data.frame(
      idc = extrema_idx,
      ppm = pp[extrema_idx],
      Int = x[extrema_idx],
      Etype = extrema_type
    )

    switch(type,
           "max" = res[res$Etype > 0, ],
           "min" = res[res$Etype < 0, ],
           "both" = res
    )
  })
}

