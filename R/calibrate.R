#' @title Chemical Shift Calibration
#'
#' @description
#' Chemical shift calibration using a reference signal for 1D NMR spectra.
#'
#' @param X Numeric matrix. NMR data matrix with spectra in rows and chemical shift positions in columns.
#' @param ppm Numeric vector. Chemical shift values corresponding to the columns of \code{X}.
#' @param type Either a character string or numeric vector. Options include:
#' \itemize{
#'   \item{\code{"tsp"}: calibrates to 0 ppm using the highest peak in the interval \code{[-0.2, 0.2]} ppm (TSP signal).}
#'   \item{\code{"glucose"}: calibrates to 5.23 ppm using glucose doublet (5.15–5.30 ppm).}
#'   \item{\code{"alanine"}: calibrates to 1.48 ppm using alanine doublet (1.40–1.53 ppm).}
#'   \item{Numeric vector of length 2: custom ppm range. Spectrum is shifted so the maximum peak in this interval aligns with \code{mean(type)}.}
#' }
#'
#' @return A numeric matrix of the same dimensions as \code{X}, with spectra calibrated to the chosen reference.
#'
#' @details
#' Calibration is typically part of the technical quality control pipeline and is commonly done using a reference signal such as TSP at 0 ppm.
#'
#' @references
#' Dona, A.C., \emph{et al.} (2014). Precision high-throughput proton NMR spectroscopy of human urine, serum, and plasma for large-scale metabolic phenotyping. \emph{Analytical Chemistry}, 86(19), 9887–9894.
#'
#' @examples
#' data(covid_raw)
#' X <- covid_raw$X
#' ppm <- covid_raw$ppm
#' X_tsp <- calibrate(X, ppm, type = "tsp")
#' X_glu <- calibrate(X, ppm, type = "glucose")
#' X_custom <- calibrate(X, ppm, type = c(1.9, 2.1))
#'
#' @seealso \code{\link{.calibrate_doub}}
#'
#' @export
calibrate <- function(X, ppm, type = c("tsp", "glucose", "alanine")) {

  if (missing(X) || missing(ppm)) {
    stop("Please provide both 'X' and 'ppm'.")
  }

  if (!is.matrix(X) || !is.numeric(ppm)) {
    stop("'X' must be a numeric matrix and 'ppm' must be a numeric vector.")
  }

  if (!.check_X_ppm(X, ppm)) {
    stop("The number of columns in 'X' must match the length of 'ppm'.")
  }

  if (is.character(type) && (length(type) != 1 || is.na(type))) {
    stop("If 'type' is a character, it must be one of 'tsp', 'glucose', or 'alanine'.")
  }

  if (is.numeric(type) && length(type) != 2) {
    stop("If 'type' is numeric, it must be a vector of length 2 indicating a ppm range.")
  }

  if (identical(type, "tsp")) {
    type <- c(-0.2, 0.2)
  }

  rnam <- rownames(X)
  cnam <- colnames(X)

  if (is.numeric(type)) {
    idx <- get_idx(type, ppm)
    if (length(idx) == 0) stop("No ppm values found in specified calibration range.")

    zeroPpm <- which.min(abs(ppm[idx] - mean(type)))
    Xc <- t(vapply(seq_len(nrow(X)), function(i) {
      iCorr <- zeroPpm - which.max(X[i, idx])
      if (iCorr < 0) {
        x <- c(X[i, -seq_len(abs(iCorr))], rep(0, abs(iCorr)))
      } else if (iCorr > 0) {
        x <- c(rep(0, abs(iCorr)), X[i, ])
      } else {
        x <- X[i, ]
      }
      x[seq_len(length(ppm))]
    }, FUN.VALUE = numeric(length(ppm))))
  } else if (type == "glucose") {
    Xc <- .calibrate_doub(X, ppm, "glu")
  } else if (type == "alanine") {
    Xc <- .calibrate_doub(X, ppm, "ala")
  } else {
    stop("Unsupported type argument. Must be one of 'tsp', 'glucose', 'alanine', or numeric vector of length 2.")
  }

  rownames(Xc) <- rnam
  colnames(Xc) <- cnam
  return(Xc)
}
