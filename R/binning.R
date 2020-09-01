#' @title NMR data binning
#' Equidistant binning of spectra based on either a specified bin witdth or the number of desired bins.
#' @export
#' @param X num matrix, NMR data with spectra in rows
#' @param ppm num array, chemical shift positions, length matches to columns in X
#' @param width num, bin size in ppm or NULL in case \code{npoints} is specified
#' @param npoints num, desired number of bins per spectrum or NULL in case  \code{width} is specified
#' @details Specify either \code{width} or \code{npoints} argument - if both are fiven, \code{npoints} is used. Input argument \code{ppm} can be omitted if chemical shift infomration is encoded in the column names of the NMR matrix \code{X}.
#' @return Numeric matrix whith spectra in rows and chemical shift bins in columns.
#' @import stats approxfun
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @family NMR

binning <- function(X, ppm, width = 0.001, npoints = NULL) {


  if( is.null(ppm) && ( is.matrix(X) | is.data.frame(X) ) && !is.null(colnames(X)) ){ ppm=as.numeric(colnames(X)); } else{
    if(!.check_X_ppm(X, ppm)) stop('Non-matching dimensions X matrix and ppm vector or missing values in ppm.')
  }


  if(is.vector(X)) {X=t(X)}

  if (ppm[1] < ppm[length(ppm)]) {
    stop("Let's stick to the convention: Revert the order of ppm and NMR matrix, so that the chemical shift decreases with increasing indices!")
  }

  if (is.null(width) & is.null(npoints)) {
    stop("Pleas specify spectral bin width or the desired number of data points per spectrum.\n")
  }
  if (is.null(width) & is.null(npoints)) {
    stop('Define bin widht or desired number of bins.')
  }
  if (!is.null(width)) {
    if (width <= abs(diff(ppm[seq_len(2)]))) {
      stop("Bin width equals or is smaller than the difference of neighbouring ppm points.")
    }
    ppm_bin <- seq(max(ppm), min(ppm), by = - width)
  }
  if (!is.null(npoints)) {
    if (npoints >= length(ppm)) {
      stop("Input variable npoints cannot be larger or equal than length of ppm vector.")
    }
    ppm_bin <- seq(max(ppm), min(ppm), length.out = npoints)
  }

  quot <- ceiling(width/(ppm[1] - ppm[2]))
  rr <- ceiling(length(ppm)/quot)
  x <- rep(seq_len(200), each = rr)
  x <- x[seq_len(length(ppm))]
  Xb <- t(apply(X, 1, function(x, ppmt = ppm_bin, ppm_fres = ppm) {
    sInter <- approxfun(ppm_fres, x)
    sInter(ppmt)
  }))
  colnames(Xb) <- ppm_bin
  rownames(Xb) <- rownames(X)
  return(Xb)
}
