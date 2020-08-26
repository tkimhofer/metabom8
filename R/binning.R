#' Equidistant binning of spectra
#' @export
#' @param X NMR matrix (spectra in rows)
#' @param ppm num ppm vector
#' @param width num, bin size
#' @param npoints num, desired number of bins per spectrum
#' @return Use either width or npoints argument. Function returns is a matrix (rows = spectra, columns=bins).
#' @import stats approxfun
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}

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
