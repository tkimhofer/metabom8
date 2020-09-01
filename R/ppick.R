#' @title Find local extrema in NMR spectra (peak picking)
#' @param X num matrix, NMR matrix with spectra in rows
#' @param ppm num vec, chemical shift matching to X
#' @param fil_p, SG filter order
#' @param fil_n num, SG filter length (must be odd)
#' @param type str, return local extrema: maxima (\code{max}), minima (\code{min}), or \code{both}
#' @details Function first smoothes the spectrum to reduce extrema related to noise fluctuations using a Savitzki-Golay (SG) filter. Local extrema are identified according to sign change of adjacent points.
#' @return List of data frames for each spectrum. Each data frame has four columns: index in X col/ppm, ppm value, signal intensity and local extrema type (+1 for local min, -1 for local max). Number of rows match to all local extrema found within a spectrum.
#' @author \email{torben.kimhofer@@gmail.com}
#' @import signal sgolayfilt
#' @seealso \code{\link[signal]{sgolayfilt}}
#' @family NMR
#' @examples
#' load(covid)
#' peaks=ppick(X, ppm)
#' # Plotting first spectrum in reagion 3.0 to 3.15 ppm
#' spec(X[1,], ppm, shift=c(4.2, 5.3))
#' points(speaks[[1]]$ppm, speaks[[1]]$Int, col=factor(speaks[[1]]$type))
ppick=function(X, ppm, fil_p=3, fil_n=5, type='max'){
  X<-.dimX(X)
  if(!.check_X_ppm(X, ppm)){stop('Dimensions of X and ppm don\'t match or contain non-numeric values.')}
  in_red=ncol(X)-1
  apply(X, 1, function(x, ic=in_red, pp=ppm){
    x[is.na(x)]=0
    # get noise level (3*sd)
    # idx=get.idx(c(9.8, 10), pp)
    # if(length(idx)>100){
    #   offset=mean(x[idx])
    #   noi=(sd(x[idx]-offset)*200)+offset
    # }else{noi=median(x)}
    x_smooth=sgolayfilt(x, p = fil_p, n = fil_n)
    xd_sign=sign(diff(x_smooth))
    ids=which(xd_sign[-ic]!=xd_sign[-1] )+1
    idx=seq(ids)
    out=data.frame(idc=ids[idx], ppm=pp[ids[idx]], Int=x[ids[idx]], Etype=xd_sign[ids[idx]])
    rownames(out)=NULL

    switch(type,
           'max'={out[out$Etype<0,]},
           'min'={out[out$Etype>0,]},
           'both'={out}
           )

  })

}
