#' Align an NMR spectra to a reference spectrum (cross-correlation)
#' @export
#' @param seg matrix or data frame, NMR submatrix containing segments (spectra in rows)
#' @param idx_ref int, rown index of reference spectrum
#' @param clim num, upper cross-correlation threshold (see Details)
#' @param med logical, use median spectrum as reference?
#' @param norm logical, minmax normalisation for each spectrum prior to alignment? Useful for viz.
#' @return Matrix of aligned spectral segments
#' @details Spectral alignment is obtained by x-positional shifting a spectral segment. The extend of the shift is determined by maximising the cross-correlation with a refrence segment. The refernce can either be a spectrum in *seg* or the segment median spectrum (parameter *idx_ref* is ignored then). Non-matching spectral segments, defined by *clim* (cross-corrleation upper limit, similarity score) are not aligned.
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom stats ccf
#' @family NMR
#' @section

alignSegment=function(seg, idx_ref=1, clim=0.7, med=TRUE, norm=FALSE){

  if(is.null(nrow(seg)) || nrow(seg)<2) {stop('Check input dimensions.')}
  if(med[1]){seg=rbind(apply(seg, 2, median), seg); idx_ref=1}
  if(norm){seg=t(apply(seg, 1, minmax))}

  out=vapply(seq(nrow(seg)), function(i, le=ncol(seg)){

    cc_mod=ccf(seg[i,], seg[idx_ref,], type='correlation', plot = FALSE)
    acf=as.vector(cc_mod$acf)
    idx<-which.max(acf)
    cval<-acf[idx]
    xlag<-as.vector(cc_mod$lag)[idx]

    if(cval > clim ){
      if(xlag>0){
        dd=c(seg[i, -seq_len(abs(xlag))], rep(0, abs(xlag)))
      }
      if(xlag<0){
        dd=c(rep(0, abs(xlag)), seg[i,-seq((le-abs(xlag)+1),le)])
      }

      if(xlag==0){dd=seg[i,]}

      return(dd)
    } else{
      # print(i)
      # print(cval)
      seg[i,]
    }
  }, FUN.VALUE = seg[1,])

  out=t(out)
  if(med){out=out[-1,]}

  return(out)
}


