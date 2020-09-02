# reference 1D NMR to glucose @ 5.233 (d)
#' @title Calibrate 1D to glucose @ 5.233 (d)
#' @param X num matrix, NMR matrix with spectra in rows
#' @param ppm num vec, chemical shift matching to X
#' @param sg_length int, nb of data points to cacl 2nd degree estimation with savitzki-golay (higher number -> smoother)
#' @return Calibrated X matrix
#' @author Torben Kimhofer \email{torben.kimhofer@@gmail.com}
#' @import signal sgolayfilt
#' @keywords internal
.calibrate1d_gluc <- function(X, ppm, sg_length=13) {

  # pick peaks glucose reagion
  idx=get.idx(c(5.15, 5.3), ppm)

  # remove broad signals and smooth
  test <- apply(X[,idx], 1, function(x, pp=ppm[idx]) {
    xs=sgolayfilt(x - asysm(x, lambda=100), p = 3, n=sg_length)
   ppick(xs, pp, type='max')[[1]]
  })

  #test=ppick(X[,idx], ppm[idx], type='max')
  #browser()
  # remove peaks below noise
  #noi=noise.est(X, ppm)


  #noi=0

  # spec(ppm[idx], X[61,idx])
  # points(test[[61]]$ppm, test[[61]]$Int)

  # calculate J const
  s=lapply(seq(length(test)), function(i){

    ptab=test[[i]][test[[i]]$Etype>0,]
    ii=combn(seq(nrow(test[[i]])), 2)
    ii_d=apply(ii, 2, function(id, pp=test[[i]]$ppm){
      abs(diff(pp[id]))
    })

    ii_idx=which(ii_d>0.006 & ii_d<0.007) # glucose doublet has J-cons of 3.85 Hz

    if(length(ii_idx)==1){

      out=test[[i]][ii[,ii_idx],]
      out$diff=diff(out$ppm)
      out$mm=mean(out$ppm)
      return(out)
    }

    if(length(ii_idx)>1){
      if(is.matrix(ii[,ii_idx]))
      crit=apply(ii[,ii_idx], 2, function(id){
        c('max'=mean(test[[i]]$Int[id]), 'ratio'=max(test[[i]]$Int)/min(test[[i]]$Int[2]))
      })

      #doubl=order(crit[1,], decreasing = T)
      doubl=which.max(crit[1,])
      out=test[[i]][ii[,ii_idx[doubl]],]
      # abline(v=out$ppm)
      out$diff=diff(out$ppm)
      out$mm=mean(out$ppm)
      return(out)
    }else{NULL}
  })

  Xc=vapply(seq(nrow(X)), function(i){
    if(is.null(s[[i]]) && nrow(s[[i]])>2){message(paste('Could not calibrate spectrum', i))}else{
      ff=approxfun(x=ppm-(max(s[[i]]$ppm)-5.233),  y=X[i,])
      news=ff(ppm)
      return(news)
    }
  }, FUN.VALUE =X[1,])
  Xc=t(Xc)
  return(Xc)
}








