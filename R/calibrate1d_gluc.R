# reference 1D NMR to glucose @ 5.233 (d)

#' @title Calibrate 1D to glucose @ 5.233 (d)
#' @param X num matrix, NMR matrix with spectra in rows
#' @param ppm num vec, chemical shift matching to X
#' @return Calibrated X matrix
#' @author Torben Kimhofer \email{torben.kimhofer@@gmail.com}
#' @import signal sgolayfilt
#' @keywords internal
.calibrate1d_gluc <- function(X, ppm) {

  ppick=function(X, ppm, fil.n=5, fil.p=3){
    in_red=ncol(X)-1
    apply(X, 1, function(x, ic=in_red, pp=ppm){
      x[is.na(x)]=0
      # get noise level (3*sd)
      idx=get.idx(c(9.8, 10), pp)
      if(length(idx)>100){
        offset=mean(x[idx])
        noi=(sd(x[idx]-offset)*200)+offset
      }else{noi=median(x)}


      x_smooth=sgolayfilt(x, p = fil.p, n = fil.n)
      xd_sign=sign(diff(x_smooth))
      ids=which(xd_sign[-ic]!=xd_sign[-1] )+1
      idx=which(x[ids]>noi)
      out=data.frame(idc=ids[idx], ppm=pp[ids[idx]], Int=x[ids[idx]])
      rownames(out)=NULL
      out
    })
  }

  # pick peaks glucose reagion
  idx=get.idx(c(5.15, 5.3), ppm)
  test=ppick(X[,idx], ppm[idx])

  # spec(ppm[idx], X[61,idx])
  # points(test[[61]]$ppm, test[[61]]$Int)

  # calculate J const
  s=lapply(test, function(x){
    ii=combn(seq(nrow(x)), 2)
    ii_d=apply(ii, 2, function(id, pp=x$ppm){
      abs(diff(pp[id]))
    })

    ii_idx=which(ii_d>0.006 & ii_d<0.007)

    if(length(ii_idx)==1){

      out=x[ii[,ii_idx],]
      out$diff=diff(out$ppm)
      out$mm=mean(out$ppm)
      return(out)
    }

    if(length(ii_idx)>1){
      if(!is.matrix(ii[,ii_idx]))
      crit=apply(ii[,ii_idx], 2, function(id){
        c('max'=mean(x$Int[id]), 'ratio'=max(x$Int)/min(x$Int[2]))
      })

      #doubl=order(crit[1,], decreasing = T)
      doubl=which.max(crit[1,])
      out=x[ii[,ii_idx[doubl]],]
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








