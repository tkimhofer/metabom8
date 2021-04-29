#' @title Normalise by the ERETIC signal
#' @export
#' @param X num matrix, NMR data, spectra in row with ppm values as column names
#' @param integr logic, FALSE: fct returns integral of ERETIC, TRUE: fct returns spectra normalised with ERETIC integral
# @param ppm ppm, num array - chemical shift positions, length matches to columns in X
#' @return
#' \itemize{
#'   \item num matrix: X-normalised with ERETIC integral (urine - 12 ppm, plasma = 15 ppm)
#'   \item num vector: ERETIC integrals (urine - 12 ppm, plasma = 15 ppm)
#' }
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @family NMR
#' @section

normErectic <- function(X, integr=FALSE){
  ppm=as.numeric(colnames(X))
  idx=get.idx(c(10,16), ppm)
  qr=ppm[which.max(X[1,idx])]
  iid=which.min(qr-c(12, 16))
  if(iid==1){iid=get.idx(c(11.9, 12.1), ppm)}else{
    iid=get.idx(c(14.9, 15.1), ppm)
  }
  if(integr){
    qr=apply(X[,idx], 1, sum)
    return(qr)
  }else{
    Xn=t(apply(X, 1, function(x, ii=iid){
      x/sum(x[ii])
    }))
    return(Xn)
  }
}
