#' @title NMR data binning
#' Equidistant binning of spectra based on either a specified bin witdth or the number of desired bins.
#' @export
#' @param X num matrix, NMR data with spectra in rows
#' @param ppm num array, chemical shift positions, length matches to columns in X
#' @param width num, bin size in ppm or NULL in case \code{npoints} is specified
#' @param npoints num, desired number of bins per spectrum or NULL in case  \code{width} is specified
#' @details Specify either \code{width} or \code{npoints} argument - if both are fiven, \code{npoints} is used. Input argument \code{ppm} can be omitted if chemical shift infomration is encoded in the column names of the NMR matrix \code{X}. The returne matrix encodes chemical shift information as column names.
#' @return Numeric matrix whith spectra in rows and chemical shift information in columns.
#' @import stats approxfun
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
# @examples
# load(covid)
# Xb=binning(X, ppm, width=0.005)
# ppm_bin=as.numeric(colnames(Xb))
# par(mfrow=c(2,1))
# spec(X[1,], ppm, shift=c(5.15, 5.3), interactive=FALSE)
# spec(Xb[1,], ppm_bin,  shift=c(5.15, 5.3), interactive=FALSE)
# @family NMR
binning <- function (X, ppm, width = 0.01, npoints=NULL)
{
    if (is.null(ppm) && (is.matrix(X) | is.data.frame(X)) &&
        !is.null(colnames(X))) {
        ppm <- as.numeric(colnames(X))
    }
    else {
        if (!.check_X_ppm(X, ppm))
            stop("Non-matching dimensions X matrix and ppm vector or missing values in ppm.")
    }
    if (is.vector(X)) {
        X <- t(X)
    }
    if (!is.null(width) & !is.null(npoints)) {
        stop("Please specify one or the other: bin width or desired number of data points per spectrum.\n")
    }
    if (is.null(width) & is.null(npoints)) {
        stop("Define bin width in ppm or desired number of bins.")
    }
    if (!is.null(width) & is.null(npoints)) {
        if (width <= abs(diff(ppm[seq_len(2)]))) {
            stop("Bin width equals or is smaller than the difference of neighbouring ppm points.")
        }

        res=(ppm[1] - ppm[2])
        new_res=width/round(width/res)
        step=round(width/res)
        ppm_new=seq(max(ppm), min(ppm), by=-new_res)
        iid=floor(length(ppm_new)/step)
        ybin=rep(seq(iid), each=step)

        Xb <- t(apply(X, 1, function(x, ppmt = ppm_new, ppm_fres = ppm, yb=ybin) {
            sInter <- approxfun(ppm_fres, x)
            s=sInter(ppmt)

            #browser()
            out=sapply(seq(max(yb)), function(i){
                iidx=which(yb==i)
                sum(s[iidx])
            })

            return(out)

        }))

        ppm_bin=sapply(seq(max(ybin)), function(i){
            iidx=which(ybin==i)
            mean(ppm_new[iidx])
        })

        colnames(Xb) <- ppm_bin
        rownames(Xb) <- rownames(X)
        return(Xb)

    }
    if (!is.null(npoints) & is.null(width)) {
        if (npoints >= length(ppm)) {
            stop("Input variable npoints cannot be larger or equal than length of ppm vector.")
        }
        ppm_bin <- seq(max(ppm), min(ppm), length.out = npoints)

        iid=floor(length(ppm)/npoints)
        ybin=rep(seq(npoints), each=iid)


        Xb <- t(apply(X, 1, function(s, yb=ybin) {
            #browser()
            out=sapply(seq(max(yb)), function(i){
                iidx=which(yb==i)
                sum(s[iidx])
            })

            return(out)

        }))

        ppm_bin=sapply(seq(max(ybin)), function(i){
            iidx=which(ybin==i)
            mean(ppm_new[iidx])
        })

        colnames(Xb) <- ppm_bin
        rownames(Xb) <- rownames(X)
        return(Xb)

    }

    # estimate line at lower resolution than current (to not loose peaks) but at multiple of width, so that sum can be accurately estimated

    return(NULL)


}
