#' @title Chemical shift calibration
#' @export
#' @aliases calibrate
#' @description Chemical shift calibration with a reference signal.
#' @param X num matrix, NMR data matrix with rows representing spectra and columns representing chemical shift position
#' @param ppm num vector, matched to columns of X
#' @param type str or num. Str: Either 'tsp' or 'glucose' for urine or blood-derived spectra, respectively (see Details). Num: ppm range of max height signal that will be used to reference to zero
#' @details Spectral calibration to a selected chemical shift reference signal.
#' \itemize{
#' \item{\code{type='tsp'}{calibration to 0 ppm using the highest peak located in interval 0 +/- 0.20 ppm (Trimethylsilylpropanoic acid resonance)}}
#' \item{\code{type='glucose'}{calibration to 5.23 ppm using the glucose doublet located in interval 5.15 ppm to 5.30 ppm}}
#' \item{\code{type='alanine'}{calibration to 1.48 ppm using the alanine doublet located in interval 1.4 ppm to 1.53 ppm}}
#' \item{\code{type=ppm_range}{ppm_range is numeric array, calibration to mean(ppm_range) using a maximum signal located in interval ppm_range}}
#' }
#' @return num matrix X, calibrated NMR data matrix.
#' @references Dona, A.C., \emph{et al.} (2014) Precision high-throughput proton NMR spectroscopy of human urine, serum, and plasma for large-scale metabolic phenotyping. \emph{Analytical Chemistry}. 86.19. 9887-94.
#' @examples
#' data(covid_raw)
#' matspec(X, ppm, shift=c(-0.1,0.1))
#' X_tsp=calibrate(X, ppm, type='tsp')
#' matspec(X, ppm, shift=c(-0.1,0.1))
#'
#' # calibrate with glucose (blood plasma)
#' matspec(X, ppm, shift=c(5.15, 5.3))
#' X_tsp=calibrate(X, ppm, type='glucose')
#' matspec(X, ppm, shift=c(5.15, 5.3))
#'
#' @author \email{torben.kimhofer@@murdoch.edu.au}


calibrate <- function(X, ppm, type = c('tsp', 'glucose', 'alanine') ){

    if(all(is.character(type)) && (any(is.na(type[1])) || length(type)>1)){stop('Check function argument `type`!')}
    if(all(is.numeric(type)) && length(type)!=2){stop('Check function argument `type`!')}

    if (type[1] == "tsp") type=c(-0.2, 0.2)

    rnam <- rownames(X)
    cnam <- colnames(X)

    if (is.numeric(type)[1]) {
        idx <- get_idx(type, ppm)
        zeroPpm <- which.min(abs(ppm[idx]-mean(type)))
        Xc<-t(vapply(seq_len(nrow(X)), function(i){
            iCorr <- zeroPpm - which.max(X[i, idx])
            if (iCorr<0) { x <- c(X[i, -c(seq_len(abs(iCorr)))], rep(0, abs(iCorr))) }
            if (iCorr>0) { x <- c(rep(0, abs(iCorr)), X[i,]) }
            if (iCorr==0) { x <- X[i,] }
            x[seq_len(length(ppm))]

            }, FUN.VALUE=ppm))

    }

    if (type[1] == "glucose") {
        Xc <- .calibrate_doub(X, ppm, 'glu')
    }
    if (type[1] == "alanine") {
        Xc <- .calibrate_doub(X, ppm, 'ala')
    }

    rownames(Xc) <- rnam
    colnames(Xc) <- cnam

    return(Xc)
}
