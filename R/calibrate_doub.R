# Calibrate 1D NMR to doublet signal
#' @title Calibrate spectra to a doublet signal
#' @param X num matrix, NMR matrix with spectra in rows
#' @param ppm num vec, chemical shift matching to X
#' @param type string/num vec, either glu or ala for glucose and alanine calibration, respectively, or num array defining chem shift of doublet to be used for calibration
#' @param j_const num array, range of J constant values used to identify doublet signal (ppm), default is for glucose (doublet has J-cons of 3.85 Hz)
#' @param sg_length int, nb of data points to cacl 2nd degree estimation with savitzki-golay (higher number -> smoother)
#' @return Calibrated X matrix
#' @author Torben Kimhofer \email{torben.kimhofer@@gmail.com}
#' @import signal sgolayfilt
#' @import utils combn
#' @keywords internal
.calibrate_doub <- function(X, ppm, type=c('glu', 'ala'), j_const=c(0.006, 0.007), sg_length = 13) {

    if(type[1] =='glu') {
        idx <- get_idx(c(5.15, 5.3), ppm);
        cent_loc=5.233;
        j_const=c(0.006, 0.007)
    }
    if(type[1] =='ala') {
        idx <- get_idx(c(1.4, 1.56), ppm);
        cent_loc=1.48
        j_const = c(0.0115, 0.0135)
    }
    if(is.numeric(type[1])) {
        idx <- get_idx(type, ppm)
        cent_loc = mean(type)
    }

    # remove broad signals and smooth
    test <- apply(X[, idx], 1, function(x, pp = ppm[idx]) {
        xs <- sgolayfilt(x - asysm(x, lambda = 100), p = 3, n = sg_length)
        ppick(xs, pp, type = "max")[[1]]
    })

    # calculate J const
    s <- lapply(seq(length(test)), function(i) {

        ptab <- test[[i]][test[[i]]$Etype > 0, ]
        ii <- combn(seq(nrow(test[[i]])), 2)
        ii_d <- apply(ii, 2, function(id, pp = test[[i]]$ppm) {
            abs(diff(pp[id]))
        })

        ii_idx <- which(ii_d > j_const[1] & ii_d < j_const[2])

        if (length(ii_idx) == 1) {

            out <- test[[i]][ii[, ii_idx], ]
            out$diff <- diff(out$ppm)
            out$mm <- mean(out$ppm)
            return(out)
        }

        if (length(ii_idx) > 1) {
            if (is.matrix(ii[, ii_idx]))
                crit <- apply(ii[, ii_idx], 2, function(id) {
                  c(max = mean(test[[i]]$Int[id]), ratio = max(test[[i]]$Int)/min(test[[i]]$Int[2]))
                })

            # doubl=order(crit[1,], decreasing = T)
            doubl <- which.max(crit[1, ])
            out <- test[[i]][ii[, ii_idx[doubl]], ]
            # abline(v=out$ppm)
            out$diff <- diff(out$ppm)
            out$mm <- mean(out$ppm)
            return(out)
        } else {
            NULL
        }
    })

    Xc <- vapply(seq(nrow(X)), function(i) {
        if (is.null(s[[i]]) && nrow(s[[i]]) > 2) {
            message(paste("Could not calibrate spectrum", i))
        } else {
            ff <- approxfun(x = ppm - (max(s[[i]]$ppm) - cent_loc), y = X[i, ])
            news <- ff(ppm)
            return(news)
        }
    }, FUN.VALUE = X[1, ])

    Xc <- t(Xc)

    return(Xc)
}








