#' Probabilistic quotient normalisation
#' @export
#' @aliases pqn
#' @param X num matrix or data.frame, metabolite data with rows representing spectra and column spectral variables
#' @param iref num or int vector of X row indices, used to calculate the reference spectrum (see Details). Set to NULL if all samples should be use to calculate reference index
#' @param TArea logical indicating if total area normalisation should be applied first (see Details).
#' @param add_DilF char string, defining name of variable exported to global namespace and containing dilution factor information  Can be set to NULL if dilution factor should not be exported.
#' @param bin named list: ppm and binning function parameters width or npoints
#' @details It is sometimes favourable not to use all spectra to calculate a dilution reference (e.g. QC samples should generally be excluded). Therefore, a vector of indices can be specified with the parameter \code{reference.idx} and the respective spectra are used to calculate the median spectrum as a dilution reference. The parameter \code{reference.idx} can also be a single index, then the respective spectrum is used as a reference. If it is set \code{N/A}, all spectra in \code{X} are used to calculate the dilution reference spectrum. to Total area normalisation is integral part of the probabilistic quotient normalisation algorithm (see References), however, this sometimes distorts the spectra, i.e. is not suitable for normalisation. The parameter \code{TArea} can be set to TRUE or FALSE, depending if total area normalisation should be applied or not.
#' @references Dieterle, F., \emph{et al.} (2006), Probabilistic Quotient Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures. Application in 1H NMR Metabonomics, \emph{Analytical Chemistry}, 78.3, 4281-90.
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}


pqn <- function(X, iref = NULL, TArea = F, add_DilF = NULL, bin=list(ppm=NULL, width=0.05, npoints=NULL)) {

  nams <- rownames(X)

  # bin specs
  if(!is.null(bin)){

    if( !any(names(bin) %in% c('width', 'npoints'))) stop('Define width or npoints binning argument.')
    if(is.null(bin$ppm)){ppm=as.numeric(colnames(X)); } else{
      if(!.check_X_ppm(X, ppm)) stop('Non-matching dimensions X matrix and ppm vector or missing values in ppm.')
    }

    X <- binning(X, ppm, unlist(bin))
  }
  # apply total area normalisation
  if (TArea == T) {
    X <- t(apply(X, 1, function(x) x/sum(x, na.rm = T) * 100))
  }

  le.ref <- length(iref)

  if(is.null(iref)){

    ref <- apply(X, 2, median, na.rm=TRUE)

  }else{

    iref=as.integer(iref)

    if( any(is.na(iref) | is.infinite(iref)) ) stop('Index vector iref contains NAs.')

    if (le.ref == 1 ) {
      ref <- X[iref, ]
    }else{
      ref <- apply(X[iref, ], 2, median, na.rm=TRUE)
    }
  }

  ref[ref == 0] <- 1e-04

  dil_F <- 1/apply(X, 1, function(x) median(x/ref, na.rm = T))
  X_pqn <- t(apply(rbind(1:nrow(X)), 2, function(i) {
    X[i, ] * dil_F[i]
  }))
  if (!is.null(add_DilF)) {
    assign(add_DilF, dil_F, envir = .GlobalEnv)
  }
  rownames(X_pqn) <- nams
  return(X_pqn)
}
