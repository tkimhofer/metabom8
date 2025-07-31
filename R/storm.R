#' @title Subset Optimisation by Reference Matching (STORM)
#' @param X  Numeric matrix or dataframe where each row represents a NMR spectrum and each column a chemical shift variable.
#' @param ppm num array of chemical shift variables, matched to columns in X
#' @param b int, boundary value provided as positive integer in chemical shift space
#' @param q num, level of significance
#' @param idx.refSpec int, reference spectrum defined as row index for X
#' @param shift num array, min and max chemical shift value defining narrow boundary of the target signal
#' @details The STORM algorithm can be used to subselect NMR spectra based on a spectrocopic target signal to match multiplicity and chemical shift position. The output are row-indices for X, defining the optimal spectral subset that is typically used as input for STOCSY.
#' @references Posma, Joram M., et al. "Subset optimization by reference matching (STORM): an optimized statistical approach for recovery of metabolic biomarker structural information from 1H NMR spectra of biofluids." Analytical chemistry 84.24 (2012): 10694-10701.
#' @return Vector of X row indices that define the optimal subset.
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @family NMR
#' @seealso stocsy
#' @importFrom stats cor cov pt
#' @importFrom scales breaks_pretty
#' @importFrom methods hasArg is
#' @export
#' @examples
#'
#' set.seed(123)
#' n <- 100; S <- 1000
#' ppm <- seq(10, 0, length.out = S)
#' gauss <- function(x, c, w, h) h * exp(-((x - c)^2) / (2 * w^2))
#' sig1 <- gauss(ppm, 7, 0.05, 10)
#' sig2 <- gauss(ppm, 3.5, 0.1, 8)
#' sig3 <- gauss(ppm, 1, 0.07, 5)
#' spectra <- matrix(0, n, S)
#' for(i in 1:n) {
#'   spectra[i, ] <- sig1 + sig2 + rnorm(S, 0, 0.1)
#'   if(i <= 25) spectra[i, ] <- spectra[i, ] + sig3
#' }
#'
#' tt=storm(X=spectra, ppm=ppm, b=30, q=0.05, idx.refSpec=1, shift=c(0.75,1.25))
storm <- function(X, ppm, b=30, q=0.05, idx.refSpec, shift){

  # center X, define reference index and reference spectrum
  X_sc <- .scaleMatRcpp(X, seq.int(0, nrow(X) - 1), center=TRUE, scale=0)
  X <- X_sc[[1]]
  ref.idx <- get_idx(range=shift, ppm)
  ref <- X[idx.refSpec, ref.idx]

  # init
  subset <- 0
  subset1 <- seq_len(nrow(X))
  i <- 1

  # perform storm
  while(length(which(!subset %in% subset1))>0){

    # correlation based subset selection
    subset <- subset1
    Xr <- X[subset, ref.idx]
    r <- cor(t(Xr), ref)
    a <- -abs(r * sqrt((length(r)-2)/(1-r^2)))
    pval <- 2*pt(a, (length(r)-2))

    subset1 <- subset[r>0]
    pval <- pval[r>0]
    index <- order(pval)
    subset1 <- subset1[index]

    # Stocsy with driver=max intensity reference spectrum
    index <- which.max(ref)
    r <- cor(X[subset1, (ref.idx[index]-(b+1)):(ref.idx[index]+(b+1))], X[subset1,ref.idx[index]])
    co <- cov(X[subset1, (ref.idx[index]-(b+1)):(ref.idx[index]+(b+1))], X[subset1,ref.idx[index]])

    # update reference spectrum and reference index
    a <- -abs(r * sqrt((length(r)-2)/(1-r^2)))
    pval <- 2*pt(a,(length(r)-2))

    ref <- co[pval<q & r>0]
    ref.idx <- (ref.idx[index]-(b+1)):(ref.idx[index]+(b+1))
    ref.idx <- ref.idx[pval<q & r>0]
  }
  return(subset1)
}
