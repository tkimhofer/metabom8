#' @title Probabilistic Quotient Normalisation (PQN)
#'
#' @description
#' Applies probabilistic quotient normalisation (PQN) to NMR spectra, optionally preceded by total area normalisation and binning.
#'
#' @param X A numeric matrix or data.frame. Each row is a sample spectrum, and each column is a chemical shift variable.
#' @param iref Integer vector of row indices in \code{X} used to compute the reference spectrum. If \code{NULL}, all rows are used.
#' @param TArea Logical. If \code{TRUE}, total area normalisation is applied before PQN.
#' @param add_DilF Character. Name of global variable to store dilution factors. Set to \code{NULL} to disable.
#' @param bin Optional named list for binning (e.g., \code{list(ppm = ppm, width = 0.05)}). If \code{NULL}, binning is not applied.
#'
#' @details
#' PQN accounts for sample dilution effects by normalising each spectrum using the median quotient relative to a reference spectrum. This reference is typically the median spectrum of QC samples or the entire dataset.
#'
#' Total area normalisation (\code{TArea = TRUE}) may help standardise sample intensities, but can distort profiles in some cases.
#'
#' If binning is applied, the function uses a reference spectrum computed on binned data.
#'
#' @return A matrix with the same dimensions as \code{X}, containing PQN-normalised spectra.
#'
#' @references
#' Dieterle F, Ross A, Schlotterbeck G, Senn H (2006). Probabilistic Quotient Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures. \emph{Analytical Chemistry}, 78(13), 4281–4290.
#'
#' @examples
#' data(covid)
#' Xn <- pqn(X, add_DilF = "dilutionFactor")
#' plot(dilutionFactor)
#'
#' @family NMR preprocessing
#' @export
pqn <- function(X, iref = NULL, TArea = FALSE, add_DilF = "dquot", bin = NULL) {

  if (is.data.frame(X)) X <- as.matrix(X)
  nams <- rownames(X)

  if (!is.null(bin)) {
    if (!all(c("ppm", "width") %in% names(bin)) && !all(c("ppm", "npoints") %in% names(bin))) {
      stop("Argument 'bin' must be a named list with 'ppm' and either 'width' or 'npoints'.")
    }

    if (is.null(bin$ppm)) {
      bin$ppm <- as.numeric(colnames(X))
    }

    if (!.check_X_ppm(X, bin$ppm)) {
      stop("Mismatch between X and ppm in binning specification.")
    }

    X1 <- binning(X, bin$ppm, unlist(bin))
  } else {
    X1 <- X
  }

  # Total area normalisation
  if (TArea) {
    X1 <- t(apply(X1, 1, function(x) x / sum(x, na.rm = TRUE) * 100))
  }

  # Determine reference spectrum
  ref <- if (is.null(iref)) {
    apply(X1, 2, median, na.rm = TRUE)
  } else {
    iref <- as.integer(iref)
    if (any(is.na(iref) | is.infinite(iref))) stop("iref contains NA or infinite values.")
    if (length(iref) == 1) X1[iref, ] else apply(X1[iref, ], 2, median, na.rm = TRUE)
  }

  ref[ref == 0] <- 1e-4

  # Compute dilution factors
  dil_F <- 1 / apply(X1, 1, function(x) median(x / ref, na.rm = TRUE))

  # Apply PQN
  X_pqn <- sweep(X, 1, dil_F, "*")

  if (!is.null(add_DilF)) {
    assign(add_DilF, dil_F, envir = .GlobalEnv)
  }

  rownames(X_pqn) <- nams
  return(X_pqn)
}
