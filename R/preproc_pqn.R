#' @title Probabilistic Quotient Normalisation (PQN)
#'
#' @description
#' Applies probabilistic quotient normalisation (PQN) to spectra. PQN estimates a
#' sample-specific dilution factor from the median of quotients relative to a
#' reference spectrum and scales each spectrum accordingly.
#'
#' @param X Numeric matrix or data.frame. Each row is a sample spectrum and each
#'   column is a variable (e.g. chemical shift point or bin).
#' @param ref_index Integer vector of row indices used to compute the reference
#'   spectrum. If \code{NULL}, all rows are used.
#' @param total_area Logical. If \code{TRUE}, total area normalisation is applied
#'   to the working copy before estimating the PQN dilution factors. See Notes.
#' @param bin Optional named list controlling binning for reference estimation,
#'   e.g. \code{list(ppm = ppm, width = 0.05)} or \code{list(ppm = ppm, npoints = 400)}.
#'   If \code{NULL}, no binning is applied.
#' @param iref Deprecated. Use \code{ref_index}.
#' @param TArea Deprecated. Use \code{total_area}.
#'
#' @details
#' \strong{Mechanics.} Let \eqn{x_i} be spectrum \eqn{i} and \eqn{r} a reference spectrum.
#' PQN computes quotients \eqn{q_{ij} = x_{ij} / r_j} and defines the dilution factor as
#' \eqn{d_i = 1 / \mathrm{median}_j(q_{ij})}. The PQN-normalised spectrum is
#' \eqn{x_i^{(PQN)} = d_i \, x_i}.
#'
#' The reference spectrum \eqn{r} is typically the median spectrum across all samples
#' or across QC samples (\code{ref_index}).
#'
#' If \code{bin} is provided, \eqn{r} and dilution factors are computed on binned spectra,
#' but applied to the original spectra.
#'
#' Dilution factors are stored in \code{attr(X, "m8_pqn")$dilution_factor}.
#'
#' @section Notes on total area normalisation:
#' Total area normalisation prior to PQN is \emph{usually not recommended}. Total area
#' scaling removes global intensity differences by enforcing equal total signal per sample.
#' PQN is itself a global scaling method intended to estimate dilution. Applying both can
#' substantially change results because PQN no longer estimates dilution alone, but also
#' compensates compositional distortions introduced by total area scaling.
#'
#' Situations where \code{total_area = TRUE} can be defensible include:
#' \itemize{
#'   \item when spectra have large, non-dilution-related amplitude differences caused by
#'         acquisition artefacts (receiver gain / baseline offset) and you explicitly want
#'         to stabilise the reference estimation step;
#'   \item when the measured total signal is expected to be constant by design (e.g. strictly
#'         controlled sample mass/volume and stable overall metabolite pool), and the main
#'         goal is to reduce technical scaling variation before PQN.
#' }
#' In most metabolomics settings, prefer PQN without total area scaling.
#'
#'#' @section On spectral alignment and binning:
#' PQN assumes that corresponding variables represent the same chemical signal
#' across spectra. If spectra are not well aligned, small peak shifts can inflate
#' the variability of pointwise quotients \eqn{x_{ij} / r_j}, leading to unstable
#' dilution factor estimates.
#'
#' In such cases, slight binning (e.g. narrow fixed-width bins) prior to reference
#' estimation is recommended. Binning reduces sensitivity to minor misalignments
#' by aggregating neighbouring variables. However, excessive binning may obscure
#' narrow signals and should be avoided.
#'
#' Alternatively, prior spectral alignment is preferable when available.
#'
#' @return Numeric matrix of PQN-normalised spectra.
#'
#' @references
#' Dieterle F, Ross A, Schlotterbeck G, Senn H (2006).
#' Probabilistic Quotient Normalization as Robust Method to Account for Dilution of
#' Complex Biological Mixtures. \emph{Analytical Chemistry}, 78(13), 4281–4290.
#'
#'
#' @examples
#' set.seed(1)
#' ppm <- seq(0, 10, length.out = 1000)
#' ref <- dnorm(ppm, 3, 0.15) + dnorm(ppm, 6, 0.20) + dnorm(ppm, 7.5, 0.18)
#' dil <- c(1, 0.8, 0.6, 0.4, 0.2)            # true dilution factors
#' X <- t(sapply(dil, function(d) d * ref + rnorm(length(ref), 0, 0.005)))
#' plot_spec(X, ppm)
#'
#' Xn <- pqn(X, ref_index=1)
#' dil_est <- attr(Xn, "m8_pqn")$dilution_factor
#'
#' cbind(true = dil, estimated = dil_est)
#'
#' @family preprocessing
#' @export
pqn <- function(X,
                ref_index = NULL,
                total_area = FALSE,
                bin = NULL,
                iref = NULL,
                TArea = NULL) {

  if (!is.null(iref)) {
    warning("`iref` is deprecated; use `ref_index` instead.", call. = FALSE)
    ref_index <- iref
  }
  if (!is.null(TArea)) {
    warning("`TArea` is deprecated; use `total_area` instead.", call. = FALSE)
    total_area <- TArea
  }

  inp <- X
  u   <- .m8_unpack_dat(X)

  X0  <- u$X
  X   <- .dimX(X0)
  ppm <- u$ppm
  meta <- u$meta
  nams <- rownames(X)

  if (!is.null(bin)) {

    if (!is.list(bin) || is.null(names(bin)))
      stop("Argument 'bin' must be a named list.", call. = FALSE)

    if (!("ppm" %in% names(bin))) {

      if (!is.null(ppm)) {
        bin$ppm <- ppm
      } else if (!is.null(colnames(X))) {
        ppm_infer <- as.numeric(colnames(X))
        if (anyNA(ppm_infer))
          stop("ppm could not be inferred from colnames(X).", call. = FALSE)
        bin$ppm <- ppm_infer

      } else {
        stop("Binning requested but no ppm supplied.", call. = FALSE)
      }
    }

    if (!.check_X_ppm(X, bin$ppm))
      stop("Mismatch between X and ppm in binning specification.", call. = FALSE)

    has_width <- "width" %in% names(bin) && !is.null(bin$width)
    has_npts  <- "npoints" %in% names(bin) && !is.null(bin$npoints)

    if (!has_width && !has_npts)
      stop("Argument 'bin' must include either 'width' or 'npoints'.", call. = FALSE)

    X1 <- binning(
      X = X,
      ppm = bin$ppm,
      width = if (has_width) bin$width else NULL,
      npoints = if (has_npts) bin$npoints else NULL
    )
    ppm_est <- attr(X1, "m8_axis")$ppm
  } else {
    X1 <- X
    ppm_est <- ppm
  }

  if (isTRUE(total_area)) {
    rs <- rowSums(X1, na.rm = TRUE)
    rs[rs == 0 | !is.finite(rs)] <- NA_real_
    X1 <- X1 / rs
  }

  if (is.null(ref_index)) {
    ref <- apply(X1, 2, median, na.rm = TRUE)
  } else {
    ref_index <- as.integer(ref_index)
    if (any(!is.finite(ref_index)))
      stop("ref_index contains non-finite values.", call. = FALSE)
    if (any(ref_index < 1L | ref_index > nrow(X1)))
      stop("ref_index contains out-of-range indices.", call. = FALSE)

    if (length(ref_index) == 1L) {
      ref <- X1[ref_index, ]
    } else {
      ref <- apply(X1[ref_index, , drop = FALSE], 2, median, na.rm = TRUE)
    }
  }

  ref[!is.finite(ref) | ref == 0] <- 1e-4

  ref_safe <- ref
  ref_safe[!is.finite(ref_safe) | ref_safe == 0] <- 1e-12

  Q <- sweep(X1, 2, ref_safe, "/")
  Q[!is.finite(Q)] <- NA_real_
  quot_med <- apply(Q, 1, stats::median, na.rm = TRUE)

  dil_F <- 1 / quot_med
  dil_F[!is.finite(dil_F)] <- 1

  X_pqn <- sweep(X, 1, dil_F, "*")
  rownames(X_pqn) <- nams
  colnames(X_pqn) <- colnames(X)

  X_pqn <- .m8_copy_attrs(X0, X_pqn)
  attr(X_pqn, "m8_axis") <- list(ppm = ppm)

  attr(X_pqn, "m8_pqn") <- list(
    dilution_factor = quot_med,
    ref_index = ref_index,
    total_area = isTRUE(total_area),
    bin = bin,
    ppm_est = ppm_est
  )

  X_pqn <- .m8_stamp(
    X_pqn,
    step = "pqn",
    params = list(
      ref_index = ref_index,
      total_area = isTRUE(total_area),
      bin = bin
    ),
    notes = "Probabilistic quotient normalisation applied."
  )

  if (is.list(inp) && !is.null(inp$X)) {
    inp$X <- X_pqn
    inp$ppm <- ppm
    inp$meta <- meta
    return(inp)
  }

  X_pqn
}
