#' @title Excise Chemical Shift Regions from 1D NMR Spectra
#'
#' @description
#' Removes specified chemical shift regions from 1D \eqn{^1}H NMR spectra.
#' By default, commonly excluded metabolomics regions are removed
#' (upfield/downfield noise, water, urea).
#'
#' @param x Numeric matrix or vector. Spectra in rows and variables in columns.
#' @param ppm Numeric vector of chemical shift positions (ppm).
#'   If omitted, \code{ppm} is inferred from \code{colnames(X)}.
#' @param regions Named list of numeric vectors (length 2), specifying ppm
#'   regions to remove. Each element must define lower and upper bounds.
#'
#' @details
#' Default regions removed (ppm):
#' \itemize{
#'   \item Upfield noise: \code{[min(ppm), 0.25]}
#'   \item Residual water: \code{[4.5, 5.2]}
#'   \item Urea region: \code{[5.5, 6.0]}
#'   \item Downfield noise: \code{[9.7, max(ppm)]}
#' }
#'
#' Removed regions are recorded in \code{attr(X, "m8_prep")} if present.
#' Updated ppm values are stored in column names and in
#' \code{attr(X, "m8_axis")$ppm}.
#'
#' @return Numeric matrix with specified chemical shift regions removed.
#' @examples
#' set.seed(1)
#' ppm <- seq(0, 10, length.out = 1000)
#' X <- matrix(rnorm(100 * length(ppm)), nrow = 100)
#'
#' Xe <- excise(X, ppm)
#'
#' dim(Xe)
#' names(attributes(Xe))
#' @export
excise <- function(x, ppm = NULL, regions = NULL) {

  inp <- x
  u   <- .m8_unpack_dat(x, ppm = ppm)
  X0  <- u$X
  X   <- .dimX(X0)
  ppm <- u$ppm

  if (!.check_X_ppm(X, ppm))
    stop("Non-matching dimensions X matrix and ppm vector or invalid ppm values.", call. = FALSE)

  ppm <- as.numeric(ppm)

  if (is.null(regions)) {
    regions <- list(
      upfield_noise  = c(min(ppm), 0.25),
      water          = c(4.5, 5.2),
      urea           = c(5.5, 6.0),
      downfield_noise = c(9.7, max(ppm))
    )
  }

  if (!is.list(regions))
    stop("`regions` must be a list of numeric vectors (length 2).", call. = FALSE)

  for (nm in names(regions)) {
    r <- regions[[nm]]
    if (!is.numeric(r) || length(r) != 2L)
      stop("Each element of `regions` must be a numeric vector of length 2.", call. = FALSE)
  }

  idx_rm <- unique(unlist(
    lapply(regions, function(r) get_idx(range(r), ppm))
  ))

  keep_idx <- setdiff(seq_along(ppm), idx_rm)

  Xc  <- X[, keep_idx, drop = FALSE]
  ppc <- ppm[keep_idx]

  colnames(Xc) <- formatC(ppc, format = "fg", digits = 10)
  rownames(Xc) <- rownames(X)

  Xc <- .m8_copy_attrs(X0, Xc)
  attr(Xc, "m8_axis") <- list(ppm = ppc)

  Xc <- .m8_stamp(
    Xc,
    step = "excise1d",
    params = list(regions = regions),
    notes = "Specified chemical shift regions removed using get_idx()."
  )

  if (is.list(inp) && !is.null(inp$X)) {
    inp$X   <- Xc
    inp$ppm <- ppc
    return(inp)
  }

  Xc
}
