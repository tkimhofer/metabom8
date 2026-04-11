#' @title Select Indices for a Chemical Shift Region
#'
#' @description
#' Returns the indices of the chemical shift vector (\code{ppm}) that fall
#' within the specified range.
#'
#' @param range Numeric vector of length 2 specifying lower and upper bounds
#' (order does not matter).
#' @param ppm Numeric vector. The full chemical shift axis (in ppm).
#'
#' @return Integer vector of indices corresponding to \code{ppm} values within the given range.
#'
#' @export
#'
#' @examples
#' data(covid_raw)
#' X <- covid_raw$X
#' ppm <- covid_raw$ppm
#' idx_tsp <- get_idx(c(-0.1, 0.1), ppm)
#' ppm[range(idx_tsp)]
#' plot(ppm[idx_tsp], X[1, idx_tsp], type = 'l')
#'
#' @family NMR
get_idx <- function(range = c(1, 5), ppm) {
  if (anyNA(ppm) || any(!is.finite(ppm)))
    stop("`ppm` contains NA or non-finite values.", call. = FALSE)

  if (length(range) != 2L)
    stop("`range` must be a numeric vector of length 2.", call. = FALSE)

  if (!is.numeric(ppm))
    stop("`ppm` must be numeric.", call. = FALSE)

  lo <- min(range)
  hi <- max(range)

  which(ppm >= lo & ppm <= hi)
}

#' @rdname get_idx
#' @export
get.idx <- function(range = c(1, 5), ppm) {
  warning("`get.idx` is deprecated and will be removed in future versions. Use `get_idx` instead.", call. = FALSE)
  get_idx(range, ppm)
}
