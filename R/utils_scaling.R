#' @title Min-Max Scaling to Arbitrary Range
#' @description
#' Rescales a numeric vector to a specified range using min-max scaling.
#' This is a generalized form of min-max normalization allowing any output range.
#'
#' @param x Numeric vector. Input values to be scaled.
#' @param ra Numeric vector of length 2. Desired output range (e.g., \code{c(5, 10)}).
#'
#' @return A numeric vector of the same length as \code{x}, scaled to the range \code{ra}.
#'
#' @details
#' The scaled values are computed as:
#' \deqn{x_{scaled} = \frac{x - \min(x)}{\max(x) - \min(x)} \cdot (r_{max} - r_{min}) + r_{min}}
#'
#' @examples
#' x <- rnorm(20)
#' plot(x, type = 'l'); abline(h = range(x), lty = 2)
#' points(scRange(x, ra = c(5, 10)), type = 'l', col = 'red');
#' abline(h = c(5, 10), col = 'red', lty = 2)
#'
#' @seealso [minmax()]
#' @export
scRange <- function(x, ra) {
  ra <- sort(ra)
  (ra[2] - ra[1]) * (x - min(x)) / (max(x) - min(x)) + ra[1]
}

#' @title Min-Max Scaling to \eqn{[0,1]}
#' @description
#' Scales a numeric vector to the range \eqn{[0,1]} using min-max normalization.
#' This is a special case of \code{\link{scRange}}.
#'
#' @param x Numeric vector. Input values to be scaled.
#' @param na.rm Logical; if `TRUE`, ignore `NA`s when computing the range.

#' @return A numeric vector of the same length as \code{x}, scaled to the range \eqn{[0,1]}.
#'
#' @details
#' The scaled values are computed as:
#' \deqn{x_{scaled} = \frac{x - \min(x)}{\max(x) - \min(x)}}
#'
#' Equivalent to \code{scRange(x, ra = c(0, 1))}.
#'
#' @examples
#' x <- rnorm(20)
#' plot(x, type = 'l'); abline(h = range(x), lty = 2)
#' points(minmax(x), type = 'l', col = 'blue')
#' abline(h = c(0, 1), col = 'blue', lty = 2)
#'
#' @seealso [scRange()] for flexible output ranges.
#' @family NMR ++
#' @export
minmax <- function(x, na.rm = FALSE) {
  stopifnot(is.numeric(x))
  r <- range(x, na.rm = na.rm)
  d <- r[2] - r[1]
  if (is.na(d) || d == 0) {
    return(ifelse(is.na(x), NA_real_, 0))
  }
  (x - r[1]) / d
}
