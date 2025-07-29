#' @title Cliff's Delta Effect Size
#'
#' @description
#' Calculates Cliff's delta, a non-parametric effect size that quantifies the degree of overlap between two distributions. The value ranges from -1 to 1, where:
#' \itemize{
#'   \item{-1: all values in \code{comp} are greater than those in \code{ref}}
#'   \item{0: complete overlap between distributions}
#'   \item{1: all values in \code{ref} are greater than those in \code{comp}}
#' }
#'
#' @param ref Numeric vector representing the reference group.
#' @param comp Numeric vector representing the comparator group.
#'
#' @return A single numeric value: Cliff's delta.
#'
#' @details
#' The effect size is calculated as the scaled difference in dominance between groups.
#' Missing or infinite values are removed with a message.
#'
#' @references
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal questions.
#' \emph{Psychological Bulletin}, 114(3), 494–509. \doi{10.1037/0033-2909.114.3.494}
#'
#' @seealso \code{\link{eruption}}, \code{\link{opls}}, \code{\link{dmodx}}
#'
#' @examples
#' ref <- rnorm(100, mean = 0)
#' comp <- rnorm(100, mean = 1)
#' es_cdelta(ref, comp)
#'
#' @export
es_cdelta <- function(ref, comp) {
  if (!is.numeric(ref) || !is.numeric(comp)) {
    stop("Both 'ref' and 'comp' must be numeric vectors.")
  }

  ref <- ref[is.finite(ref)]
  comp <- comp[is.finite(comp)]

  if (length(ref) < 5 || length(comp) < 5) {
    stop("Each group must contain at least 5 finite values.")
  }

  # Compare each pair of ref vs comp values
  dominance <- vapply(ref, function(x) {
    c(sum(x > comp), sum(x < comp))
  }, FUN.VALUE = numeric(2))

  delta <- (sum(dominance[1, ]) - sum(dominance[2, ])) / (length(ref) * length(comp))
  return(-delta)
}
