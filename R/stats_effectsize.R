#' @title Cliff's Delta Effect Size
#' @description
#' Calculates Cliff's delta (\code{Cd}) as a non-parametric effect size for
#' comparing two numeric vectors. \code{Cd} quantifies the directional difference
#' between a reference (\code{ref}) and comparison (\code{comp}) distribution
#' based on pairwise comparisons of their values.
#'
#' \code{Cd} ranges from -1 to 1, where:
#' \itemize{
#'   \item{-1: values in \code{comp} are consistently larger than those in \code{ref}}
#'   \item{0: both distributions are similar, with no systematic difference}
#'   \item{1: values in \code{ref} are consistently larger than those in \code{comp}}
#' }
#'
#' @param ref Numeric vector representing the reference group.
#' @param comp Numeric vector representing the comparator group.
#' @return A single numeric value: Cliff's delta effect size.
#' @details
#' The effect size is calculated as the scaled difference in dominance between groups.
#' Missing or infinite values are removed with a message. Synonyms are
#' @references
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal questions.
#' \emph{Psychological Bulletin}, 114(3), 494–509. \doi{10.1037/0033-2909.114.3.494}
#' @examples
#' ref <- rnorm(100, mean = 0)
#' comp <- rnorm(100, mean = 1)
#'
#' hist(ref,  col=rgb(0,0,1,0.3), breaks=50, xlim=c(-4,4))
#' hist(comp, col=rgb(1,0,0,0.3), breaks=50, add=TRUE)
#' legend("topright", legend=c("ref","comp"), fill=c("blue","red"))
#'
#' cliffs_d(ref, comp)
#' cliffs_d(comp, ref)
#'
#' comp1 <- rnorm(100, mean = 10)
#' cliffs_d(ref, comp1)
#'
#' cliffs_d(ref, ref)
#' @export
cliffs_d <- function(ref, comp) {
  if (!is.numeric(ref) || !is.numeric(comp))
    stop("'ref' and 'comp' must be numeric.", call. = FALSE)

  ref0 <- length(ref); comp0 <- length(comp)

  ref <- ref[is.finite(ref)]
  comp <- comp[is.finite(comp)]

  if (length(ref) < ref0) message("Removed non-finite values from 'ref'")
  if (length(comp) < comp0) message("Removed non-finite values from 'comp'")

  n <- length(ref); m <- length(comp)
  if (n < 5 || m < 5)
    stop("Each group must contain at least 5 finite values.", call. = FALSE)

  r <- rank(c(ref, comp), ties.method = "average")
  R1 <- sum(r[seq_len(n)])
  U  <- R1 - n*(n+1)/2

  delta <- (2*U)/(n*m) - 1
  delta
}


#' @rdname cliffs_d
#' @export
es_cdelta <- cliffs_d
