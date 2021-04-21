#' @title Cliff's delta (effect size)
#' @export
#' @description Non-parametric effect size estimate used to quantify the degree to which two variable distributions overlap. Output is scaled to range between -1 to 1, with the arithmetic sign indicating if the comparator group lies to the right (higher) or left (lower) to a the reference.
#' @param ref num, reference group
#' @param comp num, comparator group
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @references Cliff, N (1993), Dominance statistics: Ordinal analyses to answer ordinal quetsions \emph{Psychological Bulletin}, 114.(3), 494-509.<https://doi.org/10.1037%2F0033-2909.114.3.494>
#' @return This function returns Cliff's delta effect size ranging from -1 to 1
## #' @seealso \code{\link{OPLS_MetaboMate-class}} \code{\link{dmodx}}
## \code{\link{plotscores}} \code{\link{plotload}} \code{\link{specload}}
#' @examples
#' # define two distrubutions
#' a=rnorm(100, mean=0, sd=1)
#' b=rnorm(100, mean=2, sd=1)
#'
#' # Cliff's delta
#'  es_cdelta(ref=a, comp=b)
#'  es_cdelta(ref=b, comp=a)
#' @family NMR ++
#' @importFrom base sapply length
es_cdelta <- function(ref, comp) {

    if (!is.numeric(ref) | !is.numeric(comp))
        stop("Input not numeric")

    ind_ref <- is.na(ref) | is.infinite(ref)
    ind_comp <- is.na(comp) | is.infinite(comp)

    if (any(ind_ref)) {
        ref <- ref[!ind_ref]
        message("Removing NA or infinite values from reference group.")
    }
    if (any(ind_comp)) {
        comp <- comp[!ind_comp]
        message("Removing NA or infinite values from comparator group.")
    }

    if (length(ref) < 5 | length(comp) < 5)
        stop("Low number of values (< 5)")

    top_counts <- vapply(ref, function(x, y = comp) {
        names(x) <- NULL
        names(y) <- NULL
        c(length(which(x > y)), length(which(x < y)))
    }, FUN.VALUE = c(2, length(ref)))
    out <- ((sum(top_counts[1, ]) - sum(top_counts[2, ]))/(length(ref) * length(comp))) *
        (-1)
    return(out)
}

