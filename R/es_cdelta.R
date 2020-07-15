#' @title Cliff's delta (effect size)
#' @export
#' @description Non-parametric effect size estimate used to quantify the degree to which one population's variable distribution lies to the right or left of a reference group distribution.
#' @param ref num, reference group
#' @param comp num, comparator group
#' @references Cliff, N (1993), Dominance statistics: Ordinal analyses to answer ordinal quetsions \emph{Psychological Bulletin}, 114.(3), 494-509.<https://doi.org/10.1037%2F0033-2909.114.3.494>
#' @return This function returns Cliff's delta effect size ranging from -1 to 1
## #' @seealso \code{\link{OPLS_MetaboMate-class}} \code{\link{dmodx}} \code{\link{plotscores}} \code{\link{plotload}} \code{\link{specload}}
#' @examples
#' # define two distrubutions
#' ref=rnorm(100, mean=0, sd=1)
#' comp=rnorm(100, mean=2, sd=1)
#'
#' # Cliff's delta
#'  es_cdelta(comp, ref)
#'  es_cdelta(ref, comp)
#' @author Torben Kimhofer \email{torben.kimhofer@murdoch.edu.au}
#' @importFrom base sapply length
es_cdelta = function(ref, comp){

  ind_ref=is.na(ref) | is.infinite(ref)
  ind_comp=is.na(comp) | is.infinite(comp)

  if(any(ind_ref)) {ref=ref[!ind_ref]; message('Removing NA or infinite values from reference group.')}
  if(any(ind_comp)) {comp=comp[!ind_comp]; message('Removing NA or infinite values from comparator group.')}

  if(length(ref) < 5 | length(comp) < 5) stop('Low number of values (< 5)')

  top_counts=sapply(ref, function(x, y=comp){
    c(length(which(x>y)), length(which(x<y)))
  })
  out=(sum(top_counts[1,])-sum(top_counts[2,]))/(length(ref)*length(comp))
  return(out)
}
