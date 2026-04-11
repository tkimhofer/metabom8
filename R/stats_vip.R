#' Compute Variable Importance in Projection (VIP)
#'
#' Calculates VIP scores for a PLS/OPLS model using predictive
#' components only.
#'
#' @param Tx Numeric matrix (n x A). Score matrix.
#' @param W Numeric matrix (A x p). Weight matrix.
#'          Rows correspond to components, columns to variables.
#' @param C Numeric vector or matrix of length A. Y-loadings.
#'
#' @details
#' VIP is computed as:
#'
#' \deqn{
#' VIP_j = \sqrt{ p \frac{\sum_a SSY_a w_{aj}^2}{\sum_a SSY_a} }
#' }
#'
#' where
#'   SSY_a = (t_a^T t_a) * c_a^2
#'
#' By construction, mean(VIP^2) = 1.
#'
#' @return Numeric vector of length p containing VIP scores.
#'
#' @keywords internal
.vip <- function(Tx, W, C) {

  if (!is.matrix(Tx))
    stop("T must be a numeric matrix (n x A).", call. = FALSE)

  if (!is.matrix(W))
    stop("W must be a numeric matrix (A x p).", call. = FALSE)

  if (!is.numeric(Tx) || !is.numeric(W))
    stop("T and W must be numeric.", call. = FALSE)

  A <- nrow(W)
  p <- ncol(W)

  if (is.matrix(C)) {

    if (!all(dim(C) %in% c(1L, A)) || length(C) != A)
      stop("C must be length A, or a 1 x A / A x 1 matrix.", call. = FALSE)
    C <- as.numeric(C)
  } else if (is.numeric(C)) {
    if (length(C) != A)
      stop("Length of C must equal number of components (A).", call. = FALSE)
  } else {
    stop("C must be numeric.", call. = FALSE)
  }

  if (ncol(Tx) != A)
    stop("Number of components in T must match rows of W.", call. = FALSE)

  C <- as.numeric(C)

  if (length(C) != A)
    stop("Length of C must equal number of components.", call. = FALSE)

  SSY <- colSums(Tx^2) * (C^2)
  if (all(SSY == 0))
    stop("All SSY values are zero; cannot compute VIP.", call. = FALSE)
  Wnorm <- W / sqrt(rowSums(W^2))
  vip_sq <- p * colSums((Wnorm^2) * SSY) / sum(SSY)

  vip <- sqrt(vip_sq)

  vip
}
