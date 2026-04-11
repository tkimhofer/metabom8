#' @title Hotelling T^2 Statistic
#'
#' @description Computes the Hotelling T^2 statistic defining a multivariate
#' confidence region for a set of observations. The region corresponds to an
#' ellipsoid in p-dimensional space and is commonly visualised as an ellipse
#' in two-dimensional (OPLS) score plots.
#' @param X Numeric matrix with observations in rows and variables (dimensions)
#' in columns.
#' @param alpha Numeric scalar. Confidence level for the T^2 region. Default is 0.95.
#' @return A named \code{list} with elements:
#' \describe{
#'   \item{center}{Numeric vector of column means.}
#'   \item{cov}{Sample covariance matrix.}
#'   \item{c2}{Squared radius of the Hotelling T^2 region.}
#' }
#' @details
#' The Hotelling T^2 region is defined as
#' \deqn{(x - \mu)^T S^{-1} (x - \mu) \le c^2}
#' where \eqn{c^2 = (p (n - 1) / (n - p)) F_{p, n-p}(\alpha)}.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200), ncol = 2)
#' t2 <- hotellingsT2(X)
#' ell <- ellipse2d(t2)
#' plot(X, asp = 1)
#' lines(ell$x, ell$y, col = "red", lwd = 2)
#' @importFrom stats cov qf
#' @export
hotellingsT2 <- function(X, alpha = 0.95) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  center <- colMeans(X)
  S <- cov(X)

  c2 <- (p * (n - 1) / (n - p)) * qf(alpha, p, n - p)

  list(center = center, cov = S, c2 = c2)
}


#' @title Calculate 2D Hotelling T^2 Ellipse
#' @description Generates coordinates of a two-dimensional ellipse corresponding
#' to a Hotelling T^2 region projected onto selected dimensions.
#' @param obj A \code{list} as returned by \code{hotellingsT2()}.
#' @param dims Integer vector of length 2 specifying which dimensions to project
#' onto. Default is \code{c(1, 2)}.
#' @param npoints Integer. Number of points used to approximate the ellipse.
#' Default is 100.
#' @return A \code{data.frame} with columns \code{x} and \code{y} containing
#' ellipse coordinates.
#' @details
#' The ellipse is obtained by projecting the multivariate T^2 ellipsoid onto
#' the specified dimensions. The scaling is derived from the Hotelling T^2
#' statistic and accounts for sample size and dimensionality.
#' @examples
#' set.seed(1)
#' X <- cbind(rnorm(100), rnorm(100) + 0.5)
#'
#' t2 <- hotellingsT2(X)
#' ell <- ellipse2d(t2)
#'
#' plot(X[,1], X[,2], asp = 1, pch = 16,
#'      xlab = "Score 1", ylab = "Score 2")
#' lines(ell$x, ell$y, col = "red", lwd = 2)
#' @importFrom ellipse ellipse
#' @export
ellipse2d <- function(obj, dims = c(1,2), npoints = 100) {
  S <- obj$cov[dims, dims]
  center <- obj$center[dims]

  ell <- ellipse::ellipse(S, centre = center, t = sqrt(obj$c2), npoints = npoints)
  as.data.frame(ell)
}
