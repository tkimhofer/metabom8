#' @title Hotelling's T² Ellipse in 2D
#'
#' @description Computes the Hotelling's T² ellipse for a bivariate dataset, commonly used to visualise multivariate confidence intervals in score plots.
#'
#' @param x Numeric vector. First dimension (e.g., PC1 or \code{t[1]}).
#' @param y Numeric vector. Second dimension (e.g., PC2 or \code{t[2]}).
#' @param alpha Numeric. Confidence level for the ellipse. Default is 0.95.
#'
#' @return A \code{data.frame} with two columns \code{V1} and \code{V2} containing coordinates of the T² ellipse.
#'
#' @references
#' Geladi, P. and Kowalski, B.R. (1986). Partial least squares regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1–17.
#'
#' @author Torben Kimhofer
#'
#' @importFrom stats cov
#' @importFrom ellipse ellipse
#' @keywords internal
#' @family NMR ++
.hotellingsT2 <- function(x, y, alpha = 0.95) {
  if (!is.numeric(x) || !is.numeric(y)) stop("Input vectors must be numeric.")
  if (length(x) != length(y)) stop("x and y must be the same length.")
  if (length(x) < 3) stop("At least 3 observations are required to compute the ellipse.")

  XY <- cbind(x, y)
  cov_matrix <- cov(XY, use = "complete.obs")
  center <- colMeans(XY, na.rm = TRUE)

  ell <- ellipse::ellipse(cov_matrix, centre = center, level = alpha)
  colnames(ell) <- c("V1", "V2")

  return(as.data.frame(ell))
}
