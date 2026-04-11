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



#' @title Applies a preprocessing strategy to a numeric matrix.
#' @param preproc_strategy list with elements 'center' (logical) and 'scale' (character).
#' @param X numeric matrix - processing is done column-wise.
#' @return A list containing the processed matrix and parameters.
#' @examples
#' X <- matrix(c(10,10, 0,0, 0, 10), ncol=3)
#' autoscale <- uv_scaling(center=TRUE)
#' X_scaled <- prep_X(autoscale, X)
#' str(X_scaled)
#' @export
prep_X <- function(preproc_strategy, X){

            msd <- .sdRcpp(X)

            Xcs <- .scaleMatRcpp(
              X,
              0:(nrow(X)-1),
              center = preproc_strategy$center,
              scale_type =
                switch(preproc_strategy$scale,
                       "None"=0,
                       "UV"=1,
                       "Pareto"=2)
            )[[1]]

            preproc_strategy$X_mean <- msd[[2]]
            preproc_strategy$X_sd   <- msd[[1]]

            list(X=Xcs, prep=preproc_strategy)
          }


#' @title No Scaling
#' This function defines a preprocessing strategy that is applied via
#' \code{\link{prep_X}}.
#' @param center Logical. If TRUE, variables are mean-centered before scaling.
#' @return A \code{list} with elements:
#' \describe{
#'   \item{X}{Numeric matrix containing the scaled data.}
#'   \item{prep}{List describing the preprocessing, including centering and scaling
#'   parameters (\code{center}, \code{scale}, \code{X_mean}, \code{X_sd}).}
#' }#' @details Leaves variables unscaled. Optional centering.
#' @examples
#' just_centering <- unscaled(center=TRUE)
#' X <- matrix(c(10,10, 0,0, 0, 10), ncol=3)
#' X_centered <- prep_X(just_centering, X)
#' str(X_centered)
#' X_centered$X
#' @family scaling_strategies
#' @export
unscaled <- function(center = FALSE)
  list(center=center, scale="None")

#' @title Unit Variance Scaling
#' This function defines a preprocessing strategy that is applied via
#' \code{\link{prep_X}}.
#' @param center Logical. If TRUE, variables are mean-centered before scaling.
#' @return A \code{list} with elements:
#' \describe{
#'   \item{X}{Numeric matrix containing the scaled data.}
#'   \item{prep}{List describing the preprocessing, including centering and scaling
#'   parameters (\code{center}, \code{scale}, \code{X_mean}, \code{X_sd}).}
#' }
#' @details
#' Centers variables and scales each feature to unit variance.
#' UV scaling divides each variable by its standard deviation.
#' This is the default scaling in many multivariate methods such as PCA and PLS.
#' @family scaling_strategies
#' @examples
#' autoscale <- uv_scaling(center=TRUE)
#' X <- matrix(c(10,10, 0,0, 0, 10), ncol=3)
#' X_scaled <- prep_X(autoscale, X)
#' str(X_scaled)
#' X_scaled$X
#' @export
uv_scaling <- function(center = TRUE)
  list(center=center, scale="UV")


#' @title Pareto Scaling
#' Leaves variables unscaled. Optional centering.
#' @param center Logical. If TRUE, variables are mean-centered before scaling.
#' @return A \code{list} with elements:
#' \describe{
#'   \item{X}{Numeric matrix containing the scaled data.}
#'   \item{prep}{List describing the preprocessing, including centering and scaling
#'   parameters (\code{center}, \code{scale}, \code{X_mean}, \code{X_sd}).}
#' }
#' @details
#' Scales variables by the square root of their standard deviation.
#' @family scaling_strategies
#' @examples
#' paritalUV <- pareto_scaling(center=TRUE)
#' X <- matrix(c(10,10, 0,0, 0, 10, 0, 1000), ncol=4)
#' X_scaled <- prep_X(paritalUV, X)
#' str(X_scaled)
#' X_scaled$X
#' @export
pareto_scaling <- function(center = FALSE)
  list(center=center, scale="Pareto")
