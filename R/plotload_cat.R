#' #' @title Plot Loadings from PCA or O-PLS Models
#' #'
#' #' @description
#' #' Generates a scatter plot of variable loadings from PCA or O-PLS models using \pkg{ggplot2}.
#' #' The function supports flexible annotation of points using color, shape, and labels.
#' #' O-PLS plots display the predictive component on the x-axis and an orthogonal component on the y-axis.
#' #'
#' #' @param obj An object of class \code{"PCA_metabom8"} or \code{"OPLS_metabom8"}.
#' #' Typically created by the \code{\link{pca}} or \code{\link{opls}} functions.
#' #'
#' #' @param pc A numeric vector specifying the components to be plotted.
#' #' For PCA models, use \code{c(1, 2)} for PC1 vs PC2.
#' #' For O-PLS models, use a character vector like \code{c("1", "o1")}, where \code{"o1"} refers to the first orthogonal component.
#' #'
#' #' @param an Optional named list for point annotations with up to three elements:
#' #' \describe{
#' #'   \item{\code{colour}}{A vector indicating groupings used to color the points.}
#' #'   \item{\code{shape}}{A factor or grouping variable for point shapes.}
#' #'   \item{\code{label}}{A character vector of labels for each point (e.g., variable names).}
#' #' }
#' #'
#' #' @param title Optional character string to be used as the plot title.
#' #'
#' #' @param legend Logical; if \code{TRUE}, displays the plot legend. Set to \code{FALSE} to suppress the legend, useful when plotting many groups.
#' #'
#' #' @param ... Additional arguments passed to \code{\link[ggplot2]{geom_point}}.
#' #'
#' #' @return A \code{\link[ggplot2]{ggplot}} object representing the loadings scatter plot.
#' #' This object can be further customized using \pkg{ggplot2} syntax.
#' #'
#' #' @details
#' #' This function is intended for use with PCA or O-PLS models generated using the \pkg{metabom8} package.
#' #' For PCA, the x- and y-axes correspond to the selected principal components.
#' #' For O-PLS, the x-axis always corresponds to the predictive component, while the y-axis displays the specified orthogonal component.
#' #'
#' #' @seealso \code{\link{pca}}, \code{\link{opls}}, \code{\link[ggplot2]{ggplot}}, \code{\link{plotscores}}
#' #'
#' #' @references
#' #' Geladi, P., & Kowalski, B. R. (1986). Partial least squares regression: A tutorial. \emph{Analytica Chimica Acta}, 185, 1–17.
#' #'
#' #' @importFrom ggplot2 ggplot aes_string geom_hline geom_vline geom_point labs theme_bw guides theme
#' #' @importFrom ggrepel geom_label_repel
#' #' @importFrom methods is
#' #'
#' #' @family Visualization
#' #' @examples
#' #' set.seed(123)
#' #' X <- matrix(rnorm(200), 10, 20)
#' #' grp <- rep(c("A", "B"), each = 5)
#' #' X[grp == "B", 1:5] <- X[grp == "B", 1:5] + 2
#' #' X[grp == "B", 5] <- X[grp == "B", 4] * 0.9 + rnorm(5, 0, 0.1)
#' #' colnames(X) <- paste0("metabolite", 1:20)
#' #'
#' #' mod = opls(X, grp)
#' #' plotload_cat(mod)
#' plotload_cat <- function(obj, pc, an, title = NULL, legend = TRUE, ...) {
#'
#'   if (!inherits(obj, c("PCA_metabom8",'PLS_metabom8', "OPLS_metabom8"))) {
#'     stop("`obj` must be of class 'PCA_metabom8', 'PLS_metabom8' or 'OPLS_metabom8'.")
#'   }
#'
#'   # Default PC
#'   if (missing(pc)) {
#'     pc <- if (inherits(obj, "PCA_metabom8") || inherits(obj, "PLS_metabom8")) c(1, 2) else c("1", "o1")
#'   }
#'
#'   # Default annotations
#'   if (missing(an)) {
#'     if (ncol(obj@X) < 150) {
#'       an <- list(colour = "All", shape = "All", label = colnames(obj@X))
#'     } else {
#'       stop("Annotation list 'an' is missing and too many features to set a default.")
#'     }
#'   }
#'
#'   res <- .viz_df_helper(obj, pc, an, type = "p")
#'   if (is.null(res$an_le) || is.na(res$an_le)) {
#'     stop("Check helper function `.viz_df_helper`.")
#'   }
#'
#'   melted <- res$df
#'
#'   d <- ggplot(data = melted, aes(x = !!sym(colnames(melted)[1]), y = !!sym(colnames(melted)[2]))) +
#'     geom_hline(yintercept = 0, colour = "gray70") +
#'     geom_vline(xintercept = 0, colour = "gray70")
#'
#'   if (res$an_le == 1) {
#'     d <- d + aes(colour = !!sym(names(an)[1])) +
#'       geom_point() +
#'       labs(colour = names(an)[1], title = title)
#'   }
#'
#'   if (res$an_le == 2) {
#'     d <- d + aes(colour = !!sym(colnames(melted)[3]), shape = !!sym(colnames(melted)[4])) +
#'       geom_point() +
#'       labs(colour = names(an)[1], shape = names(an)[2], title = title)
#'   }
#'
#'   if (res$an_le == 3) {
#'     d <- d + aes(colour = !!sym(colnames(melted)[3]), shape = !!sym(colnames(melted)[4])) +
#'       geom_point() +
#'       geom_label_repel(aes(label = !!sym(colnames(melted)[5]))) +
#'       labs(colour = names(an)[1], shape = names(an)[2], title = title)
#'   }
#'
#'   d <- d +
#'     labs(
#'       x = colnames(melted)[1],
#'       y = colnames(melted)[2]
#'     ) +
#'     theme_bw()
#'
#'   if (!legend) {
#'     d <- d + guides(colour = FALSE, shape = FALSE)
#'   }
#'
#'   return(d)
#' }
