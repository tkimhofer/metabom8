#' @title Plot PCA, PLS or OPLS Model Scores
#' @description Generate score plots for PCA, PLS, or OPLS models with optional annotation, QC highlighting, and cross-validation scores.
#'
#' @param obj A fitted model object of class \code{PCA_metabom8}, \code{OPLS_metabom8}, \code{PCA_MetaboMate}, \code{PLS_MetaboMate}, or \code{OPLS_MetaboMate}.
#' @param pc Numeric or character vector of length 2. Principal components or score components to plot. Defaults: \code{c(1,2)} for PCA/PLS, \code{c("1", "o1")} for OPLS.
#' @param an Optional list of up to 3 elements specifying annotation: \code{list(color, shape, label)}.
#' @param title Optional character. Plot title.
#' @param qc Optional integer vector of row indices indicating QC samples to highlight.
#' @param legend Character. Legend position inside plot (default), or set to \code{NA} to suppress.
#' @param cv Logical. If \code{TRUE}, cross-validated scores are shown when available (OPLS only).
#' @param ... Additional arguments passed to \code{\link[ggplot2]{scale_colour_gradientn}}.
#'
#' @return A \code{ggplot2} object.
#'
#' @references
#' Trygg J., Wold S. (2002) Orthogonal projections to latent structures (O-PLS). \emph{Journal of Chemometrics}, 16(3):119–128.
#' Hotelling H. (1931) The generalization of Student’s ratio. \emph{Annals of Mathematical Statistics}, 2:360–378.
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_vline geom_hline geom_polygon labs theme_bw scale_colour_gradientn scale_x_continuous scale_y_continuous guides unit element_text theme
#' @importFrom rlang sym
#' @importFrom ggrepel geom_text_repel
#' @importFrom colorRamps matlab.like2
#' @importFrom ellipse ellipse
#' @importFrom stats cov
#' @importFrom scales breaks_pretty
#' @export
#' @examples
#' data("covid")
#' y = an$type
#' mod = opls(X, y)
#' plotscores(mod, an=list(Group=y, Hospital=an$hospital))
plotscores <- function(obj, pc, an, title = "", qc, legend = "in", cv = TRUE, ...) {
  # Default PC selection based on model class
  if (missing(pc)) {
    if (grepl("PCA", class(obj)[1])) pc <- c(1, 2)
    if (grepl("PLS", class(obj)[1])) pc <- c(1, 2)
    if (grepl("OPLS", class(obj)[1])) pc <- c("1", "o1")
  }

  # Determine score type (raw or cross-validated)
  etype <- if (cv && grepl("OPLS", class(obj)[1])) "t_cv" else "t"

  res <- .viz_df_helper(obj, pc, an, type = etype)
  sc <- res$df
  type <- if (ncol(sc) > 2 && is.numeric(sc[, 3])) "continuous" else "categorical"

  if (is.null(res$an_le) || is.na(res$an_le)) stop("Check helper function .viz_df_helper")

  # Hotelling's T2 ellipse
  df <- .hotellingsT2(sc[, 1], sc[, 2])

  g <- ggplot() +
    geom_polygon(data = df, aes(x = !!sym("V1"), y = !!sym("V2")),
                 fill = NA, colour = "black", linetype = 3) +
    geom_hline(yintercept = 0, colour = "gray70") +
    geom_vline(xintercept = 0, colour = "gray70") +
    theme_bw() +
    labs(caption = expression("Dashed line: Hotelling's T"^2 * " ellipse (α = 0.95)"))

  # Add points and optional annotation
  aes_args <- aes(x = !!sym(colnames(sc)[1]), y = !!sym(colnames(sc)[2]))
  if (res$an_le >= 1) {
    aes_args <- modifyList(aes_args, aes(colour = !!sym(colnames(sc)[3])))
  }
  if (res$an_le >= 2) {
    aes_args <- modifyList(aes_args, aes(shape = !!sym(colnames(sc)[4])))
  }

  g <- g + geom_point(data = sc, mapping = aes_args, alpha = 0.9)

  if (type == "continuous") {
    g <- g + scale_colour_gradientn(colours = matlab.like2(10), breaks = breaks_pretty(), ...)
  }

  if (res$an_le == 3) {
    g <- g + geom_text_repel(data = sc, aes_string(label = colnames(sc)[5]))
  }

  # Add QC points
  if (!missing(qc) && length(qc) > 0) {
    df_qc <- data.frame(x = sc[qc, 1], y = sc[qc, 2])
    g <- g + geom_point(data = df_qc, aes(x = x, y = y), colour = "black", shape = 21)
  }

  # Axis labeling per model
  g <- switch(class(obj)[1],
              "PCA_metabom8" = {
                g + scale_x_continuous(name = paste0("PC ", pc[1], " (", round(obj@Parameters$R2[pc[1]] * 100, 1), "%)")) +
                  scale_y_continuous(name = paste0("PC ", pc[2], " (", round(obj@Parameters$R2[pc[2]] * 100, 1), "%)"))
              },
              "OPLS_metabom8" = {
                axis_type <- if (etype == "t_cv") ",[cv]" else ""
                g + scale_x_continuous(name = bquote(t[pred]*.(axis_type))) +
                  scale_y_continuous(name = bquote(t[orth]*.(axis_type))) +
                  # labs(title = paste0("OPLS-", obj@type, ": 1+", obj@nPC - 1, " comp.",
                  #                     " (R2X=", round(obj@summary$R2X[obj@nPC - 1], 2),
                  #                     ", Q2=", round(obj@summary$Q2[obj@nPC - 1], 2),
                  #                     if (!is.null(obj@summary$AUROC[obj@nPC - 1])) paste0(", AUROC=", round(obj@summary$AUROC[obj@nPC - 1], 2)) else "",
                  #                     ")")) +
                  labs(title = paste0("OPLS-", obj@type, ": 1+", obj@nPC - 1, " comp.",
                                      " (R2X=", round(obj@summary$R2X[obj@nPC - 1], 2),
                                      if (!is.null(obj@summary$AUROC[obj@nPC - 1])) paste0(", AUROC=", round(obj@summary$AUROC[obj@nPC - 1], 2)) else "",
                                      ")")) +
                  theme(plot.title = element_text(size = 10))
              },
              g
  )

  # Legend positioning
  if (legend == "in") {
    g <- g + theme(legend.position = c(1.01, 1.02), legend.justification = c(1, 1))
  } else if (is.na(legend)) {
    g <- g + theme(legend.position = "none")
  }

  return(g)
}
