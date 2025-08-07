#' @title Eruption Plot for OPLS-DA Model Interpretation
#'
#' @description
#' Generates an eruption plot to visualise variable importance in an OPLS-DA model.
#' Variables are plotted by their model loading (from predictive or orthogonal component),
#' effect size (Cliff's delta), and adjusted p-value from the Kruskal-Wallis rank sum test.
#'
#' @param mod An OPLS model object of class \code{OPLS_metabom8}, generated via \code{\link{opls}}.
#' @param pc Integer or character. Component to plot: use \code{1} for predictive,
#'        or \code{'o1'}, \code{'o2'}, etc., for orthogonal components.
#' @param p_adj Character string. Method for p-value adjustment. See \code{\link[stats]{p.adjust}}.
#' @param invert_es Logical. If \code{TRUE}, swaps the reference group for Cliff's delta calculation.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{data}}{A \code{data.frame} with model loadings, Cliff's delta, and (adjusted) p-values.}
#'   \item{\code{plot}}{A \code{ggplot} object for visualisation.}
#' }
#'
#' @details
#' Designed for binary classification OPLS-DA models. Effect size is calculated with Cliff's delta,
#' and variable-level p-values are obtained via Kruskal-Wallis tests. Color indicates the negative log10
#' of (adjusted) p-values. This plot is designed for a low nb of features, ie. not full res NMR data
#'
#' @references
#' Torben Kimhofer
#' Eriksson, L. et al. (2008). CV-ANOVA for significance testing of PLS and OPLS models.
#' \emph{J. Chemometrics}, 22, 594–600.
#'
#' @importFrom ggplot2 ggplot aes_string geom_point scale_colour_gradientn
#'   theme_bw theme element_blank element_text labs scale_x_continuous
#' @importFrom ggrepel geom_label_repel
#' @importFrom colorRamps matlab.like2
#' @importFrom stats p.adjust kruskal.test
#' @seealso \code{\link{opls}}
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(200), 10, 20)
#' grp <- rep(c("A", "B"), each = 5)
#' X[grp == "B", 1:5] <- X[grp == "B", 1:5] + 2
#' X[grp == "B", 5] <- X[grp == "B", 4] * 0.9 + rnorm(5, 0, 0.1)
#' colnames(X) <- paste0("metabolite", 1:20)
#'
#' mod = opls(X,grp)
#' er = eruption(mod)
#' head(er[[1]]); plot(er[[2]])
eruption <- function(mod, pc = 1, p_adj = "BH", invert_es = FALSE) {

  if (!inherits(mod, "OPLS_metabom8")) {
    stop("Function requires an OPLS_metabom8 object.")
  }

  if (length(unique(mod@Y$ori)) != 2) {
    stop("Eruption plot is defined only for two-level classification.")
  }

  if (is.na(pc) || length(pc) != 1 || !is.character(pc) && !is.numeric(pc)) {
    stop("Argument 'pc' must be a single integer or string (e.g., 'o1').")
  }

  if (!is.character(p_adj) || length(p_adj) != 1) {
    p_adj <- "none"
  }

  # Select loading
  if (grepl("o", pc)) {
    pc_num <- as.numeric(gsub("o", "", pc))
    if (is.na(pc_num) || pc_num > nrow(mod@p_orth)) {
      stop("Invalid orthogonal component. Check `pc` argument.")
    }
    ddl <- data.frame(p1 = mod@p_orth[pc_num, ], id = colnames(mod@X))
  } else {
    ddl <- data.frame(p1 = mod@p_pred[1, ], id = colnames(mod@X))
  }

  # Set comparison group
  Y <- mod@Y$ori
  levels_Y <- unique(Y)

  comp_group <- if (invert_es) levels_Y[2] else levels_Y[1]
  message("Using ", if (invert_es) levels_Y[1] else levels_Y[2],
          " as reference group for Cliff's delta.")

  # Univariate stats
  idx <- which(Y == comp_group)
  uni_stats <- t(apply(mod@X, 2, function(x) {

    c(es_cdelta(x[idx], x[-idx]), kruskal.test(x, Y)$p.value)
  }))

  colnames(uni_stats) <- c("Cd", "pval")
  ddl <- cbind(ddl, uni_stats)
  ddl$p1 <- abs(ddl$p1)

  if (p_adj %in% p.adjust.methods) {
    ddl$pval_adjusted <- p.adjust(ddl$pval, method = p_adj)
    ddl$pval_transformed <- abs(log10(ddl$pval_adjusted))
  } else {
    ddl$pval_transformed <- abs(log10(ddl$pval))
  }

  # Plot
  gl2 <- ggplot(ddl, aes(x = !!sym("Cd"), y = !!sym("p1"), colour = !!sym("pval_transformed"))) +
    geom_label_repel(aes(label = !!sym("id")), colour = "black", min.segment.length = 0.001) +
    geom_point(size = 3, shape = 1) +
    geom_point(size = 3, shape = 16, alpha = 0.3) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_colour_gradientn(colours = matlab.like2(10)) +
    labs(
      x = "Cliff's Delta",
      y = paste0("|p_", pc, "|"),
      colour = if (p_adj %in% p.adjust.methods)
        expression("| log10( p value"["adj"] ~ ") |") else
          "| log10( p value ) |"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.tag = element_text(face = "bold", size = 25),
      legend.position = "top",
      legend.direction = "horizontal"
    )

  ddl$pval_transformed <- NULL
  list(data = ddl, plot = gl2)
}
