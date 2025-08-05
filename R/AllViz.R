#' @title Plot a Single NMR Spectrum
#'
#' @description
#' Plots a 1D NMR spectrum using either base R or interactive Plotly graphics.
#'
#' @param x Numeric vector. NMR spectrum intensity values.
#' @param ppm Numeric vector. Chemical shift positions.
#' @param shift Numeric length-2 vector. Chemical shift range to plot (e.g., \code{c(0, 11)}).
#' @param interactive Logical. Use plotly for interactive output? Defaults to TRUE.
#' @param name Character. Label for trace (Plotly only).
#' @param mode Character. Plotly trace mode: 'lines', 'markers', or 'lines+markers'.
#' @param base Logical. If TRUE, use base R plot (overrides \code{interactive}).
#' @param add Logical. If TRUE, add to an existing base R plot.
#' @param ... Additional arguments passed to base \code{plot()} or \code{points()}.
#'
#' @return If \code{interactive = TRUE}, returns a \code{plotly} object. Otherwise, returns \code{NULL}.
#'
#' @importFrom graphics plot points
#' @importFrom plotly plot_ly add_trace layout
#' @importFrom magrittr %>%
#' @examples
#' data(covid)
#' p <- spec(X[1, ], ppm, shift = c(0, 11), interactive = TRUE)
#' if (interactive()) p
#'
#' spec(X[2, ], ppm, shift = c(0, 11), interactive = FALSE)
#' spec(X[3, ], ppm, shift = c(0, 11), interactive = FALSE, add = TRUE, col = 'red')
#' @export
spec <- function(x, ppm, shift = c(0, 11), interactive = TRUE, name = "A", mode = "lines",
                 base = FALSE, add = FALSE, ...) {

  if (!is.null(ncol(x))) stop("x must be a numeric vector (single spectrum).")
  if (length(x) != length(ppm)) stop("x and ppm must have equal length.")

  idx <- get_idx(shift, ppm)
  ppm_sub <- ppm[idx]
  x_sub <- x[idx]

  if (interactive && !base) {
    plotly::plot_ly(x = ~ppm_sub, y = ~x_sub, type = "scatter", mode = mode, name = name,
                    hovertemplate = "%{x} ppm<extra></extra>") %>%
      plotly::layout(
        xaxis = list(title = "ppm", autorange = "reversed"),
        yaxis = list(title = "Intensity")
      )
  } else {
    if (add) {
      graphics::points(ppm_sub, x_sub, type = "l", ...)
    } else {
      graphics::plot(ppm_sub, x_sub, type = "l", xlim = rev(range(ppm_sub)),
                     xlab = "ppm", ylab = "Intensity", ...)
    }
    invisible(NULL)
  }
}



#' @title Plot Overlayed NMR Spectra
#'
#' @description
#' Plots multiple 1D NMR spectra using either base R or interactive Plotly graphics.
#'
#' @param X Numeric matrix. NMR data matrix with spectra in rows.
#' @param ppm Numeric vector. Chemical shift axis (must match \code{ncol(X)}).
#' @param shift Numeric vector of length 2. Chemical shift window to plot (e.g., \code{c(0, 9.5)}).
#' @param interactive Logical. Use plotly for interactive plotting (default: \code{TRUE}).
#' @param ... Additional arguments passed to \code{matplot()} for non-interactive plots.
#'
#' @return A \code{plotly} object if \code{interactive = TRUE}; otherwise, \code{NULL}.
#'
#' @importFrom plotly plot_ly add_lines layout
#' @importFrom graphics matplot
#' @importFrom reshape2 melt
#' @importFrom grDevices colorRampPalette
#' @examples
#' data(covid)
#' matspec(X[1:3,], ppm, interactive = TRUE)
#' matspec(X[1:3,], ppm, interactive = FALSE)
#' @seealso \code{\link{spec}}, \code{\link{plot}}
#' @family NMR
#' @export
matspec <- function(X, ppm, shift = c(0, 9.5), interactive = TRUE, ...) {
  if (!.check_X_ppm(X, ppm)) {
    stop("X and ppm dimensions do not match or ppm contains NA/Inf values.")
  }

  idx <- get_idx(shift, ppm)
  ppm_sub <- ppm[idx]
  X_sub <- X[, idx, drop = FALSE]

  if (interactive) {
    df <- reshape2::melt(X_sub)
    colnames(df) <- c("spectrum", "index", "intensity")
    df$ppm <- rep(ppm_sub, each = nrow(X_sub))

    col_fun <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
    col_vec <- col_fun(nrow(X_sub))
    df$color <- rep(col_vec, times = length(idx))

    p <- plotly::plot_ly(data = df, x = ~ppm, y = ~intensity, color = ~I(color),
                         type = "scatter", mode = "lines", name = ~spectrum,
                         hovertemplate = "%{x} ppm<extra></extra>") %>%
      plotly::layout(
        xaxis = list(title = "ppm", autorange = "reversed"),
        yaxis = list(title = "Intensity")
      )

    return(p)
  }

  graphics::matplot(ppm_sub, t(X_sub), type = "l", xlab = "ppm", ylab = "Intensity",
                    xlim = rev(range(ppm_sub)), ...)
  invisible(NULL)
}



#' @title Overlay NMR Spectra Using ggplot2
#'
#' @description
#' Overlay multiple NMR spectra in a specified chemical shift region using `ggplot2`.
#' Faceting, coloring, and line types can be controlled via an annotation list. Useful
#' for QC and visual inspection of subsets or trends in spectral data.
#'
#' @param X Numeric matrix. NMR data with spectra in rows.
#' @param ppm Numeric vector. Chemical shift values corresponding to columns of `X`.
#' @param shift Numeric vector of length 2. Region in ppm to plot (recommended to be small for performance).
#' @param an Named list of 1–3 elements for grouping: 1st for facetting, 2nd for color, 3rd for line type.
#' @param alp Numeric. Alpha transparency of lines (0 to 1).
#' @param size Numeric. Line width (e.g., 0.5).
#' @param title Character. Plot title.
#' @param ... Additional arguments passed to `facet_grid`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' The first element in `an` defines the facet grouping and must always be specified. Additional
#' elements define coloring and line type. Use short ppm ranges (e.g., 5.15–4.6) to improve performance.
#'
#' @examples
#' data(covid)
#' specOverlay(X, ppm, shift = c(5.15, 4.6), an = list(Group = meta$Group))
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_line aes scale_x_reverse ggtitle xlab facet_grid theme_bw theme element_text scale_y_continuous
#' @importFrom colorRamps matlab.like2
#' @importFrom scales breaks_pretty
#' @importFrom stats as.formula
#' @export
specOverlay <- function(X, ppm, shift = c(0, 0.1), an = list("Group"), alp = 0.7, size = 0.5, title = "", ...) {
  if (!.check_X_ppm(X, ppm)) {
    stop("Non-matching dimensions or invalid values in 'ppm'.")
  }

  idx <- get_idx(shift, ppm)
  specs <- X[, idx, drop = FALSE]
  colnames(specs) <- paste0("ppm_", idx)

  if (is.null(names(an))) {
    names(an) <- paste0("an", seq_along(an))
  }
  names(an) <- make.names(names(an), unique = TRUE)

  df <- data.frame(do.call(cbind.data.frame, an), ID = seq_len(nrow(specs)), alp, specs)
  colnames(df)[seq_along(an)] <- names(an)
  df <- reshape2::melt(df, id.vars = c("alp", "ID", names(an)))
  df$variable <- ppm[as.numeric(gsub("ppm_", "", df$variable))]

  g <- ggplot2::ggplot() +
    ggplot2::scale_x_reverse(name = expression(delta ^ 1 * H ~ "(ppm)"),
                             breaks = seq(min(shift), max(shift), length.out = 5))

  # Add traces based on number of grouping variables
  switch(as.character(length(an)),
         "1" = {
           g <- g + ggplot2::geom_line(data = df, aes(x = !!sym("variable"), y = !!sym("value"), group = !!sym("ID")),
                                       alpha = alp, linewidth = size)
         },
         "2" = {
           col_var <- names(an)[2]
           g <- g + ggplot2::geom_line(data = df, aes(x = !!sym("variable"), y = !!sym("value"), group = !!sym("ID"), colour = !!sym(col_var)),
                                       alpha = alp, linewidth = size)
           if (!is.factor(an[[2]]) && !is.character(an[[2]])) {
             g <- g + ggplot2::scale_colour_gradientn(colors = colorRamps::matlab.like2(100))
           }
         },
         "3" = {
           df[[names(an)[3]]] <- factor(df[[names(an)[3]]])
           g <- g + ggplot2::geom_line(data = df, aes(x = !!sym("variable"), y = !!sym("value"), group = !!sym("ID"),
                                                             colour = !!sym(names(an)[2]), linetype = !!sym(names(an)[3])),
                                       alpha = alp, linewidth = size)
           if (!is.factor(an[[2]]) && !is.character(an[[2]])) {
             g <- g + ggplot2::scale_colour_gradientn(colors = colorRamps::matlab.like2(100))
           }
         })

  g <- g +
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(), name = "Intensity") +
    ggplot2::ggtitle(title) +
    ggplot2::facet_grid(stats::as.formula(paste0(names(an)[1], " ~ .")), ...) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(g)
}


#' @title Overlay Loadings with NMR Spectra
#'
#' @description
#' This function overlays loadings from a PCA or OPLS model on top of NMR spectra using `ggplot2`.
#' It supports two types of loadings visualization: statistical reconstruction or backscaling.
#'
#' @param mod A model object of class `PCA_metabom8` or `OPLS_metabom8`.
#' @param shift Numeric vector (length 2), chemical shift region to visualize (in ppm).
#' @param pc Integer. Component to visualize (for OPLS, set to 1).
#' @param type Character. Either `"Statistical reconstruction"` or `"Backscaled"` (case-insensitive).
#' @param an List of up to three grouping variables (for facet, color, linetype). First must be defined.
#' @param alp Numeric (0–1). Alpha level for spectra.
#' @param title Character. Plot title.
#' @param size Numeric. Line width for plotted spectra.
#' @param r_scale Logical. If `TRUE`, correlation color scale is fixed to [0, 1]. Only used in statistical reconstruction.
#'
#' @details
#' - **Statistical reconstruction**: Correlates predictive scores with spectral variables.
#' - **Backscaled**: Multiplies model loadings by feature SDs.
#'
#' @return A `ggplot` object.
#'
#' @seealso \code{\link{specOverlay}}, \code{\link{pca}}, \code{\link{opls}}
#'
#' @examples
#' data(covid)
#' model <- pca(X)
#' specload(model, shift = c(1, 2), an = list(meta$Group), pc = 1, alp = 0.8)
#'
#' @importFrom ggplot2 ggplot geom_line aes_string scale_x_reverse scale_y_continuous
#' @importFrom ggplot2 facet_grid ggtitle xlab ylab theme_bw theme element_text
#' @importFrom colorRamps matlab.like2
#' @importFrom reshape2 melt
#' @importFrom stats as.formula
#' @export
specload <- function(mod, shift = c(0, 10), an, alp = 0.3, size = 0.5, pc = 1,
                     type = "Backscaled", title = "", r_scale = FALSE) {

  if (!inherits(mod, c("OPLS_metabom8", "PCA_metabom8"))) {
    stop("Requires 'PCA_metabom8' or 'OPLS_metabom8' object.")
  }

  X <- mod@X
  ppm <- as.numeric(colnames(X))
  idx <- get_idx(shift, ppm)
  if (length(idx) < 3) stop("Insufficient shift region selected.")

  type <- if (grepl("st|recon", type, ignore.case = TRUE)) {
    "Statistical reconstruction"
  } else {
    "Backscaled"
  }

  if (type == "Statistical reconstruction") {
    df_l <- .load_stat_reconstr_nmr(mod, pc, X, idx, ppm)
    y <- df_l$cov
    cols <- df_l$cor
    raCol <- if (r_scale) c(0, 1) else c(0, max(abs(cols)))
  } else {
    df_l <- .load_backscaled_nmr(mod, pc, idx, ppm)
    y <- df_l$p_bs
    cols <- abs(df_l$p_abs)
  }

  specs <- X[, idx, drop = FALSE]
  limY <- range(specs)
  colnames(specs) <- paste0("ppm_", ppm[idx])

  an <- .check_an_viz(an, mod)
  df <- data.frame(ID = seq_len(nrow(specs)), do.call(cbind.data.frame, an), specs)
  df <- melt(df, id.vars = c("ID", names(an)))
  df$variable <- as.numeric(gsub("^\\.", "-", gsub("ppm_", "", df$variable)))
  colnames(df)[match(names(an)[1], colnames(df))] <- "facet"

  # Convert y (loading) to display range within intensity
  cv1 <- (minmax(y) * (limY[2]/3)) + limY[2] * 0.67
  if (max(cv1) > limY[2]) {
    cv1 <- cv1 - abs(max(cv1 - limY[2]))
  }

  # Create loadings dataframe
  fac_lev <- unique(df$facet)
  df1 <- do.call(rbind, lapply(seq_along(fac_lev), function(i) {
    data.frame(alp, ID = nrow(X) + i, facet = fac_lev[i], Group = "load",
               ppm = ppm[idx], Intensity = cv1, load = cols)
  }))

  g <- ggplot() +
    geom_line(data = df1, aes(x=!!sym("ppm"), y=!!sym("Intensity"), color = !!sym("load"), group = !!sym("ID")),
              linewidth = 0.8) +
    geom_line(data = df, aes(x=!!sym("variable"), y=!!sym("value"), group = !!sym("ID")),
              alpha = alp, linewidth = size) +
    facet_grid(facet ~ .) +
    scale_x_reverse(breaks = round(seq(shift[1], shift[2], length.out = 10), 3)) +
    scale_y_continuous(limits = limY) +
    ggtitle(title) +
    xlab(expression(delta ^ 1 * H ~ "(ppm)")) +
    ylab("Intensity (AU)") +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  g <- g + scale_colour_gradientn(colors = matlab.like2(10),
                                  name = if (type == "Statistical reconstruction") {
                                    "cor(t, x)"
                                  } else {
                                    expression(abs * w[pred ~ "," ~ sc])
                                  },
                                  limits = if (type == "Statistical reconstruction") raCol else NULL)

  return(g)
}


#' @title Visualize PCA or OPLS Loadings for NMR Data
#'
#' @description
#' Overlay PCA or OPLS loadings on the ppm axis for NMR data. Visualizations can be based on statistical reconstruction
#' or backscaling, providing insights into spectral regions that contribute most to variation or separation.
#'
#' @param mod A PCA or OPLS model object generated via the \emph{metabom8} package.
#' @param shift Numeric vector of length 2. Chemical shift (ppm) range to display.
#' @param pc Integer. Index of component to visualize (use 1 for OPLS).
#' @param type Character. Either `"Statistical reconstruction"` or `"Backscaled"` (case-insensitive).
#' @param title Optional plot title.
#' @param r_scale Logical. If `TRUE`, correlation color gradient is fixed between 0 and 1 (only for statistical reconstruction).
#'
#' @details
#' For OPLS:
#' - **Statistical reconstruction** visualizes correlation (`r`) and covariance between scores and predictors.
#' - **Backscaled** multiplies loadings by feature standard deviations, highlighting consistent spectral influence.
#'
#' For PCA, only statistical reconstruction is available.
#'
#' @return A `ggplot2` object.
#'
#' @references
#' Cloarec, O., et al. (2005). *Anal. Chem.* 77(2), 517–526.
#'
#' @seealso \code{\link{specload}}, \code{\link{pca}}, \code{\link{opls}}
#'
#' @examples
#' data(covid)
#' model <- pca(X)
#' plotload(model, pc = 1)
#'
#' @importFrom ggplot2 ggplot geom_line aes_string scale_x_reverse
#'   ggtitle xlab ylab theme_bw scale_colour_gradientn labs
#' @importFrom colorRamps matlab.like2
#' @importFrom scales breaks_pretty
#' @export
plotload <- function(mod, shift = c(0, 10), pc = 1,
                     type = "Backscaled", title = NULL, r_scale = FALSE) {

  if (!inherits(mod, c("OPLS_metabom8", "PLS_metabom8", "PCA_metabom8"))) {
    stop("Input must be a metabom8 PCA, PLS or OPLS model.")
  }

  type <- if (grepl("st|recon", type, ignore.case = TRUE)) {
    "Statistical reconstruction"
  } else {
    "Backscaled"
  }

  X <- mod@X
  ppm <- as.numeric(colnames(X))
  idx <- get_idx(shift, ppm)

  if (length(idx) < 3) stop("Insufficient shift range selected.")

  if (type == "Statistical reconstruction") {
    df <- .load_stat_reconstr_nmr(mod, pc, X, idx, ppm)
    raCol <- if (r_scale) c(0, 1) else c(0, max(abs(df$cor)))

    g <- ggplot(df, aes(x=!!sym("ppm"), y=!!sym("cov"), colour = !!sym("cor"))) +
      geom_line() +
      scale_x_reverse(breaks = breaks_pretty(n = 15)) +
      scale_colour_gradientn(colors = matlab.like2(10), limits = raCol, name = "r") +
      labs(
        title = title,
        x = expression(delta ^ 1 * H ~ "(ppm)"),
        y = "cov(t, x)",
        caption = paste(gsub("_metabom8", "", class(mod)[1]), "-", mod@type, "component", pc)
      ) +
      theme_bw()
  } else {
    df <- .load_backscaled_nmr(mod, pc, idx, ppm)

    g <- ggplot(df, aes(x=!!sym("ppm"), y=!!sym("p_bs"), colour = !!sym("p_abs"))) +
      geom_line() +
      scale_x_reverse(breaks = breaks_pretty(n = 15)) +
      scale_colour_gradientn(colors = matlab.like2(10), name = expression("|p[sc]|")) +
      labs(
        title = title,
        x = expression(delta ^ 1 * H ~ "(ppm)"),
        y = expression(p * "*" * sigma[x]),
        caption = paste(gsub("_metabom8", "", class(mod)[1]), "-", mod@type, "component", pc)
      ) +
      theme_bw()
  }

  return(g)
}


#' @title Distance to the Model in X-Space (DModX)
#'
#' @description
#' Calculates the orthogonal distance of each observation to the OPLS model in X-space. The DModX is used for identifying outliers.
#'
#' @param mod An OPLS model object of class \code{OPLS_metabom8}.
#' @param plot Logical. If \code{TRUE}, a plot of DModX values with a cutoff line based on a t-test is shown.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{ID}{Sample index}
#'   \item{DmodX}{Distance to the model in X-space}
#'   \item{passedT.test}{Logical. TRUE if within the 95\\% confidence interval}
#' }
#'
#' @details
#' DModX is calculated as the scaled root-mean-squared residual. An approximate upper 95\\% confidence limit is drawn using a one-sided t-test.
#'
#' @references
#' Bylesjö, M. et al. (2006). *J. Chemometrics*, 20, 341–351.
#' Wold, S. (1976). *Pattern Recognition*, 8, 127–139.
#'
#' @seealso \code{\link{opls}}, \code{\link{plotload}}
#'
#' @importFrom stats t.test sd
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_hline
#'   scale_y_continuous xlab theme_bw theme element_blank element_text
#'   scale_colour_gradientn labs
#' @importFrom colorRamps matlab.like2
#'
#' @examples
#' data(covid)
#' model <- opls(X, Y = an$type)
#' dmx <- dmodx(model)
#'
#' @export
dmodx <- function(mod, plot = TRUE) {
  if (!inherits(mod, "OPLS_metabom8")) {
    stop("Please provide an OPLS_metabom8 object.")
  }

  E <- mod@X_res
  N <- nrow(E)
  K <- ncol(E)
  A <- ncol(mod@t_pred)
  A0 <- if (isTRUE(mod@Parameters$center)) 1 else 0

  ss_res <- rowSums(E^2)
  dmodX <- sqrt(ss_res / (K - A)) / sqrt(sum(ss_res) / ((N - A - A0) * (K - A)))

  # 95% CI via one-sided t-test
  tt <- t.test(dmodX, alternative = "less")
  ci95 <- tt$conf.int[2] + 2 * sd(dmodX)

  df <- data.frame(
    ID = seq_len(N),
    DmodX = dmodX,
    passedT.test = dmodX < ci95
  )

  if (plot) {
    df$Group <- mod@Y$ori

    g <- ggplot(df, aes(x = ID, y = DmodX, colour = Group)) +
      geom_segment(aes(xend = ID, yend = min(dmodX) - 0.1), colour = "grey60", linewidth = 0.1) +
      geom_point() +
      geom_hline(yintercept = ci95, linetype = 2, colour = "black") +
      scale_y_continuous(
        limits = c(min(dmodX) - 0.1, max(c(dmodX, ci95)) + 0.2),
        expand = c(0, 0),
        name = "DModX"
      ) +
      xlab("Sample index") +
      theme_bw() +
      theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(colour = "black")
      ) +
      labs(caption = "Dashed line indicates upper limit of 95% CI")

    if (!is.null(mod@type) && mod@type == "R") {
      g <- g + scale_colour_gradientn(colours = matlab.like2(10), name = expression(t[pred]))
    }

    plot(g)
  }

  return(df)
}

