#' @title Distance to the Model in X-Space (DModX)
#'
#' @description
#' Calculates the orthogonal distance of each observation to an OPLS model in X-space.
#' DModX can be used for identifying outliers.
#'
#' @param mod An OPLS model object of class \code{m8_model} (engine = "opls").
#' @param plot Logical. If TRUE, a plot of DModX values with an approximate cutoff is shown.
#'
#' @return A data frame with columns \code{ID}, \code{DmodX}, and \code{passed}.
#'
#' @details
#' DModX is computed from the X-residual matrix as a scaled RMSE of residuals.
#' The empirical cutoff (dashed line in plot) uses a t-based confidence interval
#' and assumes approximate normality of DModX values. This assumption may not be
#' satisfied in all datasets, so the resulting threshold should be regarded
#' as a pragmatic heuristic for outlier detection.
#'
#' @family model_validation
#' @importFrom stats t.test sd
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_hline
#'   scale_y_continuous xlab theme_bw theme element_blank element_text
#'   scale_colour_gradientn labs
#' @importFrom colorRamps matlab.like2
#' @importFrom rlang .data
#' @examples
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- uv_scaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' dX <- dmodx(model)
#' print(dX[[1]])
#' df <- dX[[2]]; head(df)
#'
#' @export
dmodx <- function(mod, plot = TRUE) {

  if (!inherits(mod, "m8_model") || !identical(mod@engine, "opls")) {
    stop("Please provide an m8_model object with engine = 'opls'.", call. = FALSE)
  }

  E <- xres(mod)

  if (is.null(E)) {
    stop("Could not find X residual matrix.", call. = FALSE)
  }

  N <- nrow(E)
  K <- ncol(E)
  A <- 1L

  A0 <- if (isTRUE(mod@prep$center)) 1L else 0L

  if ((K - A) <= 0) stop("Invalid dimensions for DModX: need K - A > 0.", call. = FALSE)
  if ((N - A - A0) <= 0) stop("Invalid dimensions for DModX: need N - A - A0 > 0.", call. = FALSE)

  ss_res <- rowSums(E^2)
  denom_global <- sqrt(sum(ss_res) / ((N - A - A0) * (K - A)))
  dmodX <- sqrt(ss_res / (K - A)) / denom_global

  tt <- t.test(dmodX, alternative = "less")
  ci95 <- tt$conf.int[2] + 2 * sd(dmodX)

  df <- data.frame(
    ID     = seq_len(N),
    DmodX  = dmodX,
    passed = dmodX < ci95
  )

  if (isTRUE(plot)) {

    df$outlier <- df$DmodX > ci95
    df$label <- as.character(df$ID)

    grp <- NULL
    if (!is.null(mod@fit$Y)) {
      grp <- mod@fit$Y
      if (mod@ctrl$type=='DA'){
        grp <- as.factor(grp)
      }
      }
    if (!is.null(grp) && length(grp) == nrow(df)) {
      df$Group <- grp
    } else {
      df$Group <- factor("all")
    }


    g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$ID, y = .data$DmodX)) +
      ggplot2::geom_segment(
        ggplot2::aes(xend = .data$ID, yend = 0),
        colour = "grey85", linewidth = 0.2
      ) +
      ggplot2::geom_hline(yintercept = ci95, linetype = 2, colour = "black") +
      ggplot2::scale_y_continuous(
        limits = c(0, max(df$DmodX, ci95) * 1.05),
        expand = c(0, 0),
        name = "DModX"
      ) +
      ggplot2::xlab("Sample index") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(colour = "black")
      ) +
      ggplot2::scale_x_continuous(
        breaks = scales::breaks_pretty(n = 10)
      ) +
      ggplot2::labs(caption = "Dashed line indicates one-sided empirical cutoff (t-based, 95% limit).")

    g <- g +
      ggplot2::geom_point(
        data = df[!df$outlier, , drop = FALSE],
        ggplot2::aes(colour = .data$Group),
        size = 1.7, alpha = 0.9
      )

    g <- g +
      ggplot2::geom_point(
        data = df[df$outlier, , drop = FALSE],
        colour = "red3",
        size = 3.0, alpha = 1
      )

    if (any(df$outlier)) {
      g <- g +
        ggrepel::geom_text_repel(
          data = df[df$outlier, , drop = FALSE],
          ggplot2::aes(label = .data$label),
          colour = "red3",
          size = 3,
          min.segment.length = 0,
          box.padding = 0.25,
          point.padding = 0.2,
          max.overlaps = Inf
        )
    }

    if (!is.null(mod@ctrl$type) && identical(mod@ctrl$type, "R") && !is.factor(df$Group)) {
      g <- g + scale_colour_gradientn(colours = matlab.like2(10), name = "Y")
    }

    plot(g)
  }

  df
}
