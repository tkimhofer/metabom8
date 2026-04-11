#' Prepare Spectral Data for Plotting
#'
#' Internal helper that validates input, enforces ppm ordering,
#' subsets the selected shift region, and returns matrix-form data.
#'
#' @param x Numeric vector or matrix of spectra.
#' @param ppm Numeric chemical shift vector.
#' @param shift Length-2 numeric vector specifying ppm window.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{ppm}{Subsetted chemical shift axis}
#'     \item{X}{Matrix of subsetted spectra}
#'   }
#'
#' @keywords internal
.prep_spec <- function(x, ppm, shift) {

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  } else if (!is.matrix(x)) {
    stop("x must be vector or matrix", call. = FALSE)
  }

  if (ncol(x) != length(ppm))
    stop("ppm mismatch", call. = FALSE)

  if (ppm[1] < ppm[length(ppm)]) {
    ord <- order(ppm, decreasing = TRUE)
    ppm <- ppm[ord]
    x   <- x[, ord, drop = FALSE]
  }

  idx <- get_idx(shift, ppm)

  list(
    ppm = ppm[idx],
    X   = x[, idx, drop = FALSE]
  )
}


#' Draw NMR Spectra Using Base Graphics
#'
#' Internal backend renderer for base R plotting.
#'
#' @param dat List returned by `.prep_spec()`.
#' @param add Logical. Add to existing plot.
#' @param ... Additional graphical parameters.
#'
#' @keywords internal
.draw_base <- function(dat, add, ...) {

  if (add) {
    graphics::matlines(dat$ppm, t(dat$X), ...)
  } else {
    graphics::matplot(dat$ppm, t(dat$X),
                      type="l",
                      xlim = rev(range(dat$ppm)),
                      xlab="ppm", ylab="Intensity", ...)
  }

  invisible(NULL)
}


#' Draw NMR Spectra Using Plotly
#'
#' Internal backend renderer for interactive plotting.
#'
#' @param dat List returned by `.prep_spec()`.
#' @param ... Additional arguments passed to `plotly`.
#'
#' @keywords internal
.draw_plotly <- function(dat, col = NULL, ...) {

  X   <- dat$X
  ppm <- dat$ppm
  n   <- nrow(X)

  labels <- rownames(dat$meta)
  if (is.null(labels))
    labels <- paste0("spec_", seq_len(n))
    labels <- factor(labels, levels = labels)

  df <- data.frame(
    ppm = rep(ppm, times = n),
    intensity = as.vector(t(X)),
    spec = rep(seq_len(n), each = length(ppm)),
    label     = rep(labels, each = length(ppm))
  )

  p <- plotly::plot_ly(type = "scattergl", mode = "lines")

  # -------- numeric colouring --------
  if (!is.null(col) && is.numeric(col)) {

    df$z <- rep(col, each = length(ppm))

    p <- plotly::add_trace(
      p,
      data  = df,
      x     = ~ppm,
      y     = ~intensity,
      split = ~label,
      text  = ~label,
      hovertemplate = "%{text}<extra></extra>",
      line  = list(
        color = ~z,
        colorscale = "Viridis",
        showscale  = TRUE,
        width = 1
      ),
      hoverinfo = "none"
    )

    # -------- categorical colouring --------
  } else if (!is.null(col)) {

    df$grp <- rep(as.factor(col), each = length(ppm))

    p <- plotly::add_trace(
      p,
      data  = df,
      x     = ~ppm,
      y     = ~intensity,
      split = ~label,
      color = ~grp,
      text  = ~label,
      hovertemplate = "%{text}<extra></extra>",
      line  = list(width = 1),
      hoverinfo = "none"
    )

    # -------- default --------
  } else {

    p <- plotly::add_trace(
      p,
      data  = df,
      x     = ~ppm,
      y     = ~intensity,
      split = ~label,
      text  = ~label,
      hovertemplate = "%{text}<extra></extra>",
      line  = list(width = 1),
      hoverinfo = "none"
    )
  }

  plotly::layout(
    p,
    xaxis = list(autorange = "reversed"),
    yaxis = list(title = "Intensity")
  )
}

#' Draw NMR Spectra Using ggplot2
#'
#' Internal backend renderer using `ggplot2`. Converts spectra
#' to long format prior to plotting.
#'
#' @param dat List returned by `.prep_spec()`.
#' @param ... Additional arguments for `ggplot2`.
#' @importFrom rlang .data
#' @keywords internal
.draw_ggplot2 <- function(dat, ...) {

  df <- reshape2::melt(dat$X)
  df$ppm <- rep(dat$ppm, each = nrow(dat$X))

  g <- ggplot2::ggplot(df,
                  ggplot2::aes(.data$ppm, .data$value, group = .data$Var1)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_reverse()+
    ggplot2::ylab("Intensity") +
    ggplot2::theme_bw()

  return(g)
}

#' Plot 1D NMR Spectra
#'
#' Plot one or multiple 1D ¹H NMR spectra using different rendering backends.
#'
#' @param x Numeric vector (single spectrum) or numeric matrix (spectra in rows).
#' @param ppm Numeric vector of chemical shift values. Must match `ncol(x)`.
#' @param shift Numeric vector of length 2. Chemical shift window to display
#'   (e.g., `c(0, 10)`).
#' @param backend Character. Rendering backend. One of `"plotly"`, `"base"`,
#'   or `"ggplot"`. Defaults to `"plotly"`.
#' @param add Logical. If `TRUE` and `backend = "base"`, add spectra to an
#'   existing plot.
#' @param ... Additional arguments passed to the selected backend.
#'
#' @details
#' The function accepts both single spectra and matrices of spectra. Input is
#' internally normalized to matrix form, subset to the selected ppm region,
#' and then rendered using the chosen backend.
#'
#' For large NMR datasets (e.g., >500 spectra × >10k ppm), the `"base"` and
#' `"plotly"` backends are substantially more memory-efficient than `"ggplot"`,
#' which requires reshaping to long format.
#'
#' Chemical shift axes are automatically displayed in decreasing order
#' (NMR convention).
#'
#' @return
#' - `"plotly"`: a `plotly` object.
#' - `"ggplot"`: a `ggplot2` object.
#' - `"base"`: `NULL` (invisibly).
#' @examples
#' data(hiit_raw)
#' plot_spec(hiit_raw)
#' plot_spec(hiit_raw, shift=c(-0.05,0.05), backend='base')
#' plot_spec(hiit_raw, shift=c(-0.05,0.05), backend='ggplot2')
#' @export
plot_spec <- function(x,
                      ppm=NULL,
                      shift = c(0,10),
                      backend = c("plotly","base","ggplot2"),
                      add = FALSE,
                      ...) {

  u   <- .m8_unpack_dat(x, ppm = ppm)
  X   <- .dimX(u$X)
  ppm <- u$ppm
  meta <- u$meta

  backend <- match.arg(backend)

  dat <- .prep_spec(X, ppm, shift)
  dat$meta <- meta

  switch(
    backend,
    base   = .draw_base(dat, add, ...),
    plotly = .draw_plotly(dat, ...),
    ggplot2 = .draw_ggplot2(dat, ...)
  )
}


#' @rdname plot_spec
#' @export
spec <- function(...) {
  warning("spec() is deprecated; use plot_spec() instead.", call. = FALSE)
  plot_spec(..., backend = "plotly")
}

#' @rdname plot_spec
#' @export
matspec <- function(...) {
  warning("matspec() is deprecated; use plot_spec() instead.", call.=FALSE)
  plot_spec(..., backend = "plotly")
}


