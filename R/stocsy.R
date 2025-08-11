#' @title Statistical Total Correlation Spectroscopy (STOCSY)
#'
#' @description
#' Performs STOCSY analysis on 1D NMR spectra. The function calculates correlation and covariance of all variables in `X` against a given driver variable (either a ppm value or external vector).
#'
#' @param X Numeric matrix or data frame where rows are spectra and columns are chemical shift variables.
#' @param ppm Numeric vector of chemical shift values corresponding to columns in `X`.
#' @param driver Numeric scalar or vector. A ppm value (for internal STOCSY) or external numeric vector of same length as number of rows in `X`.
#' @param plotting Logical. If TRUE, display a STOCSY plot.
#' @param title Character string for the plot title.
#'
#' @return An S4 object of class `stocsy1d_metabom8` containing correlation, covariance, driver, and metadata.
#'
#' @examples
#' data(covid)
#' X <- covid$X
#' ppm <- covid$ppm
#'
#' stcy_glucose <- stocsy(X, ppm, driver = 5.233)
#' plotStocsy(stcy_glucose, shift = c(5.15, 5.30), title = "Alpha-anomeric proton of glucose (doublet)")
#'
#' @importFrom ggplot2 ggplot aes_string geom_line geom_vline scale_x_reverse scale_colour_gradientn labs theme_bw theme element_text
#' @importFrom colorRamps matlab.like2
#' @importFrom scales breaks_pretty
#' @importFrom methods hasArg is
#' @export
stocsy <- function(X, ppm, driver, plotting = TRUE, title = NULL) {
  if (!hasArg(ppm)) {
    ppm <- as.numeric(colnames(X))
    if (!any(!is.na(ppm))) stop("Provide ppm argument")
  }

  extD <- FALSE
  if (length(driver) > 1) {
    if (length(driver) != nrow(X)) stop("External driver does not match to X.")
    idx <- which(!is.na(driver))
    driver <- driver[idx]
    X <- X[idx, , drop = FALSE]
    extD <- TRUE
  }

  if (!is.matrix(X) && !is.data.frame(X)) stop("X must be a numeric matrix or data frame")
  if (!is.numeric(driver) || (!extD && (driver < min(ppm) || driver > max(ppm)))) {
    stop("STOCSY driver must be numeric and within ppm range")
  }

  # Determine driver signal
  l <- if (extD) driver else X[, which.min(abs(ppm - driver))]

  # Compute correlation and covariance
  cc <- apply(X, 2, function(x) cor(x, l))
  cv <- apply(X, 2, function(x) cov(x, l))

  stoc_mod <- new("stocsy1d_metabom8", version = "0.9", X = X, ppm = ppm,
                  driver = driver, r = cc, cov = cv)

  if (plotting) {
    df <- data.frame(cc = cc, cv = cv, ppm = ppm)
    csc_lab <- if (extD) "r(ext, X)" else paste("r(d=", driver, ", X)", sep = "")
    caption <- if (extD) {
      extD_stats <- round(summary(driver), 2)
      paste0("Sample size: n=", length(driver), "\nExternal driver: Median=",
             extD_stats[3], " Range=", extD_stats[1], "-", extD_stats[6])
    } else {
      paste0("Sample size: n=", nrow(X))
    }

    g1 <- ggplot(df, aes_string(x = "ppm", y = "cv", colour = "abs(cc)")) +
      geom_line() +
      scale_x_reverse(breaks = breaks_pretty()) +
      scale_colour_gradientn(colours = matlab.like2(10), limits = c(0, 1), name = csc_lab) +
      labs(x = expression(delta^1 * H ~ "(ppm)"),
           y = gsub("^r", "cov", csc_lab),
           title = title, caption = caption) +
      theme_bw() +
      theme(axis.text = element_text(colour = "black"),
            panel.grid.minor.x = element_blank())

    if (!extD) {
      g1 <- g1 + geom_vline(xintercept = driver, linetype = 2, col = "grey")
    }

    plot(g1)
  }

  return(stoc_mod)
}


#' @title Plot Statistical Total Correlation Spectroscopy (STOCSY)
#'
#' @description
#' Generates a plot of the results from a STOCSY analysis.
#'
#' @param stoc_mod An object of class `stocsy1d_metabom8` returned by the `stocsy()` function.
#' @param shift Numeric vector of length 2 specifying the chemical shift range (ppm) to be plotted.
#' @param title Character string for the plot title.
#'
#' @return A `ggplot2` object displaying the STOCSY spectrum within the specified shift range.
#'
#' @examples
#' data(covid)
#' X <- covid$X
#' ppm <- covid$ppm
#' sy_glucose <- stocsy(X, ppm, driver = 5.233)
#' plotStocsy(sy_glucose, shift = c(5.15, 5.30), title = "d = 5.233")
#'
#' @importFrom ggplot2 ggplot aes_string geom_line geom_vline scale_x_reverse scale_colour_gradientn labs theme_bw theme element_text ggtitle
#' @importFrom colorRamps matlab.like2
#' @importFrom scales breaks_pretty
#' @export
plotStocsy <- function(stoc_mod, shift = c(0, 10), title = NULL) {
  if (!inherits(stoc_mod, "stocsy1d_metabom8")) stop("Input must be a 'stocsy1d_metabom8' object")
  if (!is.numeric(shift) || length(shift) != 2) stop("'shift' must be a numeric vector of length 2")

  ds <- data.frame(r = stoc_mod@r, cov = stoc_mod@cov, ppm = stoc_mod@ppm)
  extD <- length(stoc_mod@driver) > 1

  # Clamp shift range to available ppm range
  shift[1] <- max(shift[1], min(ds$ppm))
  shift[2] <- min(shift[2], max(ds$ppm))
  shift <- sort(shift)

  idx <- which(ds$ppm >= shift[1] & ds$ppm <= shift[2])
  if (length(idx) == 0) stop("No ppm values in specified shift range")
  ds <- ds[idx, ]

  csc_lab <- if (extD) "r(ext, X)" else paste("r(d=", stoc_mod@driver, ", X)", sep = "")
  caption <- if (extD) {
    extD_stats <- round(summary(stoc_mod@driver), 2)
    paste0("Sample size: n=", length(stoc_mod@driver),
           "\nExternal driver: Median=", extD_stats[3],
           " Range=", extD_stats[1], "-", extD_stats[6])
  } else {
    paste0("Sample size: n=", nrow(stoc_mod@X))
  }

  ds$r_abs <- abs(ds$r)

  g1 <- ggplot(ds, aes(x = !!sym("ppm"), y = !!sym("cov"), colour = !!sym("r_abs"))) +
    geom_line() +
    scale_x_reverse(breaks = breaks_pretty()) +
    scale_colour_gradientn(colours = matlab.like2(10), limits = c(0, 1), name = csc_lab) +
    labs(x = expression(delta^1 * H ~ "(ppm)"),
         y = gsub("^r", "cov", csc_lab),
         title = title, caption = caption) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black")) +
    ggtitle(title)

  if (!extD && shift[1] < stoc_mod@driver && shift[2] > stoc_mod@driver) {
    g1 <- g1 + geom_vline(xintercept = stoc_mod@driver, linetype = 2, col = "grey")
  }

  return(g1)
}
