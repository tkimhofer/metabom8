#' Statistical Total Correlation Spectroscopy (STOCSY)
#'
#' Performs STOCSY analysis on 1D NMR spectra. Computes correlation and covariance
#' of all variables in \code{X} against a driver signal (internal ppm position or
#' an external numeric vector).
#'
#' @param X Numeric matrix or data frame. Rows are spectra, columns are variables.
#' @param ppm Optional numeric vector of chemical shift values (length = ncol(X)).
#'   If missing or NULL, will try to read from \code{colnames(X)}.
#' @param driver Numeric scalar (ppm value, internal driver) or numeric vector
#'   (external driver, length = nrow(X)).
#' @param plotting Logical. If TRUE, plot the STOCSY spectrum.
#' @param title Optional character plot title.
#'
#' @return An object of class \code{m8_stocsy1d} (S3), a named list with entries:
#'   \code{version, X, ppm, driver, r, cov, extD}.
#'
#' @examples
#' # st <- stocsy(X, ppm, driver = 5.233, plotting = FALSE)
#' # plotStocsy(st, shift = c(5.15, 5.30))
#'
#' @family structural_annotation
#' @importFrom ggplot2 ggplot geom_line geom_vline scale_x_reverse
#'   scale_colour_gradientn labs theme_bw theme element_text
#' @importFrom colorRamps matlab.like2
#' @importFrom scales breaks_pretty
#' @examples
#' data(covid)
#' cs = 5.233
#' stocsy(covid$X, driver=cs, plotting = TRUE)
#' @export
stocsy <- function(X, ppm = NULL, driver, plotting = TRUE, title = NULL) {

  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a numeric matrix or data frame.", call. = FALSE)
  }
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.", call. = FALSE)

  if (missing(ppm) || is.null(ppm)) {
    ppm_try <- as.numeric(colnames(X))
    if (is.null(ppm_try) || all(is.na(ppm_try))) {
      stop("Provide 'ppm' or ensure colnames(X) are numeric ppm values.", call. = FALSE)
    }
    ppm <- ppm_try
  }
  if (!is.numeric(ppm) || length(ppm) != ncol(X)) {
    stop("'ppm' must be numeric and have length ncol(X).", call. = FALSE)
  }

  if (!is.numeric(driver)) stop("'driver' must be numeric.", call. = FALSE)
  extD <- length(driver) > 1L

  if (extD) {
    if (length(driver) != nrow(X)) stop("External driver must have length nrow(X).", call. = FALSE)
    idx <- which(!is.na(driver))
    if (length(idx) < 2L) stop("External driver has <2 non-NA entries.", call. = FALSE)
    driver_use <- driver[idx]
    X_use <- X[idx, , drop = FALSE]
  } else {
    if (driver < min(ppm) || driver > max(ppm)) {
      stop("Internal STOCSY driver must lie within the ppm range.", call. = FALSE)
    }
    driver_use <- driver
    X_use <- X
  }

  l <- if (extD) {
    as.numeric(driver_use)
  } else {
    X_use[, which.min(abs(ppm - driver_use))]
  }

  cc <- apply(X_use, 2, function(x) cor(x, l, use = "pairwise.complete.obs"))
  cv <- apply(X_use, 2, function(x) cov(x, l, use = "pairwise.complete.obs"))

  stoc_mod <- list(
    # version = "0.9",
    X       = X_use,
    ppm     = as.numeric(ppm),
    driver  = driver_use,
    r       = as.numeric(cc),
    cov     = as.numeric(cv),
    extD    = extD
  )

  if (isTRUE(plotting)) {
    print(plotStocsy(stoc_mod, shift = range(ppm, na.rm = TRUE), title = title))
  }

  stoc_mod
}


#' Plot STOCSY result
#'
#' @description
#' Generates a STOCSY plot (covariance trace coloured by absolute correlation).
#'
#' @param stoc_mod An object of class \code{m8_stocsy1d} returned by \code{stocsy()}.
#' @param shift Numeric vector of length 2 specifying the chemical shift range (ppm).
#' @param title Optional character plot title.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' # st <- stocsy(X, ppm, driver = 5.233, plotting = FALSE)
#' # plotStocsy(st, shift = c(5.15, 5.30), title = "Glucose")
#'
#' @importFrom ggplot2 ggplot geom_line geom_vline scale_x_reverse
#'   scale_colour_gradientn labs theme_bw theme element_text ggtitle
#' @importFrom colorRamps matlab.like2
#' @importFrom scales breaks_pretty
#' @importFrom rlang .data
#' @examples
#' data(covid)
#' cs = 5.233 # anomeric H of gluc
#' s1 <- stocsy(covid$X, driver=cs, plotting = FALSE)
#' plotStocsy(s1)
#'
#' @export
plotStocsy <- function(stoc_mod, shift = c(0, 10), title = NULL) {

  if (!is.list(stoc_mod)) {
    stop("Input must be a list instance returned by `stocsy()`.")
  }
  if (!is.numeric(shift) || length(shift) != 2L) {
    stop("'shift' must be a numeric vector of length 2.")
  }

  ds <- data.frame(
    r   = stoc_mod$r,
    cov = stoc_mod$cov,
    ppm = stoc_mod$ppm
  )

  # Clamp shift range to available ppm range
  shift[1] <- max(shift[1], min(ds$ppm, na.rm = TRUE))
  shift[2] <- min(shift[2], max(ds$ppm, na.rm = TRUE))
  shift <- sort(shift)

  idx <- which(ds$ppm >= shift[1] & ds$ppm <= shift[2])
  if (length(idx) == 0L) stop("No ppm values in specified shift range.")
  ds <- ds[idx, , drop = FALSE]

  extD <- isTRUE(stoc_mod$extD)
  csc_lab <- if (extD) {
    "r(ext, X)"
  } else {
    paste0("r(d=", signif(stoc_mod$driver, 6), ", X)")
  }

  caption <- if (extD) {
    extD_stats <- round(summary(stoc_mod$driver), 2)
    paste0("Sample size: n=", length(stoc_mod$driver),
           "\nExternal driver: Median=", extD_stats[3],
           " Range=", extD_stats[1], "-", extD_stats[6])
  } else {
    paste0("Sample size: n=", nrow(stoc_mod$X))
  }

  ds$r_abs <- abs(ds$r)

  g1 <- ggplot(ds, aes(x = .data$ppm, y = .data$cov, colour = .data$r_abs)) +
    geom_line() +
    scale_x_reverse(breaks = breaks_pretty()) +
    scale_colour_gradientn(colours = matlab.like2(10), limits = c(0, 1), name = csc_lab) +
    labs(x = expression(delta^1 * H ~ "(ppm)"),
         y = gsub("^r", "cov", csc_lab),
         title = title, caption = caption) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black")) +
    ggtitle(title)

  # Add internal driver vertical line if it lies within the shown window
  if (!extD) {
    d <- as.numeric(stoc_mod$driver)
    if (shift[1] <= d && d <= shift[2]) {
      g1 <- g1 + geom_vline(xintercept = d, linetype = 2, col = "grey")
    }
  }

  g1
}


#' @export
plot.m8_stocsy1d <- function(x, shift = c(0, 10), title = NULL, ...) {
  plotStocsy(x, shift = shift, title = title)
}

#' @export
print.m8_stocsy1d <- function(x, ...) {
  cat("m8_stocsy1d <", if (isTRUE(x$extD)) "external" else "internal", " driver>\n", sep = "")
  cat("Samples   :", nrow(x$X), "\n")
  cat("Variables :", ncol(x$X), "\n")
  invisible(x)
}
