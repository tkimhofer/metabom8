#' @title Principal Component Analysis (PCA)
#'
#' @description
#' Performs Principal Component Analysis on metabolomics datasets. Supports NIPALS via Rcpp or alternative PCA algorithms from the \pkg{pcaMethods} package.
#'
#' @param X A numeric matrix or \code{data.frame}, where rows represent observations and columns represent metabolic features.
#' @param pc Integer. The number of principal components to compute. If \code{pc} exceeds the allowed maximum, it will be set to the smaller of \code{nrow(X)} or \code{ncol(X)}.
#' @param scale Character. Scaling method to apply before PCA. One of \code{"None"}, \code{"UV"} (unit variance), or \code{"Pareto"}.
#' @param center Logical. Whether to mean-center the variables. Defaults to \code{TRUE}.
#' @param method Character. Algorithm to use for PCA. Defaults to \code{"nipals"} (Rcpp implementation). Alternatives include: \code{"svd"}, \code{"ppca"}, \code{"svdImpute"}, \code{"robustPca"}, and \code{"nlpca"}, all from the \pkg{pcaMethods} package.
#'
#' @details
#' The default \code{"nipals"} method uses an efficient NIPALS implementation in C++ via Rcpp for large datasets. If another method is specified, the function delegates to \code{pca()} from the \pkg{pcaMethods} package. Scaling and centering are performed manually before PCA, regardless of method.
#'
#' @return
#' An S4 object of class \code{"PCA_metabom8"} containing:
#' \itemize{
#'   \item \code{t}: Principal component scores (observations × components).
#'   \item \code{p}: Loadings matrix.
#'   \item \code{X_mean}: Feature means (for re-scaling).
#'   \item \code{X_sd}: Feature standard deviations (for re-scaling).
#'   \item \code{Parameters}: A list of metadata including centering/scaling choices, explained variance, and number of components.
#'   \item \code{X}: The original input matrix.
#' }
#'
#' @references
#' Geladi, P., & Kowalski, B. R. (1986). Partial least squares and regression: a tutorial. \emph{Analytica Chimica Acta}, \strong{185}, 1–17.
#'
#' @examples
#' data(covid)
#' X <- covid$X
#' an <- covid$an
#'
#' model <- pca(X, pc = 2)
#' plotscores(model, an = list(Class = an$type, Clinic = an$hospital, id = 1:nrow(an)), pc = c(1, 2))
#'
#' @seealso \code{\link[pcaMethods]{pca}}, \code{\link{plotscores}}
#'
#' @importFrom pcaMethods pca
#' @export
pca <- function(X, pc = 2, scale = "UV", center = TRUE, method = "nipals") {

  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.logical(center)) stop("`center` must be logical (TRUE or FALSE).")
  if (!is.character(scale) || !(scale %in% c("None", "UV", "Pareto")))
    stop("`scale` must be one of 'None', 'UV', or 'Pareto'.")
  if (!is.character(method)) stop("`method` must be a character string.")

  pc_max <- min(ncol(X), nrow(X))
  if (pc > pc_max) {
    warning("Too many components requested; resetting to maximum allowed (", pc_max, ").")
    pc <- pc_max
  }

  sc_num <- switch(scale,
                   "None" = 0,
                   "UV" = 1,
                   "Pareto" = 2)

  # Internal checks and preprocessing
  x_check <- .checkXclassNas(X)
  msd_x <- .sdRcpp(X)
  XcsTot <- .scaleMatRcpp(X, 0:(nrow(X) - 1), center = TRUE, scale_type = sc_num)[[1]]

  tssx <- .tssRcpp(XcsTot)

  if (tolower(method) == "nipals") {
    res <- vector("list", pc)
    for (i in seq_len(pc)) {
      res[[i]] <- if (i == 1) .nipPcaCompRcpp(XcsTot) else .nipPcaCompRcpp(X = res[[i - 1]][[1]])
    }
    Tpc <- vapply(res, "[[", FUN.VALUE = X[, 1], 2)
    Ppc <- vapply(res, "[[", FUN.VALUE = X[1, ], 3)

    ss_comp <- vapply(seq_len(ncol(Tpc)), function(i) {
      .tssRcpp(Tpc[, i] %o% Ppc[, i]) / tssx
    }, numeric(1))

    pars <- list(center = center, scale = scale, nPC = pc, R2 = ss_comp, tssx = tssx)
    mod_pca <- new("PCA_metabom8",
                   type = "NIPALS",
                   t = Tpc,
                   p = Ppc,
                   nPC = pc,
                   X_mean = msd_x[[1]],
                   X_sd = msd_x[[2]],
                   Parameters = pars,
                   X = X)

  } else {
    mod <- pcaMethods::pca(X, nPcs = pc, scale = "none", center = FALSE, method = method)
    r2 <- mod@R2cum
    r2[-1] <- diff(mod@R2cum)

    pars <- list(center = center, scale = scale, nPC = pc, R2 = r2, tssx = tssx)
    mod_pca <- new("PCA_metabom8",
                   type = method,
                   t = mod@scores,
                   p = mod@loadings,
                   nPC = pc,
                   X_mean = msd_x[[1]],
                   X_sd = msd_x[[2]],
                   Parameters = pars,
                   X = X)
  }

  return(mod_pca)
}
