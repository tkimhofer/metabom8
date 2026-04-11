#' Principal Component Analysis (PCA)
#'
#' @description
#' Fits an unsupervised Principal Component Analysis (PCA) model
#' to a spectral data matrix. Model-internal centering and scaling
#' are controlled via a \code{ScalingStrategy} object and do not
#' modify the input matrix.
#'
#' The default backend uses an internal NIPALS implementation.
#' Alternatively, PCA algorithms from the \pkg{pcaMethods} package
#' (e.g. \code{"ppca"}, \code{"svd"}, \code{"robustPca"}) can be used.
#'
#' @param X Numeric matrix (rows = samples, columns = variables).
#' @param scaling A \code{ScalingStrategy} object defining model-internal
#'   centering and scaling (e.g. \code{uv_scaling(center = TRUE)}).
#' @param ncomp Integer. Number of principal components to compute.
#' @param method Character. PCA backend. Use \code{"nipals"} for the
#'   internal implementation or any method supported by
#'   \code{pcaMethods::pca()}.
#'
#' @details
#' PCA decomposes \eqn{X} into orthogonal score and loading matrices:
#' \deqn{X \approx T P}
#' where:
#' \itemize{
#'   \item \eqn{T} contains the principal component scores
#'   \item \eqn{P} contains the loadings
#' }
#'
#' The number of components is fixed by \code{ncomp}.
#' Unlike supervised models, PCA does not use cross-validation
#' or stopping rules.
#'
#' Scaling and centering are applied internally during model fitting.
#' The original input matrix is not modified.
#'
#' @return
#' An object of class \code{m8_model} with \code{engine = "pca"}.
#' The object contains:
#' \itemize{
#'   \item \code{fit$t}: Score matrix (samples × components)
#'   \item \code{fit$p}: Loading matrix (components × variables)
#'   \item \code{ctrl}: Model control information (variance explained, scaling settings)
#'   \item \code{provenance}: Attributes inherited from the input matrix
#' }
#'
#' @seealso \code{\link{uv_scaling}}
#'
#' @family modelling
#' @examples
#' data(covid)
#'
#' uv <- uv_scaling(center=TRUE)
#' model <- pca(X=covid$X, scaling=uv, ncomp=2)
#'
#' show(model)
#' summary(model)
#'
#' Tx <- scores(model)
#' Px <- loadings(model)
#'
#' t2 <- hotellingsT2(Tx)
#' ell <-ellipse2d(t2)
#'
#' # scores plot
#' plot(Tx, asp = 1,
#'   col = as.factor(covid$an$type),
#'   xlim = range(c(Tx[1,], ell$x)),
#'   ylim = range(c(Tx[2,], ell$y))
#'  )
#' lines(ell$x, ell$y, col = "grey", lty=2)
#' @export
pca <- function(X, scaling, ncomp, method = "nipals") {

  cl <- match.call()

  engine <- .pcaEngine(X, scaling, ncomp, method=method)

  mod <- methods::new("m8_model",

      engine = "pca",
      fit = engine$full,
      prep = scaling,
      cv = NULL,
      ctrl = engine$ctrl,
      dims = engine$dims,

      provenance = attributes(X),
      # data_ref = "list",
      call = cl,
      session = list(
        engine_version = "ppca lib",
        R.version = getRversion(),
        pkg = list(
          m8 = as.character(utils::packageVersion("metabom8")),
          Rcpp = as.character(utils::packageVersion("Rcpp"))
        )
      )
  )

  methods::validObject(mod)
  mod
}

.pcaEngine <- function(X, scaling, ncomp = 2L, method) {

  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X)) stop("`X` must be a numeric matrix.", call. = FALSE)
  storage.mode(X) <- "double"

  if (!is.numeric(ncomp) || length(ncomp) != 1L || is.na(ncomp) || ncomp < 1)
    stop("`ncomp` must be a positive integer.", call. = FALSE)
  ncomp <- as.integer(ncomp)

  sc_num <- switch(as.character(scaling$scale),
                   "None"   = 0L,
                   "UV"     = 1L,
                   "Pareto" = 2L,
                   stop("Unsupported scaling@scale. Use 'None', 'UV', or 'Pareto'.", call. = FALSE)
  )

  pc_max <- min(nrow(X)-1, ncol(X))
  if (ncomp > pc_max) {
    warning("Too many components requested; resetting to maximum allowed (", pc_max, ").")
    ncomp <- pc_max
  }

  msd_x  <- .sdRcpp(X)
  XcsTot <- .scaleMatRcpp(
    X,
    0:(nrow(X) - 1),
    center     = scaling$center,
    scale_type = sc_num
  )[[1]]

  tssx <- .tssRcpp(XcsTot)

  n <- nrow(XcsTot)
  p <- ncol(XcsTot)


  method_lc <- tolower(method)

  if (identical(method_lc, "nipals")) {

    comps <- vector("list", ncomp)
    for (i in seq_len(ncomp)) {
      comps[[i]] <- if (i == 1L) .nipPcaCompRcpp(XcsTot) else .nipPcaCompRcpp(X = comps[[i - 1L]][[1]])
    }

    Tpc <- vapply(comps, function(z) as.numeric(z[[2]]), FUN.VALUE = numeric(n))
    Tpc <- matrix(Tpc, nrow = n, ncol = ncomp)

    Ppc <- vapply(comps, function(z) as.numeric(z[[3]]), FUN.VALUE = numeric(p))
    Ppc <- t(matrix(Ppc, nrow = p, ncol = ncomp))  # A x p

    r2x_comp <- vapply(seq_len(ncomp), function(i) {
      Xi <- Tpc[, i] %o% Ppc[i, ]
      .tssRcpp(Xi) / tssx
    }, numeric(1))

    full <- list(
      t = Tpc,
      p = Ppc,
      X_mean = msd_x[[2]],
      X_sd   = msd_x[[1]]
    )

    ctrl <- list(
      method = "nipals",
      type = 'unsupervised',
      ncomp_selected  = ncomp,
      ncomp_tested  = ncomp,
      r2x_comp = r2x_comp,
      tssx = tssx,
      center = scaling$center,
      scale  = as.character(scaling$scale)
    )

  } else {

    mod <- pcaMethods::pca(
      XcsTot,
      nPcs   = ncomp,
      scale  = "none",
      center = FALSE,
      method = method
    )

    Tpc <- as.matrix(mod@scores)

    Ppc <- t(as.matrix(mod@loadings))

    r2x_comp <- NA_real_
    if (!is.null(mod@R2cum) && length(mod@R2cum) >= 2L) {
      r2 <- mod@R2cum
      r2[-1L] <- diff(r2)
      r2x_comp <- r2[seq_len(min(length(r2), ncomp))]
    }

    full <- list(
      t = Tpc,
      p = Ppc,
      X_mean = msd_x[[2]],
      X_sd   = msd_x[[1]]
    )

    ctrl <- list(
      method = method,
      type = 'unsupervised',
      ncomp_selected  = ncomp,
      ncomp_tested  = ncomp,
      r2x_comp = r2x_comp,
      tssx = tssx,
      center = scaling$center,
      scale  = as.character(scaling$scale),
      type = "unsupervised"
    )
  }

  list(full = full, cv = NULL, ctrl = ctrl,  dims = list(n=n, p=p))
}

