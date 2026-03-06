#' Fit a Partial Least Squares (PLS) model
#'
#' Fits a supervised Partial Least Squares (PLS) model using a NIPALS-based
#' algorithm with optional cross-validation and automatic component selection.
#'
#' @param X Numeric matrix of predictors (rows = samples, columns = variables).
#' @param Y Numeric matrix or factor vector of responses.
#' @param scaling A scaling strategy object (e.g., \code{UVScaling(center = TRUE)}),
#'   specifying model-internal centering and/or scaling applied during fitting.
#'   This does not modify the original spectral matrix.
#' @param validation_strategy A cross-validation strategy object defining how
#'   resampling is performed (e.g., k-fold, Monte Carlo).
#' @param maxPCo Integer. Maximum number of predictive components to evaluate.
#'
#' @details
#' Model components are estimated sequentially using a NIPALS-based algorithm.
#' For each component, cross-validated performance metrics (e.g., Q², R²,
#' classification AUC) are computed according to the supplied
#' \code{validation_strategy}. Component extraction stops when the
#' \code{stopRule} indicates overfitting or when \code{maxPCo} is reached.
#'
#' Scaling specified via \code{scaling} is applied internally during model
#' fitting and does not alter the input matrix \code{X}. Spectral preprocessing
#' steps (e.g., alignment, baseline correction) should be performed prior to
#' model fitting.
#'
#' The returned model object stores:
#' \itemize{
#'   \item Fitted component models
#'   \item Cross-validation results
#'   \item Performance metrics (R², Q², AUC)
#'   \item Model control parameters
#'   \item Input data provenance metadata
#'   \item Session information for reproducibility
#' }
#'
#' @return
#' An object of class \code{m8_model} containing the fitted PLS model,
#' cross-validation results, and performance statistics.
#'
#' @seealso \code{\link{opls}}, \code{\link{UVScaling}}
#'
#' @family modelling
#' @examples
#' # example code
#'
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- UVScaling(center=TRUE)
#' model <-pls(X=covid$X, Y=covid$an$type, scaling, cv)
#' show(model)
#' summary(model)
#'
#' Tp <- scores(model)
#' Pp <- loadings(model)
#' @export
pls <- function(X, Y, scaling, validation_strategy, maxPCo=5) {

  cl <- match.call()

  engine <- .plsEngine(X, Y, scaling, validation_strategy,  maxPCo=maxPCo)

  mod <- methods::new("m8_model",

      engine = "pls",
      fit = engine$full,
      prep = scaling,
      cv = engine$cv,
      ctrl = engine$ctrl,
      dims = engine$dims,

      provenance = attributes(X),
      # data_ref = "list",
      call = cl,
      session = list(
        engine_version = "nipals_v1",
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


.plsEngine <- function(X, Y, scaling, validation_strategy, maxPCo){

  inputs <- .prepareInputs(X, Y, scaling@center, scaling@scale)
  type <- inputs$type
  is_multi_Y <- grepl('mY', type)

  cv <- instantiate(validation_strategy, inputs$Y)

  preppedX <- prep(scaling, inputs$X)
  XcsTot <- preppedX$X
  YcsTot <- .scaleMatRcpp(inputs$Y, 0:(nrow(inputs$Y) - 1), center = scaling@center, scale_type = inputs$scale_code)[[1]]
  tssx <- .tssRcpp(XcsTot)
  tssy <- .tssRcpp(YcsTot) / ncol(YcsTot)

  Ycs_fold <- lapply(cv@train, function(idc)
    .scaleMatRcpp(inputs$Y, idc - 1, scaling@center, inputs$scale_code)[[1]]
  )

  tt <- NULL
  full_mod <- list()

  n  <- nrow(XcsTot)
  q <- ncol(inputs$Y)
  k <- length(cv@train)

  acc <- list(
    sum_test  = matrix(0, n, q),
    n_test    = integer(n),
    sum_train = matrix(0, n, q),
    n_train   = integer(n)
  )


  r2_comp <- r2x_comp <- q2_comp <- aucs_tr <- aucs_te <- array()
  nc <- 1
  overfitted <- FALSE

  while (!overfitted) {

    res <- if (nc == 1) {
      .plsComponentCv(inputs$X, cv.set=cv@train, Ycs_fold = Ycs_fold,  nc=nc, mod.cv = NULL, acc=acc) # , it_max=800, eps=1e-6
    } else {
      .plsComponentCv(NULL,cv.set=cv@train, Ycs_fold = Ycs_fold,  nc=nc, mod.cv = mod_cv, acc=acc)
    }

    mod_cv <- res$mod.cv
    acc <- res$acc

    preds_test <- acc$sum_test
    preds_test[acc$n_test > 0,] <- preds_test[acc$n_test > 0,] / acc$n_test[acc$n_test > 0]

    preds_train <- acc$sum_train
    preds_train[acc$n_train > 0,] <- preds_train[acc$n_train > 0,] / acc$n_train[acc$n_train > 0]

    perf <- .evalComponentPerformance(
      cv_inst      = cv,
      type         = type,
      is_multi_Y   = is_multi_Y,
      Y            = Y,
      preds_train  = preds_train,
      preds_test   = preds_test,
      class_memb   = if (is.data.frame(inputs$Ydummy)) inputs$Ydummy$Level[inputs$Ydummy$Column] else NULL,
      YcsTot       = YcsTot,
      tssy         = tssy
    )

    q2_comp[nc] <- perf$q2
    r2_comp[nc] <- perf$r2
    aucs_te[nc] <- perf$aucs_te
    aucs_tr[nc] <- perf$aucs_tr


    if (nc == 1) {
      full_mod[[1]] <- .nipPlsCompRcpp(XcsTot, YcsTot, it_max=800, eps=1e-6)
    } else {
      Xinp <- full_mod[[nc - 1]]$x_res
      full_mod[[nc]] <- .nipPlsCompRcpp(Xinp, YcsTot, it_max=800, eps=1e-6)
    }

    full_mod[[nc]]$t_x

    full_mod[[nc]]$t_x <- .ensure_matrix(full_mod[[nc]]$t_x)
    full_mod[[nc]]$p_x <- .ensure_matrix(full_mod[[nc]]$p_x, ncol=ncol(XcsTot))

    r2x_comp[nc] <- .r2(XcsTot, full_mod[[nc]]$t_x %*% full_mod[[nc]]$p_x, NULL)

    fit_eval <- .evalFit(type, q2_comp, aucs_te, maxPCo)
    overfitted <- fit_eval$stop

    if (!overfitted) {
      nc <- nc + 1
      }
    else{


      if (nc == 1){
        message(sprintf("The PLS-%s model does not gerenalise: %s -> no stable predictive structure detected.", type, fit_eval$reason))
        full <-  list("not fitted as model does not generalise")
      }else{
        if (nc == 2){
          message(sprintf("A PLS-%s model with 1 component was fitted.", type))
        } else if (nc > 2){
          message(sprintf("A PLS-%s model with %d components was fitted.", type, nc-1))
        }

        full <- list(comp=full_mod[seq_len(nc-1)])

        full$X_mean <- preppedX$prep@X_mean
        full$X_sd <- preppedX$prep@X_sd
        full$X_prepped <- preppedX$X
        full$Y <- inputs$Y
      }

      break
    }
  }

  list(
    full = full,
    cv = cv,
    ctrl = list(
      ncomp_selected  = nc - 1,
      ncomp_tested  = nc,
      stop_reason = fit_eval$reason,
      stop_metric = fit_eval$metric,
      stop_delta  = fit_eval$delta,

      r2x_comp = r2x_comp,

      q2 = q2_comp,
      r2 = r2_comp,

      aucs_tr = aucs_tr,
      aucs_te = aucs_te,

      tssx = tssx,
      tssy = tssy,

      type = type,
      is_multy_Y = is_multi_Y,
      maxPCo = maxPCo
    ),
    dims = list(n=n, p=ncol(X))
  )

}

.plsComponentCv <- function(X, cv.set, Ycs_fold, nc, mod.cv, acc) {

  n <- nrow(Ycs_fold[[1]])

  if (nc == 1) {
    q <- ncol(Ycs_fold[[1]])
    p <- ncol(X)
    mod.cv <- vector("list", length(cv.set))

    for (i in seq_along(cv.set)) {
      mod.cv[[i]] <- list(
        t_xp = matrix(0, n, 1),
        y_pred_train = matrix(0, n, q),
        y_pred_test = matrix(0, n, q),
        x_res = matrix(0, n, p)
      )
    }

  } else {

    for (i in seq_along(mod.cv)) {
      mod.cv[[i]]$t_xp <- cbind(mod.cv[[i]]$t_xp, 0)
    }

  }

  for (k in seq_along(cv.set)) {

    idc <- cv.set[[k]]
    not_idc <- setdiff(seq_len(n), unique(idc))

    Xcs <- if (nc == 1)
      .scaleMatRcpp(X, idc - 1, center = TRUE, scale_type = 1)[[1]]
    else
      mod.cv[[k]]$x_res

    Ycs <- Ycs_fold[[k]]

    comp <- .nipPlsCompRcpp(
      X = Xcs[idc, , drop = FALSE],
      Y = Ycs[idc, , drop = FALSE],
      it_max = 800,
      eps = 1e-8
    )

    p <- length(comp$w_x)
    comp$w_x <- .ensure_matrix(comp$w_x, ncol=p)
    comp$p_x <- .ensure_matrix(comp$p_x, ncol=p)
    comp$p_y <- .ensure_matrix(comp$p_y)
    comp$b   <- as.numeric(comp$b)

    Xte <- as.matrix(Xcs[not_idc, , drop = FALSE])
    storage.mode(Xte) <- "double"

    t_pred <- Xte %*% t(comp$w_x)
    if (nrow(comp$p_y) == 1)
      comp$p_y <- t(comp$p_y)
    y_pred <- comp$b * t_pred %*% t(comp$p_y)

    mod.cv[[k]]$t_xp[not_idc, nc] <- t_pred
    mod.cv[[k]]$y_pred_train[idc, ] <- comp$y_pred
    mod.cv[[k]]$y_pred_test[not_idc, ] <- y_pred

    mod.cv[[k]]$x_res[idc, ] <- comp$x_res
    mod.cv[[k]]$x_res[not_idc, ] <- Xte - t_pred %*% comp$p_x

    acc$sum_test[not_idc, ] <- acc$sum_test[not_idc, ] + y_pred
    acc$n_test[not_idc] <- acc$n_test[not_idc] + 1

    uid <- unique(idc)
    agg_pred <- rowsum(comp$y_pred, idc)
    idx <- as.integer(rownames(agg_pred))

    acc$sum_train[idx, ] <- acc$sum_train[idx, ] + agg_pred
    acc$n_train[idx] <- acc$n_train[idx] + tabulate(idc, nbins = nrow(acc$sum_train))[idx]
  }

  return(list(mod.cv = mod.cv, acc = acc))
}
