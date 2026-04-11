#' Fit an Orthogonal Partial Least Squares (O-PLS) model
#'
#' Fits a supervised Orthogonal Partial Least Squares (O-PLS) model
#' using a NIPALS-based algorithm with optional cross-validation and
#' automatic component selection.
#'
#' O-PLS decomposes the predictor matrix into:
#' \itemize{
#'   \item One predictive component capturing variation correlated with \code{Y}
#'   \item Orthogonal components capturing structured variation in \code{X}
#'         unrelated to \code{Y}
#' }
#'
#' @param X Numeric matrix of predictors (rows = samples, columns = variables).
#' @param Y Numeric matrix or factor vector of responses.
#' @param scaling A scaling strategy object (e.g., \code{uv_scaling(center = TRUE)}),
#'   specifying model-internal centering and/or scaling applied during fitting.
#'   This does not modify the original spectral matrix.
#' @param validation_strategy A cross-validation strategy object defining how
#'   resampling is performed (e.g., k-fold, Monte Carlo).
#'
#' @details
#' Predictive and orthogonal components are estimated sequentially.
#' Cross-validated performance metrics (e.g., Q², R², classification AUC)
#' are computed for each model configuration according to the supplied
#' \code{validation_strategy}.
#'
#' The model extracts a single predictive component and iteratively
#' adds orthogonal components until the \code{stopRule} indicates
#' overfitting or \code{maxPCo} is reached.
#'
#' Scaling specified via \code{scaling} is applied internally during model
#' fitting and does not alter the input matrix \code{X}. Spectral preprocessing
#' (e.g., alignment or baseline correction) should be performed prior to
#' model fitting.
#'
#' The returned model object stores:
#' \itemize{
#'   \item Predictive and orthogonal component models
#'   \item Cross-validation results
#'   \item Performance metrics (R², Q², AUC)
#'   \item Model control parameters
#'   \item Input data provenance metadata
#'   \item Session information for reproducibility
#' }
#'
#' @return
#' An object of class \code{m8_model} containing the fitted O-PLS model,
#' cross-validation results, and performance statistics.
#'
#' @seealso \code{\link{pls}}, \code{\link{uv_scaling}}
#'
#' @family modelling
#' @examples
#' data(covid)
#'
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- uv_scaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#'
#' show(model)
#' summary(model)
#'
#' # scores
#' Tp <- scores(model)
#' To <- scores(model, orth=TRUE)
#'
#' t2 <- hotellingsT2(cbind(Tp, To))
#' ell <-ellipse2d(t2)
#'
#' plot(Tp, To, asp = 1,
#'   col = as.factor(covid$an$type),
#'   xlim = range(c(Tp, ell$x)),
#'   ylim = range(c(To, ell$y))
#'  )
#' lines(ell$x, ell$y, col = "grey", lty=2)
#'
#' # loadings & vip's
#' Pp <- loadings(model)
#' Po <- loadings(model, orth=TRUE)
#' vips <- vip(model)
#'
#' x=covid$ppm
#' y = Pp * apply(covid$X, 2, sd)
#'
#' palette <-  colorRampPalette(c("blue", "cyan", "yellow", "red"))(100)
#' idx <- cut(vips, breaks = 100, labels = FALSE)
#' plot(x, y, type = "n", xlim = rev(range(x)), xlab='ppm', ylab='t_pred_sc')
#'
#' for (i in seq_len(length(x) - 1)) {
#'   segments(x[i], y[i], x[i+1], y[i+1], col = palette[idx[i]], lwd = 2)
#' }
#'
#' @export
opls <- function(X, Y, scaling, validation_strategy) {

  cl <- match.call()

  engine <- .oplsEngine(X, Y, scaling, validation_strategy, maxPCo=5)

  obj <- methods::new("m8_model",

      engine = "opls",
      fit = engine$full,
      prep = scaling,
      cv = engine$cv,
      ctrl = engine$ctrl,
      dims = engine$dims,

      provenance = attributes(X),
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

  methods::validObject(obj)
  obj

}

.oplsEngine <- function(X, Y, scaling, validation_strategy, maxPCo, cv=NULL) {

  #### checking scaling
  .arg_check_scaling(scaling)

  inputs <- .prepareInputs(X, Y, scaling$center, scaling$scale)
  type <- inputs$type
  is_multi_Y <- grepl('mY', type)
  if (is_multi_Y) {stop("Multi-response Y is currently not supported", call. = FALSE)}
  n  <- nrow(inputs$Y)
  p  <- ncol(inputs$X)
  mo <- maxPCo

  cv <- .arg_check_cv(cv_pars=validation_strategy, model_type=type, n=n, Y_prepped=inputs$Y)

  preppedX <- prep_X(scaling, inputs$X)
  XcsTot <- preppedX$X
  YcsTot <- .scaleMatRcpp(inputs$Y, 0:(nrow(inputs$Y) - 1), center = scaling$center, scale_type = inputs$scale_code)[[1]]

  tssx <- .tssRcpp(XcsTot)
  tssy <- .tssRcpp(YcsTot) / ncol(YcsTot)

  Ycs_fold <- lapply(cv$train, function(idc)
    .scaleMatRcpp(inputs$Y, idc - 1, scaling$center, inputs$scale_code)[[1]]
  )

  tt         <- NULL
  opls_filt  <- NULL

  t_orth <- matrix(0, n, mo)
  p_orth <- matrix(0, mo, p)
  w_orth <- matrix(0, mo, p)

  q <- ncol(inputs$Y)
  k <- length(cv$train)

  acc <- list(
    sum_test  = matrix(0, n, q),
    n_test    = integer(n),
    sum_train = matrix(0, n, q),
    n_train   = integer(n)
  )


  r2_comp <- r2x_comp <- q2_comp <- aucs_tr <- aucs_te <- array()
  nc <- 1
  overfitted <- FALSE

  while(!overfitted){

    res <- if (nc == 1)
      .oplsComponentCv(
        X        = inputs$X,
        cv.set   = cv$train,
        Ycs_fold = Ycs_fold,
        nc       = nc,
        mod.cv   = NULL,
        acc      = acc
      )
    else
      .oplsComponentCv(
        X        = NULL,
        cv.set   = cv$train,
        Ycs_fold = Ycs_fold,
        nc       = nc,
        mod.cv   = tt,
        acc      = acc
      )

    tt  <- res$mod.cv
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

    fit_eval <- .evalFit(type, q2_comp, aucs_te, maxPCo)
    overfitted <- fit_eval$stop


    if (nc >= 2) {
      if (nc == 2) {
        opls_filt <- .nipOplsRcpp(XcsTot, YcsTot)
      } else {
        opls_filt <- .nipOplsRcpp(opls_filt$X_res, YcsTot)
      }
      t_orth[,nc-1] <- opls_filt$t_o
      p_orth[nc-1,] <- opls_filt$p_o
      w_orth[nc-1,] <- t(opls_filt$w_o)

      pred_comp <- .nipPlsCompRcpp(opls_filt$X_res, YcsTot, it_max = 800, eps = 1e-8)

      Tx <- .ensure_matrix(pred_comp$t_x, ncol = 1L)
      Px <- .ensure_matrix(pred_comp$p_x, ncol=length(pred_comp$p_x))

      X_str <- Tx %*% Px

      r2x_comp[nc - 1] <- .r2(opls_filt$X_res, X_str, NULL)
    }

    if (overfitted) {
      if(nc >1){
        message(sprintf("An O-PLS-%s model with 1 predictive and %d orthogonal components was fitted.", type, nc - 1))
        pred_comp$t_o_cv <- matrix(
          .extMeanCvFeat(tt, "t_xo")[, seq_len(nc-1), drop = FALSE],
          ncol = nc - 1
        )

        pred_comp$t_cv <- matrix(
          .extMeanCvFeat(tt, "t_xp")[, nc-1, drop = FALSE],
          ncol = 1
        )
        pred_comp$t_o <- t_orth
      } else{
        pred_comp <- "not fitted as model does not generalise"
        message(sprintf("The O-PLS-%s model does not gerenalise: %s -> no stable predictive structure detected.", type, fit_eval$reason))
      }

      break

    }else{
      nc <- nc + 1
    }
  }

  nc_orth <- max(0L, nc - 1L)

  t_orth_out <- if (nc_orth > 0L) t_orth[, seq_len(nc_orth), drop = FALSE] else
    matrix(numeric(0), nrow = n, ncol = 0)

  p_orth_out <- if (nc_orth > 0L) p_orth[seq_len(nc_orth), , drop = FALSE] else
    matrix(numeric(0), nrow = 0, ncol = p)

  w_orth_out <- if (nc_orth > 0L) w_orth[seq_len(nc_orth), , drop = FALSE] else
    matrix(numeric(0), nrow = 0, ncol = p)

  list(

    full = list(
      pred_comp = pred_comp,
      t_orth = t_orth_out,
      p_orth = p_orth_out,
      w_orth = w_orth_out,
      X_mean = preppedX$prep$X_mean,
      X_sd = preppedX$prep$X_sd,
      X_prepped = preppedX$X,
      Y = inputs$Y
    ),
    cv = cv,
    ctrl = list(
      nc_pred = if(nc > 1) 1 else 0,
      nc_orth = nc_orth,
      r2x_comp = r2x_comp,

      q2 = q2_comp,
      r2 = r2_comp,

      aucs_tr = aucs_tr,
      aucs_te = aucs_te,

      tssx = tssx,
      tssy = tssy,

      type = type,
      is_multy_Y = is_multi_Y,
      maxPCo = maxPCo,

      ncomp_selected = (if(nc > 1) 1 else 0) + nc_orth,  # orth comps
      ncomp_tested   = nc + 1,  # total orth stages evaluated
      stop_reason = fit_eval$reason,
      stop_metric = fit_eval$metric,
      stop_delta  = fit_eval$delta

    ),
    dims = list(n=n, p=p)
  )
}

.evalComponentPerformance <- function(
    cv_inst,
    type,
    is_multi_Y,
    Y,
    preds_train,
    preds_test,
    class_memb,
    YcsTot,
    tssy
){

  if (grepl('DA', type)) {

    if (is_multi_Y) {

      colnames(preds_test)  <- class_memb
      colnames(preds_train) <- class_memb

      auc_te <- multiclass.roc(
        response  = factor(Y),
        predictor = preds_test,
        quiet     = TRUE
      )$auc

      auc_tr <- multiclass.roc(
        response  = factor(Y),
        predictor = preds_train,
        quiet     = TRUE
      )$auc

    } else {

      auc_te <- roc(
        response  = as.vector(Y),
        predictor = as.vector(preds_test),
        quiet     = TRUE
      )$auc

      auc_tr <- roc(
        response  = as.vector(Y),
        predictor = as.vector(preds_train),
        quiet     = TRUE
      )$auc
    }

    return(list(
      q2      = NA,
      r2      = NA,
      aucs_te = auc_te,
      aucs_tr = auc_tr
    ))
  }

  list(
    q2      = .r2(YcsTot, preds_test, tssy),
    r2      = .r2(YcsTot, preds_train, NULL),
    aucs_te = NA,
    aucs_tr = NA
  )
}

.prepareInputs <- function(X, Y, center, scale) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.logical(center)) stop("Check center parameter argument!", call. = FALSE)

  sc_num <- switch(scale,
                   "None" = 0,
                   "UV" = 1,
                   "Pareto" = 2,
                   stop("Check scale parameter argument!", call. = FALSE))

  x_check <- .checkXclassNas(X)
  y_check <- .checkYclassNas(Y)
  .checkDimXY(X, y_check[[1]])

  list(X = X, Y = y_check[[1]], Ydummy = y_check[[2]],
       type = y_check[[3]], scale_code = sc_num)
}
