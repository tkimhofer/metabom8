#' #' @title Orthogonal Partial Least Squares (O-PLS)
#' #'
#' #' @description
#' #' Fits an Orthogonal Projections to Latent Structures (O-PLS) model for regression or classification.
#' #' The number of orthogonal components is automatically selected using internal cross-validation.
#' #' To avoid underfitting, components are added incrementally, while overfitting is prevented by requiring
#' that each new component improves predictive performance beyond a defined threshold (ΔQ² / ΔAUROC > 0.05).
#' This ensures the model captures relevant structure without modelling noise or irrelevant variation.
#'
#' @param X Numeric matrix or data frame of predictors. Rows are samples; columns are features.
#' @param Y Response variable. Vector or matrix. For classification, should be a factor or class labels; for regression, a numeric vector or matrix.
#' @param center Logical. Should features be mean-centered? Default is TRUE.
#' @param scale Character. Scaling method. Supported: "None", "UV" (unit variance), or "Pareto".
#' @param cv Named list specifying cross-validation settings: method ("k-fold", "k-fold_stratified"), split (training set proportion), and k (number of folds).
#' @param maxPCo Integer. Maximum number of orthogonal components to consider. Default is 5.
#' @param plotting Logical. If TRUE, shows model summary plot.
#'
#' @return An S4 object of class OPLS_metabom8.
#'
#' @references
#' Trygg J, Wold S. (2002). Orthogonal projections to latent structures (O-PLS).
#' Journal of Chemometrics, 16(3), 119–128.
#'
#' @author Torben Kimhofer
#' @importFrom stats cov sd var
#' @importFrom graphics plot
#' @importFrom methods new
#' @importFrom ggplot2 ggplot aes scale_fill_manual theme_bw labs
#' @importFrom pROC roc multiclass.roc
#' @export
#'
#' @examples
#' data("covid")
#' mod <- opls(X, an$type)
#' summary(mod)
opls <- function(X, Y, center = TRUE, scale = "UV",
                 cv = list(method = "k-fold_stratified", k = 7, split = 2/3),
                 maxPCo = 5, plotting = TRUE) {

  inputs <- .prepareInputs(X, Y, center, scale)
  type <- inputs$type
  is_multi_Y <- grepl('mY', type)
  class_memb <- if (is_multi_Y) sub("^Numeric\\.", "", names(apply(inputs$Ydummy[,-1], 2, which.max))) else NULL

  msd_y <- .sdRcpp(inputs$Y)
  msd_x <- .sdRcpp(inputs$X)
  XcsTot <- .scaleMatRcpp(inputs$X, 0:(nrow(inputs$X) - 1), center = TRUE, scale_type = inputs$scale_code)[[1]]
  YcsTot <- .scaleMatRcpp(inputs$Y, 0:(nrow(inputs$Y) - 1), center = TRUE, scale_type = inputs$scale_code)[[1]]
  tssx <- .tssRcpp(XcsTot)
  tssy <- .tssRcpp(YcsTot) / ncol(YcsTot)

  cv <- .generateCv(inputs$Y, type, cv)
  r2_comp <- r2x_comp <- q2_comp <- aucs_tr <- aucs_te <- array()
  nc <- 1
  overfitted <- FALSE

  while (!overfitted) {
    tt <- if (nc == 1) .oplsComponentCv(inputs$X, inputs$Y, cv$cv_sets, nc, mod.cv = NULL)
    else .oplsComponentCv(X = NA, inputs$Y, cv$cv_sets, nc, mod.cv = tt)

    preds_test <- .extMeanCvFeat(tt, feat = 'y_pred_test', cv_type = cv$method, model_type = type)
    preds_train <- .extMeanCvFeat(tt, feat = 'y_pred_train', cv_type = cv$method, model_type = type)

    perf <- .evalComponentPerformance(strsplit(cv$method, '_')[[1]][1], type, is_multi_Y, Y, preds_train, preds_test, class_memb, YcsTot, tssy)
    q2_comp[nc] <- perf$q2
    r2_comp[nc] <- perf$r2
    aucs_te[nc] <- perf$aucs_te
    aucs_tr[nc] <- perf$aucs_tr

    overfitted <- .evalFit(type, q2_comp, aucs_te, maxPCo+1)

    if (nc >= 2) {
      if (nc == 2) {
        opls_filt <- .nipOplsRcpp(XcsTot, YcsTot)
        t_orth <- opls_filt$t_o
        p_orth <- opls_filt$p_o
        w_orth <- t(opls_filt$w_o)
      } else {
        opls_filt <- .nipOplsRcpp(opls_filt$X_res, YcsTot)
        t_orth <- cbind(t_orth, opls_filt$t_o)
        p_orth <- rbind(p_orth, opls_filt$p_o)
        w_orth <- rbind(w_orth, t(opls_filt$w_o))
      }
      pred_comp <- .nipPlsCompRcpp(opls_filt$X_res, YcsTot, it_max = 800, eps = 1e-8)
      r2x_comp[nc - 1] <- .r2(opls_filt$X_res, pred_comp$t_x %*% pred_comp$p_x, NULL)
    }

    if (overfitted) {
      message(sprintf("An O-PLS-%s model with 1 predictive and %d orthogonal components was fitted.", type, nc - 1))
      pred_comp$t_o_cv <- matrix(.extMeanCvFeat(tt, 't_xo', cv_type = cv$method, model_type = type)[,seq_len(nc-1), drop=FALSE], ncol = nc - 1)
      pred_comp$t_cv <- matrix(.extMeanCvFeat(tt, 't_xp', cv_type = cv$method, model_type = type)[,nc-1, drop=FALSE], ncol = 1)
      pred_comp$t_o <- t_orth
    }

    nc <- nc + 1
  }

  .finaliseOplsModel(
    pred_comp = pred_comp,
    t_orth = t_orth,
    p_orth = p_orth,
    w_orth = w_orth,
    XcsTot = XcsTot,
    YcsTot = YcsTot,
    msd_x = msd_x,
    msd_y = msd_y,
    r2x_comp = r2x_comp,
    r2_comp = r2_comp,
    q2_comp = q2_comp,
    aucs_tr = aucs_tr,
    aucs_te = aucs_te,
    type = type,
    nc = nc,
    maxPCo = maxPCo,
    cv = cv,
    plotting = plotting,
    original_X = X,
    original_Y = Y,
    Y_dummy_cnam = colnames(inputs$Y)
  )
}


# Helper: Prepare and validate input
.prepareInputs <- function(X, Y, center, scale) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.logical(center)) stop("Check center parameter argument!")

  sc_num <- switch(scale,
                   "None" = 0,
                   "UV" = 1,
                   "Pareto" = 2,
                   stop("Check scale parameter argument!"))

  x_check <- .checkXclassNas(X)
  y_check <- .checkYclassNas(Y)
  .checkDimXY(X, y_check[[1]])

  list(X = X, Y = y_check[[1]], Ydummy = y_check[[2]],
       type = y_check[[3]], scale_code = sc_num)
}

# Helper: Generate cross-validation folds
.generateCv <- function(Y, type, cv) {
  cv$cv_sets <- .cvSetsMethod(Y, type, method = cv$method, k = cv$k, split = cv$split)
  cv$k <- length(cv$cv_sets)
  cv
}

# Helper: Evaluate component performance
.evalComponentPerformance <- function(cv_method, type, is_multi_Y, Y, preds_train, preds_test, class_memb, YcsTot, tssy) {
  if (cv_method == 'MC') {
    if (grepl('DA', type)) {
      if (is_multi_Y) {
        colnames(preds_test) <- colnames(preds_train) <- class_memb
        auc_te <- multiclass.roc(response = factor(Y), predictor = preds_test, quiet = TRUE)$auc
        auc_tr <- multiclass.roc(response = factor(Y), predictor = preds_train, quiet = TRUE)$auc
      } else {
        auc_te <- roc(response = Y, predictor = as.vector(preds_test), quiet = TRUE)$auc
        auc_tr <- roc(response = Y, predictor = as.vector(preds_train), quiet = TRUE)$auc
      }
      list(q2 = NA, r2 = NA, aucs_te = auc_te, aucs_tr = auc_tr)
    } else {
      list(q2 = .r2(YcsTot, preds_test, tssy), r2 = .r2(YcsTot, preds_train, NULL), aucs_te = NA, aucs_tr = NA)
    }
  } else {
    if (grepl('DA', type)) {
      if (is_multi_Y) {
        colnames(preds_test) <- colnames(preds_train) <- class_memb
        auc_te <- multiclass.roc(response = Y, predictor = preds_test, quiet = TRUE)$auc
        auc_tr <- multiclass.roc(response = Y, predictor = preds_train, quiet = TRUE)$auc
      } else {
        auc_te <- roc(response = Y, predictor = as.vector(preds_test), quiet = TRUE)$auc
        auc_tr <- roc(response = Y, predictor = as.vector(preds_train), quiet = TRUE)$auc
      }
      list(q2 = NA, r2 = NA, aucs_te = auc_te, aucs_tr = auc_tr)
    } else {
      list(q2 = .r2(YcsTot, preds_test, tssy), r2 = .r2(YcsTot, as.vector(preds_train), NULL), aucs_te = NA, aucs_tr = NA)
    }
  }
}

.finaliseOplsModel <- function(pred_comp, t_orth, p_orth, w_orth,
                               XcsTot, YcsTot, msd_x, msd_y,
                               r2x_comp, r2_comp, q2_comp, aucs_tr, aucs_te,
                               type, nc, maxPCo, cv, plotting,
                               original_X, original_Y, Y_dummy_cnam) {

  m_summary <- .orthModelCompSummary(
    type = type,
    r2x_comp = c(r2x_comp, NA),
    r2_comp = r2_comp,
    q2_comp = q2_comp,
    aucs_tr = aucs_tr,
    aucs_te = aucs_te,
    cv = cv
  )
  if (plotting) plot(m_summary[[2]])

  # Final matrices
  Xorth <- t_orth %*% p_orth
  Xpred <- pred_comp$t_x %*% pred_comp$p_x
  E <- XcsTot - (Xpred + Xorth)
  Y_res <- YcsTot - (pred_comp$t_y %*% t(pred_comp$p_y))

  # Parameters
  pars <- list(
    center = TRUE,
    scale = "UV",
    nPredC = 1,
    nOrthC = nc - 1,
    maxPCo = maxPCo,
    cv = cv,
    tssx = .tssRcpp(XcsTot),
    tssy = .tssRcpp(YcsTot) / ncol(YcsTot)
  )

  # Return S4 object
  new("OPLS_metabom8",
      type        = type,
      t_pred      = pred_comp$t_x,
      p_pred      = pred_comp$p_x,
      w_pred      = pred_comp$w_x,
      betas_pred  = drop(pred_comp$b),
      Qpc         = pred_comp$w_y,
      t_pred_cv   = pred_comp$t_cv,
      t_orth_cv   = pred_comp$t_o_cv,
      t_orth      = t_orth,
      w_orth      = w_orth,
      p_orth      = p_orth,
      nPC         = nc - 1,
      summary     = m_summary[[1]],
      X_orth      = Xorth,
      Y_res       = Y_res,
      X_res       = E,
      X_mean      = msd_x[[2]],
      X_sd        = msd_x[[1]],
      Y_mean      = msd_y[[2]],
      Y_sd        = msd_y[[1]],
      Parameters  = pars,
      X           = original_X,
      X_scaled    = XcsTot,
      Y           = list(ori = original_Y, dummy = YcsTot, cname=Y_dummy_cnam)
  )
}

