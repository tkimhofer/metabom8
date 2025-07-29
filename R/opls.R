#' @title Orthogonal Partial Least Squares (O-PLS)
#'
#' @description
#' Fits an O-PLS model for regression or classification (discriminant analysis).
#' The number of orthogonal components is determined via cross-validation and
#' overfitting criteria.
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
#' plot(mod)
opls <- function(X, Y, center = TRUE, scale = "UV",
                 cv = list(method = "k-fold_stratified", k = 7, split = 2/3),
                 maxPCo = 5, plotting = TRUE) {

  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.logical(center)) stop("Check center parameter argument!")

  sc_num <- switch(scale,
                   "None" = 0,
                   "UV" = 1,
                   "Pareto" = 2,
                   stop("Check scale parameter argument!"))

  x_check <- .checkXclassNas(X)
  y_check <- .checkYclassNas(Y)
  xy_check <- .checkDimXY(X, y_check[[1]])
  type <- y_check[[3]]

  msd_y <- .sdRcpp(y_check[[1]])
  msd_x <- .sdRcpp(X)
  XcsTot <- .scaleMatRcpp(X, 0:(nrow(X) - 1), center = TRUE, scale_type = sc_num)[[1]]
  YcsTot <- .scaleMatRcpp(y_check[[1]], 0:(nrow(y_check[[1]]) - 1), center = TRUE, scale_type = sc_num)[[1]]

  tssx <- .tssRcpp(XcsTot)
  tssy <- .tssRcpp(YcsTot) / ncol(YcsTot)

  cv$cv_sets <- .cvSetsMethod(Y = y_check[[1]], type = type,
                              method = cv$method, k = cv$k, split = cv$split)
  cv$k <- length(cv$cv_sets)

  r2_comp <- r2x_comp <- q2_comp <- aucs_tr <- aucs_te <- array()
  nc <- 1
  overfitted <- FALSE

  while (!overfitted) {
    if (nc == 1) {
      tt <- .oplsComponentCv(X, Y = y_check[[1]], cv$cv_sets, nc, mod.cv = NULL)
    } else {
      tt <- .oplsComponentCv(X = NA, Y = y_check[[1]], cv$cv_sets, nc, mod.cv = tt)
    }

    preds_test <- .extMeanCvFeat(tt, feat = 'y_pred_test', cv_type = cv$method, model_type = type)
    preds_train <- .extMeanCvFeat(tt, feat = 'y_pred_train', cv_type = cv$method, model_type = type)

    cv_method <- strsplit(cv$method, '_')[[1]][1]
    is_DA <- grepl('DA', type)
    is_multi_Y <- grepl('mY', type)

    if (cv_method == 'MC') {
      if (is_DA) {
        if (is_multi_Y) {
          mod <- multiclass.roc(response = factor(Y), predictor = preds_test[1, , ], quiet=TRUE)
          aucs_te[nc] <- mod$auc
          mod <- multiclass.roc(response = factor(Y), predictor = preds_train[1, , ], quiet=TRUE)
          aucs_tr[nc] <- mod$auc
        } else {
          mod <- roc(response = Y, predictor = preds_test[1, ], quiet=TRUE)
          aucs_te[nc] <- mod$auc
          mod <- roc(response = Y, predictor = preds_train[1, ], quiet=TRUE)
          aucs_tr[nc] <- mod$auc
        }
      } else {
        r2_comp[nc] <- .r2(YcsTot, preds_train[1, , drop = FALSE], NULL)
        q2_comp[nc] <- .r2(YcsTot, preds_test[1, , drop = FALSE], tssy)
      }
    } else if (cv_method == 'k-fold') {
      if (is_DA) {
        if (is_multi_Y) {
          mod <- multiclass.roc(response = Y, predictor = apply(preds_test, 2, as.numeric), quiet=TRUE)
          aucs_te[nc] <- mod$auc
          mod <- multiclass.roc(response = Y, predictor = preds_train[1, , ], quiet=TRUE)
          aucs_tr[nc] <- mod$auc
        } else {
          mod <- roc(response = Y, predictor = preds_test[,1], quiet=TRUE)
          aucs_te[nc] <- mod$auc
          mod <- roc(response = Y, predictor = preds_train[1, ], quiet=TRUE)
          aucs_tr[nc] <- mod$auc
        }
      } else {
        r2_comp[nc] <- .r2(YcsTot, preds_train[1, , drop = FALSE], NULL)
        q2_comp[nc] <- .r2(YcsTot, preds_test[, 1], tssy)
      }
    }

    overfitted <- .evalFit(type, q2_comp, aucs_te, maxPCo)

    if (nc >= 2) {
      if (nc == 2) {
        opls_filt <- .nipOplsRcpp(X = XcsTot, Y = YcsTot)
        t_orth <- opls_filt$t_o
        p_orth <- opls_filt$p_o
        w_orth <- t(opls_filt$w_o)
      } else {
        opls_filt <- .nipOplsRcpp(X = opls_filt$X_res, Y = YcsTot)
        t_orth <- cbind(t_orth, opls_filt$t_o)
        p_orth <- rbind(p_orth, opls_filt$p_o)
        w_orth <- rbind(w_orth, t(opls_filt$w_o))
      }
      pred_comp <- .nipPlsCompRcpp(X = opls_filt$X_res, Y = YcsTot)
      r2x_comp[nc - 1] <- .r2(opls_filt$X_res, pred_comp$t_x %*% pred_comp$p_x, NULL)
    }

    if (overfitted) {
      message(sprintf("An O-PLS-%s model with 1 predictive and %d orthogonal components was fitted.",
                      type, nc - 1))
      pred_comp$t_cv <- matrix(.extMeanCvFeat(tt, 't_xp', cv_type = cv$method, model_type = type), ncol = 1)
      pred_comp$t_o_cv <- matrix(t(.extMeanCvFeat(tt, 't_xo', cv_type = cv$method, model_type = type)[, nc - 1]),
                                 ncol = nc - 1)
      pred_comp$t_o <- t_orth
    }
    nc <- nc + 1
  }

  m_summary <- .orthModelCompSummary(type, c(r2x_comp, NA), r2_comp, q2_comp, aucs_te, aucs_tr, cv)
  if (plotting) plot(m_summary[[2]])

  Xorth <- t_orth %*% p_orth
  Xpred <- pred_comp$t_x %*% pred_comp$p_x
  E <- XcsTot - (Xpred + Xorth)
  Y_res <- YcsTot - (pred_comp$t_y %*% t(pred_comp$p_y))

  pars <- list(center = center, scale = scale, nPredC = 1, nOrthC = nc - 1,
               maxPCo = maxPCo, cv = cv, tssx = tssx, tssy = tssy)

  new("OPLS_metabom8",
      type = type,
      t_pred = pred_comp$t_x,
      p_pred = pred_comp$p_x,
      w_pred = pred_comp$w_x,
      betas_pred = drop(pred_comp$b),
      Qpc = pred_comp$w_y,
      t_pred_cv = pred_comp$t_cv,
      t_orth_cv = pred_comp$t_o_cv,
      t_orth = t_orth,
      w_orth = w_orth,
      p_orth = p_orth,
      nPC = nc - 1,
      summary = m_summary[[1]],
      X_orth = Xorth,
      Y_res = Y_res,
      X_res = E,
      X_mean = msd_x[[2]],
      X_sd = msd_x[[1]],
      Y_mean = msd_y[[2]],
      Y_sd = msd_y[[1]],
      Parameters = pars,
      X = X,
      X_scaled = XcsTot,
      Y = list(ori = Y, dummy = y_check[[1]])
  )
}
