#' @title Partial Least Squares (PLS)
#' @description Fit a Partial Least Squares (PLS) model for regression (R) or discriminant analysis (DA), with automated cross-validated component selection.
#'
#' @param X A numeric matrix or data frame. Each row represents an observation, and each column a metabolic variable.
#' @param Y Response vector or matrix. Must match the number of rows in \code{X}.
#' @param center Logical. Should data be mean-centered? Default is \code{TRUE}.
#' @param scale Character. Scaling method: \code{"None"}, \code{"UV"} (unit variance), or \code{"Pareto"}.
#' @param cv Named list specifying cross-validation settings:
#'   \describe{
#'     \item{\code{method}}{Cross-validation type: \code{"k-fold"}, \code{"k-fold_stratified"}, \code{"MC"}, or \code{"MC_balanced"}.}
#'     \item{\code{split}}{Fraction of observations used for training (used in Monte Carlo CV).}
#'     \item{\code{k}}{Number of folds or repetitions.}
#'   }
#' @param maxPCo Integer. Maximum number of orthogonal components to test.
#' @param plotting Logical. If \code{TRUE}, model summary (e.g. R2X, Q2, AUROC) is plotted. Default is \code{TRUE}.
#'
#' @details
#' Cross-validation is used to select the optimal number of predictive components based on Q2 or AUROC. The method supports both regression and classification with binary or multi-class responses. Model interpretability is often best with pairwise class comparisons.
#'
#' @return An object of class \code{PLS_metabom8}, an S4 class with scores, loadings, predictions, and validation statistics.
#'
#' @references
#' Geladi, P. & Kowalski, B.R. (1986). Partial least squares regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1–17.
#'
#' @examples
#' data(covid)
#' model <- pls(X, Y = an$type)
#' plotscores(model, an = list(Class = an$type, Clinic = an$hospital, id = 1:nrow(an)), pc = c(1, 2))
#'
#' @importFrom graphics plot
#' @importFrom methods new
#' @importFrom stats cov sd var
#' @importFrom ggplot2 ggplot aes labs theme_bw
#' @importFrom pROC roc multiclass.roc
#' @export

pls <- function(X, Y,
                center = TRUE,
                scale = "UV",
                cv = list(method = "k-fold_stratified", k = 7, split = 2 / 3),
                maxPCo = 5,
                plotting = TRUE) {

  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.logical(center)) stop("Invalid 'center' parameter.")

  sc_num <- switch(scale,
                   "None" = 0,
                   "UV" = 1,
                   "Pareto" = 2,
                   NULL
  )
  if (is.null(sc_num)) stop("Invalid 'scale' parameter. Use 'None', 'UV', or 'Pareto'.")

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

  cv$cv_sets <- .cvSetsMethod(Y = y_check[[1]], type = type, method = cv$method, k = cv$k, split = cv$split)
  cv$k <- length(cv$cv_sets)

  r2_comp <- r2x_comp <- q2_comp <- aucs_tr <- aucs_te <- array()
  nc <- 1
  overfitted <- FALSE

  while (overfitted == FALSE) {
    if (nc == 1) {
      tt <- .plsComponentCv(X, Y = y_check[[1]], cv$cv_sets, nc, mod.cv = NULL)
    } else {
      tt <- .plsComponentCv(X = NA, Y = y_check[[1]], cv$cv_sets, nc, mod.cv = tt)
    }

    # preds_test <- .extMeanCvFeat_pls(cv_obj = tt, feat = "y_pred_test", cv_type = cv$method, model_type = type)
    # preds_train <- .extMeanCvFeat_pls(tt, feat = "y_pred_train", cv_type = cv$method, model_type = type)
    #
    preds_test <- .extMeanCvFeat(cv_obj = tt, feat = "y_pred_test", cv_type = cv$method, model_type = type)
    preds_train <- .extMeanCvFeat(tt, feat = "y_pred_train", cv_type = cv$method, model_type = type)

    switch(strsplit(cv$method, '_')[[1]][1],
           'MC' = {
             if (grepl('DA', type)) {
               if (grepl('mY', type)) {
                 pred_mean <- preds_test[1, , ]
                 colnames(pred_mean) <- colnames(y_check[[1]])
                 aucs_te[nc] <- multiclass.roc(response = factor(Y), predictor = pred_mean, quiet=TRUE)$auc

                 pred_tr_mean <- preds_train[1, , ]
                 colnames(pred_tr_mean) <- colnames(y_check[[1]])
                 aucs_tr[nc] <- multiclass.roc(response = factor(Y), predictor = pred_tr_mean, quiet=TRUE)$auc
               } else {
                 aucs_te[nc] <- roc(response = Y, predictor = preds_test[1, ], quiet=TRUE)$auc
                 aucs_tr[nc] <- roc(response = Y, predictor = preds_train[1, ], quiet=TRUE)$auc
               }
             } else {
               r2_comp[nc] <- .r2(YcsTot, preds_train[1, ], NULL)
               q2_comp[nc] <- .r2(YcsTot, preds_test[1, ], tssy)
             }
           },
           'k-fold' = {
             if (grepl('DA', type)) {
               if (grepl('mY', type)) {
                 colnames(preds_test) <- colnames(y_check[[1]])
                 aucs_te[nc] <- multiclass.roc(response = Y, predictor = apply(preds_test, 2, as.numeric), quiet=TRUE)$auc
                 preds_te <- preds_train
                 colnames(preds_te) <- colnames(y_check[[1]])
                 aucs_tr[nc] <- multiclass.roc(response = Y, predictor = preds_te, quiet=TRUE)$auc
               } else {
                 aucs_te[nc] <- roc(response = Y, predictor = preds_test[, 1], quiet=TRUE)$auc
                 aucs_tr[nc] <- roc(response = Y, predictor = preds_train[1, ], quiet=TRUE)$auc
               }
             } else {
               r2_comp[nc] <- .r2(YcsTot, preds_train[1, ], NULL)
               q2_comp[nc] <- .r2(YcsTot, preds_test[, 1], tssy)
             }
           }
    )

    # if (nc==1){ X_pred = XcsTot}
    #
    #
    # pred_comp <- .nipPlsCompRcpp(X = opls_filt$X_res, Y = YcsTot)
    # r2x_comp[nc - 1] <- .r2(opls_filt$X_res, pred_comp$t_x %*% pred_comp$p_x, NULL)

    overfitted <- .evalFit(type, q2_comp, aucs_te, maxPCo)
    nc <- nc + 1
  }

  message(sprintf("A PLS-%s model with %d components was fitted.", type, nc - 1))
  mod_summary <- .orthModelCompSummary(type = type, r2x_comp = c(r2x_comp, NA), r2_comp = r2_comp, q2_comp = q2_comp, aucs_te = aucs_te, aucs_tr = aucs_tr, cv = cv)
  if (plotting) plot(mod_summary[[2]] + labs(x = 'Predictive Components'))

  final_model <- .nipPlsCompRcpp(X = XcsTot, Y = YcsTot)
  Xpred <- final_model$t_x %*% final_model$p_x
  E <- XcsTot - Xpred
  Y_res <- YcsTot - (final_model$t_y %*% t(final_model$p_y))

  pars <- list(
    center = center,
    scale = scale,
    nPredC = 1,
    nOrthC = nc - 1,
    maxPCo = maxPCo,
    cv = cv,
    tssx = tssx,
    tssy = tssy
  )

  new("PLS_metabom8",
      type = type,
      t_pred = final_model$t_x,
      p_pred = final_model$p_x,
      w_pred = final_model$w_x,
      betas_pred = drop(final_model$b),
      Qpc = final_model$w_y,
      t_pred_cv = final_model$t_x,
      nPC = nc - 1,
      summary = mod_summary[[1]],
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
