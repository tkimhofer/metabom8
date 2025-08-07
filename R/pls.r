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

  prep <- .prepareInputsPLS(X, Y, center, scale)
  type <- prep$y_check[[3]]

  fit <- .runPlsCVLoop(
    Xcs = prep$Xcs, Ycs = prep$Ycs, Y_ = prep$y_check[[1]], Y_vec=Y,
    type = type, cv = cv, maxPCo = maxPCo
  )

  mod_summary <- .plsModelCompSummary(
    type = type,
    r2x_comp = cumsum(fit$r2x_comp),
    r2_comp = fit$r2_comp,
    q2_comp = fit$q2_comp,
    aucs_te = fit$aucs_te,
    aucs_tr = fit$aucs_tr,
    cv = cv
  )
  if (plotting) plot(mod_summary[[2]])

  pars <- list(
    center = center,
    scale = scale,
    nPredC = 1,
    nOrthC = fit$n_components,
    maxPCo = maxPCo,
    cv = cv,
    tssx = prep$tssx,
    tssy = prep$tssy,
    r2x_comp = fit$r2x_comp
  )

  .finalisePlsModel(fit$full_mod, prep$Xcs, prep$Ycs, prep$msd_x, prep$msd_y, type, pars, mod_summary[[1]], Y)
}


.prepareInputsPLS <- function(X, Y, center, scale) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.logical(center)) stop("Invalid 'center' parameter.")

  sc_num <- switch(scale,
                   "None" = 0,
                   "UV" = 1,
                   "Pareto" = 2,
                   stop("Invalid 'scale' parameter.")
  )

  y_check <- .checkYclassNas(Y)
  x_check <- .checkXclassNas(X)
  .checkDimXY(X, y_check[[1]])

  Xcs <- .scaleMatRcpp(X, 0:(nrow(X) - 1), center = TRUE, scale_type = sc_num)[[1]]
  Ycs <- .scaleMatRcpp(y_check[[1]], 0:(nrow(y_check[[1]]) - 1), center = TRUE, scale_type = sc_num)[[1]]

  Y <- as.matrix(Y)[, 1, drop = FALSE]

  list(
    X = X, Y = Y,
    Xcs = Xcs, Ycs = Ycs,
    tssx = .tssRcpp(Xcs),
    tssy = .tssRcpp(Ycs) / ncol(Ycs),
    msd_x = .sdRcpp(X),
    msd_y = .sdRcpp(y_check[[1]]),
    y_check = y_check
  )
}

.runPlsCVLoop <- function(Xcs, Ycs, Y_, Y_vec, type, cv, maxPCo) {
  cv_sets <- .cvSetsMethod(Y = Y_, type = type, method = cv$method, k = cv$k, split = cv$split)
  preds <- list()
  full_mod <- list()
  r2_comp <- r2x_comp <- q2_comp <- aucs_tr <- aucs_te <- array()
  overfitted <- FALSE
  nc <- 1

  while (!overfitted) {
    mod_cv <- if (nc == 1) {
      .plsComponentCv(Xcs, Y = Y_, cv_sets, nc, mod.cv = NULL) # , it_max=800, eps=1e-6
    } else {
      .plsComponentCv(NA, Y = Y_, cv_sets, nc, mod.cv = mod_cv)
    }

    preds_test <- .extMeanCvFeat_alt(mod_cv, feat = "y_pred_test", cv_type = cv$method, model_type = type)
    preds_train <- .extMeanCvFeat_alt(mod_cv, feat = "y_pred_train", cv_type = cv$method, model_type = type)

    perf <- .evaluatePlsPerformance(type, Y_, Y_vec, Ycs, preds_test, preds_train, cv$method, nc)
    aucs_te[nc] <- perf$auc_te
    aucs_tr[nc] <- perf$auc_tr
    r2_comp[nc] <- perf$r2
    q2_comp[nc] <- perf$q2

    if (nc == 1) {
      full_mod[[1]] <- .nipPlsCompRcpp(Xcs, Ycs, it_max=800, eps=1e-6)
    } else {
      Xinp <- full_mod[[nc - 1]]$x_res
      full_mod[[nc]] <- .nipPlsCompRcpp(Xinp, Ycs, it_max=800, eps=1e-6)
    }

    r2x_comp[nc] <- .r2(Xcs, full_mod[[nc]]$t_x %*% full_mod[[nc]]$p_x, NULL)
    overfitted <- .evalFit(type, q2_comp, aucs_te, maxPCo)
    if (!overfitted) nc <- nc + 1
  }

  list(
    full_mod = full_mod[seq_len(nc)],
    r2_comp = r2_comp,
    q2_comp = q2_comp,
    r2x_comp = r2x_comp,
    aucs_tr = aucs_tr,
    aucs_te = aucs_te,
    n_components = nc
  )
}

.evaluatePlsPerformance <- function(type, Y_, Y_vec, Ycs, preds_test, preds_train, cv_method, nc) {
  is_multi_Y <- grepl('mY', type)
  is_DA <- grepl('DA', type)

  result <- list(r2 = NA, q2 = NA, auc_te = NA, auc_tr = NA)

  if (grepl('MC', cv_method) || grepl('k-fold', cv_method)) {
    if (is_DA) {
      if (is_multi_Y) {
        colnames(preds_test) <- colnames(Y_)
        colnames(preds_train) <- colnames(Y_)
        result$auc_te <- multiclass.roc(response = Y_vec, predictor = preds_test, quiet = TRUE)$auc
        result$auc_tr <- multiclass.roc(response = Y_vec, predictor = preds_train, quiet = TRUE)$auc
      } else {
        if (length(dim(preds_test)) == 2) {pte <- preds_test[, 1]} else {pte <- preds_test}
        result$auc_te <- roc(response = Y_vec, predictor = pte, quiet = TRUE)$auc
        if (length(dim(preds_train)) == 2) {ptr <- preds_train[, 1]} else {ptr <- preds_train}
        result$auc_tr <- roc(response = Y_vec, predictor = ptr, quiet = TRUE)$auc
      }
    } else {
      result$r2 <- .r2(Ycs, preds_train, NULL)
      result$q2 <- .r2(Ycs, preds_test, .tssRcpp(Ycs) / ncol(Ycs))
    }
  }
  result
}

.finalisePlsModel <- function(mod_list, XcsTot, YcsTot, msd_x, msd_y, type, pars, summary, Y) {
  t_x <- do.call(cbind, lapply(mod_list, `[[`, "t_x"))
  p_x <- do.call(rbind, lapply(mod_list, `[[`, "p_x"))
  w_x <- do.call(rbind, lapply(mod_list, `[[`, "w_x"))
  b   <- do.call(rbind, lapply(mod_list, `[[`, "b"))
  w_y <- do.call(rbind, lapply(mod_list, `[[`, "w_y"))

  Xpred <- t_x %*% p_x
  E <- XcsTot - Xpred
  Y_res <- YcsTot - (mod_list[[length(mod_list)]]$t_y %*% t(mod_list[[length(mod_list)]]$p_y))

  new("PLS_metabom8",
      type = type,
      t_pred = t_x,
      p_pred = p_x,
      w_pred = w_x,
      betas_pred = drop(b),
      Qpc = w_y,
      t_pred_cv = t_x,
      nPC = length(mod_list),
      summary = summary,
      Y_res = Y_res,
      X_res = E,
      X_mean = msd_x[[2]],
      X_sd = msd_x[[1]],
      Y_mean = msd_y[[2]],
      Y_sd = msd_y[[1]],
      Parameters = pars,
      X = XcsTot,
      X_scaled = XcsTot,
      Y = list(ori = Y, dummy = YcsTot)
  )
}


# pls <- function(X, Y,
#                 center = TRUE,
#                 scale = "UV",
#                 cv = list(method = "k-fold_stratified", k = 7, split = 2 / 3),
#                 maxPCo = 5,
#                 plotting = TRUE) {
#
#   if (is.data.frame(X)) X <- as.matrix(X)
#   if (!is.logical(center)) stop("Invalid 'center' parameter.")
#
#   sc_num <- switch(scale,
#                    "None" = 0,
#                    "UV" = 1,
#                    "Pareto" = 2,
#                    NULL
#   )
#   if (is.null(sc_num)) stop("Invalid 'scale' parameter. Use 'None', 'UV', or 'Pareto'.")
#
#   x_check <- .checkXclassNas(X)
#   y_check <- .checkYclassNas(Y)
#   xy_check <- .checkDimXY(X, y_check[[1]])
#   type <- y_check[[3]]
#
#   msd_y <- .sdRcpp(y_check[[1]])
#   msd_x <- .sdRcpp(X)
#   XcsTot <- .scaleMatRcpp(X, 0:(nrow(X) - 1), center = TRUE, scale_type = sc_num)[[1]]
#   YcsTot <- .scaleMatRcpp(y_check[[1]], 0:(nrow(y_check[[1]]) - 1), center = TRUE, scale_type = sc_num)[[1]]
#
#   tssx <- .tssRcpp(XcsTot)
#   tssy <- .tssRcpp(YcsTot) / ncol(YcsTot)
#
#   cv$cv_sets <- .cvSetsMethod(Y = y_check[[1]], type = type, method = cv$method, k = cv$k, split = cv$split)
#   # cv$k <- length(cv$cv_sets)
#
#   r2_comp <- r2x_comp <- q2_comp <- aucs_tr <- aucs_te <- array()
#   nc <- 1
#   overfitted <- FALSE
#   full_mod = list()
#
#   while (overfitted == FALSE) {
#
#     if (nc == 1) {
#       tt <- .plsComponentCv(X, Y = y_check[[1]], cv$cv_sets, nc, mod.cv = NULL)
#     } else {
#       tt <- .plsComponentCv(X = NA, Y = y_check[[1]], cv$cv_sets, nc, mod.cv = tt)
#     }
#
#     # preds_test <- .extMeanCvFeat_pls(cv_obj = tt, feat = "y_pred_test", cv_type = cv$method, model_type = type)
#     # preds_train <- .extMeanCvFeat_pls(tt, feat = "y_pred_train", cv_type = cv$method, model_type = type)
#     # in-sample predictions
#     preds_test <- .extMeanCvFeat_alt(cv_obj = tt, feat = "y_pred_test", cv_type = cv$method, model_type = type)
#     # out-of-sample predictions
#     preds_train <- .extMeanCvFeat_alt(tt, feat = "y_pred_train", cv_type = cv$method, model_type = type)
#
#     if(length(preds_test)!=length(y_check[[1]])){
#       error('outcome of cv routine: length of Y-predictions of in-sample and/or out-of-sample set is/are uneuql to length of Y')
#     }
#
#     if(length(preds_test)!=length(preds_train)){
#       error('outcome of cv routine: in-sample and out-of-sample sets have uneuql length of Y-predictions')
#     }
#
#     switch(strsplit(cv$method, '_')[[1]][1],
#            'MC' = {
#              if (grepl('DA', type)) {
#                if (grepl('mY', type)) {
#                  pred_mean <- preds_test
#                  colnames(pred_mean) <- colnames(y_check[[1]])
#                  aucs_te[nc] <- multiclass.roc(response = Y, predictor = pred_mean, quiet=TRUE)$auc
#
#                  pred_tr_mean <- preds_train
#                  colnames(pred_tr_mean) <- colnames(y_check[[1]])
#                  aucs_tr[nc] <- multiclass.roc(response = Y, predictor = pred_tr_mean, quiet=TRUE)$auc
#                } else {
#                  aucs_te[nc] <- roc(response = Y, predictor = preds_test, quiet=TRUE)$auc
#                  aucs_tr[nc] <- roc(response = Y, predictor = preds_train, quiet=TRUE)$auc
#                }
#              } else {
#                r2_comp[nc] <- .r2(YcsTot, preds_train, NULL)
#                q2_comp[nc] <- .r2(YcsTot, preds_test, tssy)
#              }
#            },
#            'k-fold' = {
#              if (grepl('DA', type)) {
#                if (grepl('mY', type)) {
#                  colnames(preds_test) <- colnames(y_check[[1]])
#                  aucs_te[nc] <- multiclass.roc(response = Y, predictor = preds_test, quiet=TRUE)$auc
#                  preds_te <- preds_train
#                  colnames(preds_te) <- colnames(y_check[[1]])
#                  aucs_tr[nc] <- multiclass.roc(response = Y, predictor = preds_te, quiet=TRUE)$auc
#                } else {
#                  aucs_te[nc] <- roc(response = Y, predictor = preds_test[, 1], quiet=TRUE)$auc
#                  aucs_tr[nc] <- roc(response = Y, predictor = preds_train, quiet=TRUE)$auc
#                }
#              } else {
#                r2_comp[nc] <- .r2(YcsTot, preds_train, NULL)
#                q2_comp[nc] <- .r2(YcsTot, preds_test, tssy)
#              }
#            }
#     )
#
#
#     overfitted <- .evalFit(type, q2_comp, aucs_te, maxPCo)
#
#     if (nc == 1){
#       Xinp <- XcsTot
#       full_mod[[1]] = .nipPlsCompRcpp(X = Xinp, Y = YcsTot)
#       r2xc <- .r2(Xinp, (full_mod[[1]]$t_x %*% full_mod[[1]]$p_x), NULL)
#     }else{
#       c_idx = length(full_mod)+1
#       Xinp <- full_mod[[length(full_mod)]]$x_res
#       full_mod[[c_idx]] = .nipPlsCompRcpp(X = Xinp, Y = YcsTot)
#       r2xc <- .r2(XcsTot, (full_mod[[c_idx]]$t_x %*% full_mod[[c_idx]]$p_x), NULL)
#     }
#     r2x_comp[nc] <- r2xc
#
#     if (!overfitted) {nc <- nc + 1}
#   }
#
#   message(sprintf("A PLS-%s model with %d components was fitted.", type, nc))
#   full_mod = full_mod[1:nc]
#
#
#   mod_summary <- .plsModelCompSummary(type = type, r2x_comp = cumsum(r2x_comp), r2_comp = r2_comp, q2_comp = q2_comp, aucs_te = aucs_te, aucs_tr = aucs_tr, cv = cv)
#   if (plotting) plot(mod_summary[[2]])# + labs(x = 'Component(s)'))
#
#   # final_model <- .nipPlsCompRcpp(X = XcsTot, Y = YcsTot)
#   t_x = do.call(cbind, lapply(full_mod, function(x){x$t_x}))
#   p_x = do.call(rbind, lapply(full_mod, function(x){x$p_x}))
#   w_x = do.call(rbind, lapply(full_mod, function(x){x$w_x}))
#   b = do.call(rbind, lapply(full_mod, function(x){x$b}))
#   w_y = do.call(rbind, lapply(full_mod, function(x){x$w_y}))
#
#   Xpred <- t_x %*% p_x
#   E <- XcsTot - Xpred
#   Y_res <- YcsTot - (full_mod[[length(full_mod)]]$t_y %*% t(full_mod[[length(full_mod)]]$p_y))
#
#   pars <- list(
#     center = center,
#     scale = scale,
#     nPredC = 1,
#     nOrthC = nc,
#     maxPCo = maxPCo,
#     cv = cv,
#     tssx = tssx,
#     tssy = tssy,
#     r2x_comp = r2x_comp
#   )
#
#   new("PLS_metabom8",
#       type = type,
#       t_pred = t_x,
#       p_pred = p_x,
#       w_pred = w_x,
#       betas_pred = drop(b),
#       Qpc = w_y,
#       t_pred_cv = t_x,
#       nPC = nc,
#       summary = mod_summary[[1]],
#       Y_res = Y_res,
#       X_res = E,
#       X_mean = msd_x[[2]],
#       X_sd = msd_x[[1]],
#       Y_mean = msd_y[[2]],
#       Y_sd = msd_y[[1]],
#       Parameters = pars,
#       X = X,
#       X_scaled = XcsTot,
#       Y = list(ori = Y, dummy = y_check[[1]])
#   )
# }
