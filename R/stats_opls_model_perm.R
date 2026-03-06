#' @title OPLS Model Validation via Y-Permutation
#'
#' @description
#' Performs Y-permutation tests to assess the robustness of OPLS models by comparing model performance statistics on real vs. permuted response labels.
#'
#' @param smod An OPLS model object of class \code{OPLS_metabom8}.
#' @param n Integer. Number of permutations to perform.
#' @param plot Logical. If \code{TRUE}, generates a visual summary of permutation statistics.
#' @param mc Logical. If \code{TRUE}, enables multicore processing (currently not implemented).
#'
#' @return A \code{data.frame} with model metrics (e.g., R2, Q2, AUROC) from both permuted and original models.
#'
#' @details
#' Each permutation shuffles the response labels and fits a new OPLS model. The function captures model statistics (R2, Q2, AUC) to compare against the non-permuted model. This helps determine whether the original model performance is better than expected by chance.
#'
#' @references
#' Wiklund, S. et al. (2008). A Tool for Improved Validation of OPLS Models. *Journal of Chemometrics*, 22(11–12), 594–600.
#'
#' @family model_validation
#' @importFrom ggplot2 ggplot aes_string theme_bw labs geom_point element_rect
#' @importFrom reshape2 melt
#' @importFrom stats cor median
#' @importFrom parallel detectCores
#' @importFrom progress progress_bar
#' @importFrom rlang .data
# examples
# data(covid)
# X <- covid$X
# an <- covid$an
#
# model <- opls(X, Y = an$type)
# perm_results <- opls_perm(model, n = 10)
opls_perm <- function(smod, n = 10, plot = TRUE, mc = FALSE) {
  if (!inherits(smod, "m8_model") || smod@engine != "opls") {
    stop("Input must be an m8_model instance with opls engine.")
  }

  nc_o <- smod@ctrl$ncomp_selected - 1
  type <- smod@ctrl$type

  # full X/Y
  Y <- smod@fit$Y
  XcsTot <- smod@fit$X_prepped

  cv <- smod@cv

  pb <- progress::progress_bar$new(
    format = "  Permutations [:bar] :percent ETA: :eta",
    total = n,
    clear = FALSE,
    width = 60
  )

  perms <- lapply(seq_len(n), function(i) {
    pb$tick()

    Y_perm <- Y[sample(seq_len(nrow(Y))), ]
    stats  <- .permYmod(XcsTot, Y_perm, cv, type, nc_o)

    list(
      ind = stats,
      r   = mean(abs(diag(cor(Y_perm, Y))))
    )
  })

  out_df <- as.data.frame(do.call(rbind, lapply(perms, function(p){
    ind_vec <- unlist(p$ind)
    r_val <- unlist(p$r)
    c(ind_vec, r = p$r)
  } )))
  colnames(out_df) <- c("q2_comp", "r2_comp", "aucs_te", "aucs_tr", "r")
  out_df$model <- paste0("perm_", seq_len(nrow(out_df)))

  idx <- smod@ctrl$nc_orth

  add <- as.data.frame(list(q2_comp = smod@ctrl$q2[idx], r2_comp = smod@ctrl$r2[idx], aucs_te = smod@ctrl$aucs_te[idx], aucs_tr = smod@ctrl$aucs_tr[idx], r=1, model='non-permuted'))
  out_df <- rbind(out_df, add)

  idx_rm <- vapply(out_df, function(x) all(is.na(x)), logical(1))
  out_df <- out_df[, !idx_rm, drop = FALSE]
  out_df$r_abs <- abs(as.numeric(out_df$r))

  if (plot) {

    signed_log1p <- function(x) sign(x) * log1p(abs(x))
    out_df$q2_plot <- signed_log1p(out_df$q2_comp)

    if (grepl('DA', type)){
      dd <- reshape2::melt(out_df[,c('aucs_tr', 'aucs_te', 'r_abs', 'model')], id.vars = c("model", "r_abs"))
    }else{
      dd <- reshape2::melt(out_df[, c('r2_comp', 'q2_plot', 'r_abs', 'model')], id.vars = c("model", "r_abs"))
    }

    dd$r_abs <- as.numeric(dd$r_abs)
    dd$value <- as.numeric(dd$value)
    dd$variable <- as.character(dd$variable)

    map_parameter <- c(aucs_tr = "AUROC_train", aucs_te = "AUROC_test", r2_comp = "R2", q2_comp = "Q2", q2_plot = "Q2_s_log1p")
    dd$variable <- map_parameter[match(dd$variable, names(map_parameter))]

    y_lwer <- min(dd$value, na.rm = TRUE)

    g <- ggplot(dd, aes(x = .data$r_abs, y = .data$value, colour = .data$variable)) +
      geom_point() +
      theme_bw() +
      labs(
        x = "|cor(Y_perm, Y)|",
        y = NULL,
        colour = "Metric",
        caption = paste0("O-PLS-", type, ": ", n, " Y-permutations.")
      ) +
      scale_y_continuous(limits = c(y_lwer, 1)) +
      theme(
        legend.position = c(0.95, 0.05),  # x, y in npc coordinates
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = "white", colour = "black")
      )

    plot(g)
  }

  .perm_test_from_table(out_df)
}


#' Permutation-test summary from an opls_perm out_df table
#'
#' @description
#' Takes the `out_df` you showed (perm rows + one "non-permuted" row) and returns
#' observed values, permutation distributions, and permutation p-values.
#'
#' @param out_df data.frame with columns like q2_comp, r2_comp, aucs_te, aucs_tr,
#'   model (with one row "non-permuted"), and optionally r/r_abs.
#' @param observed_label character. Label used in `model` for the observed row.
#' @param alternative character. "greater" (default) tests obs > perm; "less" tests obs < perm.
#' @param add_one logical. If TRUE, uses (1 + count)/(B + 1) p-value correction.
#' @param na_rm logical. Drop NA values before computing p-values.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{observed}: named numeric vector of observed metrics
#'   \item \code{perm}: list of numeric vectors for each metric
#'   \item \code{p_value}: named numeric vector of permutation p-values
#'   \item \code{B}: number of permutations used
#' }
.perm_test_from_table <- function(out_df,
                                 observed_label = "non-permuted",
                                 alternative = c("greater", "less"),
                                 add_one = TRUE,
                                 na_rm = TRUE) {
  alternative <- match.arg(alternative)

  if (!is.data.frame(out_df)) stop("out_df must be a data.frame.")
  if (!("model" %in% names(out_df))) stop("out_df must have a 'model' column.")

  obs_idx <- which(out_df$model == observed_label)
  if (length(obs_idx) != 1L) {
    stop("Expected exactly one observed row with model == '", observed_label, "'. Found: ", length(obs_idx))
  }

  perm_df <- out_df[out_df$model != observed_label, , drop = FALSE]
  B <- nrow(perm_df)
  if (B < 1L) stop("No permutation rows found (model != '", observed_label, "').")

  drop_cols <- c("model", "r_abs")  # keep r_abs out unless you want it as a metric
  metric_cols <- setdiff(names(out_df), drop_cols)
  metric_cols <- metric_cols[vapply(out_df[metric_cols], is.numeric, logical(1))]

  if (length(metric_cols) == 0L) stop("No numeric metric columns found in out_df.")

  obs <- as.list(out_df[obs_idx, metric_cols, drop = FALSE])
  obs <- unlist(obs, use.names = TRUE)

  perm <- lapply(metric_cols, function(col) perm_df[[col]])
  names(perm) <- metric_cols

  perm_p <- function(perm_vals, obs_val) {
    if (na_rm) perm_vals <- perm_vals[!is.na(perm_vals)]
    if (is.na(obs_val)) return(NA_real_)
    B_eff <- length(perm_vals)
    if (B_eff == 0L) return(NA_real_)

    if (alternative == "greater") {
      k <- sum(perm_vals >= obs_val, na.rm = TRUE)
    } else {
      k <- sum(perm_vals <= obs_val, na.rm = TRUE)
    }

    if (add_one) (k + 1) / (B_eff + 1) else k / B_eff
  }

  pvals <- vapply(metric_cols, function(col) perm_p(perm[[col]], obs[[col]]), numeric(1))
  names(pvals) <- metric_cols

  list(
    observed = obs,
    perm = perm,
    p_value = pvals,
    B = B
  )
}


#' @title OPLS Y-permutation Modeling
#' @description
#' Performs OPLS modeling on permuted Y to estimate cross-validated model performance under the null hypothesis.
#' Extracts predictive performance (R2, Q2, AUROC) from cross-validation for either regression or classification.
#'
#' @param XcsTot Numeric matrix. Input data matrix (samples x features).
#' @param Y Numeric matrix. Response variable (numeric or dummy-coded).
#' @param cv List. Cross-validation parameters, including \code{method} and \code{cv_sets}.
#' @param type Character. Model type: \code{"DA"}, \code{"R"}, possibly with \code{"-mY"} suffix for multi-column Y.
#' @param nc_o Integer. Number of orthogonal components to be fitted.
#'
#' @return List with elements:
#' \itemize{
#'   \item \code{r2_comp}: Numeric, R2 statistic for test data.
#'   \item \code{q2_comp}: Numeric, Q2 statistic from cross-validation.
#'   \item \code{aucs_tr}: Numeric, AUROC on training data (for classification only).
#'   \item \code{aucs_te}: Numeric, AUROC on test data (for classification only).
#' }
#'
#' @details
#' When \code{nc_o > 1}, the function fits additional orthogonal components sequentially, updating the model object.
#' Classification performance is computed using AUROC (via \code{pROC}). Regression performance uses R2 and Q2.
#'
#' @importFrom pROC roc multiclass.roc
#' @keywords internal
.permYmod <- function(XcsTot, Y, cv, type, nc_o) {

  if(is.vector(Y) && is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
  }

  Ycs_fold <- lapply(cv@train, function(idc)
    .scaleMatRcpp(Y, idc - 1, 0, 0)[[1]]
  )

  n  <- nrow(XcsTot)
  p  <- ncol(XcsTot)
  q <- ncol(Y)

  acc <- list(
    sum_test  = matrix(0, n, q),
    n_test    = integer(n),
    sum_train = matrix(0, n, q),
    n_train   = integer(n)
  )

  for (i in seq_len(nc_o)) {
    res <- if (nc_o == 1)
      .oplsComponentCv(
        X        = XcsTot,
        cv.set   = cv@train,
        Ycs_fold = Ycs_fold,
        nc       = nc_o,
        mod.cv   = NULL,
        acc      = acc
      )
    else
      .oplsComponentCv(
        X        = NULL,
        cv.set   = cv@train,
        Ycs_fold = Ycs_fold,
        nc       = nc_o,
        mod.cv   = tt,
        acc      = acc
      )
    tt  <- res$mod.cv
  }

  tt  <- res$mod.cv
  acc <- res$acc

  preds_test <- acc$sum_test
  preds_test[acc$n_test > 0,] <- preds_test[acc$n_test > 0,] / acc$n_test[acc$n_test > 0]

  preds_train <- acc$sum_train
  preds_train[acc$n_train > 0,] <- preds_train[acc$n_train > 0,] / acc$n_train[acc$n_train > 0]

  tssy <- .tssRcpp(Y) / ncol(Y)

  perf <- .evalComponentPerformance(
    cv_inst      = cv,
    type         = type,
    is_multi_Y   = FALSE,
    Y            = Y,
    preds_train  = preds_train,
    preds_test   = preds_test,
    class_memb   = NULL,
    YcsTot       = Y,
    tssy         = tssy
  )

  list(q2 = perf$q2, r2 = perf$r2, aucs_te = perf$aucs_te, aucs_tr = perf$aucs_tr)

}

