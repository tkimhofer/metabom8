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
#' @importFrom ggplot2 ggplot aes_string theme_bw labs geom_point
#' @importFrom reshape2 melt
#' @importFrom stats cor median
#' @importFrom parallel detectCores
#' @export
#' @examples
#' data(covid)
#' model <- opls(X, Y = an$type)
#' perm_results <- opls_perm(model, n = 10)
opls_perm <- function(smod, n = 10, plot = TRUE, mc = FALSE) {
  if (!inherits(smod, "OPLS_metabom8")) {
    stop("Input must be an OPLS_metabom8 object.")
  }

  nc_o <- nrow(smod@summary) - 1
  Xs <- smod@X
  Y <- smod@Y$dummy
  cv <- list(cv_sets = smod@Parameters$cv$cv_sets,
             k = length(smod@Parameters$cv$cv_sets),
             method = smod@Parameters$cv$method)
  type <- smod@type

  perms <- lapply(seq_len(n), function(i) {
    Y_perm <- matrix(sample(Y[, 1]), ncol = 1)
    stats <- .permYmod(Xs, Y_perm, cv, type, nc_o)
    list(ind = stats, r = cor(Y[, 1], Y_perm[, 1])[1])
  })

  out_df <- as.data.frame(do.call(rbind, lapply(perms, function(p){
    ind_vec <- unlist(p$ind)
    r_val <- unlist(p$r)
    c(ind_vec, r = p$r)
  } )))
  colnames(out_df) <- c("r2_comp", "q2_comp", "aucs_tr", "aucs_te", "r")
  out_df$model <- paste0("perm_", seq_len(nrow(out_df)))

  real_stats <- .permYmod(Xs, Y, cv, type, nc_o)
  real_row <- as.data.frame(t(c(real_stats, r = 1)))
  real_row$model <- "Non-permuted"
  colnames(real_row) <- colnames(out_df)

  out_df <- rbind(out_df, unlist(real_row))

  # Clean up NA-only columns
  idx_rm <- vapply(out_df, function(x) all(is.na(x)), logical(1))
  out_df <- out_df[, !idx_rm, drop = FALSE]
  out_df$r_abs <- abs(as.numeric(out_df$r))

  if (plot) {


    if (smod@type=='DA'){
      dd <- reshape2::melt(out_df[,c('aucs_te', 'r_abs', 'model')], id.vars = c("model", "r_abs"))
    }
    else{
      dd <- reshape2::melt(out_df[, c('q2_comp', 'r_abs', 'model')], id.vars = c("model", "r_abs"))
    }

    dd$r <- as.numeric(dd$r)
    dd$value <- as.numeric(dd$value)
    dd$variable <- as.character(dd$variable)

    map_parameter <- c(aucs_tr = "AUROC_train", aucs_te = "AUROC_test", r2_comp = "R2", q2_comp = "Q2")
    dd$variable <- map_parameter[match(dd$variable, names(map_parameter))]
    ylab <- if (any(grepl("AU", dd$variable))) "AUROC" else expression(R^2 ~ "/" ~ Q^2)

    g <- ggplot(dd, aes(x = !!sym("r_abs"), y = !!sym("value"), colour = !!sym("variable"))) +
      geom_point() +
      theme_bw() +
      labs(x = "|r|", y = ylab, colour = "Metric",
           caption = paste0("O-PLS-", type, ": ", n, " Y-permutations."))

    plot(g)
  }

  return(out_df)
}
