#' @title PLS Cross-Validation Component Builder
#' @description Performs Partial Least Squares (PLS) modelling for each cross-validation fold and collates output.
#'
#' @param X Numeric matrix. Predictor matrix (e.g., preprocessed metabolomics data).
#' @param Y Numeric matrix. Response matrix. Should be dummy-coded for classification tasks.
#' @param cv.set List of integer vectors. Each element contains row indices for training in the corresponding fold.
#' @param nc Integer. Number of predictive components to compute.
#' @param mod.cv List or NULL. Existing model CV structure for extending previous components, or NULL for first component.
#'
#' @return A named list of matrices containing predicted values, latent scores, and residuals per fold.
#'
#' @keywords internal
.plsComponentCv <- function(X, Y, cv.set, nc, mod.cv) {
  out <- lapply(seq_along(cv.set), function(k) {
    idc <- cv.set[[k]]  # training indices

    if (nc == 1) {
      Xcs <- .scaleMatRcpp(X, idc - 1, center = TRUE, scale_type = 1)[[1]]
    } else {
      Xcs <- mod.cv[[k]]$x_res
    }

    Ycs <- .scaleMatRcpp(Y, idc - 1, center = TRUE, scale_type = 1)[[1]]

    pred_comp <- .nipPlsCompRcpp(X = Xcs[idc, , drop = FALSE], Y = Ycs[idc, , drop = FALSE], it_max=800, eps=1e-8)
    Xte <- Xcs[-idc, , drop = FALSE]
    if (!is.matrix(Xte)) Xte <- matrix(Xte, nrow = 1)

    pls_pred <- .plsPredRcpp(pls_mod = pred_comp, Xnew = Xte)

    if (nc == 1) {
      mod.cv <- list(
        t_xp = matrix(NA, nrow = nrow(Y), ncol = 1),
        y_pred_train = matrix(NA, nrow = nrow(Y), ncol = ncol(Y)),
        y_pred_test = matrix(NA, nrow = nrow(Y), ncol = ncol(Y)),
        x_res = matrix(NA, nrow = nrow(Y), ncol = ncol(Xcs))
      )

      mod.cv$t_xp[-idc, nc] <- pls_pred$t_pred
      mod.cv$y_pred_train[idc, ] <- pred_comp$y_pred
      mod.cv$y_pred_test[-idc, ] <- pls_pred$y_pred
      mod.cv$x_res[idc, ] <- pred_comp$x_res
      mod.cv$x_res[-idc, ] <- pls_pred$Xres

      return(mod.cv)
    } else {
      t_xp <- matrix(NA, nrow = nrow(Y), ncol = 1)
      t_xp[-idc, 1] <- pls_pred$t_pred
      mod.cv[[k]]$t_xp <- cbind(mod.cv[[k]]$t_xp, t_xp)

      mod.cv[[k]]$y_pred_test <- {
        tmp <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
        tmp[-idc, ] <- pls_pred$y_pred
        tmp
      }

      mod.cv[[k]]$y_pred_train <- {
        tmp <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
        tmp[idc, ] <- pred_comp$y_pred
        tmp
      }

      mod.cv[[k]]$x_res[idc, ] <- pred_comp$x_res
      mod.cv[[k]]$x_res[-idc, ] <- pls_pred$Xres

      return(mod.cv[[k]])
    }
  })

  return(out)
}
