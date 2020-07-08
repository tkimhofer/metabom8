#' OPLS model predictions
#' @export
#' @description Calculation of OPLS model predictions using new data
#' @aliases pred.opls
#' @param opls_model OPLS model (regression of discriminant analysis) of class \code{OPLS_MetaboMate}.
#' @param newdata NMR data matrix or dataframe with rows representing spectra and identical features in columns as data matrix used to calculate original OPLS model.
#' @param idx_scale int vector, row-indices of newdata used to subselect samples to determine scale and center pars. Recommded: set to NULL: use center and scaling parameters from opsl training data
#' @return Returned is a list with the following elements:
#' \describe{
#' \item{Y_predicted}{Class or numeric outcome predictions for discriminant analysis or regression, repspectively.}
#' \item{t_pred}{Predicted OPLS model scores for predictive component(s).}
#' \item{t_orth}{Predicted OPLS model scores for orthogonal component(s).}
#' \item{t_orth_pca}{Scores of a PCA model (first component) calculated using all predicted OPLS orthogonal component scores - only done when there are more than one orthogonal components in \code{opls_model}.}
#' }
#' @details Class predictions for discriminant analysis are not adjusted for unbalanced sample sizes and therefore, predictions can be biased towards the group with the largest number of samples. The list element \code{t_orth_pca} represent scores of the first principal component of a PCA model caclulated with all orthogonal components, therefore, summarises all orthogonal components into a single one. This can only be done if there are more than one orthogonal components in \code{opls_modelel}, otherwise this list element is \code{NULL}.
#' @references Trygg J. and Wold, S. (2002) Orthogonal projections to latent structures (O-PLS). \emph{Journal of Chemometrics}, 16.3, 119-128.
#' @references Geladi, P and Kowalski, B.R. (1986), Partial least squares and regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1-17.
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @seealso \code{\link{opls}}
#' @export
predict_opls <- function(opls_model, newdata, idx_scale=NULL) {
  if (class(opls_model)[1] != "OPLS_MetaboMate") {
    cat("Error: Model input does not belong to class OPLS_Torben!\n")
    return(NULL)
  }
  # In case of one sample scenario: define X as column vector
  if (is.null(ncol(newdata))) {
    X <- rbind(newdata)
  } else {
    X <- newdata
  }

  if(!is.null(idx_scale)){ # use provided scaling vector
    X=scale_rcpp(X, idx_scale, center = opls_model@Parameters$center, scale_type = opls_model@Parameters$scale)
    #X=MetaboMate:::center_scale(ds, idx_scale, center = cent_scale$center, scale = cent_scale$scale )
  }else{ # use model paramters for scaling

    if(all(!is.null(opls_model@X_mean)) && all(!is.null(opls_model@X_sd))){
      X=t(sapply(seq(ncol(X)), function(i){
        (X[,1] - opls_model@X_mean[i]) / opls_model@X_sd[i]
      } ) )
    }

  }

  # center and scale X<-scale(newdata, center=opls_model@Xcenter, scale=opls_model@Xscale) iteratively remove all orthogonal components from
  # prediction data set
  e_new_orth <- X
  t_orth <- matrix(NA, nrow = nrow(X), ncol = opls_model@nPC)
  for (i in 1:opls_model@nPC) {
    t_orth[, i] <- e_new_orth %*% t(t(opls_model@w_orth[i, ]))/drop(crossprod(t(t(opls_model@w_orth[i, ]))))
    e_new_orth <- e_new_orth - (cbind(t_orth[, i]) %*% t(opls_model@p_orth[i, ]))
  }
  if (opls_model@nPC > 1) {
    pc.orth <- pca(t_orth, pc = 1, method = "ppca", scale = "UV")
    t_orth_pca <- pc.orth@t[, 1]
  } else {
    t_orth_pca <- NULL
  }
  # calc predictive component score predictions and residuals
  t_pred <- e_new_orth %*% t(opls_model@w_pred)
  E_h <- e_new_orth - (t_pred %*% opls_model@p_pred)
  betas <- opls_model@betas_pred
  q_h <- opls_model@Qpc
  res <- matrix(NA, nrow = nrow(X), ncol = ncol(opls_model@t_pred))
  for (i in 1:ncol(opls_model@t_pred)) {
    opts <- t(cbind(betas[i]) %*% t_pred[, i]) %*% rbind(q_h[, i])
    res[, i] <- apply(opts, 1, sum)
  }
  totalPrediction <- apply(res, 1, sum)  # sum over all components
  Y_predicted <- (totalPrediction * opls_model@Yscale) + opls_model@Ycenter
  if (opls_model@type == "DA") {
    levs <- opls_model@Yout
    Y_predicted <- levs$Original[apply(sapply(levs$Numeric, function(x, y = Y_predicted) {
      abs(x - y)
    }), 1, which.min)]
  }
  out <- list("Y_predicted" <- Y_predicted, "t_pred" <- t_pred, "t_orth" <- t_orth, "t_orth_pca" <- t_orth_pca)
  return(out)
}
