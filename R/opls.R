#' @title Orthogonal-Partial Least Squares (O-PLS) modelling
#' @export
#' @description This function is used to fit  Orthogonal-Partial Least Squares (O-PLS) models for regression (R) or discriminant analysis (DA). In the latter case the outcome can have two or more levels.
#' @param X Numeric input matrix or dataframe (usually measurements obtained through NMR spectroscopy or mass spectrometry) with each row representing an observation and each column a metabolic feature.
#' @param Y Response vector or matrix with same length or number of columns than rows in X, respectively.
#' @param t_pred Parameter specifying the maximum number of predictive components (needed only for multi-factor Y)
#' @param center Logical value (TRUE or FALSE) indicating if features should be mean centered.
#' @param scale Desired scaling method (currently only no or unit variance scaling (UV) implemented).
#' @param cv List of cross-validation paramters to derive the optimal number of components: method: 'k-fold', 'k-fold_stratified', 'MC', 'MC_stratified' (see Details), training: fraction of observations used for model training, k=fold paramter). The number of cross-validation sets. This depends on the number of observations in X but typically takes a value between 3 and 9.
#' @param plotting Logical value (TRUE or FALSE) indicating if model parameters (R2X, Q2, etc) should be visualised once the model is trained.
#' @param maxPCo The maximum number of orthogonal components (in case stop criteria fail).
#' @details Models are fully statistically validated, currently only k-fold cross validation (CV) and class-balanced k-fold cross validation is implemented. Further extensions, e.g. Monte-Carlo CV, are work in progress. Although the algorithm accepts three and more levels as Y, model interpretation is more straightforward for pairwise group comparisons.
#' @references Trygg J. and Wold, S. (2002) Orthogonal projections to latent structures (O-PLS). \emph{Journal of Chemometrics}, 16.3, 119-128.
#' @references Geladi, P and Kowalski, B.R. (1986), Partial least squares and regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1-17.
#' @return This function returns an \code{\link{OplsMate}} S4 object.
## #' @seealso \code{\link{OPLS_MetaboMate-class}} \code{\link{dmodx}} \code{\link{plotscores}} \code{\link{plotload}} \code{\link{specload}}
#' @author Torben Kimhofer \email{torben.kimhofer@murdoch.edu.au}
#' @examples
#' data(iris)
#' smod=opls(X=iris[,seq_len(4)], Y=iris$Species)
#' @importFrom graphics plot
#' @importFrom methods getSlots new representation setClass
#' @importFrom stats cov sd var
#' @importFrom utils methods
#' @importFrom ggplot2 ggplot aes aes_string scale_fill_manual scale_y_continuous theme_bw labs scale_x_discrete scale_alpha theme element_blank element_line element_rect geom_bar element_text
#' @importFrom pROC roc multiclass.roc


# source('/Volumes/Torben_1 1/Rproj/MetaboMate/R/cv_sets_method.R')
# source('/Volumes/Torben_1 1/Rproj/MetaboMate/R/add_fct.R')
# source('/Volumes/Torben_1 1/Rproj/MetaboMate/R/eval_fit.R')
# Rcpp::sourceCpp('/Volumes/Torben_1 1/Rproj/MetaboMate/src/sd_mat_rcpp.cpp')
# Rcpp::sourceCpp('/Volumes/Torben_1 1/Rproj/MetaboMate/src/nip_opls_mlevcomp.cpp')
# source('/Volumes/Torben_1 1/Rproj/MetaboMate/R/opls_component_cv.R')

# data(iris)
# X=as.matrix(iris[,1:4])
# Y=cbind(as.character(iris[,5]))
# idx=which(Y %in% c('versicolor',  'virginica' ))
# X=X[idx,]
# Y=Y[idx,]
#
# center = T
# scale = "UV"
# plotting = T
# maxPCo = 5
# cv=list(method='MC_balanced', split=2/3, k=10)

# load('/Volumes/Torben_1 1/Rproj/KMeyer/NMR_cleanSpectra.Rdata')
# idx=which(!is.na(an$HBP30))
# X=X[idx,]
# Y=factor(an$HBP30[idx])

#mod=opls_rcpp(X,factor(iris$Species))
#plotscores(mod, an=list(class=factor(Y)), cv.scores = F)

# Load data
# data(bariatric)

# Declare variables
# X <-bariatric$X.pqn # NMR matrix
# X[,100]=999
# ppm<-bariatric$ppm     # ppm vector
# meta<-bariatric$meta   # Metadata
# an<-bariatric$an       # Sample annotation
# Y=factor(an$Timepoint)
# readBruker('/Volumes/Torben_2/Ariwave/')

# test=opls_rcpp(X, Y, t_pred = 1, center = T, scale = "UV", plotting = T, maxPCo = 5, cv=list(method='k-fold_stratified', split=2/3, k=20))


opls <- function(X, Y, t_pred = 1, center = TRUE, scale = 'UV', cv=list(method='MC_balanced', split=2/3, k=7), maxPCo = 5, plotting = TRUE) {

  {
    if(!is.logical(center)){stop('Check center parameter argument!')}
    sc_num<-switch(scale,
                   'None'={0},
                   'UV'={1},
                   'Pareto'={2})
    if(is.null(sc_num)){stop('Check scale parameter argument!')}

    x_check<-.checkXclassNas(X)
    y_check<-.checkYclassNas(Y)
    xy_check<-.checkDimXY(X, y_check[[1]])
    type<-y_check[[3]]

    msd_y<-.sdMatRcpp(y_check[[1]]);
    msd_x<-.sdMatRcpp(X); # returns list with elements mean and sd
    XcsTot<-.scaleMatRcpp(X, 0:(nrow(X)-1), center=TRUE, scale_type = sc_num)[[1]]
    YcsTot<-.scaleMatRcpp(y_check[[1]], 0:(nrow(y_check[[1]])-1), center=TRUE, scale_type = sc_num)[[1]]

    tssx<-.tssRcpp(XcsTot)
    tssy<-.tssRcpp(YcsTot)/ncol(YcsTot)
  }

  cv$cv_sets<-.cvSetsMethod(Y=y_check[[1]], type = y_check[[3]], method = cv[[1]], k = cv[[3]],  split=cv[[2]])

  r2_comp <- r2x_comp <- q2_comp <- aucs <- r2x_comp_cvr <- array()
  nc <- 1
  enough <- FALSE

  while (enough == FALSE) {
    if(nc == 1){
      tt<-.oplsComponentCv(X, Y=y_check[[1]], cv$cv_sets, nc,  mod.cv=NULL)
    }else{
      tt<-.oplsComponentCv(X=NA, Y=y_check[[1]], cv$cv_sets, nc,  mod.cv=tt)
    }
    # extract cv stats
    preds_test<-.extMeanCvFeat(tt, 'y_pred_test')
    #r2_comp_cvr<-.extMeanCvFeat(cv_obj = tt, feat = 'r2x_pred_comp_cv')
    #r2x_comp_cvr[nc]<-r2_comp_cvr[4,1]
    q2_comp[nc] <-.r2(YcsTot, preds_test[4,], tssy)
    if (type == "DA") {
      if (ncol(YcsTot) == 1) {
        mod <- roc(response = Y, predictor = preds_test[4,], quiet = TRUE)
        aucs[nc] <- mod$auc
      } else {
        mod <- multiclass.roc(response = Y, predictor = preds_test[4,])
        aucs[nc] <- mod$auc
      }
    } else {
      aucs[nc] <- 0
    }
    enough<-.evalFit(type, q2_comp, aucs, maxPCo)
    if (nc >= 2) {
      if(nc==2){
        opls_filt <- .nipOplsRcpp(X = XcsTot, Y = YcsTot)
        t_orth<-opls_filt$t_o
        p_orth<-opls_filt$p_o
      }else{
        opls_filt <- .nipOplsRcpp(X = opls_filt$X_res, Y = YcsTot)
        t_orth<-cbind(t_orth, opls_filt$t_o)
        p_orth<-rbind(p_orth, opls_filt$p_o)
      }
      pred_comp <- .nipPlsCompRcpp(X = opls_filt$X_res, Y = YcsTot)
      r2x_comp[nc-1]<-.r2(opls_filt$X_res, pred_comp$t_x %*% pred_comp$p_x, tssx)
      r2x_comp[nc-1]<-.r2(opls_filt$X_res, pred_comp$t_x %*% pred_comp$p_x, tssx)
    }
    if (enough == TRUE) {
      message(paste0("An O-PLS-", type, " model with 1 predictive and ", nc - 1, " orthogonal components was fitted."))
      pred_comp$t_cv<-t(t(.extMeanCvFeat(tt, 't_xp')[4,]))
      pred_comp$t_o<-t_orth
      pred_comp$t_o_cv<-.extMeanCvFeat(tt, 't_xo')[,-nc]
      if(!is.matrix(pred_comp$t_o_cv)){pred_comp$t_o_cv=t(t(pred_comp$t_o_cv))}
    }
    nc <- nc + 1
  }

  xfilt_obsolete <- .nipOplsRcpp(X = opls_filt$X_res, Y = YcsTot)
  pred_comp_obsolete <- .nipPlsCompRcpp(X = xfilt_obsolete$X_res, Y = YcsTot)
  r2x_comp[nc-1]=.r2(xfilt_obsolete$X_res, pred_comp_obsolete$t_x %*% pred_comp_obsolete$p_x, tssx)
  m_summary=.orthModelCompSummary(type, r2x_comp, r2_comp, q2_comp, aucs)
  if (plotting == TRUE) { plot(m_summary[[2]]) }

  Xorth <- t_orth %*% p_orth
  Xpred <- pred_comp$t_x %*% pred_comp$p_x
  E <- XcsTot - (Xpred + Xorth)  # this is for calculation of DModX
  Y_res = YcsTot -  (pred_comp$t_y %*% pred_comp$p_y)

  pars <- list(
    'center'=center,
    'scale'=scale,
    'nPredC'=1,
    'nOrthC'=nc-1,
    'maxPCo'=maxPCo,
    'cv'=cv,
    'tssx'=tssx,
    'tssy'=tssy
  )

  mod_opls <- new("OPLS_MetaboMate",
                  type = type,
                  t_pred = pred_comp$t_x,
                  p_pred = pred_comp$p_x,
                  w_pred = pred_comp$w_x,
                  betas_pred = drop(pred_comp$b),
                  Qpc = pred_comp$w_y,
                  t_cv = pred_comp$t_cv,
                  t_orth_cv = pred_comp$t_o_cv,
                  t_orth = t_orth,
                  w_orth = p_orth,
                  nPC = nc - 1,
                  summary = m_summary[[1]],
                  X_orth = Xorth,
                  Y_res = Y_res,
                  X_res = E,
                  X_mean =  msd_x[[1]],
                  X_sd = msd_x[[2]],
                  Y_mean = msd_y[[1]],
                  Y_sd =msd_y[[2]],
                  Y_dummy = as.data.frame(y_check[[1]]),
                  Parameters = pars,
                  X = X)

  return(mod_opls)
}






























































































