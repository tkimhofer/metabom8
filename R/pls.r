#' @title Partial Least Squares (PLS)
#' @export
#' @description This function fits Partial Least Squares (O-PLS) models for regression (R) or discriminant analysis (DA). The optimal number of components is determined automatically using automated stop criteria based on statistical cross validation indices.
#' @param X Numeric input matrix or dataframe (usually measurements obtained through NMR spectroscopy or mass spectrometry) with each row representing an observation and each column a metabolic feature.
#' @param Y Response vector or matrix with same length or number of columns than rows in X, respectively. Y have multiple columns.
#' @param t_pred Parameter specifying the maximum number of predictive components (needed only for multi-factor Y)
#' @param center Logical value (TRUE or FALSE) indicating if features should be mean centered.
#' @param scale Desired scaling method (currently only no or unit variance scaling (UV) implemented).
#' @param cv Named list of cross-validation paramters to derive the optimal number of components: method, one of 'k-fold', 'k-fold_stratified', 'MC', 'MC_balanced' (see Details), split: fraction of observations used for model training, k: k-fold paramter, ie., the number of cross-validation sets. The latter depends on the number of observations in X but typically takes a value between 3 and 9.
#' @param plotting Logical value (TRUE or FALSE) indicating if model parameters (R2X, Q2, etc) should be visualised once the model is trained.
#' @param maxPCo The maximum number of orthogonal components (in case stop criteria fail).
#' @details Models are fully statistically validated, currently only k-fold cross validation (CV) and class-balanced k-fold cross validation is implemented. Further extensions, e.g. Monte-Carlo CV, are work in progress. Although the algorithm accepts three and more levels as Y, model interpretation is more straightforward for pairwise group comparisons.
#' @references Geladi, P and Kowalski, B.R. (1986), Partial least squares and regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1-17.
#' @return This function returns a PLS-metabom8 object (S4).
#' @examples
#' data(covid)
#' model=opls(X, Y=an$type)
#' plotscores(model, an=list(Class=an$type, Clinic=an$hospital, id=1:nrow(an)), pc=c(1,2))
#' @author \email{torben.kimhofer@murdoch.edu.au}
#' @importFrom graphics plot
#' @importFrom methods getSlots new representation setClass
#' @importFrom stats cov sd var
#' @importFrom utils methods
#' @importFrom ggplot2 ggplot aes aes_string scale_fill_manual scale_y_continuous theme_bw labs scale_x_discrete scale_alpha theme element_blank element_line element_rect geom_bar element_text
#' @importFrom pROC roc multiclass.roc

# data(iris)
# X=as.matrix(iris[,1:4])
# Y=cbind(as.character(iris[,5]))
# center = T
# scale = "UV"
# plotting = T
# maxPCo = 5
# cv=list(method='k-fold', split=2/3, k=3)

# load('/Volumes/Torben_1 1/Rproj/KMeyer/NMR_cleanSpectra.Rdata')
# idx=which(!is.na(an$HBP30))
# X=X[idx,]
# Y=factor(an$HBP30[idx])

#mod=pls(X,factor(iris$Species))
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


pls <- function(X, Y, center = TRUE, scale = 'UV', cv=list(method='k-fold_stratified', k=7, split=2/3), maxPCo = 5, plotting = TRUE) {

  {

    if(is.data.frame(X)) X=as.matrix(X)

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

    msd_y<-.sdRcpp(y_check[[1]])
    msd_x<-.sdRcpp(X)
    XcsTot<-.scaleMatRcpp(X, 0:(nrow(X)-1), center=TRUE, scale_type = sc_num)[[1]]
    YcsTot<-.scaleMatRcpp(y_check[[1]], 0:(nrow(y_check[[1]])-1), center=TRUE, scale_type = sc_num)[[1]]

    tssx<-.tssRcpp(XcsTot)
    tssy<-.tssRcpp(YcsTot)/ncol(YcsTot)
  }

  cv$cv_sets<-.cvSetsMethod(Y=y_check[[1]], type = type, method = cv$method, k = cv$k,  split=cv$split)
  cv$k <- length(cv$cv_sets)


  r2_comp <- r2x_comp <- q2_comp <- aucs_tr <- aucs_te <- array()
  nc <- 1
  overfitted <- FALSE

  while (overfitted == FALSE) {

    if(nc == 1){ tt<-.plsComponentCv(X, Y=y_check[[1]], cv$cv_sets, nc,  mod.cv=NULL)
    #browser()
    }else{
      tt<-.plsComponentCv(X=NA, Y=y_check[[1]], cv$cv_sets, nc,  mod.cv=tt)
    }

    # extract cv stats -> if Y is multi-column, fct .extrMeanCVFeat function is likely not working
    preds_test<-.extMeanCvFeat_pls(cv_obj = tt, feat = 'y_pred_test', cv_type = cv$method, model_type=type)
    preds_train<-.extMeanCvFeat_pls(tt, feat = 'y_pred_train', cv_type = cv$method, model_type=type)



    # distinguish between MULTI-Y/SINGLE-Y, as well as MCCV vs k-fold


    # calculate R2 and auroc for cv compounds
    switch(strsplit(cv$method, '_')[[1]][1],
           'MC'={
             if( grepl('DA', type) ) {
               if ( grepl('mY', type) ) {                           # multi Y DA using MCCV
                 pred_mean=preds_test[1,,]
                 colnames(pred_mean)=colnames(y_check[[1]])
                 mod <- suppressWarnings(multiclass.roc(response = factor(Y), predictor = pred_mean))
                 aucs_te[nc] <- mod$auc
                 pred_tr_mean=preds_train[1,,]
                 colnames(pred_tr_mean)=colnames(y_check[[1]])
                 mod <- suppressWarnings(multiclass.roc(response = factor(Y), predictor = pred_tr_mean))
                 aucs_tr[nc] <- mod$auc
               }else{                                             # single Y DA using MCCV
                 # need the same decision boundary
                 mod <- suppressWarnings(roc(response = Y, predictor = preds_test[1,], quiet = TRUE))
                 aucs_te[nc] <- mod$auc
                 mod <- suppressWarnings(roc(response = Y, predictor = preds_train[1,], quiet = TRUE))
                 aucs_tr[nc] <- mod$auc
               }
             } else{
               if (grepl('mY', type)) {                         # multi Y regression using MCCV
                 r2_comp[nc] <- .r2(YcsTot, preds_train, NULL)
                 q2_comp[nc] <- .r2(YcsTot, preds_test, tssy)
               }else{

                 # single Y regression using MCCV
                 r2_comp[nc] <- .r2(YcsTot, preds_train[1,], NULL)
                 q2_comp[nc] <- .r2(YcsTot, preds_test[1,], tssy)
               }
             }
           },
           'k-fold'={
             if( grepl('DA', type) ) {
               if ( grepl('mY', type) ) {
                 colnames(preds_test)=colnames(y_check[[1]])
                 mod <- suppressWarnings(multiclass.roc(response = Y, predictor = apply(preds_test, 2, as.numeric)))
                 aucs_te[nc] <- mod$auc
                 preds_te = preds_train # this extracts mean values
                 colnames(preds_te)=colnames(y_check[[1]])
                 mod <- suppressWarnings(multiclass.roc(response = Y, predictor = preds_te, quiet = TRUE))
                 aucs_tr[nc] <- mod$auc
               }else{
                 mod <- suppressWarnings(roc(response = Y, predictor = preds_test, quiet = TRUE))
                 aucs_te[nc] <- mod$auc
                 mod <- suppressWarnings(roc(response = Y, predictor = preds_train[1,], quiet = TRUE))
                 aucs_tr[nc] <- mod$auc
               }
             }else{
               if (grepl('MC', type)) {
                 r2_comp[nc] <- .r2(YcsTot, preds_train[1,,], NULL)
                 q2_comp[nc] <- .r2(YcsTot, preds_test[4,], tssy)
               }else{
                 r2_comp[nc] <- .r2(YcsTot, preds_train[1,], NULL)
                 q2_comp[nc] <- .r2(YcsTot, preds_test[,1], tssy)
               }
             }
           }
    )
    overfitted <- .evalFit(type, q2_comp, aucs_te, maxPCo)

    if (nc >= 2) {
      if(nc==2){
        pls_full <- .nipPlsCompRcpp(X = XcsTot, Y = YcsTot)
        r2x_comp[nc-1]<-.r2(XcsTot, pls_full$t_x %*% pls_full$p_x, NULL)
      }else{
        x_sq_bef=pls_full$x_res
        pls_full <- .nipPlsCompRcpp(X = pls_full$x_res, Y = YcsTot)
        r2x_comp[nc-1]<-.r2(x_sq_bef, pls_full$t_x %*% pls_full$p_x, NULL)
      }
    }



    if (overfitted == TRUE) {
      message(paste0("A PLS-", type, " model with ", nc - 1, " components was fitted."))
      if(grepl('MC', cv$method)){
        # this means output will be mean, sd and %cv
        if(grepl('mY', type)){
          pls_full$t_cv <- t(t(.extMeanCvFeat(tt, 't_xp', cv_type = cv$method, model_type=type)[,1:(nc-1)]))
          #pred_comp$t_o_cv <- t(t(.extMeanCvFeat(tt, 't_xo', cv_type = cv$method, model_type=type)[1,,-nc]))
        }else{
          #extract mean value
          pls_full$t_cv <- t(t(.extMeanCvFeat(tt, 't_xp', cv_type = cv$method, model_type=type)[1,]))
          #pred_comp$t_o_cv <- t(t(.extMeanCvFeat(tt, 't_xo', cv_type = cv$method, model_type=type)[1,]))
        }
      }else{
        #pred_comp$t_cv<-matrix(.extMeanCvFeat(tt, 't_xp', cv_type = cv$method, model_type=type), ncol=1)
        pls_full$t_cv <- t(t(.extMeanCvFeat(tt, 't_xp', cv_type = cv$method, model_type=type)[,1:(nc-1)]))
      }
      # pred_comp$t_cv<-t(t(.extMeanCvFeat(tt, 't_xp', cv_type = cv$method, model_type=type)[1,,]))
      #pred_comp$t_o<-t_orth
      # pred_comp$t_o_cv<-.extMeanCvFeat(tt, 't_xo', cv_type = cv$method, model_type=type)[1,,-nc]
      #if(!is.matrix(pred_comp$t_o_cv)){pred_comp$t_o_cv=t(t(pred_comp$t_o_cv))}
    }
    nc <- nc + 1
  }

  # xfilt_obsolete <- .nipOplsRcpp(X = opls_filt$X_res, Y = YcsTot)
  m_summary=.orthModelCompSummary(type = type, r2x_comp = c(r2x_comp, NA), r2_comp = r2_comp, q2_comp = q2_comp, aucs_te = aucs_te, aucs_tr=aucs_tr, cv = cv)
  if (plotting == TRUE) { suppressWarnings(print(m_summary[[2]]+labs(x='Predictive Components', titel='PLS-DA'))) }

  Xpred <- pls_full$t_x %*% pls_full$p_x
  E <- XcsTot - (Xpred)  # this is for calculation of DModX
  Y_res = YcsTot -  (pls_full$t_y %*% t(pls_full$p_y)) # this is scaled and mean centered Y residual

  pars <- list(
    'center'= center,
    'scale'= scale,
    'nPredC'= 1,
    'nOrthC'= nc-1,
    'maxPCo'= maxPCo,
    'cv'= cv,
    'tssx'= tssx,
    'tssy'= tssy
  )


  mod_pls <- new("PLS_metabom8",
                  type = type,
                  t_pred = pls_full$t_x,
                  p_pred = pls_full$p_x,
                  w_pred = pls_full$w_x,
                  betas_pred = drop(pls_full$b),
                  Qpc = pls_full$w_y,
                  t_pred_cv = pls_full$t_cv,
                  nPC = nc - 1,
                  summary = m_summary[[1]],
                  Y_res = Y_res,
                  X_res = E,
                  X_mean =  msd_x[[2]],
                  X_sd = msd_x[[1]],
                  Y_mean = msd_y[[2]],
                  Y_sd = msd_y[[1]],
                  Parameters = pars,
                  X = X,
                  X_scaled = XcsTot,
                  Y = list(ori=Y, dummy=y_check[[1]])
  )

  return(mod_pls)
}




















































































