#PLS


#' @title Performs PLS modelling for each CV set and collates output
#' @param X num matrix X: preproc NMR data
#' @param Y num matrix Y: outcome, dummy matrix in case of categorical outcome
#' @param cv.set list of k elements containing vector of X and Y row-indices representing training set for cv round k
#' @param nc int, max number of orthogonal components to fit
#' @param mod.cv cv parameters
#' @return Named list of collated OPLS data for respective component
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @keywords internal
#' @section
.plsComponentCv=function(X, Y, cv.set, nc,  mod.cv){

  out=lapply(seq_along(cv.set), function(k){

    idc=cv.set[[k]]
    if(nc==1){
      Xcs <- .scaleMatRcpp(X, idc-1, center=TRUE, scale_type = 1)[[1]] # subtract 1 since Rcpp indexing starts at zero
    }else{
      Xcs=mod.cv[[k]]$x_res
    }

    Y_scale <- .scaleMatRcpp(Y, idc-1, center=TRUE, scale_type = 1)
    Ycs <- Y_scale[[1]] # subtract 1 since Rcpp indexing starts at zero
    pred_comp <- .nipPlsCompRcpp(X = Xcs[idc, ], Y = cbind(Ycs[idc, ]))
    Xte=Xcs[-idc, ]
    if(!is.matrix(Xte)){ Xte=matrix(Xte, nrow=1)} # loocv
    pls_pred=.plsPredRcpp(pls_mod = pred_comp,  Xnew=Xte)

    # create list pls_mod, what is nc?
    if(nc==1){

      mod.cv=list(
        t_xp = matrix(NA, nrow=nrow(Y), ncol=1),
        y_pred_train = matrix(NA, nrow=nrow(Y), ncol=ncol(Y)),
        y_pred_test = matrix(NA, nrow=nrow(Y), ncol=ncol(Y)),
        x_res = matrix(NA, nrow=nrow(Y), ncol=ncol(Xcs))
      )

      mod.cv$t_xp[-idc,nc]=pls_pred$t_pred # cv scores
      mod.cv$y_pred_train[idc,]=pred_comp$y_pred
      mod.cv$y_pred_test[-idc,]=pls_pred$y_pred
      mod.cv$x_res[idc,] = pred_comp$x_res
      mod.cv$x_res[-idc,] = pls_pred$Xres

      return(mod.cv)

    }else{

      t_xp=matrix(NA, nrow=nrow(Y), ncol=1)
      t_xp[-idc,1]=pls_pred$t_pred
      mod.cv[[k]]$t_xp=cbind(mod.cv[[k]]$t_xp, t_xp)

      y_pred_add=matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
      y_pred_add[-idc,]=pls_pred$y_pred
      mod.cv[[k]]$y_pred_test=y_pred_add


      y_pred_add1=matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
      y_pred_add1[idc,]=pred_comp$y_pred
      mod.cv[[k]]$y_pred_train=y_pred_add1


      mod.cv[[k]]$x_res[idc,] = pred_comp$x_res
      mod.cv[[k]]$x_res[-idc,] = pls_pred$Xres

      #mod.cv[[k]]$r2x_pred_comp_cv = .r2(opls_filt$X_res, pred_comp$t_x %*% pred_comp$p_x, tssx)

      return(mod.cv[[k]])
    }

  })
  return(out)
}
