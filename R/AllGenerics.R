#' @title Check for missing values (PLS context)
#' @description This function checks for missing values in Y (NA, NAN, infinite) and establishes analysis type accroding to the class of the input data Y: regression (R) or discriminant analysis (DA). It also converts Y to a matrix.
#' @param Y Input data (uni or multivar) formatted as vector or matrix or data.frame
#' @return List of i) Y matrix, ii) Y levels (empty if numeric Y), iii) Y type (R or DA)
#' @keywords internal
#' @section

# check Y for regression or DA
.checkYclassNas = function(Y) {

    if (any(!is(Y) %in% c("matrix", "data.frame"))) {
        Y = matrix(Y, nrow = length(Y), ncol = 1)
    } else {
        Y = as.matrix(Y)
    }

    if (any(is.na(Y)) | any(is.nan(Y)) | any(is.infinite(Y))) {
        stop("Error: Input Y contains na/nan/inf values.")
    }

    if (is.numeric(Y)) {
        type <- "R"
        Y_out <- .createDummyY(Y)
        Y <- Y_out[[1]]
    } else {
        type <- "DA"
        Y_out <- .createDummyY(Y)
        Y <- Y_out[[1]]
        levs <- unique(apply(Y, 2, function(x) {
            length(unique(x))
        }))
        if (length(levs) == 1 & levs[1] == 1) {
            stop("Error: Input Y has only a single level.")
        }
        message("Performing discriminant analysis.")
    }

    if(ncol(Y)>1){
        type=paste0(type, '-mY')
    }


    return(c(Y_out, type))

}


#' @title Check X matrix (PLS context)
#' @description This function checks for missing values in Y (NA, NAN, infinite), numeric variable format.
#' @param X Input data (uni or multivar) formatted as vector or matrix or data.frame
#' @return NULL if all's fine and throws error otherwise
#' @keywords internal
#' @section
.checkXclassNas = function(X) {

    if (!is.matrix(X) | !is.numeric(X[1, 1])) {
        stop("Input X must be matrix of class numeric")
    }

    if (any(is.na(X)) | any(is.nan(X)) | any(is.infinite(X))) {
        stop("Error: Input X contains na/nan/inf values.")
    }

}


#' @title Check X and Y shapes (PLS context)
#' @description This function checks for matching dimensions (X - Y).
#' @param X X matrix (observations x variables)
#' @param Y Y matrix (observations x variables)
#' @return NULL if all's fine and throws error otherwise
#' @keywords internal
#' @section

.checkDimXY = function(X, Y) {

    if (nrow(X) != nrow(Y)) {
        stop("Error: Dimensions of input X and Y do not match.")
    }

    if (ncol(X) <= ncol(Y)) {
        stop("Error: Number of variables (columns) in X should be higher than in Y.")
    }

}



#' @title Criteria to evaluate number of components
#' @description The number of component is decided based on indices tthat estimate a models generalisation criteria: Q2 for regression and AUROC (calculated with predictive Y in the innter CV) for discriminant analytis.
#' @param type char, 'R' for regression or 'DA' for discriminant analysis, with prefix '-mY' if mutli-column Y
#' @param q2s num array with Q2 entries of principal components
#' @param aucs_te num array with AUROC entries of principal components
#' @param pc_max maximum number of components
#' @return bool: FALSE if another component should be fitted, TRUE if the current component is starting to overfit OR improves only little in prediction capcity when compared to the previous one.
#' @keywords internal
#' @section
#'
#' evalFit(type, q2_comp, aucs_tr aucs_te, maxPCo)
.evalFit = function(type, q2s, aucs_te, pc_max) {
    type=strsplit(type, '-')[[1]][1]

    if (type == "R") {
        cv_est = q2s
    } else {
        cv_est = aucs_te
    }
    nc = length(cv_est)
    if (nc == pc_max)
        return(TRUE)

    switch(strsplit(type, '-')[[1]][1],
           R = {
               if (any(is.na(cv_est))) {
                   stop(paste("Something went wrong: Q2 is NA"))
               }
               if (nc == 1 && cv_est < 0.05) stop(paste0("No significant component (r2=", round(cv_est, 2), ")"))
               if (nc > 1 && (diff(cv_est[(nc - 1):nc]) < 0.05 | cv_est[nc] > 0.98)) return(TRUE)
           },
           DA = {
               if (any(is.na(cv_est))) {
                   stop(paste("Something went wrong: Q2 is NA"))
               }
               if (nc == 1 && cv_est < 0.6) stop(paste0("No significant component (AUROC_cv=", round(cv_est, 2), ")"))
               if (nc > 1 && (diff(cv_est[(nc - 1):nc]) < 0.05 | cv_est[nc] > 0.97)) return(TRUE)
           })

    return(FALSE)
}






#' @title k-fold subsets for cross validation (CV)
#' @description This function ccreates a list of indices that represent k fold subsets for CV
#' @param k int k parameter
#' @param Y matrix (observations times variables)
#' @return list of Y row-indices for each fold
#' @keywords internal
#' @section
.kFold = function(k, Y){
    #browser()
    if(k > nrow(Y) || k < 2 || is.na(k) || is.infinite(k)) {message('Check cross-validation k paramter - using LOO-CV.'); k<-nrow(Y)}
    #browser()

    if (nrow(Y)/k < 2) {message('Cross-validation paramter k too low - using LOO-CV.');  sets <-sample(rep_len(seq_len(k), nrow(Y)))} else { # LOO-CV
        sets <-sample(rep_len(seq_len(k), nrow(Y)))
    }
    sets_list <- lapply(seq_len(k), function(i, idc = sets) {which(idc != i)})
    return(sets_list)
}


# use this for unit testing
# Y=cbind(c(rep('A', 50), rep('B', 50), rep('C', 6)))
# stratified=list(type='DA', Y=Y, probs=c(0, 0.5, 1))
# Y=cbind(sample(1:100))
# stratified=list(type='R', Y=Y, probs=c(0, 0.33, 0.66, 1))
# k=10
# out=sample_kFold_stratified(k, stratified)
# tt=sapply(out, function(x) length(unique(x)))
#' @title Y-stratified k-fold subsets for cross validation (CV)
#' @description This function creates a list of indices that represent k fold subsets for CV, stratified by the first column in Y. If type is regression, then Y will be categorised accroding to quantile function using probabilities igven in the thrid element of stratified list (see parameter).
#' @param k int CV parameter k
#' @param stratified list with 3 elements: 1) char: R for regression (Y is numeric), DA for discriminant analysis (Y treated as categorical), prefix '-mY' indicates multi-column Y (see Y) 2) Y matrix (observations times variables), always with ncol = 1 (stratified CV generation for multi-Y to be implemented) and 3) quantile function probabilities for stratification of numberic Y
#' @return list of Y row-indices for each fold
#' @keywords internal
#' @section
.kFoldStratified = function(k, stratified){
    if (grepl('R', stratified[[1]])) {
        Yori = stratified[[2]];
        Y = cbind(cut(Yori[, 1], breaks = as.numeric(quantile(Yori[, 1], stratified[[3]])), include.lowest = TRUE)) }else{
        Y=stratified[[2]][,1]
        }
    ct <- table(Y)
    if(min(ct)/max(ct) < 0.1){
        message(paste0("Skewed frequencies of Y classes- results are strongly biased towards higher represented classes (check Y). Use k-fold CV (non-stratified)."))
        return(NULL)}
    levs = names(ct)
    if (min(ct) <= k || min(ct)/k < 2 || is.na(k) || is.infinite(k) || k<2) {
        message(paste0("Setting CV k-fold parameter to ", min(ct) ,". Number of observations in group ", names(which.min(ct)), " is too low (n=", min(ct), ") for k=", k, "."))
        k = min(ct)
    }
    nobs <- floor(min(ct)/k)
    if(nobs == 0) stop('Number of observations in Y level too low.')
    set_lev <- lapply(levs, function(x, ks=k, y = Y, nob = nobs) {
        idx <- sample(which(y == x))
        k1 <- rep(seq_len(ks), length.out = length(idx))
        names(k1) <- c(idx)
        return(k1)
    })
    sets_list <- lapply(seq_len(k), function(i, idc = unlist(set_lev)) {
        idx <- which(idc != i)
        as.numeric(names(idc[idx]))
    })
    return(sets_list)
}




#' @title k-fold Monte Carlo Cross Validation (MCCV) subsets
#' @description This function creates a list of indices that represent k fold subsets for MCCV
#' @param k int k parameter
#' @param Y matrix (observations times variables) with ncol(Y)=1, always (column reduction for multicolumn Y  in parent function: .cvSetsMethod)!
#' @param split Ratio of the number of samples in training and test set
#' @return list of Y row-indices for each fold
#' @keywords internal
#' @section
.mc = function(k, Y, split){
    if (!is.integer(k)) {k=ceiling(k); warning(paste0("The k-fold parameter should be an integer. I'm rounding to k=", k))}
    if(k==0 || k > 1e6 || is.na(k) || is.infinite(k))  {stop('Check MC-CV k fold paramter!')}
    if (split >= 1 || split <= 0 || is.na(split) || is.infinite(split)) {stop("Check MC-CV split parameter - anything from 0.6 to 0.8 works best usually.")}
    if (split > 0.9 || split < 0.3) {message("What an unusual choice of the MC-CV split parameter (anything from 0.5 to 0.9 works best usually).")}
    sets_list <- lapply(seq_len(k), function(i, le = nrow(Y)) {sample(seq_len(le), le * split, replace = TRUE)})
    return(sets_list)
}


#' @title Monte Carlo Cross Validation (MCCV) generation of training sets indices
#' @description This function creates a list of indices that represent k fold subsets for MCCV
#' @param k int k parameter
#' @param split Fraction of samples used to generate training set for class/group-balanced Monte Carlo (MC) CV
#' @param stratified list with 3 elements: 1) char: R for regression (Y is numeric), DA for discriminant analysis (Y treated as categorical), prefix '-mY' indicates multi-column Y (see Y) 2) Y matrix (observations times variables), always with ncol = 1 (stratified CV generation for multi-Y to be implemented) and 3) quantile function probabilities for stratification of numberic Y
#' @keywords internal
#' @return list of Y row-indices for each fold
#' @section
.mcBalanced = function(k, split, stratified){
    if (!is.integer(k)) {k=ceiling(k); warning(paste0("The k-fold parameter should be an integer. I'm rounding to k=", k))}
    if(k<=2 || k > 1e6 || is.na(k) || is.infinite(k)) {stop('Check MC-CV k fold paramter!')}
    if (split >= 1 || split <= 0 || is.na(split) || is.infinite(split)) {stop("Check MC-CV split parameter - anything from 0.6 to 0.8 works best usually.")}
    if (split > 0.9 || split < 0.3) {message("What an unusual choice of the MC-CV split parameter (anything from 0.5 to 0.9 works best usually).")}
    if (grepl('R', stratified[[1]])) { Yori = stratified[[2]]; Y = cbind(cut(Yori[, 1], breaks = as.numeric(quantile(Yori[, 1], stratified[[3]])), include.lowest = TRUE))}else{
        Y=stratified[[2]][,1]
        }
    ct <- table(Y)
    if(min(ct)/max(ct) < 0.1){
        stop(paste0("Skewed frequencies of Y classes- results are strongly biased towards higher represented classes (check Y). Continuing with k-fold CV (non-stratified)."))
    }
    nobs <- floor(min(ct) * split)
    if (nobs == 0) {  stop("Insufficient nb of members of Y group. Check Y input and increase MCCV split argument.") }
    sets_list <- lapply(seq_len(k), function(i, y =  Y, nob = nobs, ctab = ct) {
        k1 <- list()
        for (j in seq_len(length(ctab))) {
            k1[[j]] <- sample(which(y == names(ctab)[j]), nob, replace = TRUE)
        }
        unlist(k1)
    })
    return(sets_list)
}

# Y=cbind(sample(1:100))
# stratified=list(type='R', Y=Y, probs=c(0, 0.33, 0.66, 1))
# out=.mcBalanced(k=2, split=2/4, stratified)
# sapply(out, length)
# sapply(out, summary)


#' @title Create Cross Validation (CV) row-indices for Y
#' @description This function creates a list of indices that represent k fold subsets for MCCV
#' @param Y matrix, if multi-column, only the first column will be used
#' @param type char Regression (R) or disciminant analysis (DA), combined with '-mY' indicating multi-column Y
#' @param method Type of cross validation: k-fold or Monte Carlo-Cross Validation (MCCV), sampling can be performed group-blind (\code{'k-fold'}, \code{'MC'}) or group-stratified (\code{'k-fold_stratified'}, \code{'MC_balanced'}).
#' @param k Number of training sets to generate.
#' @param split Fraction of samples used to generate training set for class/group-balanced Monte Carlo (MC) CV
#' @return list of Y row-indices for each CV fold
#' @keywords internal
#' @section
.cvSetsMethod <- function(Y, type, method = "k-fold_stratified", k = 7, split = 2/3) {

    if (grepl('mY', type)) {Y <- as.matrix(Y[,1], ncol = 1)} # reduce Y to first column
    if (!method %in% c("k-fold", "k-fold_stratified", "MC", "MC_balanced")) { stop(paste0("Check method argument, valid methods are: k-fold, k-fold_stratified, MC, MC_balanced"))}
    if (k%%1 ==0) {k<-as.integer(k)}
    if (grepl("MC", method) && split >= 1) { message("Paramter mc.split can only take values < 1 (ratio training samples / test samples)")}

    switch(method,
           `k-fold_stratified` = { sets_list<-.kFoldStratified(k, stratified=list(type, Y, probs=c(0, 0.33, 0.66, 1)))},
           `k-fold` = { sets_list<-.kFold(k, Y) },
           MC = { sets_list<-.mc(k, Y, split) },
           MC_balanced = { sets_list<-.mcBalanced(k, split, stratified=list(type, Y, probs=c(0, 0.33, 0.66, 1))) })

    if(is.null(sets_list)){stop('Check cross validation input arguments!')}
    return(sets_list)
}



#' @title Create dummy matrix and check numeric variable for OPLS analysis
#' @param Y Dependent variable (numeric or factor for regression or discriminat analysis, resepctively)
#' @return List of two: Dummy matrix and data frame mapping Y-levels to numeric representations.
#' @aliases create_dummy_Y
#' @author Torben Kimhofer \email{MetaboMate@@tkimhofer.com}
#' @keywords internal
#' @section
.createDummyY <- function(Y) {
    if (!is.numeric(Y)) {
        Y_levels <- unique(Y)
        if (length(Y_levels) == 2) {
            Y_new <- cbind(as.numeric(as.factor(Y)))
        } else {
            Y_new <- matrix(-1, nrow = length(Y), ncol = length(Y_levels))
            for (i in seq_len(length(Y_levels))) {
                Y_new[which(Y == Y_levels[i]), i] <- 1
            }
            colnames(Y_new) <- Y_levels
        }
        Y_levs <- unique(data.frame(Original = Y, Numeric = Y_new, stringsAsFactors = FALSE))
        # return(list(Y_new, Y_levs))
        return(list(Y_new, Y_levs))
    } else {
        return(list(cbind(Y), data.frame()))
    }
}

# #' @title Summarise model indices using tabular and graphic modes
# #' @param type char: regression (\code{R}) or discriminant analysis (\code{DA})
# #' @param r2x_comp num array: r2x
# #' @param r2_comp num array: r2y
# #' @param q2_comp num array: q2
# #' @param  aucs num array: auroc
# #' @param  cv list of cross-validion paramters (see opls function)
# #' @return List of two: 1. summary data.frame and 2. ggplot2 figure
# #' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
# #' @section
# .orthModelCompSummary = function(type, r2x_comp, r2_comp, q2_comp, aucs, cv) {
#
#
#     switch(type, DA = {
#         model_summary <- data.frame(PC_pred = 1, PC_orth = seq(q2_comp), R2X = round(r2x_comp, 2), R2Y = round(r2_comp, 2), Q2 = round(q2_comp, 2), AUROC = round(aucs, 2))
#     }, R = {
#         model_summary <- data.frame(PC_pred = 1, PC_orth = seq(q2_comp), R2X = round(r2x_comp, 2), R2Y = round(r2_comp, 2), Q2 = round(q2_comp, 2))
#     }, )
#
#     model_summary$PC_pred = 1
#     model_summary$PC_orth = seq_len(nrow(model_summary))
#
#     mm <- melt(model_summary, id.vars = c("PC_orth", "PC_pred"))
#     mm$PC <- paste0(mm$PC_pred, "+", mm$PC_orth)
#
#     mm$alpha1 <- 1
#     mm$alpha1[mm$PC_orth == max(mm$PC_orth)] <- 0.7
#
#     g <- ggplot(mm, aes_string("PC", " value", fill = "variable")) + geom_bar(stat = "identity", position = "dodge", colour = NA, aes_string(alpha = "alpha1")) + scale_y_continuous(limits = c(min(c(0,
#         max(model_summary$Q2 - 0.02, -0.05))), 1), breaks = breaks_pretty(), expand = c(0, 0)) + scale_fill_manual(values = c(R2X = "lightgreen", R2Y = "lightblue", Q2 = "red",
#         AUROC = "black"), labels = c(expression(R^2 * X), expression(R^{
#         2
#     } * Y), expression(Q^2), expression(AUROC["cv"])), name = "") + scale_alpha(guide = FALSE, limits = c(0, 1)) + labs(x = "Predictive + Orthogonal Component(s)", y = "", title = paste("O-PLS-",
#         type, "  - Component Summary", sep = ""), caption = paste0("\nCross validation: ", cv$method, " (k=", cv$k, ", ratio test/training set=", round(cv$split, 2), ")")) + theme_bw() +
#         theme(legend.text.align = 0, panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "black", size = 0.15), panel.grid.minor = element_blank(),
#             panel.border = element_blank(), axis.line.x = element_line(color = "black", size = 0.55), axis.line.y = element_blank(), axis.ticks = element_blank(), legend.key = element_rect(colour = "white"),
#             text = element_text(family = "Helvetica"))
#
#     return(list(model_summary, g))
#
# }



#' @title Performs OPLS modelling for each CV set and collates output
#' @param X num matrix X: preproc NMR data
#' @param Y num matrix Y: outcome, dummy matrix in case of categorical outcome
#' @param cv.set list of k elements containing vector of X and Y row-indices representing training set for cv round k
#' @param nc int, max number of orthogonal components to fit
#' @param mod.cv cv parameters
#' @return Named list of collated OPLS data for respective component
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @keywords internal
#' @section


# TODO: add scaling and centering
.oplsComponentCv=function(X, Y, cv.set, nc,  mod.cv){

    out=lapply(seq_along(cv.set), function(k){

        idc=cv.set[[k]]
        if(nc==1){
            Xcs <- .scaleMatRcpp(X, idc-1, center=TRUE, scale_type = 1)[[1]] # subtract 1 since Rcpp indexing starts at zero
        }else{
            Xcs=mod.cv[[k]]$x_res
        }

        Y_scale <- .scaleMatRcpp(Y, idc-1, center=TRUE, scale_type = 1)
        Ycs <- Y_scale[[1]] # subtract 1 since Rcpp indexing starts at zero
        opls_filt <- .nipOplsRcpp(X = Xcs[idc, ], Y = cbind(Ycs[idc, ]))
        pred_comp <- .nipPlsCompRcpp(X = opls_filt$X_res, Y = cbind(Ycs[idc, ]))

        #browser()
        Xte=Xcs[-idc, ]
        if(!is.matrix(Xte)){ Xte=matrix(Xte, nrow=1)}
        opls_pred=.oplsPredRcpp(opls_mod = opls_filt, pred_mod=pred_comp,  Xnew=Xte)


        # create list opls_mod:
        # with w_xo, p_xo, (matrix with columns, rows equating to the number of orthogonal components, resp) and
        # w_xp, p_xp, p_yp (a single one for predictive component)
        if(nc==1){
            mod.cv=list(
                # w_xo = matrix(NA, nrow=ncol(X), ncol=1),
                # p_xo = matrix(NA, nrow=1, ncol=ncol(X)),
                t_xo = matrix(NA, nrow=nrow(Y), ncol=1),
                # w_xp = matrix(NA, nrow=ncol(X), ncol=1),
                # p_xp = matrix(NA, nrow=1, ncol=ncol(X)),
                # p_yp = matrix(NA, nrow=ncol(Y), ncol=1),
                t_xp = matrix(NA, nrow=nrow(Y), ncol=1),
                y_pred_train = matrix(NA, nrow=nrow(Y), ncol=ncol(Y)),
                y_pred_test = matrix(NA, nrow=nrow(Y), ncol=ncol(Y)),
                x_res = matrix(NA, nrow=nrow(Y), ncol=ncol(Xcs))
                #r2x_pred_comp_cv = array()
            )
            # mod.cv$w_xo[,nc]=opls_filt$w_o
            # mod.cv$p_xo[nc,]=opls_filt$p_o
            # mod.cv$t_xo[idc, nc]=opls_filt$t_o
            mod.cv$t_xo[-idc, nc]=opls_pred$t_xo_new
            # mod.cv$w_xp[,nc]=pred_comp$w_x
            # mod.cv$p_xp[nc,]=pred_comp$p_x
            #mod.cv$t_xp[idc,nc]=pred_comp$t_p
            mod.cv$t_xp[-idc,nc]=opls_pred$t_pred # cv scores
            #mod.cv$p_py[idc,nc]=pred_comp$p_y
            #mod.cv$p_py[,nc]=opls_pred$y_pred
            mod.cv$y_pred_train[idc,]=pred_comp$y_pred
            mod.cv$y_pred_test[-idc,]=opls_pred$y_pred #t(apply(opls_pred$y_pred, 1, function(x, me= Y_scale$mean, ssd=Y_scale$sd) { (x+me) * ssd}))
            mod.cv$x_res[idc,] = opls_filt$X_res
            mod.cv$x_res[-idc,] = opls_pred$Xres
           # mod.cv$r2x_pred_comp_cv = .r2(opls_filt$X_res, pred_comp$t_x %*% pred_comp$p_x, tssx)
            return(mod.cv)
        }else{

            #browser()
            # mod.cv$w_xo=cbind(mod.cv$w_xo, opls_filt$w_o)
            # mod.cv$p_xo=rbind(mod.cv$p_xo, opls_filt$p_o)
            #
            mod.cv[[k]]$t_xo=cbind(mod.cv[[k]]$t_xo, matrix(NA, nrow=nrow(Y), ncol=1))
            # mod.cv$t_xo[idc, nc]=opls_filt$t_o
            mod.cv[[k]]$t_xo[-idc, nc]=opls_pred$t_xo_new
            #
            # mod.cv$w_xp[,nc]=pred_comp$w_p
            # mod.cv$p_xp[nc,]=pred_comp$p_p

            #mod.cv$t_xp[idc,nc]=pred_comp$t_p
            #browser()
            mod.cv[[k]]$t_xp=matrix(NA, nrow=nrow(Y), ncol=1) # predictive component will alwas be overwritten (set to NA, since MC methods will generated non-equal indies in each round)
            mod.cv[[k]]$t_xp[-idc,1]=opls_pred$t_pred

            #mod.cv$p_py[idc,nc]=pred_comp$p_y
            #mod.cv$p_py[-idc,nc]=opls_pred$y_pred

            y_pred_add=matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
            y_pred_add[-idc,]=opls_pred$y_pred #t(apply(opls_pred$y_pred, 1, function(x, me= Y_scale$mean, ssd=Y_scale$sd) { (x+me) * ssd}))
            mod.cv[[k]]$y_pred_test=y_pred_add

            y_pred_add1=matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
            y_pred_add1[idc,]=pred_comp$y_pred
            mod.cv[[k]]$y_pred_train=y_pred_add1


            mod.cv[[k]]$x_res[idc,] = opls_filt$X_res
            mod.cv[[k]]$x_res[-idc,] = opls_pred$Xres

            #mod.cv[[k]]$r2x_pred_comp_cv = .r2(opls_filt$X_res, pred_comp$t_x %*% pred_comp$p_x, tssx)

            return(mod.cv[[k]])
        }

    })
    return(out)
}


#' @title Extract data structures from output of function .oplsComponentCv
#' @param cv_obj named list, output of .oplsComponentCv (for resp OPLS component)
#' @param feat char, feature name in output of tt (names(tt))
#' @param cv_type char indicating cross validation type (k-fold, k-fold_stratified, MC, MC_balanced)
#' @param model_type char indicating single or multi-column Y (grep mY)
#' @return num array or matrix of resp component feature. For MCCV and single-column (multi-column) Y this is a 2D (3D) matrix where D1: mean, sd, or length, D2: nrow(Y) and D3: ncol(Y). For k-fold CV single-column (multi-column) Y output is 1D (2D) with D1: nrow(Y) and D2: ncol(Y)
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom abind abind
#' @importFrom stats sd
#' @keywords internal
#' @section
#'
# this has to be adjusted from multi-column Y
.extMeanCvFeat=function(cv_obj, feat='t_xp', cv_type, model_type){

    if(!feat %in% names(cv_obj[[1]])) stop('Feature name not in feature list')
    inter <- lapply(cv_obj, '[[', feat)


    if( grepl('MC', cv_type) | (grepl('k-fold', cv_type) & grepl('train', feat )) ) {
        if(grepl('mY', model_type)){
            f_mat=abind(inter, along=3);
            fout=apply(f_mat, c(1,2), function(x){ # 3d matrix: d1: mean, sd, or length, d2: nrow(Y), d3: ncol(Y)
                le=length(which(!is.na(x)))/length(x)
                xmean=sum(x, na.rm=TRUE)/le # calc mean
                xsd=sd(x, na.rm=TRUE) # calc sd
                c('mean'=xmean, 'sd'=xsd, '%cv'=le)
            })
        }else{
            # vector
            f_mat=do.call(cbind, inter);
            fout=apply(f_mat, 1, function(x){
                idx=which(!is.na(x))
                c('mean'=mean(x[idx], na.rm = TRUE), 'sd'=sd(x[idx], na.rm = TRUE), '%cv'=length(idx)/length(x))
            })
        }
    }else{
        if(grepl('mY', model_type) | !grepl('mY', model_type) ){
            f_mat=abind(inter, along=3);
            fout=apply(f_mat, c(1,2), function(x){
                idx=which(!is.na(x))
                x[idx]
            })
        }
        # else{
        #     # vector
        #     f_mat=abind(inter, along=3);
        #     f_mat=do.call(cbind, inter);
        #     fout=apply(f_mat, 1, function(x){
        #         idx=which(!is.na(x))
        #         x[idx]
        #     })
        # }
    }

    return(fout)
}

#' @title R2 and Q2
#' @param Y num matrix Y: OPLS outcome, dummy matrix in case it is categorical
#' @param Yhat num matrix, predicted Y
#' @param ytss num, total sum of squares of Y
#' @return R2=1-(PRESS/TSS)
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @keywords internal
#' @section
.r2 <- function(Y, Yhat, ytss) {


    press <- sum(as.vector((Y - Yhat)^2), na.rm=TRUE) / ncol(Y)

    if(is.null(ytss)){  ytss <- sum(as.vector(Y^2), na.rm=TRUE) / ncol(Y) }

    1-(press / ytss)

}



#' @title Summarise component indices
#' @param type char, OPLS type: DA (discriminant analysis) or R (regression), with prefix '-mY' if mutli-column Y
#' @param r2x_comp num array, r2x value for each component
#' @param r2_comp num array, r2y value for each component
#' @param q2_comp num array, q2 value for each component
#' @param aucs_tr num array, cv training set aucs value for each component (defined only for DA)
#' @param aucs_te num array, cv test set aucs value for each component (defined only for DA)
#' @param cv list, see cv argument in function opls
#' @return list: 1. data.frame, summary table 2. ggplot2, viz of summary tbl
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom scales breaks_pretty
#' @keywords internal
#' @section
.orthModelCompSummary=function(type, r2x_comp, r2_comp, q2_comp, aucs_tr, aucs_te, cv){
    type=strsplit(type, '-')[[1]][1]

    switch(type,
           'DA' = {model_summary <- data.frame(PC_pred=1, PC_orth=seq(q2_comp), R2X = round(r2x_comp, 2), AUROC = round(aucs_tr, 2), AUROC_CV = round(aucs_te, 2))},
           'R' =  {model_summary <- data.frame(PC_pred=1, PC_orth=seq(q2_comp), R2X = round(r2x_comp, 2), R2Y = round(r2_comp, 2), Q2 = round(q2_comp, 2))}
    )
    model_summary$PC_pred=1
    model_summary$PC_orth=seq_len(nrow(model_summary))

    mm <- melt(model_summary, id.vars = c("PC_orth", "PC_pred"))
    mm$PC <- paste0(mm$PC_pred, '+', mm$PC_orth)

    mm$alph <- 1
    mm$alph[mm$PC_orth==max(mm$PC_orth)] <- 0.7

    idx <- which(mm$value < (-0.015) & mm$variable == 'Q2')
    mm$value[idx] <- (-0.01)

    g <- ggplot(mm, aes_string("PC", " value", fill = "variable")) +
        geom_bar(stat = "identity", position = "dodge", colour = NA, aes_string(alpha = "alph")) +
        scale_fill_manual(values = c(R2X = "lightgreen", R2Y = "lightblue", Q2 = "red", AUROC = "black", AUROC_CV = "red"), labels = c(expression(R^2 * X), expression(R^{2} * Y), expression(Q^2), expression(AUROC[cv])), name = "") +
        scale_alpha(guide = FALSE, limits = c(0, 1)) +
        labs(x = "Predictive + Orthogonal Component(s)",
             y = "",
             title = paste("O-PLS-", type, "  - Component Summary", sep = ""),
             caption=paste0('\nCross validation: ', cv$method, ' (k=', cv$k, ', ratio test/training set=', round(cv$split, 2), ')')) +
        theme_bw() +
        theme(
            legend.text.align = 0,
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "black", size = 0.15),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line.x = element_line(color = "black", size = 0.55),
            axis.line.y = element_blank(),
            axis.ticks = element_blank(),
            legend.key = element_rect(colour = "white"),
            text = element_text(family = "Helvetica"))

    if( any(mm$value < 0, na.rm = TRUE) ) { g=g+scale_y_continuous(limits = c(-0.015, 1), breaks = breaks_pretty(), expand = c(0,0)) } else{
        g=g+scale_y_continuous(limits = c(0, 1), breaks = breaks_pretty(), expand = c(0,0))
    }

    return(list(model_summary, g))

}


#' Check dimension, NA and infinity
#' @param X NMR matrix with spectra represented in rows.
#' @param ppm ppm vector.
#' @return  Logical indicating if X and ppm match and absence of NA and inf in ppm
#' @keywords internal
#' @section
.check_X_ppm <- function(X, ppm){

    if( any( is.na(ppm) | is.infinite(ppm) ) ) return(FALSE)

    if( ncol(X) != length(ppm) ) return(FALSE)

    return(TRUE)

}

#' Find indices in ppm vector for respective chemical shift range
#' @export
#' @param range num, range chemical shift (in ppm vector)
#' @param ppm num, ppm vector
#' @export
#' @aliases get.idx
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @section
get.idx <- function(range = c(1, 5), ppm) {
    range <- sort(range, decreasing = TRUE)
    which(ppm <= range[1] & ppm >= range[2])
}



#' Min-max scaling
#' @export
#' @param x Numeric vector to be scaled.
#' @return Scaled x vector.
#' @details Data are scaled to range between zero and one:
#' \deqn{X_{scaled}=\frac{x-x_{min}}{x_{max}-x_{min}}}
#' @usage minmax(x)
#' @author Torben Kimhofer \email{torben.kimhofer@@gmail.com}
#' @section
minmax <- function(x) {
    (x - min(x))/(max(x) - min(x))
}

#' Calculating full width at half max
#' @export
#' @param X NMR matrix with rows representing spectra.
#' @param ppm ppm vector with its length equals to \code{nrow(X)}.
#' @param shift Signal shift to calculate line width
#' @description Calculating full width at half maximum (FWHM, aka line width). This function simply returns the ppm difference where peak line crosses half of the peak height. It requires one signal across all spectra within ppm ranges specified in \code{shift}.
#' @return Array of line widths in ppm. To convert from ppm to Hertz (Hz), multiply values with the spectrometer frequency (column \code{a_SF01} in \code{meta} data frame).
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @section
lw <- function(X, ppm, shift = c(-0.01, 0.01)) {
    idx <- get.idx(shift, ppm)
    fwhm <- apply(X[, idx], 1, function(x, pp = ppm[idx]) {
        x <- x + abs(min(x))  # no negative values
        baseline <- min(x)
        height <- max(x) - baseline
        hw <- baseline + (0.5 * height)
        f <- approxfun(pp, x)
        x_new <- seq(-0.01, 0.01, by = 1e-05)
        y_new <- f(x_new)
        diff(x_new[range(which(y_new > hw))])
    })
    return(fwhm)
}


#' @title Estimation of noise level of 1D proton spectrum
#' @export
#' @param NMR Input matrix with rows representing spectra.
#' @param ppm ppm vector with its length equals to \code{ncol(X)}.
#' @param where Signal free region across all NMR spectra (see Details).
#' @details Estimation of noise level in NMR spectra. This is useful for quality control checks (e.g., before and after spectral normalisation). Noise estimation requires a signal-free ppm region across all spectra, usually this is at the extreme ends or the spectrum. This function requires a minimum number of 50 data points to robustly estimate noise levels.
#' @return Returned is a vector of noise levels for each spectrum.
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom ptw asysm
#' @section
noise.est <- function(NMR, ppm, where = c(14.6, 14.7)) {
    # set ppm range where noise should be estimated, i.e., no signals
    idx <- get.idx(where, ppm)
    if (length(idx) < 50) {
        stop("No or too few data points for noise estimation!")
    }
    noise <- apply(NMR[, idx], 1, function(x) {
        x_driftcorrected <- x - asysm(x, lambda = 1e+10)
        noi <- (max(x_driftcorrected) - min(x_driftcorrected))
        return(noi)
    })
    return(noise)
}



#' @title Prep data.frame for plotting functions
# @importFrom ptw asysm
#' @keywords internal
#' @section
.viz_df_helper=function(obj, pc, an, type='p'){

    an=.check_an_viz(an, obj)
    # #browser()
    # if (missing(an)) {
    #     if(class(obj)[1]=='OPLS_metabom8') {an <- list(Y=obj@Y$ori)} else{ an <- list(NA)}
    #     }
    #
    # if(length(an)>3){message('an is a list of maximum three elements, see documentation for more infos.')}
    #
    # an=an[!sapply(an, is.null)]
    #
    # # create df for ggplot
    # if(is.null(names(an))){
    #     names(an)=paste0('Var', seq(length(an)))
    # }
    #
    # idx_dup = duplicated(names(an))
    # names(an)[idx_dup]=paste0(names(an)[idx_dup], seq(length(which(idx_dup))))
    #
    #
    # idx=is.na(names(an))
    # if(any(idx)){
    #     names(an)[idx]=paste0('Var', seq(length(idx)))
    # }
    # get orth comp
    idx_orth=grepl('o', pc)
    if(any(idx_orth)){
        pc1=as.numeric(gsub('o', '', pc))
        if(any(is.na(pc1)) || any(is.infinite(pc1)) || any(pc1[idx_orth] > nrow(obj@p_orth))){stop('Check pc argument and see help section.')}
    }else{pc1=pc}

    com=list()


    if(type=='p'){
        # extract data from obj
        switch(class(obj)[1],
               "PCA_metabom8"={
                   for(i in seq_len(length(pc))){
                       com[[i]]=obj@p[,pc[i]]
                   }
               },
               "OPLS_metabom8"={
                   for(i in seq_len(length(pc))){
                       if( grepl('o', pc[i]) ){ com[[i]]=obj@p_orth[pc1[1],]}else{com[[i]]=obj@p_pred[1,]}
                   }
               }
        )

        melted=as.data.frame(com)
        colnames(melted)=paste0('pc_', pc)

    }



    if(type=='t'){
        # extract data from obj
        switch(class(obj)[1],
               "PCA_metabom8"={
                   for(i in seq_len(length(pc))){
                       com[[i]]=obj@t[,pc[i]]
                   }
               },
               "OPLS_metabom8"={
                   for(i in seq_len(length(pc))){
                       if( grepl('o', pc[i]) ){ com[[i]]=obj@t_orth[,pc1[1]]}else{com[[i]]=obj@t_pred[,1]}
                   }

               }
        )
        melted=as.data.frame(com)
        colnames(melted)=paste0('pc_', pc)
    }


    if(type=='t_cv'){
        # extract data from obj
        switch(class(obj)[1],
               "PCA_metabom8"={
                   stop('t_cv not defined in PCA context')
               },
               "OPLS_metabom8"={
                   for(i in seq_len(length(pc))){
                       if( grepl('o', pc[i]) ){ com[[i]]=obj@t_orth_cv[,pc1[1]]}else{com[[i]]=obj@t_pred_cv[,1]}
                   }
               }
        )
        melted=as.data.frame(com)
        colnames(melted)=paste0('pc_', pc)
    }

    #browser()
    # add annotation data
    an_full=lapply(an, function(x, le_m=nrow(melted)){
        #browser()
        le_x=length(x)
        if(le_x==1){return(rep(x, le_m))}
        if(le_x==le_m){return(x)}
        if(le_x >1 & le_x < le_m){stop('Check list elements of an argument for length discrepancies.')}
    })

    an_df=as.data.frame(an_full, row.names = NULL)
    colnames(an_df)=names(an)
    melted=cbind(melted, an_df)
    #browser()
    return(list(df=melted, an_le=ncol(an_df)))

}



.check_pc_model <- function(pc, mod, le=1, type='p') {

    if(is.na(pc) || is.infinite(pc) || length(pc)>le) stop('Check pc argument.')

    if(class(mod)[1]=='PCA_metabom8'){
        if(!is.numeric(pc)){stop('PC value is not numeric.')}
        if(max(pc)>nrow(mod@p)){stop('PC value exceeds number of principal components.')}
    }

    if(class(mod)[1]=='OPLS_metabom8'){


        idx_orth=grepl('o', pc)
        if(any(idx_orth)){
            pc1=as.numeric(gsub('o', '', pc))
            if(any(is.na(pc1)) || any(is.infinite(pc1)) || any(pc1[idx_orth] > nrow(mod@o_orth))){stop('Check pc argument and see help section.')}

        }else{
            ddl=data.frame(x=mod@p_pred[1,], id=colnames(mod@X))
        }


    }

}




#' @title Backscaling for NMR data
#' @keywords internal
#' @section
.load_backscaled_nmr=function(mod, pc, idx, ppm){
    p_mod=.viz_df_helper(mod, pc, an=NA, type='p')
    # backscaling p
    p_bs <- p_mod$df[,1] * mod@X_sd
    p_abs <- minmax(abs(p_mod$df[,1]))
    df <- data.frame(p_bs, p_abs, ppm)
    df <- df[idx, ]

    return(df)
}







#' @title cor / cov model scores and NMR data
#' @keywords internal
#' @section
.load_stat_reconstr_nmr=function(mod, pc, X, idx, ppm){
    t_mod <- .viz_df_helper(mod, pc, an=NA, type='t')
    cc <- cor(t_mod$df[,1], X)[1, ]
    cv <- cov(t_mod$df[,1], X)[1, ]
    df <- data.frame(cor = abs(cc), cov = cv, ppm = ppm)
    df <- df[idx, ]

    return(df)
}



# check input arguments for plotting pca / opls scores
.check_an_viz=function(an, obj){
    #browser()
    if (missing(an)) {
        if(class(obj)[1]=='OPLS_metabom8') {an <- list(Y=obj@Y$ori)} else{ an <- list('All samples')}
    }

    if(length(an)>3){message('an is a list of maximum three elements, see documentation for more infos.')}
    an=an[!vapply(an, is.null, FUN.VALUE = TRUE)]

    # create df for ggplot
    if(is.null(names(an))){
        names(an)=paste0('Var', seq(length(an)))
    }

    idx_dup = duplicated(names(an))
    names(an)[idx_dup]=paste0(names(an)[idx_dup], seq(length(which(idx_dup))))


    idx=is.na(names(an))
    if(any(idx)){
        names(an)[idx]=paste0('Var', seq(length(idx)))
    }

    return(an)

}





# OPLS Y-oermutations
.permYmod=function(Xs, Y, cv, type, nc_o){

    # remove orthogonal component(s)
    for(i in seq(nc_o)){
        if(nc_o == 1){ tt<-.oplsComponentCv(Xs, Y=Y, cv$cv_sets, nc_o[i],  mod.cv=NULL)
        }else{
            tt<-.oplsComponentCv(X=NA, Y=Y, cv$cv_sets, nc_o[i],  mod.cv=tt)
        }
    }

    # extract cv stats -> if Y is multi-column, fct .extrMeanCVFeat function is likely not working
    preds_test<-.extMeanCvFeat(cv_obj = tt, feat = 'y_pred_test', cv_type = cv$method, model_type=type)
    preds_train<-.extMeanCvFeat(tt, feat = 'y_pred_train', cv_type = cv$method, model_type=type)

    r2_comp <- q2_comp <- aucs_tr <- aucs_te <- array()
    nc=1

    tssy<-.tssRcpp(Y)/ncol(Y)

    # calculate R2 and auroc for cv compounds
    switch(strsplit(cv$method, '_')[[1]][1],
           'MC'={
               if( grepl('DA', type) ) {
                   if ( grepl('mY', type) ) {                           # multi Y DA using MCCV
                       pred_mean=preds_test[1,,]
                       colnames(pred_mean)=colnames(Y)
                       mod <- multiclass.roc(response = factor(Y), predictor = pred_mean)
                       aucs_te[nc] <- mod$auc
                       pred_tr_mean=preds_train[1,,]
                       colnames(pred_tr_mean)=colnames(Y)
                       mod <- multiclass.roc(response = factor(Y), predictor = pred_tr_mean)
                       aucs_tr[nc] <- mod$auc
                   }else{                                             # single Y DA using MCCV
                       # need the same decision boundary
                       mod <- roc(response = Y, predictor = preds_test[1,], quiet = TRUE)
                       aucs_te[nc] <- mod$auc
                       mod <- roc(response = Y, predictor = preds_train[1,], quiet = TRUE)
                       aucs_tr[nc] <- mod$auc
                   }
               } else{
                   if (grepl('mY', type)) {                         # multi Y regression using MCCV
                       r2_comp[nc] <- .r2(Y, preds_test[1,,], NULL)
                       q2_comp[nc] <- .r2(Y, preds_test[1,,], tssy)
                   }else{                                           # single Y regression using MCCV
                       r2_comp[nc] <- .r2(Y, preds_test[1,], NULL)
                       q2_comp[nc] <- .r2(Y, preds_test[1,], tssy)
                   }
               }
           },
           'k-fold'={
               if( grepl('DA', type) ) {
                   if ( grepl('mY', type) ) {
                       colnames(preds_test)=colnames(Y)
                       mod <- multiclass.roc(response = Y, predictor = apply(preds_test, 2, as.numeric))
                       aucs_te[nc] <- mod$auc
                       preds_te = preds_train[1,,] # this extracts mean values
                       colnames(preds_te)=colnames(Y)
                       mod <- multiclass.roc(response = Y, predictor = preds_te, quiet = TRUE)
                       aucs_tr[nc] <- mod$auc
                   }else{
                       mod <- roc(response = as.vector(Y), predictor = as.vector(preds_test), quiet = TRUE)
                       aucs_te[nc] <- mod$auc
                       mod <- roc(response = as.vector(Y), predictor = preds_train[1,], quiet = TRUE)
                       aucs_tr[nc] <- mod$auc
                   }
               }else{

                   if (grepl('MC', type)) {
                       r2_comp[nc] <- .r2(Y, preds_test[1,,], NULL)
                       q2_comp[nc] <- .r2(Y, preds_test[4,], tssy)
                   }else{
                       r2_comp[nc] <- .r2(Y, preds_test[1,,], NULL)
                       q2_comp[nc] <- .r2(Y, preds_test[4,], tssy)
                   }
               }
           }
    )

    return(list( r2_comp, q2_comp, aucs_tr, aucs_te))

}


















































































