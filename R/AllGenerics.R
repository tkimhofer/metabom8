#' @title Check for missing values (PLS context)
#' @description This function checks for missing values in Y (NA, NAN, infinite) and establishes analysis type accroding to the class of the input data Y: regression (R) or discriminant analysis (DA).
#' @param Y Input data (uni or multivar) formatted as vector or matrix or data.frame
#' @return List of i) Y matrix, ii) Y levels (empty if numeric Y), iii) Y type (R or DA)
#' @section

# check Y for regression or DA
.checkYclassNas = function(Y) {

    if (!class(Y) %in% c("matrix", "data.frame")) {
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

    return(c(Y_out, type))

}


#' @title Check X matrix (PLS context)
#' @description This function checks for missing values in Y (NA, NAN, infinite), numeric variable format.
#' @param X Input data (uni or multivar) formatted as vector or matrix or data.frame
#' @return NULL if all's fine and throws error otherwise
#' @section
.checkXclassNas = function(X) {

    if (!class(X) %in% c("matrix") | !is.numeric(X[1, 1])) {
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
#' @param type char, 'R' for regression or 'DA' for discriminant analysis
#' @param q2s num array with Q2 entries of principal components
#' @param aucs num array with AUROC entries of principal components
#' @param pc_max maximum number of components
#' @return bool: FALSE if another component should be fitted, TRUE if the current component is starting to overfit OR improves only little in prediction capcity when compared to the previous one.
#' @section
.evalFit = function(type, q2s, aucs, pc_max) {

    if (type == "R") {
        cv_est = q2s
    } else {
        cv_est = aucs
    }
    nc = length(cv_est)

    if (nc == pc_max)
        return(TRUE)

    switch(type, R = {

        if (any(is.na(cv_est))) {
            stop(paste("Something went wrong: Q2 is NA"))
        }

        if (nc == 1 && cv_est < 0.05) stop(paste0("No significant component (r2=", round(cv_est, 2), ")"))

        if (nc > 1 && (diff(cv_est[(nc - 1):nc]) < 0.05 | cv_est > 0.98)) return(TRUE)

    }, DA = {

        if (any(is.na(cv_est))) {
            stop(paste("Something went wrong: Q2 is NA"))
        }

        if (nc == 1 && cv_est < 0.6) stop(paste0("No significant component (AUROC_cv=", round(cv_est, 2), ")"))


        if (nc > 1 && (diff(cv_est[(nc - 1):nc]) < 0.05 | cv_est > 0.97)) return(TRUE)

    })


    return(FALSE)


}






#' @title k-fold subsets for cross validation (CV)
#' @description This function ccreates a list of indices that represent k fold subsets for CV
#' @param k int k parameter
#' @param Y matrix (observations times variables)
#' @return list of Y row-indices for each fold
#' @section
.kFold = function(k, Y){
    if (nrow(Y)/k < 2) {sets <- seq_len(nrow(Y))} else { # LOO-CV
        sets <- sample.int(k, nrow(Y), replace = TRUE)
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
#' @description This function ccreates a list of indices that represent k fold subsets for CV, stratified by the first column in Y
#' @param k int k parameter
#' @param stratified list with 3 elements: 1) type char: R for regression (Y is numeric) or DA for discriminant analysis (Y treated as categorical), 2) Y matrix (observations times variables) and 3) quantile function probabilities for stratification of numberic Y
#' @return list of Y row-indices for each fold
#' @section
.kFoldStratified = function(k, stratified){
    if (stratified[[1]] == "R") { Yori = stratified[[2]]; Y = cbind(cut(Yori[, 1], breaks = as.numeric(quantile(Yori[, 1], stratified[[3]])), include.lowest = TRUE))}
    ct <- table(Y)
    if(min(ct)/max(ct) < 0.1){
        message(paste0("Skewed frequencies of Y classes- results are strongly biased towards higher represented classes (check Y). Continuing with k-fold CV (non-stratified)."))
        return(NULL)}
    levs = names(ct)
    if (min(ct) <= k || min(ct)/k < 2) {
        message(paste0("Setting CV k-fold parameter to ", min(ct) ,". Number of observations in group ", names(which.min(ct)), " is too low (n=", min(ct), ") for k=", k, "."))
        k = min(ct)
    }
    nobs <- floor(min(ct)/k)

    set_lev <- lapply(levs, function(x, ks=k, y = Y[, 1], nob = nobs) {
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
#' @param Y matrix (observations times variables)
#' @param split Ratio of the number of samples in training and test set
#' @return list of Y row-indices for each fold
#' @section
.mc = function(k, Y, split){
    if (!is.integer(k)) {k=ceiling(k); warning(paste0("The k-fold parameter should be an integer. I'm rounding to k=", k))}
    if(k==0 || k > 1e6) {stop('Check MC-CV k fold paramter!')}
    if (split == 1 || split == 0) {stop("The MC-CV split parameter can't take a value of zero or one. Anything from 0.6 to 0.8 works best usually.")}
    if (split > 0.9 || split < 0.3) {message("What an unusual choice of the MC-CV split parameter (anything from 0.5 to 0.9 works best usually).")} else { # LOO-CV
        sets <- sample.int(k, nrow(Y), replace = TRUE)
    }
    sets_list <- lapply(seq_len(k), function(i, le = nrow(Y)) {sample(seq_len(le), le * split, replace = TRUE)})
    return(sets_list)
}

# out=.mc(k=10, Y=cbind(1:100), split=1)
# out=.mc(k=10, Y=cbind(1:100), split=0)
# out=.mc(k=10, Y=cbind(1:100), split=NA)
#
# out=.mc(k=1, Y=cbind(1:100), split=0.2)
# out=.mc(k=1, Y=cbind(1:100), split=0.99)
#
# sapply(out, length)


#' @title k-fold Monte Carlo Cross Validation (MCCV) subsets
#' @description This function creates a list of indices that represent k fold subsets for MCCV
#' @param k int k parameter
#' @param split Fraction of samples used to generate training set for class/group-balanced Monte Carlo (MC) CV
#' @param stratified list with 3 elements: 1) type char: R for regression (Y is numeric) or DA for discriminant analysis (Y treated as categorical), 2) Y matrix (observations times variables) and 3) quantile function probabilities for stratification of numberic Y
#' @return list of Y row-indices for each fold
#' @section
.mcBalanced = function(k, split, stratified){
    if (!is.integer(k)) {k=ceiling(k); warning(paste0("The k-fold parameter should be an integer. I'm rounding to k=", k))}
    if(k==0 || k > 1e6) {stop('Check MC-CV k fold paramter!')}
    if (split == 1 || split == 0) {stop("The MC-CV split parameter can't take a value of zero or one. Anything from 0.6 to 0.8 works best usually.")}
    if (split > 0.9 || split < 0.3) {message("What an unusual choice of the MC-CV split parameter (anything from 0.5 to 0.9 works best usually).")} else { # LOO-CV
        if (stratified[[1]] == "R") { Yori = stratified[[2]]; Y = cbind(cut(Yori[, 1], breaks = as.numeric(quantile(Yori[, 1], stratified[[3]])), include.lowest = TRUE))}
        ct <- table(Y)
        if(min(ct)/max(ct) < 0.1){
            stop(paste0("Skewed frequencies of Y classes- results are strongly biased towards higher represented classes (check Y). Continuing with k-fold CV (non-stratified)."))
        }
        levs = names(ct)
        nobs <- floor(min(ct) * split)
        if (nobs == 0) {
            stop("Insufficient nb of members of Y group. Check Y input and increase MCCV split argument.")
        }
        sets_list <- lapply(seq_len(k), function(i, y = Y[, 1], nob = nobs, ctab = ct) {
            k1 <- list()
            for (j in seq_len(length(ctab))) {
                k1[[j]] <- sample(which(y == names(ctab)[j]), nob, replace = TRUE)
            }
            unlist(k1)
        })


    }
    sets_list <- lapply(seq_len(k), function(i, le = nrow(Y)) {sample(seq_len(le), le * split, replace = TRUE)})
    return(sets_list)
}

# Y=cbind(sample(1:100))
# stratified=list(type='R', Y=Y, probs=c(0, 0.33, 0.66, 1))
# out=.mcBalanced(k=2, split=2/4, stratified)
# sapply(out, length)
# sapply(out, summary)


#' @title Create Cross Validation (CV) row-indices for Y
#' @description This function creates a list of indices that represent k fold subsets for MCCV
#' @param Y matrix, CV row indices to be generated for
#' @param type char Regression (R) or disciminant analysis (DA)
#' @param method Type of cross validation: k-fold or Monte Carlo-Cross Validation (MCCV), sampling can be performed group-blind (\code{'k-fold'}, \code{'MC'}) or group-stratified (\code{'k-fold_stratified'}, \code{'MC_balanced'}).
#' @param k Number of training sets to generate.
#' @param split Fraction of samples used to generate training set for class/group-balanced Monte Carlo (MC) CV
#' @return list of Y row-indices for each CV fold
#' @section
.cvSetsMethod <- function(Y, type = c("R", "DA"), method = "k-fold_stratified", k = 7, split = 2/3) {

    if (is.null(ncol(Y))) {Y <- as.matrix(Y, ncol = 1)}
    if (!method %in% c("k-fold", "k-fold_stratified", "MC", "MC_balanced")) {
        stop(paste0("Check method argument, valid methods are: k-fold, k-fold_stratified, MC, MC_balanced"))
    }
    if (grepl("MC", method) & split >= 1) {
        message("Paramter mc.split can only take values < 1 (ratio training samples / test samples)")
    }


    switch(method, `k-fold_stratified` = {
        sets_list<-.kFoldStratified(k, stratified=list(type, Y, probs=c(0, 0.33, 0.66, 1)))
    }, `k-fold` = {
        sets_list<-.kFold(k, Y)
    }, MC = {
        sets_list<-.mc(k, Y, split)
    }, MC_balanced = {
        sets_list<-.mcBalanced(k, split, stratified=list(type, Y, probs=c(0, 0.33, 0.66, 1)))
    })

    if(is.null(sets_list)){stop('Check cross validation input arguments!')}

    return(sets_list)
}

# Y=cbind(sample(1:100))
# .cvSetsMethod(Y, type='R', method="k-fold_stratified",  k = 7, split = 2/3)
# .cvSetsMethod(Y, type='R', method="k-fold",  k = 7, split = 2/3)
# .cvSetsMethod(Y, type='R', method="MC",  k = 7, split = 2/3)
# .cvSetsMethod(Y, type='R', method="MC_balanced",  k = 7, split = 2/3)




#' @title Create dummy matrix and check numeric variable for OPLS analysis
#' @param Y Dependent variable (numeric or factor for regression or discriminat analysis, resepctively)
#' @return List of two: Dummy matrix and data frame mapping Y-levels to numeric representations.
#' @aliases create_dummy_Y
#' @author Torben Kimhofer \email{MetaboMate@@tkimhofer.com}
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



#' @title Summarise model indices using tabular and graphic modes
#' @param type char: regression (\code{R}) or discriminant analysis (\code{DA})
#' @param r2x_comp num array: r2x
#' @param r2_comp num array: r2y
#' @param q2_comp num array: q2
#' @param  aucs num array: auroc
#' @param  cv list of cross-validion paramters (see opls function)
#' @return List of two: 1. summary data.frame and 2. ggplot2 figure
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @section
.orthModelCompSummary = function(type, r2x_comp, r2_comp, q2_comp, aucs, cv) {


    switch(type, DA = {
        model_summary <- data.frame(PC_pred = 1, PC_orth = seq(q2_comp), R2X = round(r2x_comp, 2), R2Y = round(r2_comp, 2), Q2 = round(q2_comp, 2), AUROC = round(aucs, 2))
    }, R = {
        model_summary <- data.frame(PC_pred = 1, PC_orth = seq(q2_comp), R2X = round(r2x_comp, 2), R2Y = round(r2_comp, 2), Q2 = round(q2_comp, 2))
    }, )

    model_summary$PC_pred = 1
    model_summary$PC_orth = seq_len(nrow(model_summary))

    mm <- melt(model_summary, id.vars = c("PC_orth", "PC_pred"))
    mm$PC <- paste0(mm$PC_pred, "+", mm$PC_orth)

    mm$alpha1 <- 1
    mm$alpha1[mm$PC_orth == max(mm$PC_orth)] <- 0.7

    g <- ggplot(mm, aes_string("PC", " value", fill = "variable")) + geom_bar(stat = "identity", position = "dodge", colour = NA, aes(alpha = mm$alpha1)) + scale_y_continuous(limits = c(min(c(0,
        max(model_summary$Q2 - 0.02, -0.05))), 1), breaks = breaks_pretty(), expand = c(0, 0)) + scale_fill_manual(values = c(R2X = "lightgreen", R2Y = "lightblue", Q2 = "red",
        AUROC = "black"), labels = c(expression(R^2 * X), expression(R^{
        2
    } * Y), expression(Q^2), expression(AUROC["cv"])), name = "") + scale_alpha(guide = FALSE, limits = c(0, 1)) + labs(x = "Predictive + Orthogonal Component(s)", y = "", title = paste("O-PLS-",
        type, "  - Component Summary", sep = ""), caption = paste0("\nCross validation: ", cv$method, " (k=", cv$k, ", ratio test/training set=", round(cv$split, 2), ")")) + theme_bw() +
        theme(legend.text.align = 0, panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "black", size = 0.15), panel.grid.minor = element_blank(),
            panel.border = element_blank(), axis.line.x = element_line(color = "black", size = 0.55), axis.line.y = element_blank(), axis.ticks = element_blank(), legend.key = element_rect(colour = "white"),
            text = element_text(family = "Helvetica"))

    return(list(model_summary, g))

}

