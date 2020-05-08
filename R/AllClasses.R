#' @title An S4 class to represent an OPLS model constructed with MetaboMate
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}

# define slots for OPLS_Torben object
setClass("OplsMate", representation(type = "character", t_pred = "matrix", p_pred = "matrix", w_pred = "matrix", betas_pred = "numeric", Qpc = "matrix", t_cv = "matrix", t_orth_cv = "matrix", 
    t_orth = "matrix", p_orth = "matrix", w_orth = "matrix", nPC = "numeric", summary = "data.frame", X_orth = "matrix", Y_res = "matrix", X_mean = "numeric", X_sd = "numeric", 
    Y_mean = "numeric", Y_sd = "numeric", Y_dummy = "data.frame", Parameters = "list", X = "matrix"))
