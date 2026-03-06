#####################
### PREPROCESSING ###
#####################

#' @rdname prep
setMethod("prep", "m8_preprocess",
          function(object, X){

            msd <- .sdRcpp(X)

            Xcs <- .scaleMatRcpp(
              X,
              0:(nrow(X)-1),
              center = object@center,
              scale_type =
                switch(object@scale,
                       "None"=0,
                       "UV"=1,
                       "Pareto"=2)
            )[[1]]

            object@X_mean <- msd[[2]]
            object@X_sd   <- msd[[1]]

            list(X=Xcs, prep=object)
          })


#########################################################
### STATISTICAL VALIDATION - RESAMPLING INSTANTIATION ###
#########################################################

# seed from random number generator (rng) does not exist on new R session,
# that instantiation will fail (needs seed values). To avoid this,
# rng is set up by calling a sampling function runif.
#' @importFrom stats runif
.ensureRNG <- function() {
  if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)
}

setMethod(".rngPreflight",
          "ResamplingStrategy",
          function(object) .ensureRNG())


#' @rdname instantiate
setMethod("instantiate",
          "KFold",
          function(object, Y){

            .rngPreflight(object)
            seed <- .Random.seed

            idx <- .kFold(object@k, Y)

            methods::new("ResamplingInstance",
                train    = idx,
                strategy = object,
                n        = nrow(Y),
                seed     = seed)
          })

#' @rdname instantiate
setMethod("instantiate",
          "StratifiedKFold",
          function(object, Y){

            .rngPreflight(object)
            seed <- .Random.seed

            idx <- .kFoldStratified(
              object@k,
              stratified = list(object@type,
                                as.matrix(Y),
                                object@probs)
            )

            methods::new("ResamplingInstance",
                train    = idx,
                strategy = object,
                n        = nrow(Y),
                seed     = seed)
          })

#' @rdname instantiate
setMethod("instantiate",
          "MonteCarlo",
          function(object, Y){

            .rngPreflight(object)
            seed <- .Random.seed

            idx <- .mc(object@k,
                       as.matrix(Y),
                       object@split)

            methods::new("ResamplingInstance",
                train    = idx,
                strategy = object,
                n        = nrow(Y),
                seed     = seed)
          })

#' @rdname instantiate
setMethod("instantiate",
          "BalancedMonteCarlo",
          function(object, Y){

            .rngPreflight(object)
            seed <- .Random.seed

            idx <- .mcBalanced(
              object@k,
              object@split,
              stratified = list(object@type,
                                as.matrix(Y),
                                object@probs)
            )

            methods::new("ResamplingInstance",
                train    = idx,
                strategy = object,
                n        = nrow(Y),
                seed     = seed)
          })

#' @rdname instantiate
setMethod("instantiate",
          "BalancedBootstrap",
          function(object, Y){

            .rngPreflight(object)
            seed <- .Random.seed

            idx <- .balanced_bootstrap_mc(
              object@k,
              object@split,
              stratified = list(object@type,
                                as.matrix(Y),
                                object@probs)
            )

            methods::new("ResamplingInstance",
                train    = idx,
                strategy = object,
                n        = nrow(Y),
                seed     = seed)
          })

#' @rdname train_idc
setMethod("train_idc", "ResamplingInstance", function(object) {
  object@train
})

######################
### MODEL FEATURES ###
######################

#' @rdname vip
setMethod("vip", "m8_model",
          function(object){

            if (!object@engine %in% c("opls", "pls"))
              stop("Method not available for this model")

            if (object@engine == "opls") {

              if(object@ctrl$ncomp_selected == 0)
                stop('No model components found')

              W <- t(as.matrix(object@fit$pred_comp$w_x))
              Tx <- as.matrix(object@fit$pred_comp$t_x)
              C <- object@fit$pred_comp$p_y

            } else if (object@engine == "pls") {

              W <- do.call(rbind, lapply(object@fit$comp, `[[`, "w_x"))
              Tx <- do.call(cbind, lapply(object@fit$comp, `[[`, "t_x"))
              C <- do.call(cbind, lapply(object@fit$comp, `[[`, "p_y"))

            }

            vip <- .ensure_matrix(.vip(Tx, W, C), ncol=ncol(W))

            return(vip)
          }
)

#' @rdname fitted
setMethod('fitted', "m8_model",
          function(object){
            if(object@engine == "opls"){

              if(object@ctrl$ncomp_selected == 0)
                stop('No model components found')

              return(object@fit$pred_comp$y_pred)
            }
          }
          )

#' @rdname weights
#' @param orth Logical indicating whether orthogonal scores should be returned
#'   (only applicable for OPLS models).
setMethod("weights","m8_model",
          function(object, orth = FALSE){


            if (!object@engine %in% c('opls', 'pls'))
              stop("Method not available for this model")


            if(isTRUE(orth) && object@engine != "opls")
              stop("Orthogonal weights not available for this model")


            key <- paste(if(isTRUE(orth) || is.numeric(orth)) "orth" else "pred",
                         "fit",
                         sep="_")

            if(object@engine == "opls"){

              if(object@ctrl$ncomp_selected == 0)
                stop('No model components found')

              map <- list(
                pred_fit = .ensure_matrix(object@fit$pred_comp$w_x, ncol=object@dims$p),
                orth_fit = object@fit$w_orth
              )} else if(object@engine == "pls"){
                map <- list(
                  pred_fit = .ensure_matrix(do.call(rbind, lapply(object@fit$comp, `[[`, "w_x")), ncol=object@dims$p)
                )
              }

            if(!key %in% names(map))
              stop("Requested weight type not available.")

            out <- map[[key]]

            if(is.numeric(orth)){

              if(orth > nrow(out))
                stop("Orthogonal component index exceeds available components")

              out <- out[orth, ,drop = FALSE]
            }


            return(out)
          }
    )


#' Model loadings
#' @param x An object of class \code{m8_model}.
#' @param orth Logical indicating whether orthogonal scores should be returned
#'   (only applicable for OPLS models).
#' @param ... Additional arguments (currently ignored).
#' @return Numeric vector or matrix containing loadings.
#' @examples
#' data(covid)
#' cv <- balanced_mc(k=5, split=2/3)
#' scaling <- UVScaling(center=TRUE)
#' model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#' show(model)
#' P <- loadings(model)
#' Po <- loadings(model, orth = TRUE)
#' dim(P)
#' dim(Po) == dim(P)
#' @export
setMethod("loadings","m8_model",
          function(x, orth = FALSE, ...){

            object <- x

            if (!object@engine %in% c('opls', 'pls', 'pca'))
              stop("Object needs to be of class m8_model.")

            if(isTRUE(orth) && object@engine != "opls")
              stop("Orthogonal loadings not available for this model")

            key <- paste(if(isTRUE(orth) || is.numeric(orth)) "orth" else "pred",
                         "fit",
                         sep="_")

            if(object@engine == "opls"){

              if(object@ctrl$ncomp_selected == 0)
                stop('No model components found')

              map <- list(
                pred_fit = .ensure_matrix(object@fit$pred_comp$p_x, ncol=length(object@fit$pred_comp$p_x)),
                orth_fit = object@fit$p_orth
              )} else if(object@engine == "pls"){
                map <- list(
                  pred_fit = do.call(rbind, lapply(object@fit$comp, `[[`, "p_x"))
                )
              }else if(object@engine == "pca"){
                map <- list(
                  pred_fit = object@fit$p
                )
              }

              out <- map[[key]]

              if(is.numeric(orth)){

                if(orth > nrow(out))
                  stop("Orthogonal component index exceeds available components")

                out <- out[orth, ,drop = FALSE]
              }


            return(out)


          })

#' @rdname scores
setMethod("scores","m8_model",
          function(object, orth = FALSE, cv=FALSE){

            if (!object@engine %in% c('opls', 'pls', 'pca'))
              stop("Method not available for this model")


            needs_orth <- isTRUE(orth) || is.numeric(orth)

            if(needs_orth && object@engine != "opls")
              stop("Orthogonal scores not available for this model")

            if(needs_orth && isTRUE(cv) && object@engine == "opls")
              stop("Cross-validated orthogonal scores not available for this model")

            if(isTRUE(cv) && object@engine %in% c("pls", "pca"))
              stop("Cross-validated scores not available for this model")

            key <- paste(if(isTRUE(orth) || is.numeric(orth)) "orth" else "pred",
                         if(cv)"cv"   else "fit",
                         sep="_")

            if(object@engine == "opls"){

              if(object@ctrl$ncomp_selected == 0)
                stop('No model components found')

              map <- list(
                pred_fit = .ensure_matrix(object@fit$pred_comp$t_x, ncol=1),
                orth_fit = object@fit$t_orth,
                pred_cv  = object@fit$pred_comp$t_cv,
                orth_cv  = object@fit$pred_comp$t_o_cv
              )} else if(object@engine == "pls"){
                map <- list(
                  pred_fit = do.call(cbind, lapply(object@fit$comp, `[[`, "t_x"))
                  # pred_cv  = object@fit$pred_comp$t_cv
                )
              }else if(object@engine == "pca"){
                map <- list(
                  pred_fit = object@fit$t
                )
              }


              out <- map[[key]]

              if(is.numeric(orth)){

                if(object@engine != "opls")
                  stop("Orthogonal components only available for OPLS")

                if(orth > ncol(out))
                  stop("Orthogonal component index exceeds available components")

                out <- out[, orth, drop = FALSE]
              }


           return(out)
          })

#' @rdname xres
setMethod("xres", "m8_model",
          function(object) {

            if (!identical(object@engine, "opls")) {
              stop("xres() is only available for OPLS models (engine = 'opls').")
            }

            if(object@ctrl$ncomp_selected == 0)
              stop('No model components found')

            Xcs <- object@fit$X_prepped
            if (is.null(Xcs)) stop("Missing fit$X_prepped.")

            To <- object@fit$t_orth
            Po <- object@fit$p_orth
            Tp <- .ensure_matrix(object@fit$pred_comp$t_x, ncol=1)
            Pp <- .ensure_matrix(object@fit$pred_comp$p_x, ncol=ncol(Po))

            if (is.null(To) || is.null(Po)) stop("Missing orthogonal components fit$t_o and/or fit$p_o.")
            if (is.null(Tp) || is.null(Pp)) stop("Missing predictive components fit$t_p and/or fit$p_p.")

            Xhat_o <- To %*% Po
            Xhat_p <- Tp %*% Pp

            Xcs - Xhat_p - Xhat_o
          }
)


############################
### MODEL SUMMARY / SHOW ###
############################

#' @describeIn m8_model-class Summarise model performance and component selection.
#' @param object An object of class \code{m8_model}.
#' @export
setMethod("summary", "m8_model",
          function(object) {

            eng  <- object@engine
            ctrl <- object@ctrl

            if (eng == "pls") {

              n_sel <- as.integer(ctrl$ncomp_selected)
              n_tst <- as.integer(ctrl$ncomp_tested)

              # Fallback if older objects exist
              if (is.na(n_sel) || length(n_sel) == 0L) n_sel <- as.integer(ctrl$nc)
              if (is.na(n_tst) || length(n_tst) == 0L) n_tst <- length(ctrl$q2)

              # Guard
              n_sel <- max(0L, n_sel)
              n_tst <- max(n_sel, n_tst, 0L)

              df <- data.frame(comp = seq_len(n_tst))

              if (identical(ctrl$type, "R")) {
                df$R2 <- ctrl$r2[seq_len(min(length(ctrl$r2), n_tst))]
                df$Q2 <- ctrl$q2[seq_len(min(length(ctrl$q2), n_tst))]
              } else {
                df$AUCt <- ctrl$aucs_tr[seq_len(min(length(ctrl$aucs_tr), n_tst))]
                df$AUCv <- ctrl$aucs_te[seq_len(min(length(ctrl$aucs_te), n_tst))]
              }

              # Mark which are actually selected (kept)
              df$selected <- df$comp <= n_sel

              # Optional: human-readable label
              df$fit <- ifelse(df$selected, "fitted", "not fitted")

              df$stop_reason <- NA_character_
              df$stop_reason[n_tst] <- ctrl$stop_reason
              df$stop_metric <- NA_real_
              df$stop_metric[n_tst] <- ctrl$stop_metric
              df$stop_delta <- NA_real_
              df$stop_delta[n_tst] <- ctrl$stop_delta


              return(methods::new("m8_modelSummary",
                         perf   = df,
                         engine = eng,
                         y_type = ctrl$type
              ))
            }

            if (eng == "opls") {

              n_sel <- as.integer(ctrl$ncomp_selected)-1 # combined 1/0 pred + x orth
              n_tst <- as.integer(ctrl$ncomp_tested)-1 # combined 1 pred + x orth

              if (is.na(n_tst) || length(n_tst) == 0L) n_tst <- length(ctrl$q2)
              if (is.na(n_sel) || length(n_sel) == 0L) n_sel <- n_tst

              n_sel <- max(0L, n_sel)
              n_tst <- max(n_sel, n_tst, 0L)

              df <- data.frame(
                comp_total = seq_len(n_tst)+1,
                comp_pred  = 1,
                comp_orth  = seq_len(n_tst)
              )

              if (identical(ctrl$type, "R")) {
                df$R2 <- ctrl$r2[seq_len(min(length(ctrl$r2), n_tst))]
                df$Q2 <- ctrl$q2[seq_len(min(length(ctrl$q2), n_tst))]
              } else {
                df$AUCt <- ctrl$aucs_tr[seq_len(min(length(ctrl$aucs_tr), n_tst))]
                df$AUCv <- ctrl$aucs_te[seq_len(min(length(ctrl$aucs_te), n_tst))]
              }

              df$R2X <- NA_real_
              if (!is.null(ctrl$r2x_comp) && length(ctrl$r2x_comp) > 0L) {
                idx <- seq_len(min(n_tst, length(ctrl$r2x_comp)))
                df$R2X[idx] <- ctrl$r2x_comp[seq_len(length(idx))]
              }

              df$selected <- df$comp_total <= n_sel+1
              df$fit <- ifelse(df$selected, "fitted", "not fitted")

              df$stop_reason <- NA_character_
              df$stop_metric <- NA_real_
              df$stop_delta  <- NA_real_

              if (!is.null(ctrl$stop_reason) && n_tst >= 1L) {
                df$stop_reason[n_tst] <- ctrl$stop_reason
                df$stop_metric[n_tst] <- ctrl$stop_metric
                df$stop_delta[n_tst]  <- ctrl$stop_delta
              }

              return(methods::new("m8_modelSummary",
                         perf   = df,
                         engine = eng,
                         y_type = ctrl$type
              ))
            }

            if (eng == "pca") {

              expl <- ctrl$r2x_comp
              if (is.null(expl)) {
                stop("PCA summary requires ctrl$r2x_comp (proportion variance explained per component).")
              }

              nc <- length(expl)
              df <- data.frame(
                comp   = seq_len(nc),
                var    = as.numeric(expl),
                cumvar = cumsum(as.numeric(expl))
              )


              return(methods::new("m8_modelSummary",
                         perf   = df,
                         engine = eng,
                         y_type = "unsupervised"
              ))
            }

            methods::new("m8_modelSummary",
                perf   = data.frame(),
                engine = eng,
                y_type = if (!is.null(ctrl$type)) ctrl$type else NA_character_
            )
          }
)



#' @describeIn m8_model-class Show a compact model header.
#' @param object An object of class \code{m8_model}.
#' @export
setMethod("show", "m8_model",
          function(object) {

            cat("\n")
            cat("m8_model <", object@engine, ">\n", sep = "")
            cat(strrep("-", 40), "\n", sep = "")

            if (!is.null(object@dims)) {
              cat("Dimensions  :", object@dims$n, "samples x",
                  object@dims$p, "variables\n")
            }

            y_type <- tryCatch(object@ctrl$type, error = function(e) NULL)
            if (!is.null(y_type)) {
              if (object@engine == "pca") {
                cat("Mode        : unsupervised\n")
              } else {
                cat("Mode        :", if (y_type == "R") "regression" else "classification", "\n")
              }
            }

            prep <- tryCatch(object@prep, error = function(e) NULL)
            if (!is.null(prep)) {
              cat("Preprocess  :",
                  ifelse(prep@center, "center", "no-center"),
                  "|",
                  prep@scale, "\n")
            }

            ctrl <- object@ctrl

            n_selected <- tryCatch(ctrl$ncomp_selected, error = function(e) NULL)
            n_tested   <- tryCatch(ctrl$ncomp_tested,   error = function(e) NULL)
            if (is.null(n_selected)) {
              n_selected <- tryCatch(ctrl$nc, error = function(e) NULL)
            }

            if (!is.null(n_selected)) {
              cat("Components  :", n_selected)
              if (!is.null(n_tested) && n_tested > n_selected) {
                cat(" (", n_tested, " tested)", sep = "")
              }
              cat("\n")
            }

            if (!is.null(object@cv)) {
              cv_method <- tryCatch(class(object@cv@strategy)[1], error = function(e) NULL)
              cv_k <- tryCatch(object@cv@strategy@k, error = function(e) NULL)
              if (!is.null(cv_method)) {
                line <- paste0("Validation  : ", cv_method)

                if (!is.null(cv_k)) {
                  line <- paste0(line, " (k = ", cv_k, ")")
                }

                cat(line, "\n")
              }
            }

            stop_reason <- tryCatch(ctrl$stop_reason, error = function(e) NULL)
            if (!is.null(stop_reason)) {
              cat("Stop rule   :", stop_reason, "\n")
            }

            cat(strrep("-", 40), "\n", sep = "")
            cat("Use summary() for performance metrics.\n\n")
          }
)


