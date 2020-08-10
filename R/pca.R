#' Perform Principal Component Analysis
#' @export
#' @description This function is used to perform Principal Component Analysis (PCA).
#' @param X Numeric input matrix with each row representing an observation and each column a metabolic feature.
#' @param pc Desired number of principal components.
#' @param scale Desired scaling method: \code{None}, \code{UV} (unit variance) or \code{Pareto} (Pareto scaling).
#' @param method Algorithm for computing PCA. NIPALS is standard and usually fine, see Details for other methods.
#' @details Other PCA methods: These rely on implementation inf the pcaMethods R package and include 'svd', 'ppca', 'svdImpute', 'robustPca', 'nlpca'. For complete list of available methods see pca function documentation of the \code{pcaMethods} package.
#' @param center Logical indicating if data should be mean centered.
#' @references Geladi, P and Kowalski, B.R. (1986), Partial least squares and regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1-17.
#' @return This function returns a \emph{PCA_MetaboMate} S4 object.
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom pcaMethods pca

pca <- function(X, pc = 2, scale = "UV", center = T, method = "nipals") {

  {

    if(is.data.frame(X)) X=as.matrix(X)

    pc_max=min(c(ncol(X), min(nrow(X))))
    if(pc>=pc_max){ message(paste('Too many number of components, setting pc to', pc_max)); pc=pc_max }

    if(!is.logical(center)){stop('Check center parameter argument!')}
    sc_num<-switch(scale,
                   'None'={0},
                   'UV'={1},
                   'Pareto'={2})
    if(is.null(sc_num)){stop('Check scale parameter argument!')}

    x_check<-.checkXclassNas(X)
    msd_x<-.sdRcpp(X); # returns list with elements mean and sd
    XcsTot<-.scaleMatRcpp(X, 0:(nrow(X)-1), center=TRUE, scale_type = sc_num)[[1]]
  }

  if (method == "nipals") {
    res <- list()
    for (i in 1:pc) {
      if (i == 1) {
        res[[i]] <- .nipPcaCompRcpp(XcsTot)
      } else (res[[i]] <- .nipPcaCompRcpp(X = res[[i - 1]][[1]]))
    }
    Tpc <- sapply(res, "[[", 2)
    Ppc <- sapply(res, "[[", 3)
    # total.var<-sum(diag(cov(X))) #Calculate total variance in
    tssx<-.tssRcpp(XcsTot)
    #Calculate proportion of variance explained and cumulative
    ss_comp <- rep(NA, ncol(Tpc))
    for (i in 1:ncol(Tpc)) {
      ss_comp[i] <- .tssRcpp(Tpc[, i] %o% Ppc[, i]) / tssx
    }

    pars <- list(
      'center'= center,
      'scale'= scale,
      'nPC'= pc,
      'R2' = ss_comp,
      'tssx'= tssx
    )

    mod_pca <- new("PCA_metabom8", type ='NIPALS', t = Tpc, p = Ppc, nPC = pc, X_mean=msd_x[[1]], X_sd=msd_x[[2]], Parameters=pars,  X=X )
  } else {

    mod <- pcaMethods::pca(X, nPcs = pc, scale = "none", center = F, method = method)
    r2 <- mod@R2cum
    for (i in 1:pc) {
      if (i == 1) {
        next
      } else {
        r2[i] <- r2[i] - cumsum(r2[1:(i - 1)])
      }
    }

    pars=list(
      'center'= center,
      'scale'= scale,
      'nPC'= pc,
      'R2' = r2,
      'tssx'= tssx
    )

    mod_pca <-  new("PCA_metabom8", type =method, t = mod@scores, p = mod@loadings, nPC = pc, X_mean=msd_x[[1]], X_sd=msd_x[[2]], Parameters=pars,  X=X)
  }
  return(mod_pca)
}



