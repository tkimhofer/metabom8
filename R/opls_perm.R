#' @title OPLS model Y-permutations
#' @export
#' @description Model validation using Y-permutations
#' @param smod OPLS model of the package \emph{metabom8}
#' @param n, num Number of permutations
#' @param plot logical, indicating if results should be visualised
#' @param mc logical, indicating if tasked should be parallelised using multiple cores
#' @return  data.frame with perutation indices
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom ggplot2 ggplot aes_string theme_bw labs
#' @importFrom reshape2 melt
#' @importFrom plyr ddply geom_point
#' @importFrom stats median
#' @importFrom parallel detectCores
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @family NMR ++
#' @examples
#' data(covid)
#' model=opls(X, Y=an$type)
#' perm=opls_perm(model)
#' @section

opls_perm <- function(smod, n = 10, plot = TRUE, mc = FALSE) {
    
    
    # how many orth components
    nc_o <- nrow(smod@summary) - 1
    
    
    Xs <- smod@X  # metabom8:::.scaleMatRcpp(mod@X, seq(nrow(mod@X))-1, center=mod@Parameters$center, scale_type=sc_num)$X_prep
    Y <- smod@Y$dummy
    
    cv <- list()
    cv$cv_sets <- smod@Parameters$cv$cv_sets
    cv$k <- length(cv$cv_sets)
    cv$method <- smod@Parameters$cv$method
    
    type <- smod@type
    
    # browser() if(mc==TRUE){ numCores <- detectCores()-2 out=mclapply(seq(n),
    # do_perm, Y, cv, type, nc_o, mc.cores=numCores) }else{
    out <- lapply(seq(n), function(i) {
        Y_perm <- t(t(sample(Y[, 1])))
        ind <- .permYmod(Xs, Y_perm, cv, type, nc_o)
        list(ind = ind, r = cor(Y, Y_perm)[1, 1])
    })
    # }
    
    out_df <- as.data.frame(t(vapply(out, unlist, FUN.VALUE = rep(1, 5))))
    out_df <- data.frame(apply(out_df, 2, function(x) as.numeric(unlist(x))), stringsAsFactors = FALSE)
    colnames(out_df) <- c("r2_comp", "q2_comp", "aucs_tr", "aucs_te", "r")
    out_df$model <- paste0("perm_", seq_len(nrow(out_df)))
    mod <- .permYmod(Xs, Y, cv, type, nc_o)
    add <- as.data.frame(t(c(unlist(mod), 1)))
    add$model <- "Non-permuted"
    colnames(add) <- colnames(out_df)
    
    out_df <- rbind(out_df, add)
    
    idx_rm <- apply(out_df, 2, function(x) {
        all(is.na(x))
    })
    
    out_df <- out_df[, !idx_rm]
    
    
    
    
    if (plot == TRUE) {
        
        # x.val=out_df$aucs_te[out_df$model!='Non-permuted'] dn=pnorm(seq(0,1, by=0.01),
        # mean(x.val), sd(x.val))
        
        
        dd <- melt(out_df, id.vars = c("model", "r"))
        dd$variable <- as.character(dd$variable)
        
        map_parameter <- c(aucs_tr = "Train", aucs_te = "Test", r2_comp = "R2", q2_comp = "Q2")
        ylab <- ifelse(any(grepl("aucs", dd$variable)), "AUROC", expression(R^2 ~ 
            "/" ~ Q^2))
        
        dd$variable <- map_parameter[match(dd$variable, names(map_parameter))]
        
        out_df$model <- gsub("_.*", "", out_df$model)
        
        add_summary <- data.frame(t(apply(out_df[out_df$model == "perm", -ncol(out_df)], 
            2, median)), stringsAsFactors = FALSE)
        add_summary$model <- "perm"
        
        
        
        add_summary <- rbind(add_summary, out_df[out_df$model == "Non-permuted", 
            ])
        
        g <- ggplot(dd, aes_string(x = "abs(r)", y = "value", colour = "variable")) + 
            geom_point() + theme_bw() + labs(x = "|r|", y = ylab, colour = "Sample Set", 
            caption = paste0("O-PLS-", type, ", ", n, " Y-permutations."))
        
        plot(g)
        
    }
    
    return(out_df)
    
}

# .do_perm <- function(i, Y, cv, type, nc_o){ Y_perm=t(t(sample(Y[,1])))
# ind=.permYmod(Xs, Y_perm, cv, type, nc_o) list(ind=ind, r=cor(Y, Y_perm)[1,1])
# }




























































