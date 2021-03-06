#' @title Eruption plot for OPLS-DA model visualisation

#' @export

#' @description Visualise variable importances in OPLS context, p value is generated with Kruskal Walis rank sum test, effect size measure is Cliff's delta.

#' @param mod OPLS model generated with \code{metabom8}.

#' @param pc num or char, indicating OPLS component to visualise (1 for predictive component, prefix 'o' for orthogonal components, e.g., pc='o1' for first orthogonal component)

#' @param p_adj str, p value adjustment method, see \code{?p.adjust}

#' @param invert_es logical, swap reference distribution for calculation of effect size estimate Cliff's delta

#' @return dataframe of effect size estimate (Cd), model loadings of OPLS predictive component (p_pred), p value (adjusted, if specified)

#' @seealso \code{\link{opls}}

#' @importFrom ggplot2 ggplot aes_string geom_point scale_colour_gradientn geom_label_repel theme_bw theme element_blank element_text geom_segment scale_x_continuous

#' @importFrom colorRamps matlab.like2

#' @importFrom stats p.adjust

#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}

#' @section



eruption <- function(mod, pc = 1, p_adj = "BH", invert_es = FALSE) {





    # check if mod is discriminant analysis



    if (is.na(p_adj) || is.infinite(p_adj) || length(p_adj) > 1) {

        p_adj <- "none"

    }



    if (length(unique(mod@Y$ori)) > 2) {

        stop("Eruption plot defined for two-level outcome.")

    }



    if (is.na(pc) || is.infinite(pc) || length(pc) > 1)

        stop("Check pc argument.")



    # if(class(mod)='PCA_metabom8'){ if(pc>nrow(mod@p)){stop('PC value exceeds number

    # of principal components.')} ddl=data.frame(x=mod@p[,pc], id=colnames(mod@X)) }



    if (grepl("o", pc)) {

        pc1 <- as.numeric(gsub("o", "", pc))

        if (is.na(pc1) || is.infinite(pc1) || pc1 > nrow(mod@p_orth)) {

            stop("Check pc argument and help section.")

        }

        ddl <- data.frame(x = mod@p_orth[pc1, ], id = colnames(mod@X))

    } else {

        ddl <- data.frame(x = mod@p_pred[1, ], id = colnames(mod@X))

    }



    Y <- mod@Y$ori

    if (invert_es) {

        comp <- unique(Y)[2]

        cat('Using ', unique(Y)[1], ' as reference group for Cliff\'s delta. Swap reference / comparator with invert_es=TRUE (see function help).')

    } else {

        comp <- unique(Y)[1]

        cat('Using ', unique(Y)[2], ' as reference group for Cliff\'s delta.')

    }



    uni <- t(apply(mod@X, 2, function(x, idx = which(Y == comp), y = Y) {

        c(es_cdelta(x[idx], x[-idx]), kruskal.test(x, y)$p.value)

    }))



    ddl <- cbind(ddl, uni)

    colnames(ddl) <- c("p1", "id", "Cd", "pval")

    ddl$p1 <- abs(ddl$p1)

    #browser()



    adj_log <- p_adj %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",

                            "fdr")



    if (adj_log) {

        ddl$pval_adjusted <- p.adjust(ddl$pval, method = p_adj)

        ddl$pval_transformed <- abs(log10(ddl$pval_adjusted))

    }else{

        ddl$pval_transformed <- abs(log10(ddl$pval))

    }



    gl2 <- ggplot(ddl, aes_string(x = "Cd", y = "p1", colour = "pval_transformed")) + geom_label_repel(aes_string(label = "id"),

                                                                                                       colour = "black", min.segment.length = 0.001) + geom_point(size = 3, shape = 1) +

        geom_point(size = 3, shape = 16, alpha = 0.3) + theme_bw() + scale_x_continuous(limits = c(-1,

                                                                                                   1)) + scale_colour_gradientn(colours = matlab.like2(10)) + labs(x = "Cliff's Delta",

                                                                                                                                                                   y = paste0("|p_", pc, "|")) + theme(panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold",

                                                                                                                                                                                                                                                                   size = 25), legend.position = "top", legend.direction = "horizontal")



    if (adj_log) {

        gl2 <- gl2 + labs(colour = expression("| log10( p value"["adj"] ~ ") |"))

    } else {

        gl2 <- gl2 + labs(colour = "| log10( p value ) |")

    }



    ddl=ddl[, !colnames(ddl) %in% 'pval_transformed']





    return(list(data=ddl, plot=gl2))



}
