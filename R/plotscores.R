
#' @title Plotting PCA or OPLS model scores
#' @export
#' @aliases plotscores
#' @description Function for plotting PCA, PLS or OPLS model scores.
#' @param model PCA, PLS or OPLS model of the package \emph{metabom8} or \emph{MetaboMate}.
#' @param pc Model (principal) components to be plotted on abscissa and ordinate (see Details).
#' @param an Scatter colour, shape and label specification (see Details).
#' @param title Plot title.
#' @param qc row indices of QC samples (can be left NA). If specified, QC samples will be highlighted in the plot.
#' @param legend Position of the plot legend, set to NA if legend should be outside of the plotting area.
#' @param cv.scores Logical value (TRUE or FALSE) indicating if cross-validated scores should be plotted for the predictive component(s) (only for PLS and O-PLS).
#' @param ... Additional values passed on to \code{\link[ggplot2]{scale_colour_gradientn}}.
#' @details Scores colouring is specified with the argument \code{an}, which is a list of three elements. The first list element specifies the colour (class factor required for categorical variables), the second list element specifies a point labeling (class character or factor) and the third list element specifies point shape. The Hotelling's \out{T<sup>2</sup>} ellipse is automatically included and calculated for the dimensions selected by the \code{pc} argument.
#' @references Trygg J. and Wold, S. (2002) Orthogonal projections to latent structures (O-PLS). \emph{Journal of Chemometrics}, 16.3, 119-128.
#' @references Hotelling, H. (1931) The generalization of Student’s ratio. \emph{Ann. Math. Stat.}, 2, 360-378.
#' @return This function returns a \emph{ggplot2} plot object.
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @importFrom stats cov
#' @importFrom ellipse ellipse
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot aes geom_polygon geom_hline geom_vline geom_point scale_colour_manual ggtitle labs theme labs scale_colour_gradientn scale_x_continuous unit
#' @importFrom colorRamps matlab.like2
#' @importFrom ggrepel geom_text_repel
#' @importFrom graphics plot
#' @importFrom reshape2 melt
#' @importFrom scales breaks_pretty
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @section


plotscores1 <- function(obj, pc = c(1, 2), an, title = "", qc, legend = "in", cv.scores = T, ...) {

  #if (missing(qc)) an <- list("NA")

  if(cv.scores==T){etype='t_cv'} else{etype='t'}
  res=.viz_df_helper(obj, pc, an, type=etype)
  sc=res$df

  type='categorical'
  if(is.numeric(sc[,3])){type='continuous'}


  if( is.null(res$an_le) || is.na(res$an_le) ) {stop('Check helper function .viz_df_helper')}


  # Calculate Hotellings T2 elipse
  df=.hotellingsT2(x=sc[,1], y=sc[,2])


  #browser()


  # y.labb <- diff(range(ylim)) * 0.015
  # x.labb <- diff(range(xlim)) * 0.015
  #
  #
  # df1 <- sc
  # df1[, 1] <- df1[, 1] + x.labb
  # df1[, 2] <- df1[, 2] + y.labb

  #browser()

  g <- ggplot() +
    geom_polygon(data = df, aes_string(x = "V1", y = "V2"), fill = NA, alpha = 0.4, colour='black', linetype=3) +
    geom_hline(yintercept = 0, colour = "gray70") +
    geom_vline(xintercept = 0, colour = "gray70")+
    theme_bw()+
    labs(caption = expression('Ellipse: Hotelling\'s T'^2~'(95%)'))


  if(res$an_le==1){

    switch(type,
           'categorical' = {
             g <- g +
             geom_point(data = sc, aes_string(colnames(sc)[1], colnames(sc)[2], colour = colnames(sc)[3]), alpha = 0.7, shape = 16)
           },
           'continuous' = {
             g <- g +
               geom_point(data = sc, aes_string(colnames(sc)[1], colnames(sc)[2], colour = colnames(sc)[3]), alpha = 1,  shape = 16) +
               scale_colour_gradientn(colours = matlab.like2(10), breaks = breaks_pretty(), ...)
           }
    )


    g <- g +
      labs(title=title, colour = colnames(sc)[3])

    }


  if(res$an_le==2){

    switch(type,
           'categorical' = {
             g <- g +
               geom_point(data = sc, aes_string(colnames(sc)[1], colnames(sc)[2], colour = colnames(sc)[3], shape=colnames(sc)[4]), alpha = 1)
           },
           'continuous' = {
             g <- g +
               geom_point(data = sc, aes_string(colnames(sc)[1], colnames(sc)[2], colour = colnames(sc)[3], shape=colnames(sc)[4]), alpha = 1) +
               scale_colour_gradientn(colours = matlab.like2(10), breaks = breaks_pretty(), ...)
           }
    )


    g <- g +
      labs(title=title, colour = colnames(sc)[3], shape = colnames(sc)[4])

  }



  if(res$an_le==3){

    switch(type,
           'categorical' = {
             g <- g +
               geom_point(data = sc, aes_string(colnames(sc)[1], colnames(sc)[2], colour = colnames(sc)[3], shape=colnames(sc)[4]), alpha = 1)+
               geom_text_repel(data = sc, aes_string(colnames(sc)[1], colnames(sc)[2], label = colnames(sc)[5]))
           },
           'continuous' = {
             g <- g +
               geom_point(data = sc, aes_string(colnames(sc)[1], colnames(sc)[2], colour = colnames(sc)[3], shape=colnames(sc)[4]), alpha = 1) +
               scale_colour_gradientn(colours = matlab.like2(10), breaks = breaks_pretty(), ...)+
               geom_text_repel(data = sc, aes_string(colnames(sc)[1], colnames(sc)[2], label = colnames(sc)[5]))
           }
    )

    g <- g +
      labs(title=title, colour = colnames(sc)[3], shape = colnames(sc)[4], label=colnames(sc)[5])

  }


  # prep sub-df for QC samples
  if (!missing(qc) && all(!is.na(qc)) && is.numeric(qc)) {
    df_qc <- data.frame(x = sc[qc,1], y = sc[qc,1])
    g <- g +
      geom_point(data = df.qc, aes_string("x", "y"), alpha = 0.7, colour = "black")
  }




  ###### AXIS LABELLING

  switch(class(obj),
         'PCA_MetaboMate' = {
           g <- g +
             scale_x_continuous(name = paste("PC ", pc[1], " (", round(obj@R2[pc[1]] * 100, 1), " %)", sep = "")) +
             scale_y_continuous(name = paste("PC ",pc[2], " (", round(obj@R2[pc[2]] * 100, 1), " %)", sep = ""))
         },
         'OPLS_MetaboMate' = {
           if (ncol(obj@t_orth) > 1) {
             comp <- "orthogonal components"
           } else {
             comp <- "orthogonal component"
           }
           if (cv.scores == F) {
             if (obj@type == "DA") {

               g <- g +
                 scale_x_continuous(name = expression(t[pred])) +
                 scale_y_continuous(name = expression(t[orth])) +
                 ggtitle(paste("OPLS-DA ", obj@nPC, " ", comp, " (R2X=", round(obj@summary$R2X[obj@nPC], 2), ", R2Y=", round(obj@summary$R2Y[obj@nPC], 2), ", Q2=",
                               round(obj@summary$Q2[obj@nPC], 2), ", CV-AUROC=", round(obj@summary$AUROC[obj@nPC], 2), ")", sep = "")) +
                 theme(plot.title = element_text(size = 10))


             } else {
               g <- g +
                 scale_x_continuous(name = expression(t[pred])) +
                 scale_y_continuous(name = expression(t[orth])) +
                 ggtitle(paste("OPLS-R ",obj@nPC, " ", comp, " (R2X=", round(obj@summary$R2X[obj@nPC], 2), ", R2Y=", round(obj@summary$R2Y[obj@nPC], 2), ", Q2=", round(obj@summary$Q2[obj@nPC], 2), ")", sep = "")) +
                 theme(plot.title = element_text(size = 10))
             }
           } else {
             if (obj@type == "DA") {
               g <- g +
                 scale_x_continuous(name = expression(t[pred][","][cv])) +
                 scale_y_continuous(name = expression(t[orth][","][cv])) +
                 ggtitle(paste("OPLS - ", obj@nPC, " ", comp, " (R2X=", round(obj@summary$R2X[obj@nPC], 2), ", R2Y=", round(obj@summary$R2Y[obj@nPC], 2), ", Q2=", round(obj@summary$Q2[obj@nPC], 2), ", AUROC=", round(obj@summary$AUROC[obj@nPC], 2), ")", sep = "")) +
                 theme(plot.title = element_text(size = 10))

             } else {
               g <- g +
                 scale_x_continuous(name = expression(t[pred][","][cv])) +
                 scale_y_continuous(name = expression(t[orth][","][cv])) +
                 ggtitle(paste("OPLS - ", obj@nPC, " ", comp, " (R2X=", round(obj@summary$R2X[obj@nPC], 2), ", R2Y=", round(obj@summary$R2Y[obj@nPC], 2), ", Q2=", round(obj@summary$Q2[obj@nPC], 2), ")", sep = "")) +
                 theme(plot.title = element_text(size = 10))
             }
           }
         },

         'OPLS_metabom8' = {
           if (ncol(obj@t_orth) > 1) {
             comp <- "orthogonal components"
           } else {
             comp <- "orthogonal component"
           }
           if (cv.scores == F) {
             if (obj@type == "DA") {

               g <- g +
                 scale_x_continuous(name = expression(t['pred,cv'])) +
                 scale_y_continuous(name = expression(t['orth,cv'])) +
                 ggtitle(paste("OPLS-DA: ", obj@nPC-1, " ", comp, " (R2X=", round(obj@summary$R2X[obj@nPC-1], 2), ", AUROC=", round(obj@summary$AUROC[obj@nPC-1], 2),
                               ", CV-AUROC=", round(obj@summary$AUROC_CV[obj@nPC-1], 2), ")", sep = "")) +
                 theme(plot.title = element_text(size = 10))


             } else {
               g <- g +
                 scale_x_continuous(name = expression(t[pred])) +
                 scale_y_continuous(name = expression(t[orth])) +
                 ggtitle(paste("OPLS-R ",obj@nPC-1, " ", comp, " (R2X=", round(obj@summary$R2X[obj@nPC-1], 2), ", R2Y=", round(obj@summary$R2Y[obj@nPC-1], 2), ", Q2=", round(obj@summary$Q2[obj@nPC-1], 2), ")", sep = "")) +
                 theme(plot.title = element_text(size = 10))
             }
           } else {
             if (obj@type == "DA") {
               g <- g +
                 scale_x_continuous(name = expression(t[pred][","][cv])) +
                 scale_y_continuous(name = expression(t[orth][","][cv])) +
                 ggtitle(paste("OPLS-DA: ", obj@nPC-1, " ", comp, " (R2X=", round(obj@summary$R2X[obj@nPC-1], 2), ", AUROC=", round(obj@summary$AUROC[obj@nPC-1], 2),
                               ", CV-AUROC=", round(obj@summary$AUROC_CV[obj@nPC-1], 2), ")", sep = "")) +
                 theme(plot.title = element_text(size = 10))

             } else {
               g <- g +
                 scale_x_continuous(name = expression(t[pred][","][cv])) +
                 scale_y_continuous(name = expression(t[orth][","][cv])) +
                 ggtitle(paste("OPLS-R ",obj@nPC-1, " ", comp, " (R2X=", round(obj@summary$R2X[obj@nPC-1], 2), ", R2Y=", round(obj@summary$R2Y[obj@nPC-1], 2), ", Q2=", round(obj@summary$Q2[obj@nPC-1], 2), ")", sep = "")) +
                 theme(plot.title = element_text(size = 10))
             }
           }
         },









         )
  if (legend == "in") {
    g <- g + theme(legend.position = c(1.01, 1.02), legend.justification = c(1, 1))
  }
  return(g)
}







































































































































































































































