#' @title Plot PCA or O-PLS loadings
#' @export
#' @description Visualisation of loadings using scatter plot
#' @param obj PCA or OPLS object from packages metabom8 or MetaboMate
#' @param pc, num vector descrbing which components to visualise in on abcissa and ordinate (in case of OPLS the first element in pc is ignored as abcissa always shows predicted component)
#' @param an list, three elements describing colour, point shape and point labels (all is aes context)
#' @param title Plot title
#' @param legend logical, indicating if legend should be included (legend takes up too much space when high number of categorical varaibles)
#' @param ... Optional arguments passed to \code{geom_point()}
#' @references Geladi, P and Kowalski, B.R. (1986), Partial least squares and regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1-17.
#' @return data.frame containing H.T2 ellipse cooredinates
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom ellipse ellipse
#' @importFrom ggplot2 aes_string geom_hline geom_vline geom_point labs scale_x_continuous scale_y_continuous theme_bw geom_label_repel guides theme guides
#' @importFrom stats cov
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @family dataviz
#' @section

## loadiings function
plotload_cat=function(obj, pc, an, title=NULL, legend=TRUE, ...){

  if(missing(pc)){
    if(grepl('PCA', class(obj)[1])){pc=c(1,2)}
    if(grepl('PLS', class(obj)[1])){pc=c('1', 'o1')}
  }
  if (missing(an)) {
    if(class(obj)[1]=='OPLS_metabom8') {if(ncol(obj@X)<150) an <- list('All', 'All', colnames(obj@X))}
    if(class(obj)[1]=='PCA_metabom8') {if(ncol(obj@X)<150) an <- list('All', 'All', colnames(obj@X))}
  }
  res=.viz_df_helper(obj, pc, an, type='p')
  if( is.null(res$an_le) || is.na(res$an_le) ) {stop('Check helper function .viz_df_helper')}
  melted=res$df
  if(res$an_le==1){
    ### ploting
    d=ggplot(data=melted, aes_string(x= colnames(melted)[1], y= colnames(melted)[2], colour=names(an)[1]))+
      geom_hline(yintercept = 0, colour = "gray70") +
      geom_vline(xintercept = 0, colour = "gray70") +
      geom_point(...)+
      labs(colour=names(an)[1], title = title, x=paste(colnames(melted)[1], sep="" ), y=paste(colnames(melted)[2], sep="" ))+
      theme_bw()
  }

  if(res$an_le==2){
    d=ggplot(data=melted, aes_string(x= colnames(melted)[1], y= colnames(melted)[2], colour=colnames(melted)[3], shape=colnames(melted)[4]))+
      geom_hline(yintercept = 0, colour = "gray70") +
      geom_vline(xintercept = 0, colour = "gray70") +
      geom_point(...)+
      labs(shape=names(an)[2], colour=names(an)[1], title = title, x=paste(colnames(melted)[1], sep="" ), y=paste(colnames(melted)[2], sep="" ))+
      theme_bw()
  }

  if(res$an_le==3){
    d=ggplot(data=melted, aes_string(x= colnames(melted)[1], y= colnames(melted)[2], colour=colnames(melted)[3], shape=colnames(melted)[4]))+
      geom_hline(yintercept = 0, colour = "gray70") +
      geom_vline(xintercept = 0, colour = "gray70") +
      geom_point(...)+
      labs(shape=names(an)[2], colour=names(an)[1], title = title, x=paste(colnames(melted)[1], sep="" ), y=paste(colnames(melted)[2], sep="" ))+
      geom_label_repel( aes_string(label=colnames(melted)[5]))+
      theme_bw()
  }


  if(legend==FALSE){
    d=d+guides(colour=FALSE)+theme(panel.grid.minor = element_blank())
  }


  return(d)
}



