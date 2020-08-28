# Function performs STOCSY analysis (NMR)
# Input: driver peak (ppm value), ppm vector, NMR matrix with samples in rows and variables in columns, TRUE/FALSE for stocsy graph
# output: data frame with correlation, covariance and ppm value OR if plotting=TRUE: list of length two with [[1]] data frame and [[2]] ggplot2 object
# R-code written by Torben Kimhofer, Imperial College London (15/03/2015)


#' @title Statistical Total Correlation Spectroscopy (STOCSY)
#' @param X  Numeric matrix or dataframe where each row represents a NMR spectrum and each column a chemical shift variable.
#' @param ppm num array of chemical shift variables, matched to columns in X
#' @param driver num, chemical shift indicating STOCSY driver signal
#' @param plotting logical, indicing if results should be plotted
#' @param title char, plot title
#' @param plotly logic, interactive plottin using plotly
#' @return S4 object of class stocsy1d_metabom8
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @family metID functions
#' @family nmr functions
#' @importFrom ggplot2 ggplot aes_string geom_line geom_v_line scales_x_reverse scale_colour_gradientn labs theme_bw theme element_text
#' @importFrom colorRamps matlab.like2
#' @importFrom scales sapply breaks_pretty
#' @importFrom methods hasArg
#' @export
#' @section

stocsy<-function(X, ppm, driver, plotting=TRUE, title=NULL, plotly=TRUE){

  if(!hasArg(ppm)){
    ppm<-as.numeric(colnames(X))
    if(!any(!is.na(ppm))){
      stop('Provide ppm argument')
    }
  }

  if(! 'matrix' %in% is(X) & is.numeric(X)){stop('X argument is not a numeric matrix')}
  if(!is.numeric(driver) || !( driver[1] <= max(ppm) & driver[1] >= min(ppm)) ) {stop('STOCSY driver not numeric or outside ppm range')}
  if(!is.logical(plotting)){plotting=TRUE}


  if (length(driver)==nrow(X)){
    l=driver
    ext_driver=TRUE
  }else{
    ext_driver=FALSE
    idx=which.min(abs(ppm-driver))
    l=X[,idx]
  }

  browser()
  cc=apply(X, 2, function(x, ll=l) {cor(x,ll)})
  cv=apply(X, 2, function(x, ll=l) {cov(x,ll)})
  df=data.frame(cc, cv, ppm)

  stoc_mod <- new("stocsy1d_metabom8",
                  version = "0.9",
                  X = X,
                  ppm = ppm,
                  driver = driver,
                  r = cc,
                  cov=cv
  )
  if(plotting){
    g1=ggplot(df, aes_string(x='ppm', y='cv', colour='abs((cc))'))+
      geom_line()+
      scale_x_reverse(breaks=breaks_pretty())+
      scale_colour_gradientn(colours = matlab.like2(10), limits=c(0,1), name=paste("r"))+
      labs(x=expression(delta~{}^1*H~(ppm)), y=paste("cov"), title=title)+
      theme_bw()+
      theme(axis.text=element_text(colour="black"))+
      geom_vline(xintercept = driver, linetype=2, col='grey')

    if(plotly){
      ggplotly(g1, layerData = 2)
    }
  }
  return(stoc_mod)
}



#' @title Plot Statistical Total Correlation Spectroscopy (STOCSY)
#' @param stoc_mod STOCSY object, as created with function 'stocsy' ('see ?stocsy()')
#' @param shift num array, chemical shift range defining plotting area
#' @param title char, plot title
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @family metID functions
#' @family NMR functions
#' @importFrom ggplot2 ggplot aes_string geom_line geom_v_line scales_x_reverse scale_colour_gradientn labs theme_bw theme element_text
#' @importFrom colorRamps matlab.like2
#' @export
#' @section

plotStocsy=function(stoc_mod, shift=c(0,10), title=NULL){

  if(!'stocsy1d_metabom8' %in% is(stoc_mod)) stop('Provide a STOCSY object')
  ds=data.frame(r=stoc_mod@r, cov=stoc_mod@cov, ppm=stoc_mod@ppm, stringsAsFactors = FALSE)
  if(!all(is.numeric(shift))) stop('Check shift argument')
  if(min(shift)<min(ds$ppm)){shift[which.min(shift)]=min(ds$ppm)}
  if(max(shift)>max(ds$ppm)){shift[which.max(shift)]=max(ds$ppm)}

  shift=sort(shift)
  idx=get.idx(shift, ds$ppm)
  if(length(idx)==0) stop('Check shift argument')
  ds=ds[idx,]
  g1=ggplot(ds, aes_string(x='ppm', y='cov', colour='abs(r)'))+
    geom_line()+
    scale_x_reverse(breaks=breaks_pretty())+
    scale_colour_gradientn(colours = matlab.like2(10), limits=c(0,1), name=paste("r(d=", stoc_mod@driver, ', X)', sep=''))+
    labs(x=expression(delta~{}^1*H~(ppm)), y=paste("cov(d=", stoc_mod@driver, ', X)', sep=''), title=title)+
    theme_bw()+
    theme(axis.text=element_text(colour="black"))+
    ggtitle(title)

  if(shift[1]<stoc_mod@driver & shift[2]>stoc_mod@driver){
    g1<-g1+geom_vline(xintercept =stoc_mod@driver, linetype=2, col='grey')
  }

  return(g1)

}


























































































































