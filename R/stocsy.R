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
#' @return S4 object of class stocsy1d_metabom8
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @family NMR
#' @importFrom ggplot2 ggplot aes_string geom_line geom_v_line scales_x_reverse scale_colour_gradientn labs theme_bw theme element_text
#' @importFrom colorRamps matlab.like2
#' @importFrom scales sapply breaks_pretty
#' @importFrom methods hasArg
#' @export
#'@examples
#' data(covid)
#' stcy_glucose=stocsy(X, ppm, driver=5.233)
#' plotStocsy(stcy_glucose, shift=c(5.15, 5.30), title='Alpha-anomeric proton of glucose (doublet)')
#' points(speaks[[1]]$ppm, speaks[[1]]$Int, col=factor(speaks[[1]]$type))
#' @section

stocsy<-function(X, ppm, driver, plotting=TRUE, title=NULL){

  if(!hasArg(ppm)){
    ppm<-as.numeric(colnames(X))
    if(!any(!is.na(ppm))){
      stop('Provide ppm argument')
    }
  }

  dr_le=length(driver)
  if(dr_le>1){
    if(dr_le!=nrow(X)){stop('Ext driver does not match to X.')}
    idx=which(!is.na(driver))
    driver=driver[idx]
    X=X[idx,]
    extD=TRUE
    }else{extD=FALSE}

  if(!'matrix' %in% is(X) & is.numeric(X)){stop('X argument is not a numeric matrix')}
  if(!is.numeric(driver) || (!extD & length(driver)>1 && !( driver[1] <= max(ppm) & driver[1] >= min(ppm)) )) {stop('STOCSY driver not numeric or outside ppm range')}
  if(!is.logical(plotting)){plotting=TRUE}


  if (extD){
    l=driver
  }else{
    idx=which.min(abs(ppm-driver))
    l=X[,idx]
  }

  cc=apply(X, 2, function(x, ll=l) {cor(x,ll)})
  cv=apply(X, 2, function(x, ll=l) {cov(x,ll)})

  stoc_mod <- new("stocsy1d_metabom8",
                  version = "0.9",
                  X = X,
                  ppm = ppm,
                  driver = driver,
                  r = cc,
                  cov=cv)


  if(plotting){
    df=data.frame(cc, cv, ppm)
    if(extD){
      csc_lab="r(ext, X)"
      extD_stats=round(summary(stoc_mod@driver), 2)
      cap=paste0('Sample size: n=', length(stoc_mod@driver), '\nExternal driver: Median=', extD_stats[3], ' range=', extD_stats[1], '-', extD_stats[6])
    }else{
        csc_lab=paste("r(d=", stoc_mod@driver, ', X)', sep='')
        cap=paste0('Sample size: n=', nrow(X))
        }

    g1=ggplot(df, aes_string(x='ppm', y='cv', colour='abs(cc)'))+
      geom_line()+
      scale_x_reverse(breaks=breaks_pretty())+
      scale_colour_gradientn(colours = matlab.like2(10), limits=c(0,1), name=csc_lab)+
      labs(x=expression(delta~{}^1*H~(ppm)), y=gsub('^r', 'cov', csc_lab), title=title, caption=cap)+
      theme_bw()+
      theme(axis.text=element_text(colour="black"), panel.grid.minor.x = element_blank())
    if(!extD){g1=g1+geom_vline(xintercept = driver, linetype=2, col='grey')}
    plot(g1)
  }
  return(stoc_mod)
}



#' @title Plot Statistical Total Correlation Spectroscopy (STOCSY)
#' @param stoc_mod STOCSY object, as created with function 'stocsy' ('see ?stocsy()')
#' @param shift num array, chemical shift range defining plotting area
#' @param title char, plot title
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @family NMR functions
#' @importFrom ggplot2 ggplot aes_string geom_line geom_v_line scales_x_reverse scale_colour_gradientn labs theme_bw theme element_text
#' @importFrom colorRamps matlab.like2
#' @export
#'@examples
#' data(covid)
#' stcy_glucose=stocsy(X, ppm, driver=5.233)
#' plotStocsy(stcy_glucose, shift=c(5.15, 5.30), title='Alpha-anomeric proton of glucose (doublet)')
#' points(speaks[[1]]$ppm, speaks[[1]]$Int, col=factor(speaks[[1]]$type))
#' @section

plotStocsy=function(stoc_mod, shift=c(0,10), title=NULL){
  if(!'stocsy1d_metabom8' %in% is(stoc_mod)) stop('Provide a STOCSY object')
  ds=data.frame(r=stoc_mod@r, cov=stoc_mod@cov, ppm=stoc_mod@ppm, stringsAsFactors = FALSE)
  if(!all(is.numeric(shift))) stop('Check shift argument')

  dr_le=length(stoc_mod@driver)
  if(dr_le>1){extD=TRUE}else{extD=FALSE}

  if(min(shift)<min(ds$ppm)){shift[which.min(shift)]=min(ds$ppm)}
  if(max(shift)>max(ds$ppm)){shift[which.max(shift)]=max(ds$ppm)}

  shift=sort(shift)
  idx=get.idx(shift, ds$ppm)
  if(length(idx)==0) stop('Check chemical shift argument')
  ds=ds[idx,]


  if(extD){
    csc_lab="r(ext, X)"
    extD_stats=round(summary(stoc_mod@driver), 2)
    cap=paste0('Sample size: n=', length(stoc_mod@driver), '\nExternal driver: Median=', extD_stats[3], ' Range=', extD_stats[1], '-', extD_stats[6])
  }else{
      csc_lab=paste("r(d=", stoc_mod@driver, ', X)', sep='')
      cap=paste0('Sample size: n=', nrow(stoc_mod@X))
      }


  g1=ggplot(ds, aes_string(x='ppm', y='cov', colour='abs(r)'))+
    geom_line()+
    scale_x_reverse(breaks=breaks_pretty())+
    scale_colour_gradientn(colours = matlab.like2(10), limits=c(0,1), name=csc_lab)+
    labs(x=expression(delta~{}^1*H~(ppm)), y=gsub('^r', 'cov', csc_lab), title=title, caption=cap)+
    theme_bw()+
    theme(axis.text=element_text(colour="black"))+
    ggtitle(title)

  if(!extD && (shift[1]<stoc_mod@driver & shift[2]>stoc_mod@driver)){
    g1<-g1+geom_vline(xintercept =stoc_mod@driver, linetype=2, col='grey')
  }

  return(g1)

}


























































































































