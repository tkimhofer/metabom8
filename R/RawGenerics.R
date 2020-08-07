#  apodisation functions
#' @keywords internal
#' @section
.em<-function(n, lb){
  idx<-seq(n)
  out<-exp(-( idx * lb * pi) / (length(idx)))
  minmax(out)
}

# .gmod<-function(n, a, b){
#   idx<-seq(n)
#   out<-exp(((-a*idx) - (b*idx^2))  / n)
#   minmax(out)
# }

#' @keywords internal
#' @section
.cosine<-function(n){
  idx<-seq(n)
  out=cos(pi*idx/n)
  minmax(out)
}
#' @keywords internal
#' @section
.sine<-function(n){
  idx=seq(n)
  out<-sin(pi*idx/n)
  minmax(out)
}
#' @keywords internal
#' @section
.sem<-function(n, lb=1.5){
  idx=seq(n)
  out<-sin((pi*idx)/n) * exp(-( idx * lb * pi) / n)
  minmax(out)
}


#' @section
.fidApodisationFct<-function(n, pars){

  if(is.null(pars$fun) || !pars$fun %in% c('uniform', 'exponential', 'cosine', 'sine', 'sem')){stop('Check apodisation function argument (fun).')}

  switch(pars$fun,
         'uniform'={afun<-rep(1, n)},
         'exponential'={ if( 'lb' %in% names(pars)) {afun<-.em(n, pars$lb)} else{stop('Check aposation fct arguments: Exponention function requires lb parameter')}},
         #'modGauss'={},
         #'expDampJmod'={},
         'cosine'={afun<-.cosine(n)},
         'sine'={afun<-.sine(n)},
         #'sineShift'={if( all(c('offs', 'end', 'exp') %in% names(pars))){afun<-.sineMod(n, pars$offs, pars$end, pars$exp)} else{stop('Check aposation fct arguments: sineShift function requires offs, end, exp parameters')}},
         #'triangle'={},
         'sem'={if( 'lb' %in% names(pars)) {afun<-.sem(n, pars$lb)} else{ stop('Check aposation fct arguments: SEM function requires lb parameter')}}
         )

  if(!is.null(pars$plot) && pars$plot && exists('afun')){ plot(afun, type='l', main=paste('Apodisation function:', pars$fun))}

  return(afun)

  }




#' @title Read Bruker NMR paramter files - helper function read1d
#' @param f_list list, intact files system for NMR experiments. See fct checkFiles1d
#' @param procs_exp num or char, which processing experiment should be extracted
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
# @importFrom base sapply
#' @section
#' @keywords internal
#' @section
.extract_acq_pars1d <- function(f_list) {

  out=lapply(f_list[[1]], function(fil){
    # acqus
    f_acqu=paste0(fil, .Platform$file.sep, 'acqus')
    # extract procs information for t2
    fhand <- file(f_acqu, open = "r")
    f_acqu <- readLines(fhand, n = -1, warn = FALSE)
    close(fhand)

    out=strsplit(gsub('^##\\$', '',  grep('^##\\$', f_acqu, value=T, fixed = F), fixed = F), '=')
    d_acqu_val=gsub('^ ', '', sapply(out, '[[', 2))
    names(d_acqu_val) = paste0('a_', sapply(out, '[[', 1))
    # change date
    idx=grep('date', names(d_acqu_val), ignore.case = T)
    d_acqu_val[idx]=as.character(as.POSIXct(x = '01/01/1970 00:00:00', format='%d/%m/%Y %H:%M:%S')+(as.numeric(d_acqu_val[idx])))
    return(d_acqu_val)
  })

  #out_le = sapply(out, length)

  return(out)
}






#' @title Check for intact file systems - helper function read1d
#' @param datapath char, File directory containing spectra
#' @param procs_exp num, Topspin processing experiment ID
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @keywords internal
#' @section
.check1d_files_fid <- function(datapath, n_max=10, filter=T, recursive=TRUE) {

  datapath=gsub(paste0(.Platform$file.sep, '$'), '', datapath)
#browser()
  # searches for status acquisition parameter files
  f_acqus <- list.files(path = datapath, pattern ="^acqus$",
                        all.files = FALSE, full.names = TRUE, recursive=recursive,
                        ignore.case = TRUE)

  # searches for 1r files
  f_fid <- list.files(path = datapath, pattern ="^fid$",
                      all.files = FALSE, full.names = TRUE, recursive=recursive,
                      ignore.case = TRUE)

  # check esxperiment folder ID
  id_a <- gsub(paste("^", datapath, "/|/acqus$", sep = ""), "", f_acqus)
  id_fid <- gsub(paste("^", datapath, "/|/fid$", sep = ""), "", f_fid)
  idx_a <- which(id_a %in% id_fid)
  idx_fid <- which(id_fid %in% id_a)

  if(length(idx_a) != length(id_a) | length(idx_fid) != length(id_fid)){
    if (filter == T) {
      message('Filtering experiment processing folders.')
      f_acqus <- f_acqus[idx_a]; id_a=id_a[idx_a]
      f_fid <- f_fid[idx_fid]; id_fid=id_fid[idx_fid]
    }else{
      message('File system seems to be corrupt for some experiments. Consider function argument \`filter=TRUE\`')
      return(NULL)
    }
  }


  # test if these can all be matched
  idm_fid=match(id_a, id_fid)

  if(any( is.na(idm_fid)) | any(  diff(idm_fid)>1 )){
    stop('check matching of this functions')
  }



  if(length( unique(c(length(f_acqus),length(f_fid)))) != 1) {stop('Somethings wrong after filtering!')}


  if(n_max < length(f_fid) ) {
    f_acqus=f_acqus[1:n_max];
    f_fid=f_fid[1:n_max];
    id_fid=id_fid[1:n_max];

    message('Reached n_max - not all spectra read-in.')
  }

  p_intact=gsub('/acqus$', '', f_acqus)
  exp_no=id_fid

  return(list(path=p_intact, exp_no=exp_no))
}


