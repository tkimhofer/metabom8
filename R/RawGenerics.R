#  apodisation functions
#' @keywords internal
#' @section
.em<-function(n, lb){
  idx<-seq(n)
  out<-exp(-( idx * lb * pi) / (length(idx)))
  minmax(out)
}

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
#' @keywords internal
#' @section
.expGaus_resyG<-function(n, lb=-10, gb=0.12, aq_t){
  idx=seq(n)
  b=-lb/(2*gb*aq_t)
  out<-exp((-lb*idx)-(b*idx^2))
  minmax(out)
}

#' @section
.fidApodisationFct<-function(n, pars){

  if(is.null(pars$fun) || !pars$fun %in% c('uniform', 'exponential', 'cosine', 'sine', 'sem', 'expGaus_resyG')){stop('Check apodisation function argument (fun).')}

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






# filter experiments
.filterExp_files<-function(pars, exp_type, f_list){
  idx<-match(toupper(names(exp_type)), toupper(gsub('[ap]_', '', colnames(pars))))
  if(length(idx)==0){stop('No files found that match the specified parameter specification. Check input argument exp_type.') }

  idx_na<-which(is.na(idx))
  if(length(idx_na)>0){
    message(paste('Experiment filter', names(exp_type)[idx_na]  , 'not in NMR acquisition list. Using remaining arguments to filter :', names(exp_type)[-idx_na]))
    idx=idx[-idx_na]
  }
  fmat<-vapply(seq(length(idx)), function(i){
    vars=gsub('^<|>$', '', pars[,idx[i]])
    vars %in% exp_type[[i]]
  }, FUN.VALUE = pars[,1])
  idx_filt=apply(fmat==1, 1, all)
  if(!any(idx_filt)){stop('No files found that match the specified parameter specification levels.')}
  f_list[[1]]=f_list[[1]][idx_filt]
  f_list[[2]]=f_list[[2]][idx_filt]
  pars=pars[idx_filt,]

  return(list(f_list, pars))
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

    # TODO: add visualisation of pulse sequence
    idx<-grep('..', f_acqu, fixed=TRUE)
    f_acqu[idx]<-vapply(idx, function(i){ gsub(' .*', f_acqu[i+1], f_acqu[i])}, FUN.VALUE = '')

    out=strsplit(gsub('^##\\$', '',  grep('^##\\$', f_acqu, value=TRUE, fixed = FALSE), fixed = FALSE), '=')
    d_acqu_val=gsub('^ ', '', vapply(out, '[[', 2 , FUN.VALUE = ''))
    names(d_acqu_val) = paste0('a_', vapply(out, '[[', 1, FUN.VALUE = ''))
    # change date
    idx=grep('date', names(d_acqu_val), ignore.case = TRUE)
    d_acqu_val[idx]=as.character(as.POSIXct(x = '01/01/1970 00:00:00', format='%d/%m/%Y %H:%M:%S')+(as.numeric(d_acqu_val[idx])))
    return(d_acqu_val)
  })

  #out_le = sapply(out, length)

  if(is.list(out)){
    le<-unique(vapply(out, length, FUN.VALUE = 1))
    if(length(le)==1){
      out<-do.call(rbind, out)
    }else{
      unam<-unique(names(unlist(out)))
      out<-do.call(rbind, lapply(out, '[', unam))
      colnames(out)<-unam
    }
  }

  if(nrow(out)!=length(f_list[[1]])){
    out<-t(out)
  }

  # convert to numeric where possible
  dtype_num<-apply(out, 2, function(x){
    !any(is.na(suppressWarnings(as.numeric(x))))
  })
  out<-as.data.frame(out)
  out[,dtype_num]<-apply(out[,dtype_num], 2, as.numeric)
  rownames(out)<-f_list[[2]]

  return(out)
}

#' @title Check for intact file systems - helper function read1d
#' @param datapath char, File directory containing spectra
#' @param procs_exp num, Topspin processing experiment ID
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @return List of strings, each one is a path to fid and experiment number
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @keywords internal
#' @section
.check1d_files_fid <- function(datapath, n_max=10, filter=TRUE, recursive, verbose) {

  datapath=gsub(paste0(.Platform$file.sep, '$'), '', datapath)
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
    if (filter) {
      if(verbose>1){ message('Incomplete files detected - reading experiments with intact folder structure.')}
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
  if(length( unique(c(length(f_acqus),length(f_fid)))) != 1) {stop('Something\'s wrong after filtering!')}
  if(n_max < length(f_fid) ) {
    idx_nmax=seq_len(n_max)
    f_acqus=f_acqus[idx_nmax]
    f_fid=f_fid[idx_nmax]
    id_fid=id_fid[idx_nmax]
    message('Reached n_max - not all spectra read-in.')
  }
  p_intact=gsub('/acqus$', '', f_acqus)
  exp_no=id_fid
  return(list(path=p_intact, exp_no=exp_no))
}










































