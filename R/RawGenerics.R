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

# lb to the negative of the peak width
#' @keywords internal
#' @section
.gauss<-function(n, lb=-1.2, gb=0.3,  para){
  sw_hz=para$a_SW_h
  AQ=para$a_TD/(2*sw_hz)
  idx=seq(n)-1
  e=(-pi*idx*lb) / sw_hz
  g=((lb*pi) / (2*gb*AQ)) * (idx / sw_hz)
  out<-exp(e-g^2)
  minmax(out)
}

#test=.gauss(n=length(fid), lb=-0.65, gb=0.3,sw_hz = 12019.23, para)
#plot(test, type='l')

# pars is list
#' @section
.fidApodisationFct<-function(n, pars){
  if(is.null(pars$fun) || !pars$fun %in% c('uniform', 'exponential', 'cosine', 'sine', 'sem', 'expGaus_resyG', 'gauss')){stop('Check apodisation function argument (fun).')}
  switch(pars$fun,
         'uniform'={afun<-rep(1, n)},
         'exponential'={ if( 'lb' %in% names(pars)) {afun<-.em(n, pars$lb)} else{stop('Check aposation fct arguments: Exponention function requires lb parameter')}},
         'cosine'={afun<-.cosine(n)},
         'sine'={afun<-.sine(n)},
         'sem'={if( 'lb' %in% names(pars)) {afun<-.sem(n, pars$lb)} else{ stop('Check aposation fct arguments: SEM function requires lb parameter')}},
         'gauss'={afun<-.gauss(n, pars$lb, pars$gb, para)},
         )
  if(!is.null(pars$plot) && pars$plot && exists('afun')){ plot(afun, type='l', main=paste('Apodisation function:', pars$fun))}
  return(afun)
  }


# filter experiments
#' @keywords internal
#' @section
#' @importFrom plyr ddply .
.filterExp_files<-function(pars, exp_type, f_list, n_max){


  idx<-match(toupper(names(exp_type)), toupper(gsub('[ap]_', '', colnames(pars))))
  if(length(idx)==0){stop('No parameter(s) found that that match the specification. Check for exp_type argument for typos and parameter choices in acqus and procs files.')}

  #browser()
  idx_na<-which(is.na(idx))
  if(length(idx_na)>0){
    if(length(idx_na)==length(idx)){
      stop('No matching paramter names found. Check input argument exp_type.')
    }else{
      message(paste('Experiment filter', names(exp_type)[idx_na]  , 'not in NMR acquisition list. Using remaining arguments to filter :', names(exp_type)[-idx_na]))
      idx=idx[-idx_na]
    }
  }
  fmat<-sapply(seq(length(idx)), function(i){
    vars=gsub('^<|>$', '', pars[,idx[i]])
    vars %in% exp_type[[i]]
  })
  idx_filt=apply(fmat==1, 1, all)
  if(!any(idx_filt)){stop('No files found that match the specified parameter specification levels.')}
  f_list=lapply(f_list, function(x){
    x[idx_filt]
    })



  idx_filt=which(idx_filt)
  if(length(idx_filt)>n_max){

    idx_filt=idx_filt[seq_len(n_max)]
    f_list=lapply(f_list, function(x,  idx=seq_len(n_max)){
      x[idx]
    })

  }

  pars<-pars[idx_filt,]


  # order and add rownames
  fnam <- strsplit(f_list$f_1r, .Platform$file.sep)
  idx_keep <- which((apply(do.call(rbind, fnam), 2, function(x) length(unique(x)))) > 1)

  # find order of racks
  if(length(idx_keep)>1){

    rack_<-t(vapply(seq_len(length(fnam)), function(i, iid=idx_keep[length(idx_keep)-1]){
      c(fnam[[i]][iid], pars$a_DATE[i])
    }, FUN.VALUE = c('', '')))
    colnames(rack_)<-c('a', 'b')
    rack_order_<-ddply(as.data.frame(rack_), .(a), function(x){
      mean(as.POSIXct(x$b))
    })
    rord_fac<-order(rack_order_$V1)*1e5

    fnam1 <- vapply(fnam, function(x, st = idx_keep[length(idx_keep)-1]) {
      x[st]
    }, FUN.VALUE = "")

    rord_fac=rord_fac[match(fnam1, rack_order_$a)]


    exp_ <- vapply(fnam, function(x, st = idx_keep[length(idx_keep)]) {
      x[st]
    }, FUN.VALUE = "")
    if (any(is.na(suppressWarnings(as.numeric(exp_))))) { # Is there at least one non-numeric ID?
      exp_ <- factor(exp_)
    }
    rr<-order(as.numeric(exp_)+rord_fac)

  }else{
    exp_ <- vapply(fnam, function(x, st = idx_keep) {
      x[st]
    }, FUN.VALUE = "")
    if (any(is.na(suppressWarnings(as.numeric(exp_))))) { # Is there at least one non-numeric ID?
      exp_ <- factor(exp_)
    }
    rr<-order(as.numeric(exp_))
  }

  #browser()


  # re-order f_list and pars according to rr
  pars=pars[rr,]
  f_list=lapply(f_list, function(x){
    x[rr]
  })

  pars$a_DATE=as.POSIXct(pars$a_DATE)

  return(list(f_list, pars))
}



#' @title Read Bruker NMR paramter files - helper function read1d
#' @param f_list list, intact files system for NMR experiments. See fct checkFiles1d
#' @param procs_exp num or char, which processing experiment should be extracted
#' @return data frame of spectrometer acquisition parameters
#' @author \email{torben.kimhofer@@murdoch.edu.au}
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










































