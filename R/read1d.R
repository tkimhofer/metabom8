#' @title Read-in 1D NMR spectra
#' @export
#' @param path char, File directory containing spectra
#' @param procs_exp num, Topspin processing experiment ID
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
# @importFrom base list.files readBin seq gsub
#' @importFrom stats approxfun
#' @section

# path='/Users/TKimhofer/Desktop/COVID_plasma_DIRE_IVDR02_PN_170620/'
# procs_exp = '1'
# n_max = 2
# filter = T

# read Bruker 1d new
read1d <- function(path,  procs_exp=1, n_max=1000, filter=T){

  if(grepl('^~', path, fixed=F)){
    getwd()
    path=gsub('^~',path.expand("~"), path)
  }


  # check file system intact
  f_list=.checkFiles1d(path, procs_exp, n_max, filter)

  # extract parameters from acqus and procs (for f1 and f2)

  pars <- t(.extract_pars1d ( f_list, procs_exp))

  if(is.list(pars)){
    le=unique(sapply(pars, length))
    if(length(le)==1){
      pars=do.call(rbind, pars)
    }else{
      unam<-names(unlist(pars))
      pars<-do.call(rbind, lapply(pars, '[', unam))
      colnames(pars)<-unam
    }

  }

  if(nrow(pars)!=length(f_list[[1]])){
    pars <- t(pars)
  }


  # convert to numeric where possible
  dtype_num<-apply(pars, 2, function(x){
    !any(is.na(suppressWarnings(as.numeric(x))))
  })

  pars<-as.data.frame(pars)
  pars[,dtype_num]<-apply(pars[,dtype_num], 2, as.numeric)
  rownames(pars)<-f_list[[2]]

  # chem shift
  ppm_ref <- .chemShift(swidth=pars$a_SW[1], offset=pars$p_OFFSET[1], si=pars$p_SI[1])

  # read in binary file and
  out <- sapply(1:length(f_list[[1]]), function(s, ppref=ppm_ref){

    # chem shift
    csF2_ppm <- .chemShift(swidth=pars$a_SW[s], offset=pars$p_OFFSET[s], si=pars$p_SI[s])
    #csF1_hz <- chem_shift(swidth=pars$af1_SW_h[s], offset=pars$pf1_OFFSET[s]*pars$pf1_SF[s], si=pars$pf1_SI[s])

    byteorda=c('little'=0, 'big'=1)
    #names(byteorda)[match(pars$af2_BYTORDA[s], byteorda)]


    spec <- readBin(paste0(f_list[[1]][s], .Platform$file.sep, 'pdata', .Platform$file.sep, procs_exp, .Platform$file.sep, '1r'),
                    what = "int",
                    n = pars$p_FTSIZE[s],
                    size = 4,
                    signed = T,
                    endian = names(byteorda)[match(pars$a_BYTORDA[s], byteorda)]
    ) # this is spectra
    spec <- ( spec * (2^pars$a_NC[s]) )
    nspec <- length(spec)

    f_spec=approxfun(x=csF2_ppm, y=spec)
    spec_inter=f_spec(ppref)

    return(spec_inter)

  })

  out=t(out)
  colnames(out)=ppm_ref

  fnam=strsplit(f_list[[1]], .Platform$file.sep)
  idx_rm=min(sapply(fnam, length))-1
  fnam=sapply(fnam, function(x, st=idx_rm){
    paste(x[idx_rm:length(x)], collapse=.Platform$file.sep)
  })

  rownames(out)=fnam
  rownames(pars)=fnam

  assign("X", out, envir = .GlobalEnv)
  assign("ppm", ppm_ref, envir = .GlobalEnv)
  assign("meta", pars, envir = .GlobalEnv)


}




#' @title Read Bruker NMR paramter files - helper function read1d
#' @param f_list list, intact files system for NMR experiments. See fct checkFiles1d
#' @param procs_exp num or char, which processing experiment should be extracted
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
# @importFrom base sapply
#' @section
.extract_pars1d <- function( f_list, procs_exp ) {


  out=lapply(f_list[[1]], function(fil){
    f_procs=paste0(fil, .Platform$file.sep, 'pdata', .Platform$file.sep, procs_exp, .Platform$file.sep, 'procs')
    # extract procs information for t2
    fhand <- file(f_procs, open = "r")
    f_procs <- readLines(fhand, n = -1, warn = FALSE)
    close(fhand)

    out=strsplit(gsub('^##\\$', '',  grep('^##\\$', f_procs, value=T, fixed = F), fixed = F), '=')
    d_procs_val=gsub('^ ', '', sapply(out, '[[', 2))
    names(d_procs_val) = paste0('p_', sapply(out, '[[', 1))

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

    pars = c(d_acqu_val, d_procs_val)

    return(pars)

  })

  out_le = sapply(out, length)

  if( length(unique(out_le)) > 1 ){

    cnam=unique(unlist(sapply(out, names)))
    out_df=matrix(NA, nrow = 1, ncol=length(cnam))

    out=as.data.frame(t(sapply(out, function(x, odf=out_df, cc=cnam){
      odf[1, match(names(x), cnam)]=x
      return(odf)
    })))
    colnames(out)=cnam
  }


  return(out)

}

# chem shift defined in read2d.R
# #' @title Calc chemical shift axis points - helper function read1d
# #' @param swidth num vector, sweep width extracted from procs file. See fct extract_pars1d
# #' @param offset num vector, offset value extracted from acqus file. See fct extract_pars1d
# #' @param si num, number of data points. See fct extract_pars1d
# #' @author Torben Kimhofer \email{torben.kimhofer@murdoch.edu.au}
# #' @section
# .chem_shift <- function(swidth, offset, si){
#
#   dppm <- swidth/(si - 1) # ho
#   cshift <- seq(offset, (offset - swidth), by = -dppm)
#
#   return(cshift)
#
# }
#

#' @title Check for intact file systems - helper function read1d
#' @param datapath char, File directory containing spectra
#' @param procs_exp num, Topspin processing experiment ID
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @section
.checkFiles1d <- function(datapath, procs_exp=1, n_max=10, filter=T) {

  datapath=gsub(paste0(.Platform$file.sep, '$'), '', datapath)

  # searches status processing parameter files
  f_procs <- list.files(path = datapath, pattern ="^procs$",
                        all.files = FALSE, full.names = TRUE, recursive = TRUE,
                        ignore.case = TRUE)
  f_procs=grep(paste0('pdata', .Platform$file.sep, procs_exp), f_procs, value=T)

  # searches for status acquisition parameter files
  f_acqus <- list.files(path = datapath, pattern ="^acqus$",
                        all.files = FALSE, full.names = TRUE, recursive = TRUE,
                        ignore.case = TRUE)


  # searches for 1r files
  f_1r <- list.files(path = datapath, pattern ="^1r$",
                     all.files = FALSE, full.names = TRUE, recursive = TRUE,
                     ignore.case = TRUE)


  # check file systems
  # get experiement folder id
  id_p <- gsub(paste("^", datapath, "/|/pdata/1/procs$", sep = ""), "", f_procs)
  id_a <- gsub(paste("^", datapath, "/|/acqus$", sep = ""), "", f_acqus)
  id_f1 <- gsub(paste("^", datapath, "/|/pdata/1/1r$", sep = ""), "", f_1r)
  # ensure that all three files are present

  idx_a <- which(id_a %in% id_p & id_a %in% id_f1)
  idx_p <- which(id_p %in% id_a & id_p %in% id_f1)
  idx_f1 <- which(id_f1 %in% id_p & id_f1 %in% id_a)

  if(length(idx_a)!=length(id_a) | length(idx_p)!=length(id_p)| length(idx_f1)!=length(id_f1)){
    if (filter == T) {
      message('Filtering experiment processing folders.')
      f_acqus <- f_acqus[idx_a]; id_a=id_a[idx_a]
      f_procs <- f_procs[idx_p]; id_p=id_p[idx_p]
      f_1r <- f_1r[idx_f1]; id_f1=id_f1[idx_f1]
    }else{
      message('File system seesm to be corrupt for some experiments. Consider function argument \`filter=TRUE\`')
      return(NULL)
    }
  }


  # test if these can all be matched
  idm_a=match(id_p, id_a)
  idm_f1=match(id_p, id_f1)

  if(any( is.na(idm_a) | is.na(idm_f1)) | any( diff(idm_a)>1 | diff(idm_f1)>1 )){
    stop('check matching of this functions')
  }



  if(length( unique(c(length(f_acqus), length(f_procs), length(f_1r)))) != 1) {stop('Somethings wrong after filtering!')}
  # if(length(f_acqus) > n_max) { f_procs=f_procs[1:n_max]; message('Reached n_max - not all spectra read-in.') }
  if(length(f_procs) == 0 ) { message('No spectrum found'); return(NULL) }

  if(n_max < length(f_procs) ) {
    f_procs=f_procs[1:n_max];
    f_acqus=f_acqus[1:n_max];
    f_1r=f_1r[1:n_max];

    message('Reached n_max - not all spectra read-in.') }
  if(length(f_procs) == 0) { message('No spectrum found'); return(NULL) }

  p_intact=gsub('/acqus$', '', f_acqus)
  exp_no=id_p

  return(list(path=p_intact, exp_no=exp_no))
}


.chemShift <- function(swidth, offset, si){

  dppm <- swidth/(si - 1) # ho
  cshift <- seq(offset, (offset - swidth), by = -dppm)

  return(cshift)

}

