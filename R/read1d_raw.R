#### read fids


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

# path='/Volumes/ANPC_ext1/COVID_plasma_DIRE_150720/'
# datapath=path
# procs_exp = '1'
# n_max = 20000
# filter = T

# read Bruker 1d new
read1d_raw <- function(path,  n_max=1000, filter=T){

  if(grepl('^~', path, fixed=F)){
    getwd()
    path=gsub('^~',path.expand("~"), path)
  }


  # check file system intact
  f_list=.check1d_files_fid(path, n_max, filter)

  # extract parameters from acqus and procs (for f1 and f2)

  pars <- .extract_acq_pars1d(f_list)

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
  #ppm_ref <- .chemShift(swidth=pars$a_SW[1], offset=pars$p_OFFSET[1], si=pars$a_SI[1])

  # read in binary file and
  out <- sapply(1:length(f_list[[1]]), function(s){

    # chem shift
    # csF2_ppm <- .chemShift(swidth=pars$a_SW[s], offset=pars$p_OFFSET[s], si=pars$p_SI[s])
    #csF1_hz <- chem_shift(swidth=pars$af1_SW_h[s], offset=pars$pf1_OFFSET[s]*pars$pf1_SF[s], si=pars$pf1_SI[s])

    byteorda=c('little'=0, 'big'=1)
    #names(byteorda)[match(pars$af2_BYTORDA[s], byteorda)]

    dtypa=c('int'=0, 'double'=2)

    spec <- readBin(paste0(f_list[[1]][s], .Platform$file.sep, 'fid'),
                    what = names(dtypa)[match(pars$a_DTYPA[s], dtypa)],
                    n = pars$a_TD[s],
                    #size = 4,
                    signed = T,
                    endian = names(byteorda)[match(pars$a_BYTORDA[s], byteorda)]
    ) # this is spectra
    spec <- ( spec * (2^pars$a_NC[s]) )
    nspec <- length(spec)

    # apply processing (line broadening, apodisation)

    # zero-filling
    # add =  2^32 - length(spec)
    spec_zf=c(spec, rep(0, length(spec)))

    # lb
    spec_lb=.em(spec_zf, lb=1.3, pars$a_SW_h[s])

    # plot(spec_lb[1000:30000], type='l')
    # points(spec[1000:30000], type='l', col='cyan')
    #
    # plot(spec[1000:60000], type='l', col='cyan')



    # transformation
    sp=fft(spec_lb)

    # absorption mode
    sp_ab=Re(sp)
    sp_disp=Im(sp) # dispersive mode ( used for phase correction)
    sp_mag=sp_ab+sp_disp
    sp_pow=sp_ab^2+sp_disp^2

    plot(sp_pow, type='l')
    plot(sp_ab, type='l')

    f_spec=approxfun(x=csF2_ppm, y=spec)
    spec_inter=f_spec(ppref)

    return(spec_inter)

  })

  out=t(out)
  colnames(out)=ppm_ref

  fnam=strsplit(f_list[[1]], .Platform$file.sep)
  idx_rm=min(sapply(fnam, length))
  fnam=sapply(fnam, function(x, st=idx_rm){
    paste(x[idx_rm:length(x)], collapse=.Platform$file.sep)
  })

  rownames(out)=fnam
  rownames(pars)=fnam

  assign("X", out, envir = .GlobalEnv)
  assign("ppm", ppm_ref, envir = .GlobalEnv)
  assign("meta", pars, envir = .GlobalEnv)


}



#' @title Check for intact file systems - helper function read1d
#' @param datapath char, File directory containing spectra
#' @param procs_exp num, Topspin processing experiment ID
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @section
.check1d_files_fid <- function(datapath, n_max=10, filter=T) {

  datapath=gsub(paste0(.Platform$file.sep, '$'), '', datapath)

  # searches for status acquisition parameter files
  f_acqus <- list.files(path = datapath, pattern ="^acqus$",
                        all.files = FALSE, full.names = TRUE, recursive = TRUE,
                        ignore.case = TRUE)


  # searches for 1r files
  f_fid <- list.files(path = datapath, pattern ="^fid$",
                     all.files = FALSE, full.names = TRUE, recursive = TRUE,
                     ignore.case = TRUE)


  # check file systems
  # get experiement folder id
  id_a <- gsub(paste("^", datapath, "/|/acqus$", sep = ""), "", f_acqus)
  id_fid <- gsub(paste("^", datapath, "/|/fid$", sep = ""), "", f_fid)
  # ensure that all three files are present

  idx_a <- which(id_a %in% id_fid)
  idx_fid <- which(id_fid %in% id_a)

  if(length(idx_a)!=length(id_a) | length(idx_fid)!=length(id_fid)){
    if (filter == T) {
      message('Filtering experiment processing folders.')
      f_acqus <- f_acqus[idx_a]; id_a=id_a[idx_a]
      f_fid <- f_fid[idx_fid]; id_fid=id_fid[idx_fid]
    }else{
      message('File system seesm to be corrupt for some experiments. Consider function argument \`filter=TRUE\`')
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




#' @title Read Bruker NMR paramter files - helper function read1d
#' @param f_list list, intact files system for NMR experiments. See fct checkFiles1d
#' @param procs_exp num or char, which processing experiment should be extracted
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
# @importFrom base sapply
#' @section
.extract_acq_pars1d <- function( f_list ) {


  out=lapply(f_list[[1]], function(fil){
    # f_procs=paste0(fil, .Platform$file.sep, 'pdata', .Platform$file.sep, procs_exp, .Platform$file.sep, 'procs')
    # # extract procs information for t2
    # fhand <- file(f_procs, open = "r")
    # f_procs <- readLines(fhand, n = -1, warn = FALSE)
    # close(fhand)
    #
    # out=strsplit(gsub('^##\\$', '',  grep('^##\\$', f_procs, value=T, fixed = F), fixed = F), '=')
    # d_procs_val=gsub('^ ', '', sapply(out, '[[', 2))
    # names(d_procs_val) = paste0('p_', sapply(out, '[[', 1))

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

  out_le = sapply(out, length)

  # if( length(unique(out_le)) > 1 ){
  #
  #   cnam=unique(unlist(sapply(out, names)))
  #   out_df=matrix(NA, nrow = 1, ncol=length(cnam))
  #
  #   out=as.data.frame(t(sapply(out, function(x, odf=out_df, cc=cnam){
  #     odf[1, match(names(x), cnam)]=x
  #     return(odf)
  #   })))
  #   colnames(out)=cnam
  # }
  #

  return(out)

}





















































































