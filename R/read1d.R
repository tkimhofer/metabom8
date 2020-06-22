#' @title Read-in 1D NMR spectra
#' @export
#' @param path char, File directory containing spectra
#' @param procs_exp num, Topspin processing experiment ID
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @author Torben Kimhofer \email{torben.kimhofer@murdoch.edu.au}
#' @importFrom base list.files readBin seq gsub
#' @importFrom stats approxfun
#' @section

# read Bruker 1d new
read1d <- function(path,  procs_exp=1, n_max=1000, filter=T){

  # check file system intact
  f_list=checkFiles1d(path, procs_exp, n_max, filter)

  # extract parameters from acqus and procs (for f1 and f2)
  pars <- t(extract_pars1d ( f_list ))

  # convert to numeric where possible
  dtype_num<-apply(pars, 2, function(x){
    !any(is.na(suppressWarnings(as.numeric(x))))
  })

  pars<-data.frame(pars)
  pars[,dtype_num]<-apply(pars[,dtype_num], 2, as.numeric)
  rownames(pars)<-f_list[[2]]


  # chem shift
  ppm_ref <- chem_shift(swidth=pars$af2_SW[1], offset=pars$pf2_OFFSET[1], si=pars$pf2_SI[1])

  # read in binary file and
  out <- sapply(1:length(f_list[[1]]), function(s, ppref=ppm_ref){

    # chem shift
    csF2_ppm <- chem_shift(swidth=pars$af2_SW[s], offset=pars$pf2_OFFSET[s], si=pars$pf2_SI[s])
    #csF1_hz <- chem_shift(swidth=pars$af1_SW_h[s], offset=pars$pf1_OFFSET[s]*pars$pf1_SF[s], si=pars$pf1_SI[s])

    byteorda=c('little'=0, 'big'=1)
    #names(byteorda)[match(pars$af2_BYTORDA[s], byteorda)]


    spec <- readBin(paste0(f_list[[1]][s], .Platform$file.sep, 'pdata', .Platform$file.sep, procs_exp, .Platform$file.sep, '1r'),
                    what = "int",
                    n = pars$pf2_FTSIZE[s],
                    size = 4,
                    signed = T,
                    endian = names(byteorda)[match(pars$af2_BYTORDA[s], byteorda)]
    ) # this is spectra
    spec <- ( spec * (2^pars$af2_NC[s]) )
    nspec <- length(spec)

    f_spec=approxfun(x=csF2_ppm, y=spec)
    spec_inter=f_spec(ppref)

    return(spec_inter)

  })
  out=t(out)
  colnames(out)=ppm_ref

  rnames=do.call(rbind, strsplit(f_list[[1]], .Platform$file.sep))
  rnames.idx=which(apply(rnames, 2, function(x){
    length(unique(x))>1
  }))

  if(length(rnames.idx) > 1 ){
    rnam=apply(rnames[, rnames.idx], 1, function(x){
      paste(x, collapse=.Platform$file.sep)
    })
  }else{
    rnam=rnames[, rnames.idx]
  }

  rownames(out)=rnam
  rownames(meta)=rnam

  assign("X", out, envir = .GlobalEnv)
  assign("ppm", ppm_ref, envir = .GlobalEnv)
  assign("meta", pars, envir = .GlobalEnv)


}




#' @title Read Bruker NMR paramter files - helper function read1d
#' @param f_list list, intact files system for NMR experiments. See fct checkFiles1d
#' @author Torben Kimhofer \email{torben.kimhofer@murdoch.edu.au}
#' @section
extract_pars1d <- function( f_list ) {


  out=sapply(f_list[[1]], function(fil){

    f_procs=paste0(fil, .Platform$file.sep, 'pdata/1/procs')
    # extract procs information for t2
    fhand <- file(f_procs, open = "r")
    f_procs <- readLines(fhand, n = -1, warn = FALSE)
    close(fhand)

    out=strsplit(gsub('^##\\$', '',  grep('^##\\$', f_procs, value=T, fixed = F), fixed = F), '=')
    d_procs_val=gsub('^ ', '', sapply(out, '[[', 2))
    names(d_procs_val) = paste0('pf2_', sapply(out, '[[', 1))

    # acqus
    f_acqu=paste0(fil, .Platform$file.sep, 'acqus')
    # extract procs information for t2
    fhand <- file(f_acqu, open = "r")
    f_acqu <- readLines(fhand, n = -1, warn = FALSE)
    close(fhand)

    out=strsplit(gsub('^##\\$', '',  grep('^##\\$', f_acqu, value=T, fixed = F), fixed = F), '=')
    d_acqu_val=gsub('^ ', '', sapply(out, '[[', 2))
    names(d_acqu_val) = paste0('af2_', sapply(out, '[[', 1))

    # change date
    idx=grep('date', names(d_acqu_val), ignore.case = T)
    d_acqu_val[idx]=as.character(as.POSIXct(x = '01/01/1970 00:00:00', format='%d/%m/%Y %H:%M:%S')+(as.numeric(d_acqu_val[idx])))

    pars = c(d_acqu_val, d_procs_val)

    return(pars)

  })


  return(out)

}

#' @title Calc chemical shift axis points - helper function read1d
#' @param swidth num, sweep width extracted from procs file. See fct extract_pars1d
#' @author Torben Kimhofer \email{torben.kimhofer@murdoch.edu.au}
#' @section
chem_shift <- function(swidth, offset, si){

  dppm <- swidth/(si - 1) # ho
  cshift <- seq(offset, (offset - swidth), by = -dppm)

  return(cshift)

}

#' @title Check for intact file systems - helper function read1d
#' @param datapath char, File directory containing spectra
#' @param procs_exp num, Topspin processing experiment ID
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @author Torben Kimhofer \email{torben.kimhofer@murdoch.edu.au}
#' @section
checkFiles1d <- function(datapath, procs_exp=1, n_max=10, filter=T) {

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

  idx_a <- id_a %in% id_p & id_a %in% id_f1
  idx_p <- id_p %in% id_a & id_p %in% id_f1
  idx_f1 <- id_f1 %in% id_p & id_f1 %in% id_a

  if(!any(idx_a | idx_p | idx_f1)){
    if (filter == T) {
      message('File system seesm to be corrupt for some experiments - filtering for intact file systems.')
      f_acqus <- f_acqus[idx_a]
      f_procs <- f_procs[idx_p]
      f_1r <- f_1r[idx_r]
    }else{
      message('File system seesm to be corrupt for some experiments. Consider function argument \`filter=TRUE\`')
      return(NULL)
    }
  }

  if(length(f_acqus) > n_max) { f_procs=f_procs[1:n_max]; message('Reached n_max - not all spectra read-in.') }
  if(length(f_procs) == 0 ) { message('No spectrum found'); return(NULL) }

  if(n_max < length(f_procs) ) {
    f_procs=f_procs[1:n_max];
    f_acqus=f_acqus[1:n_max];
    f_f1=f_f1[1:n_max];
    message('Reached n_max - not all spectra read-in.') }
  if(length(f_procs) == 0) { message('No spectrum found'); return(NULL) }

  p_intact=gsub('/acqus$', '', f_acqus)
  exp_no=id_p

  return(list(path=p_intact, exp_no=exp_no))
}