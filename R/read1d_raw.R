#### read fids


#' @title Read-in 1D NMR spectra
#' @export
#' @param path char, File directory containing spectra
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @param apodisation list, apodisation type and parameter
#' @param zerofil int, zero-filling, exponent summand of base 2
#' @param type char, return mode: absorption, dispersion, magnitude
#' @param exp char, experiments to read-in (e.g. <PROF_URINE_NOESY>)
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
# @importFrom base list.files readBin seq gsub
#' @importFrom stats approxfun
#' @section

# read1d(path)
#
# path='/Volumes/Torben_1 1/Rproj/AUTISM/NMR_Urine/Autism_Urine_Rack1_NH_030816/'
# datapath=path
#
# n_max = 20
# filter = T
#
# par(mfrow=c(3,1))
# ppm_ra=c(2,4)
# read1d_raw(path, n_max = 10, type='magnitude')
# matspec(X, ppm, shift=ppm_ra)
# read1d_raw(path, n_max = 10, type='dispersion')
# matspec(X, ppm, shift=ppm_ra)
# read1d_raw(path, n_max = 100, type='absorption', exp = '<PROF_URINE_NOESY>', zerofil = 1L, apodisation=list(fun='cosine', offs=1, end=0.0005, exp=3, plot=T))
# matspec(X, ppm, shift=ppm_ra)
#
# matspec(X, ppm, shift=c(3,3.1), main='None')
#

# read Bruker 1d new
read1d_raw <- function(path,  n_max=1000, filter=T, apodisation=list(fun='exponential', lb=0.2), zerofil=1L, type='absorption', exp='<PROF_URINE_NOESY>'){

  if(grepl('^~', path, fixed=F)){
    getwd()
    path=gsub('^~',path.expand("~"), path)
  }

  if(!type %in% c('absorption', 'dispersion', 'magnitude')){ type='absorption'; message('Check argument type. Returning absorption spectrum.') }

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

  idx_filt<-which(pars$a_EXP == exp)
  if(length(idx_filt)==0){stop(paste('Specified experiment (', exp, ') not in path. Select one of the following: ', unique(pars$a_EXP), sep=''))}

  pars=pars[idx_filt,]
  f_list[[1]]=f_list[[1]][idx_filt]
  f_list[[2]]=f_list[[2]][idx_filt]

  pars[,dtype_num]<-apply(pars[,dtype_num], 2, as.numeric)
  rownames(pars)<-f_list[[2]]

  # chem shift
  ppm_ref <- seq(-2, 13, by=( pars$a_SW_h[1])/pars$a_TD[1]/pars$a_SFO1[1])

  if(length(unique(pars$a_TD))>2 || length(unique(pars$a_GRPDLY))>1){stop('Number of points collected in time domain is unqual across experiments.')}
  apoFct<-.fidApodisationFct(n=(pars$a_TD[1]-(pars$a_GRPDLY[1]*2)), apodisation) # subtract group delay from TD (digital filtering artefact)



  # read in binary file and
  out <- sapply(1:length(f_list[[1]]), function(s, pref=ppm_ref, afun=apoFct, zf=zerofil){
    # chem shift
    # csF2_ppm <- .chemShift(swidth=pars$a_SW[s], offset=pars$p_OFFSET[s], si=pars$p_SI[s])
    #csF1_hz <- chem_shift(swidth=pars$af1_SW_h[s], offset=pars$pf1_OFFSET[s]*pars$pf1_SF[s], si=pars$pf1_SI[s])
    byteorda=c('little'=0, 'big'=1)
    #names(byteorda)[match(pars$af2_BYTORDA[s], byteorda)]
    dtypa=c('int'=0, 'double'=2)
    spec <- readBin(paste0(f_list[[1]][s], .Platform$file.sep, 'fid'),
                    what = names(dtypa)[match(pars$a_DTYPA[s], dtypa)],
                    n = pars$a_TD[s],
                    size = 4L,
                    endian = names(byteorda)[match(pars$a_BYTORDA[s], byteorda)]
    ) # this is spectra
    spec <- ( spec * (2^pars$a_NC[s]) )




    # sp=fft(spec)
    # Omega <- (0:(length(sp) - 1))/length(sp)
    # i <- complex(real = 0, imaginary = 1)
    # Spectrum <- sp* exp(i * pars$a_GRPDLY[s] * 2 * pi * Omega)
    # fid <- (stats::mvfft(t(Spectrum), inverse = TRUE))

    # remove group delay points
    if(pars$a_DSPFVS[s]<20){stop('Implement group delay digital filter correction ofr DSP firmware <20 (DSPFVS & DECIM')}
    fid1=spec[-(1:(pars$a_GRPDLY[s]*2))]
    # plot(spec[1:200], type='l')
    # points(c(rep(0, (pars$a_GRPDLY[s]*2)), (fid1)), type='l', col='cyan')
    # abline(v=76)
    spec_lb=fid1*afun

    if(!is.integer(zf)){stop('Zerofil argument nees to be an integer as it is summand of exponent in log2 space (ie., zerofil=1 doubles , zerofil=2 quadrupoles the number of data points.')}
    add_zeros=(2^(log2(length(spec))+zf) - length(spec_lb))
    #cat(add_zeros, '\n')
    spec_zf=c(spec_lb, rep(0, times= add_zeros))
    #cat(log2(length(spec_zf)), 'number of  data points fid after zerofilling (log2 exp)\n')

    spec_re<-spec_zf[seq(1,length(spec_zf), by=2)]
    spec_im<-spec_zf[seq(2,length(spec_zf), by=2)]
    fid=complex(real = spec_re, imaginary = spec_im)

    sp=fft(fid)

    spli=length(sp)/2
    idx=c((spli+1) : length(sp), 1:spli) # re-arrange -> center is sfo1

    sp_re=Re(sp[idx])
    # flip if neg max
    ii3=1:(length(sp_re)/3)
    if(abs(min(sp_re[ii3])) > max(sp_re[ii3])) sp_re=sp_re*(-1)
    sp_im=Im(sp[idx])
    sp_mag=sp_re+sp_im


    browser()
    # phasing


    ppmDist<-pars$a_SW[s]/length(sp_re)
    # define ppm and set zero for TSP
    idx_tsp=which.max(sp_re[0:(length(sp_re)/3)])
    ppm=c(-(idx_tsp-1):0, 1:((length(sp_re))-idx_tsp))* ppmDist

    switch(type,
           'absorption' ={ sp_out<-sp_re},
           'dispersion' ={ sp_out<-sp_im},
           'magnitude' ={sp_out<-sp_mag}
           )

    fspec=approxfun(ppm, sp_out)


    return(fspec(ppm_ref))

  })

  out=t(out)
  colnames(out)=ppm_ref

  #
  # par(mfrow=c(2,1))
  # matspec(out, ppm_ref, shift=c(2.5,2.6))
  #
  # matspec(out1, ppm_ref, shift=c(2.5,2.6))
  #


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





















































































