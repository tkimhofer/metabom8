#### read fids


#' @title Read-in 1D NMR spectra
#' @export
#' @param path char, File directory containing spectra
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @param apodisation list, apodisation type and parameter
#' @param zerofil int, zero-filling, exponent summand of base 2
#' @param return char, return mode: absorption, dispersion, magnitude
#' @param pulprog char, experiment type to read-in as defined by TopSpin pulprog (e.g. <noesygppr1d>)
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
# @importFrom base list.files readBin seq gsub
#' @importFrom stats approxfun
#' @section

# read1d(path)
# path='/Volumes/Torben Kimhofer/tutdata/'
# #path='/Volumes/Torben_1 1/BariatS/BariatS/dat/Cohort1_NMR_Urine/'
# path='/Volumes/Torben_1 1/Epi_data/airwave2/AIRWAVE_Urine_Rack01_RCM_061014/'
# path='/Users/tk2812/Box Sync/Matt_cys/NMR Raw/NMR_COPIED/Cys_Urine_Rack1_SLL_270814/'
# library(metabom8)
#
# source('R/RawGenerics.R')
#
# a=Sys.time()
 # tt=read1d_raw(path,  n_max=50, filter=T, apodisation=list(fun='exponential', lb=0.2), zerofil=1L, return='absorption', pulprog='<noesygppr1d>')
# b=Sys.time()
# as.numeric(b-a)/nrow(tt)
#save(tt, file='Unphased1d_list.Rdata')

# metabom8::matspec(X, ppm, shift=c(-0.1,0.1), main='None')

# read Bruker 1d new
read1d_raw <- function(path,  n_max=1000, filter=T, apodisation=list(fun='exponential', lb=0.2), zerofil=1L, return='absorption', pulprog='<noesygppr1d>', verbose=TRUE, recursive=TRUE) {

  path=path.expand(path)

  if(!return %in% c('absorption', 'dispersion', 'magnitude')){ type='absorption'; message('Check argument type. Returning absorption spectrum.') }

  if(verbose){message('Looking for spectral data...')}
  # check file system intact
  f_list=.check1d_files_fid(path, n_max, filter)
  if(verbose){message(paste(length(f_list[[1]]), 'experiment files found'))}

  # extract parameters from acqus and procs (for f1 and f2)
  pars <- .extract_acq_pars1d(f_list)

  if(is.list(pars)){
    le=unique(sapply(pars, length))
    if(length(le)==1){
      pars=do.call(rbind, pars)
    }else{
      unam<-unique(names(unlist(pars)))
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

  idx_filt<-which(pars$a_PULPROG == pulprog)
  if(length(idx_filt)==0){stop(paste('Specified experiment (', pulprog, ') not in path. Select one of the following: ', unique(pars$a_PULPROG), sep=''))}

  pars=pars[idx_filt,]
  f_list[[1]]=f_list[[1]][idx_filt]
  f_list[[2]]=f_list[[2]][idx_filt]

  pars[,dtype_num]<-apply(pars[,dtype_num], 2, as.numeric)
  rownames(pars)<-f_list[[2]]

  # chem shift
  ppm_ref=defineChemShiftppm(pars$a_SFO1[1], pars$a_SW_h[1], pars$a_TD[1], dref = 4.79, ref=TRUE)[,1] # 4.79: distance water to TSP
  ppm_ref=ppm_ref-ppm_ref[which.min(abs(ppm_ref-0))]

  if(length(unique(pars$a_TD))>2 || length(unique(pars$a_GRPDLY))>1){stop('Number of points collected in time domain is unqual across experiments.')}
  apoFct<-.fidApodisationFct(n=(pars$a_TD[1]-(pars$a_GRPDLY[1]*2)), apodisation) # subtract group delay from TD (digital filtering artefact)

  # read in binary file and
  out <- sapply(1:length(f_list[[1]]), function(s, pref=ppm_ref, afun=apoFct, zf=zerofil){
    #browser()
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

    # zerofill, fft
    spec_zf=zerofil(fid = spec_lb, zf = zf, le_ori = length(spec))
    sp=cplx_fft(spec_zf)[,1]
    sp_re=Re(sp)
    sp_im=Im(sp)
    sp_mag=sp_re+sp_im

    # define ppm
    ppm=defineChemShiftppm(pars$a_SFO1[s], pars$a_SW_h[s], length(sp_re), dref = 4.79, ref=FALSE) # 4.79: distance water to TSP

    #phasing
    sp_re=phaseTsp(sp_re, sp_im, ppm, seq(0, pi, by=0.01), 0, idx_tsp=get.idx(c(-0.05, 0.05), ppm)-1)[,1]
    if(abs(min(sp_re[0:(length(sp_re)/3)]))>max(sp_re[0:(length(sp_re)/3)])) {sp_re=sp_re*(-1)}


    # calibration
    #browser()
    ppm=calibTsp(sp_re, ppm)

    # browser()

    #return(list(sp_re, sp_im, ppm))
    switch(return,
           'absorption' ={ sp_out<-sp_re},
           'dispersion' ={ sp_out<-sp_im},
           'magnitude' ={sp_out<-sp_mag}
           )

    fspec=approxfun(ppm, sp_out)
    #

    spec_out=fspec(ppm_ref)

    #
    return(spec_out)

    #browser()
  })

  #return(out)

  out=t(out)
  colnames(out)=ppm_ref

  fnam=strsplit(f_list[[1]], .Platform$file.sep)
  idx_keep=which((apply(do.call(rbind, fnam), 2, function(x) length(unique(x))))>1)
  fnam=sapply(fnam, function(x, st=idx_keep){
    paste(x[idx_keep], collapse=.Platform$file.sep)
  })

  rownames(out)=fnam
  rownames(pars)=fnam

  assign("X", out, envir = .GlobalEnv)
  assign("ppm", ppm_ref, envir = .GlobalEnv)
  assign("meta", pars, envir = .GlobalEnv)


}


