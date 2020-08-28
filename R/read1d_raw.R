#### read fids


#' @title Read-in 1D NMR spectra
#' @export
#' @param path char, File directory containing spectra
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @param apodisation list, apodisation type and parameter
#' @param zerofil int, zero-filling, exponent summand of base 2
#' @param return char, return mode: absorption, dispersion, magnitude
#' @param exp_type list, specific acquisition paramters use to select spectra for read-in (e.g. experiment type or pulse sequence)
#' @param verbose num, different verbose levels: 0 (no info), 1 (overview), 2 (detailed), 3 (debugging mode detail)
#' @param recursive logic, if TRUE recursively search all subfolders of path for specified NMR files
#' @return Three objects: NMR data matrix (2D: rows=spectra, cols=chem shift variables), ppm num vector matched to NMR data columns, meta data.frame containing spectrometer metadata
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @family Import NMR data functions
# @seealso \code{\reference{read1d}}
#' @importFrom stats approxfun
#' @section

# read1d(path)
# path='/Volumes/Torben Kimhofer/tutdata/'
# #path='/Volumes/Torben_1 1/BariatS/BariatS/dat/Cohort1_NMR_Urine/'
# path='/Volumes/Torben_1 1/Epi_data/airwave2/AIRWAVE_Urine_Rack01_RCM_061014/'
# path='/Users/tk2812/Box Sync/Matt_cys/NMR Raw/NMR_COPIED/Cys_Urine_Rack1_SLL_270814/'
# library(metabom8)
# path='/Volumes/Torben_2/adjNS_LTRexp/'
#
# source('R/RawGenerics.R')
#
# a=Sys.time()
 # tt=read1d_raw(path,  n_max=50, filter=T, apodisation=list(fun='exponential', lb=0.2), zerofil=1L, return='absorption', pulprog='<cpmgesgp1d>')
# b=Sys.time()
# as.numeric(b-a)/nrow(tt)
#save(tt, file='Unphased1d_list.Rdata')

# metabom8::matspec(X, ppm, shift=c(-0.1,0.1), main='None')

# read Bruker 1d new
read1d_raw<-function(path,  n_max=1000, filter=TRUE, apodisation=list(fun='exponential', lb=0.2), zerofil=1L, return='absorption', exp_type=list(exp=c('PROF_PLASMA_CPMG128_3mm', 'PROF_PLASMA_NOESY128_3mm'), pulprog=c('noesygppr1d')),  verbose=1, recursive=TRUE) {

  path=path.expand(path)

  if(!return %in% c('absorption', 'dispersion', 'magnitude')){ return='absorption'; message('Check argument type. Returning absorption spectrum.') }

  if(verbose>1){message('Searching for spectral data...')}
  # check file system intact
  f_list<-.check1d_files_fid(path, n_max, filter, recursive, verbose)
  if(verbose>1){message('Found', paste(length(f_list[[1]]), 'experiments files in path.'))}

  if(verbose>1){message('Extracting spectrometer meta-data.')}
  pars <- .extract_acq_pars1d(f_list)

  if(verbose>1){message('Filtering for experiments using user-defined parameters (ext_type argument)')}
  exp_filt<-.filterExp_files(pars, exp_type, f_list)
  f_list=exp_filt[[1]]
  pars<-exp_filt[[2]]
  if(verbose>0){message('Reading ', length(f_list[[1]]), ' experiments.')}

  # chem shift
  if(verbose>=2){message('Defining chemical shift axis.')}
  ppm_ref=.defineChemShiftPpm(pars$a_SFO1[1], pars$a_SW_h[1], pars$a_TD[1], dref = 4.79, ref=TRUE)[,1] # 4.79: distance water to TSP
  ppm_ref=ppm_ref-ppm_ref[which.min(abs(ppm_ref-0))]

  if(length(unique(pars$a_TD))>2 || length(unique(pars$a_GRPDLY))>1){stop('Number of points collected in time domain is unqual across experiments.')}

  if(verbose>=2){message('Defining apidisation function.')}
  apoFct<-.fidApodisationFct(n=(pars$a_TD[1]-(floor(pars$a_GRPDLY[1])*2)), apodisation) # subtract group delay from TD (digital filtering artefact)

  if(verbose>1){message('Reading FIDs and transform to spectra.')}
  # read in binary file and
  out <- vapply(seq(f_list[[1]]), function(s, pref=ppm_ref, afun=apoFct, zf=zerofil){
    if(verbose>1){message(f_list[[1]][s])}
    byteorda<-c('little'=0, 'big'=1)
    dtypa<-c('int'=0, 'double'=2)
    if(verbose==3){message('Read FID')}
    fid <- readBin(paste0(f_list[[1]][s], .Platform$file.sep, 'fid'),
                    what = names(dtypa)[match(pars$a_DTYPA[s], dtypa)],
                    n = pars$a_TD[s],
                    size = 4L,
                    endian = names(byteorda)[match(pars$a_BYTORDA[s], byteorda)]
    )
    fid <- ( fid * (2^pars$a_NC[s]) )

    # remove group delay points
    if(verbose==3){message('Adjust for group delay (=dig filter)')}
    if(pars$a_DSPFVS[s]<20){stop('Implement group delay digital filter correction ofr DSP firmware <20 (DSPFVS & DECIM')}
    fid_corF<-fid[-(seq_along(floor(pars$a_GRPDLY[s])*2))]

    # apppodisatoin
    if(verbose==3){message('Apodisation fct')}
    if(length(afun)!=length(fid_corF)){message('Apd fct mismatch in length w spec')}
    spec_lb<-fid_corF*afun

    if(!is.integer(zf)){stop('Zerofil argument nees to be an integer as it is summand of exponent in log2 space (ie., zerofil=1 doubles , zerofil=2 quadruples the number of data points.')}

    # zerofill, fft
    if(verbose==3){message('Zero-fill, fft')}
    spec_zf<-.zerofil(fid = spec_lb, zf = zf, le_ori = length(fid))
    sp<-.cplxFft(spec_zf)[,1]

    sp_re<-Re(sp)
    sp_im<-Im(sp)
    sp_mag<-sp_re+sp_im

    # define ppm
    if(verbose==3){message('Establish indiv. ppm scale')}
    ppm<-.defineChemShiftPpm(pars$a_SFO1[s], pars$a_SW_h[s], length(sp_re), dref = 4.79, ref=FALSE) # 4.79: distance water to TSP

    #phasing
    if(verbose==3){message('\tPhasing')}

    if(abs(min(sp_re[0:(length(sp_re)/3)]))>max(sp_re[0:(length(sp_re)/3)])) {sp_re<-sp_re*(-1)}
    sp_re<-.phaseTsp(sp_re, sp_im, ppm, seq(0, pi, by=0.01), 0, idx_tsp=get.idx(c(-0.15, 0.15), ppm)-1)[,1]
    if(abs(min(sp_re[0:(length(sp_re)/3)]))>max(sp_re[0:(length(sp_re)/3)])) {sp_re<-sp_re*(-1)}


    # calibration
    if(verbose==3){message('Calibrate w TSP')}
    ppm<-.calibTsp(sp_re, ppm)

    switch(return,
           'absorption'={sp_out<-sp_re},
           'dispersion'={sp_out<-sp_im},
           'magnitude'={sp_out<-sp_mag}
           )

    if(verbose==3){message('Approx. spec to common ppm scale')}
    fspec<-approxfun(ppm, sp_out)
    spec_out<-fspec(pref)
    if(verbose==3){message('###')}
    return(spec_out)
  }, FUN.VALUE = ppm_ref)


  out<-t(out)
  colnames(out)<-ppm_ref

  if(verbose==3){message('Prep rownames for X and ppm')}
  fnam<-strsplit(f_list[[1]], .Platform$file.sep)
  idx_keep<-which((apply(do.call(rbind, fnam), 2, function(x) length(unique(x))))>1)
  fnam<-vapply(fnam, function(x, st=idx_keep){
    paste(x[idx_keep], collapse=.Platform$file.sep)
  }, FUN.VALUE='')

  rownames(out)<-fnam
  rownames(pars)<-fnam
  if(verbose>0){message('Returning spectral varaibles X, ppm, meta')}
  assign("X", out, envir = .GlobalEnv)
  assign("ppm", ppm_ref, envir = .GlobalEnv)
  assign("meta", pars, envir = .GlobalEnv)


}


