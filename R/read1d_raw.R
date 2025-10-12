#' @title Read raw FIDs and process to spectra
#'
#' @description
#' Reads Bruker 1D NMR FIDs, corrects the digital filter (group delay),
#' applies apodisation (windowing), optional zero-filling, FFT, phasing and
#' ppm calibration.
#'
#' Returns either absorption-, dispersion-, or magnitude-mode
#' spectra.
#'
#' Results are assigned to the global environment.
#'
#' @param path Character. Path to the root directory containing Bruker experiment folders.
#' @param exp_type Named list. Acquisition-parameter filter to select experiments
#'   (e.g., \code{list(PULPROG = "noesygppr1d")}).
#' @param apodisation Named list. Apodisation function and parameters. \code{fun}
#'   must be one of \code{"exponential"}, \code{"cosine"}, \code{"sine"}, \code{"sem"}.
#' @param zerofill Integer. Zero-filling exponent (\code{1} doubles points, \code{2} quadruples, …).
#' @param mode Character. Spectrum type to return: \code{"absorption"}, \code{"dispersion"}, or \code{"magnitude"}.
#' @param verbose Integer. Verbosity level: \code{0} = silent, \code{1} = summary (default),
#'   \code{2} = detailed, \code{3} = debug.
#' @param recursive Logical. Recursively search subdirectories for FIDs.
#' @param n_max Integer. Maximum number of experiments to process.
#' @param filter Logical. Remove experiments with incomplete file structures.
#' @param to_global Logical. If TRUE, objects are also assigned to the global
#' environment; otherwise, only an invisible list is returned.
#' @details
#' FIDs are read from Bruker acquisition folders and processed by the following pipeline:
#' \enumerate{
#'   \item Digital-filter (group-delay) correction: the initial \emph{n} complex points
#'         are invalid due to the causal DSP decimation filter and are discarded; \emph{n}
#'         equals \code{GRPDLY} when present, or is looked up from Bruker tables indexed by
#'         \code{DECIM} and \code{DSPFVS} on older systems.
#'   \item Apodisation (windowing).
#'   \item Zero-filling (optional).
#'   \item FFT to the frequency domain.
#'   \item Phase correction.
#'   \item PPM calibration (e.g., to TSP).
#' }
#' A common ppm scale is then interpolated across spectra.
#'
#' \strong{Digital filter note:} On newer systems \code{GRPDLY} is written in \code{acqus}/\code{acqu2s}
#' and should be used directly. For older data sets (\code{GRPDLY < 0} or missing), the group delay
#' is derived from \code{DECIM} and \code{DSPFVS} via an internal look-up table.
#'
#' \strong{Apodisation functions}:
#' \itemize{
#'   \item \code{"uniform"}
#'   \item \code{"cosine"},
#'   \item \code{"sine"},
#'   \item \code{"exponential"} (parameter: \code{lb})
#'   \item \code{"sem"} (sine * exponential, parameter: \code{lb})
#'   \item \code{"gauss"} (parameter: \code{lb}, \code{'gb'}, and \code{'para'})
#'   \item \code{"expGaus_resyG"} (parameter: \code{lb}, \code{'gb'}, and \code{'aq_t'})
#' }
#'
#' @return
#' A named list with three elements:
#' \describe{
#'   \item{X}{A numeric matrix of spectra (rows = samples, columns = ppm values).}
#'   \item{ppm}{A numeric vector of chemical-shift axis (ppm).}
#'   \item{meta}{A data frame of acquisition metadata, row-aligned with \code{X}.}
#' }
#'
#' If \code{to_global = TRUE}, these objects are also assigned to the global environment.
#' In that case, any existing objects with the same names will be overwritten.
#'
#' @seealso \code{\link{read1d_proc}} for importing TopSpin-processed spectra
#'
#' @examples
#' path <- system.file("extdata", package = "metabom8")
#' read1d_raw(
#'   path,
#'   exp_type    = list(PULPROG = "noesygppr1d"),
#'   apodisation = list(fun = "exponential", lb = 0.2),
#'   zerofill    = 1,
#'   n_max       = 3
#' )
#'
#' @family NMR
#' @export
read1d_raw <- function(path,
                       exp_type = list(exp = "PROF_PLASMA_CPMG128_3mm", pulprog = "noesygppr1d"),
                       apodisation = list(fun = "exponential", lb = 0.2),
                       zerofill = 1L,
                       mode = "absorption",
                       verbose = 1,
                       recursive = TRUE,
                       n_max = 1000,
                       filter = TRUE,
                       to_global = TRUE) {

  path <- path.expand(path)
  mode <- match.arg(mode, c("absorption", "dispersion", "magnitude"))

  if (verbose > 1) message("Searching for spectral data...")
  f_list <- .check1d_files_fid(path, n_max = 1e6, filter = filter, recursive = recursive, verbose = verbose)

  if (verbose > 1) message("Found ", length(f_list[[1]]), " experiment(s).")

  if (verbose > 1) message("Extracting spectrometer metadata.")
  pars <- .extract_acq_pars1d(f_list)

  if (verbose > 1) message("Filtering experiments using exp_type conditions.")
  exp_filt <- .filterExp_files(pars, exp_type, f_list, n_max)
  f_list <- exp_filt[[1]]
  pars <- exp_filt[[2]]

  if (verbose > 0) {
    if (length(f_list[[1]])==1) message("Processing 1 experiment.") else{
      message("Processing ", length(f_list[[1]]), " experiments.")
    }
  }

  if (length(unique(pars$a_TD)) > 2 || length(unique(pars$a_GRPDLY)) > 1) {
    stop("Unequal TD or GRPDLY across experiments. Cannot proceed.")
  }

  ppm_ref <- .defineChemShiftPpm(
    pars$a_SFO1[1], pars$a_SW_h[1], pars$a_TD[1],
    dref = 4.79, ref = TRUE
  )[, 1]
  ppm_ref <- ppm_ref - ppm_ref[which.min(abs(ppm_ref - 0))]

  if (verbose >= 2) message("Defining apodisation function.")

  if (!'a_GRPDLY' %in% colnames(pars)){
    if (verbose >= 1) message('Group delay (GRPDLY) value not in acqus - using look-up table.')
    pars$a_GRPDLY <- vapply(seq_len(nrow(pars)), function(i){
      get_grpdly(pars$a_DECIM[i], pars$a_DSPFVS[i])
    }, 1.)
  }

  apoFct <- .fidApodisationFct(n = (as.numeric(pars$a_TD[1]) - (floor(as.numeric(pars$a_GRPDLY[1])) * 2)), apodisation)

  if (verbose > 1) message("Reading FIDs and processing spectra.")

  out <- vapply(seq_along(f_list[[1]]), function(s) {
    if (verbose > 1) message(f_list[[1]][s])

    byteorda <- c(little = 0, big = 1)
    dtypa <- c(int = 0, double = 2)

    fid <- readBin(
      file.path(f_list[[1]][s], "fid"),
      what = names(dtypa)[match(pars$a_DTYPA[s], dtypa)],
      n = pars$a_TD[s], size = 4L,
      endian = names(byteorda)[match(pars$a_BYTORDA[s], byteorda)]
    )
    fid <- fid * 2^(-pars$a_NC[s])

    fid_corF <- fid[-(seq_len(floor(pars$a_GRPDLY[s]) * 2))]

    if (length(apoFct) != length(fid_corF)) warning("Apodisation length mismatch.")
    spec_lb <- fid_corF * apoFct

    if (zerofill != floor(zerofill)) stop("zerofill must be integer (e.g., 1, 2).")
    spec_zf <- .zerofill(fid = spec_lb, zf = zerofill, le_ori = length(fid))
    sp <- .cplxFft(spec_zf)[, 1]

    sp_re <- Re(sp)
    sp_im <- Im(sp)
    sp_mag <- sp_re + sp_im

    ppm <- .defineChemShiftPpm(pars$a_SFO1[s], pars$a_SW_h[s], length(sp_re),
                               dref = 4.79, ref = FALSE)

    # Phase
    if (abs(min(sp_re[seq_len(length(sp_re) %/% 3)])) > max(sp_re[seq_len(length(sp_re) %/% 3)])) {
      sp_re <- -sp_re
    }
    sp_re <- .phaseTsp(sp_re, sp_im, ppm, seq(0, pi, by = 0.01), 0,
                       idx_tsp = get_idx(c(-0.15, 0.15), ppm) - 1)[, 1]
    if (abs(min(sp_re[seq_len(length(sp_re) %/% 3)])) > max(sp_re[seq_len(length(sp_re) %/% 3)])) {
      sp_re <- -sp_re
    }

    ppm <- .calibTsp(sp_re, ppm)

    sp_out <- switch(mode,
                     absorption = sp_re,
                     dispersion = sp_im,
                     magnitude = sp_mag)

    fspec <- approxfun(ppm, sp_out)
    spec_out <- fspec(ppm_ref)
    return(spec_out)
  }, FUN.VALUE = ppm_ref)

  out <- t(out)
  colnames(out) <- ppm_ref

  pars$exp_path=f_list[[1]]
  fnam <- strsplit(f_list[[1]], .Platform$file.sep)
  idx_keep <- which(apply(do.call(rbind, fnam), 2, function(x) length(unique(x))) > 1)
  fnam <- vapply(fnam, function(x) paste(x[idx_keep], collapse = .Platform$file.sep), FUN.VALUE = "")
  rownames(out) <- fnam
  rownames(pars) <- fnam


  collect_out <- list(X = out, ppm = ppm_ref, meta = pars)

  if (to_global) {
    if (verbose > 0) message("Assigning objects X, ppm, meta to global environment.")
    list2env(collect_out, envir = .GlobalEnv)
  }

  invisible(collect_out)
}


# Bruker DSP (DECIM × DSPFVS → group delay in complex points)
# filter data obtained from nmrglue
# https://github.com/jjhelmus/nmrglue/blob/master/nmrglue/fileio/bruker.py
# (link last accessed 11/10/25)
bruker_dsp_table <- list(
  `10` = c(
    `2`=44.75, `3`=33.5, `4`=66.625, `6`=59.083333333333333, `8`=68.5625,
    `12`=60.375, `16`=69.53125, `24`=61.020833333333333, `32`=70.015625,
    `48`=61.34375, `64`=70.2578125, `96`=61.505208333333333, `128`=70.37890625,
    `192`=61.5859375, `256`=70.439453125, `384`=61.626302083333333,
    `512`=70.4697265625, `768`=61.646484375, `1024`=70.48486328125,
    `1536`=61.656575520833333, `2048`=70.492431640625
  ),
  `11` = c(
    `2`=46.0, `3`=36.5, `4`=48.0, `6`=50.166666666666667, `8`=53.25,
    `12`=69.5, `16`=72.25, `24`=70.166666666666667, `32`=72.75, `48`=70.5,
    `64`=73.0, `96`=70.666666666666667, `128`=72.5, `192`=71.333333333333333,
    `256`=72.25, `384`=71.666666666666667, `512`=72.125,
    `768`=71.833333333333333, `1024`=72.0625, `1536`=71.916666666666667,
    `2048`=72.03125
  ),
  `12` = c(
    `2`=46.0, `3`=36.5, `4`=48.0, `6`=50.166666666666667, `8`=53.25,
    `12`=69.5, `16`=71.625, `24`=70.166666666666667, `32`=72.125, `48`=70.5,
    `64`=72.375, `96`=70.666666666666667, `128`=72.5, `192`=71.333333333333333,
    `256`=72.25, `384`=71.666666666666667, `512`=72.125,
    `768`=71.833333333333333, `1024`=72.0625, `1536`=71.916666666666667,
    `2048`=72.03125
  ),
  `13` = c(
    `2`=2.75, `3`=2.8333333333333333, `4`=2.875, `6`=2.9166666666666667,
    `8`=2.9375, `12`=2.9583333333333333, `16`=2.96875,
    `24`=2.9791666666666667, `32`=2.984375, `48`=2.9895833333333333,
    `64`=2.9921875, `96`=2.9947916666666667
  )
)

# Helper to get the value; returns NA if not found
get_grpdly <- function(decim, dspfvs, table = bruker_dsp_table, default = NA_real_) {
  vec <- table[[as.character(dspfvs)]]
  if (is.null(vec)) return(default)
  val <- vec[[as.character(decim)]]
  if (is.null(val)) default else unname(val)
}
