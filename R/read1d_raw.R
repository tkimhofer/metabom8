#' @title Read Raw 1D NMR FIDs and Process to Spectra
#'
#' @description
#' Reads Bruker 1D NMR FIDs from directories, applies apodisation and FFT, phases and calibrates the spectrum, and returns absorption, dispersion, or magnitude mode spectra. Results are assigned to the global environment.
#'
#' @param path Character. Path to the root directory containing Bruker experiment folders.
#' @param exp_type Named list. Filter for acquisition parameters to select experiments (e.g., \code{list(exp = "noesygppr1d")}).
#' @param apodisation Named list. Specifies apodisation function and parameters. \code{fun} must be one of: \code{"exponential"}, \code{"cosine"}, \code{"sine"}, \code{"sem"}.
#' @param zerofil Integer. Zerofilling exponent (e.g., 1 doubles points, 2 quadruples).
#' @param return Character. Spectrum type to return: \code{"absorption"}, \code{"dispersion"}, or \code{"magnitude"}.
#' @param verbose Integer. Verbosity level: 0 = silent, 1 = summary (default), 2 = detailed, 3 = debug mode.
#' @param recursive Logical. Whether to recursively search subdirectories for FIDs.
#' @param n_max Integer. Maximum number of experiments to process.
#' @param filter Logical. Whether to remove experiments with incomplete file structures.
#'
#' @details
#' Reads FID data from Bruker acquisition folders and applies digital signal processing:
#' apodisation, zerofilling, FFT, phase correction, and calibration to TSP. A common ppm scale is interpolated across spectra.
#'
#' Apodisation functions:
#' \itemize{
#'   \item \code{"exponential"} (requires \code{lb})
#'   \item \code{"cosine"}, \code{"sine"}, \code{"sem"}
#' }
#'
#' @return
#' The function assigns three objects into the global environment:
#' \describe{
#'   \item{X}{Numeric matrix of spectra (rows = samples, columns = ppm values).}
#'   \item{ppm}{Numeric vector of chemical shift axis (ppm).}
#'   \item{meta}{Data.frame with spectrometer parameters, row-matched to \code{X}.}
#' }
#' Existing objects named \code{X}, \code{ppm}, or \code{meta} will be overwritten.
#'
#' @seealso \code{\link{read1d}} for TopSpin-processed spectra
#'
#' @examples
#' path <- system.file("extdata", package = "metabom8")
#' read1d_raw(path, exp_type = list(exp = "noesygppr1d"),
#'            apodisation = list(fun = "exponential", lb = 0.2), n_max = 3)
#'
#' @family NMR
#' @export
read1d_raw <- function(path,
                       exp_type = list(exp = "PROF_PLASMA_CPMG128_3mm", pulprog = "noesygppr1d"),
                       apodisation = list(fun = "exponential", lb = 0.2),
                       zerofil = 1L,
                       return = "absorption",
                       verbose = 1,
                       recursive = TRUE,
                       n_max = 1000,
                       filter = TRUE) {

  path <- path.expand(path)
  return <- match.arg(return, c("absorption", "dispersion", "magnitude"))

  if (verbose > 1) message("Searching for spectral data...")
  f_list <- .check1d_files_fid(path, n_max = 1e6, filter = filter, recursive = recursive, verbose = verbose)

  if (verbose > 1) message("Found ", length(f_list[[1]]), " experiment(s).")

  if (verbose > 1) message("Extracting spectrometer metadata.")
  pars <- .extract_acq_pars1d(f_list)

  if (verbose > 1) message("Filtering experiments using exp_type conditions.")
  exp_filt <- .filterExp_files(pars, exp_type, f_list, n_max)
  f_list <- exp_filt[[1]]
  pars <- exp_filt[[2]]

  if (verbose > 0) message("Reading ", length(f_list[[1]]), " experiments.")

  if (length(unique(pars$a_TD)) > 2 || length(unique(pars$a_GRPDLY)) > 1) {
    stop("Unequal TD or GRPDLY across experiments. Cannot proceed.")
  }

  ppm_ref <- .defineChemShiftPpm(
    pars$a_SFO1[1], pars$a_SW_h[1], pars$a_TD[1],
    dref = 4.79, ref = TRUE
  )[, 1]
  ppm_ref <- ppm_ref - ppm_ref[which.min(abs(ppm_ref - 0))]

  if (verbose >= 2) message("Defining apodisation function.")
  apoFct <- .fidApodisationFct(n = (pars$a_TD[1] - (floor(pars$a_GRPDLY[1]) * 2)), apodisation)

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

    if (!is.integer(zerofil)) stop("zerofil must be integer (e.g., 1, 2).")
    spec_zf <- .zerofil(fid = spec_lb, zf = zerofil, le_ori = length(fid))
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
                       idx_tsp = get.idx(c(-0.15, 0.15), ppm) - 1)[, 1]
    if (abs(min(sp_re[seq_len(length(sp_re) %/% 3)])) > max(sp_re[seq_len(length(sp_re) %/% 3)])) {
      sp_re <- -sp_re
    }

    ppm <- .calibTsp(sp_re, ppm)

    sp_out <- switch(return,
                     absorption = sp_re,
                     dispersion = sp_im,
                     magnitude = sp_mag)

    fspec <- approxfun(ppm, sp_out)
    spec_out <- fspec(ppm_ref)
    return(spec_out)
  }, FUN.VALUE = ppm_ref)

  out <- t(out)
  colnames(out) <- ppm_ref

  fnam <- strsplit(f_list[[1]], .Platform$file.sep)
  idx_keep <- which(apply(do.call(rbind, fnam), 2, function(x) length(unique(x))) > 1)
  fnam <- vapply(fnam, function(x) paste(x[idx_keep], collapse = .Platform$file.sep), FUN.VALUE = "")
  rownames(out) <- fnam
  rownames(pars) <- fnam

  if (verbose > 0) message("Assigning objects X, ppm, meta to global environment.")
  assign("X", out, envir = .GlobalEnv)
  assign("ppm", ppm_ref, envir = .GlobalEnv)
  assign("meta", pars, envir = .GlobalEnv)
}
