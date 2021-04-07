#### read fids


#' @title Read-in 1D NMR FIDs and process to spectra
#' @export
#' @param path char, path to file directory containing spectra
#' @param exp_type named list, filter for acquisition paramters of experiments to read-in (see Details)
#' @param apodisation named list, apodisation function and its parameters (see Details)
#' @param zerofil int, amount of zeros to append to FID given as exponent added of base 2 (see Details)
#' @param return char, return mode of spectra: absorption, dispersion or magnitude mode
#' @param verbose num, different verbose levels: 0 (no info), 1 (overview), 2 (detailed), 3 (step-by-step for debugging)
#' @param recursive logic, if TRUE recursively search all subfolders of path for specified NMR files
#' @param n_max int, maximum number of experiments to read-in
#' @param filter logic, remove experiments with incomplete file systems (TRUE is recommended)
#' @details
#' In the first step, read-in are FIDs generated with experimental condition(s) specified with the exp_type argument. This represents a list with each element representing a parameter condition, named according to spectrometer parameters listed in \emph{acqus} file. For example, to read standard 1D NMR experiments use \code{exp_type=list(exp='noesygppr1d')}. More than one argument can be provided as list element.
#'
#' The apodisation argument is a named list specifying the function name in element \emph{fun} and functions-specific paramter arguments. There are the following different apodisation functions and arguments:
#' #' @return
#' \itemize{
#'   \item exponential, arguments: lb (line broadening factor)
#'   \item cosine, no further arguments
#'   \item sine, no further arguments
#'   \item sem, combined sine-bell - exponential fct: arguments: lb (line broadening factor)
#' }
#'
#' The zerofil argument specifies the amount of zeros to append to the FID and is expressed as exponent addand in the binary numeral system: \code{2^(1+x)}, with x being the zerofil parameter argument. Hence, \code{zerofil=1} doubles the amount of data points.
#'
# @return Three objects: NMR data matrix (2D: rows=spectra, cols=chem shift
# variables), ppm num vector matched to NMR data columns, meta data.frame
# containing spectrometer metadata
#' @return
#' The function exports the following three objects into the currently active R environment (no variable assignments needed):
#' \itemize{
#'   \item X, num matrix:  NMR data, spectra in rows
#'   \item ppm, num array - chemical shift positions, length matches to columns in X
#'   \item meta, data.frame - spectrometer metadata as extracted from individual \emph{acqus} files, row-matched to X
#' }
#' Objects in the R environment with the same variable names will be overwritten.
#' @examples
#' path<-system.file('extdata/',  package = 'metabom8')
#' read1d_raw(path,  exp_type=list(exp='PROF_PLASMA_NOESY'), apodisation=list(fun='exponential', lb=0.2), n_max=3)
#' @author \email{torben.kimhofer@@murdoch.edu.au}
# @seealso \code{\reference{read1d}}
#' @importFrom stats approxfun
#' @family NMR
#' @seealso \code{\link[=read1d]{Import TopSpin processed spectra}}
#' @section

read1d_raw <- function(path, exp_type = list(exp = c("PROF_PLASMA_CPMG128_3mm"),
    pulprog = c("noesygppr1d")), apodisation = list(fun = "exponential", lb = 0.2),
    zerofil = 1L, return = "absorption", verbose = 1, recursive = TRUE, n_max = 1000,
    filter = TRUE) {

    path <- path.expand(path)

    if (!return %in% c("absorption", "dispersion", "magnitude")) {
        return <- "absorption"
        message("Check argument type. Returning absorption spectrum.")
    }

    if (verbose > 1) {
        message("Searching for spectral data...")
    }
    # check file system intact
    f_list <- .check1d_files_fid(path, n_max, filter, recursive, verbose)
    if (verbose > 1) {
        message("Found", paste(length(f_list[[1]]), "experiments files in path."))
    }

    if (verbose > 1) {
        message("Extracting spectrometer meta-data.")
    }
    pars <- .extract_acq_pars1d(f_list)

    if (verbose > 1) {
        message("Filtering for experiments using user-defined parameters (ext_type argument)")
    }
    exp_filt <- .filterExp_files(pars, exp_type, f_list)
    f_list <- exp_filt[[1]]
    pars <- exp_filt[[2]]
    if (verbose > 0) {
        message("Reading ", length(f_list[[1]]), " experiments.")
    }

    # chem shift
    if (verbose >= 2) {
        message("Defining chemical shift axis.")
    }
    ppm_ref <- .defineChemShiftPpm(pars$a_SFO1[1], pars$a_SW_h[1], pars$a_TD[1],
        dref = 4.79, ref = TRUE)[, 1]  # 4.79: distance water to TSP
    ppm_ref <- ppm_ref - ppm_ref[which.min(abs(ppm_ref - 0))]

    if (length(unique(pars$a_TD)) > 2 || length(unique(pars$a_GRPDLY)) > 1) {
        stop("Number of points collected in time domain is unqual across experiments.")
    }

    if (verbose >= 2) {
        message("Defining apidisation function.")
    }
    apoFct <- .fidApodisationFct(n = (pars$a_TD[1] - (floor(pars$a_GRPDLY[1]) * 2)),
        apodisation)  # subtract group delay from TD (digital filtering artefact)

    if (verbose > 1) {
        message("Reading FIDs and transform to spectra.")
    }
    # read in binary file and
    out <- vapply(seq(f_list[[1]]), function(s, pref = ppm_ref, afun = apoFct, zf = zerofil) {
        if (verbose > 1) {
            message(f_list[[1]][s])
        }
        byteorda <- c(little = 0, big = 1)
        dtypa <- c(int = 0, double = 2)
        if (verbose == 3) {
            message("Read FID")
        }
        fid <- readBin(paste0(f_list[[1]][s], .Platform$file.sep, "fid"), what = names(dtypa)[match(pars$a_DTYPA[s],
            dtypa)], n = pars$a_TD[s], size = 4L, endian = names(byteorda)[match(pars$a_BYTORDA[s],
            byteorda)])
        fid <- (fid * (2^pars$a_NC[s]))

        # remove group delay points
        if (verbose == 3) {
            message("Adjust for group delay (=dig filter)")
        }
        if (pars$a_DSPFVS[s] < 20) {
            stop("Implement group delay digital filter correction ofr DSP firmware <20 (DSPFVS & DECIM")
        }
        fid_corF <- fid[-(seq(floor(pars$a_GRPDLY[s]) * 2))]

        # apppodisatoin
        if (verbose == 3) {
            message("Apodisation fct")
        }
        if (length(afun) != length(fid_corF)) {
            #browser()
            message("Apd fct mismatch in length w spec")
        }
        spec_lb <- fid_corF * afun

        if (!is.integer(zf)) {
            stop("Zerofil argument nees to be an integer as it is summand of exponent in log2 space (ie., zerofil=1 doubles , zerofil=2 quadruples the number of data points.")
        }

        # zerofill, fft
        if (verbose == 3) {
            message("Zero-fill, fft")
        }
        spec_zf <- .zerofil(fid = spec_lb, zf = zf, le_ori = length(fid))
        sp <- .cplxFft(spec_zf)[, 1]

        sp_re <- Re(sp)
        sp_im <- Im(sp)
        sp_mag <- sp_re + sp_im

        # define ppm
        if (verbose == 3) {
            message("Establish indiv. ppm scale")
        }
        ppm <- .defineChemShiftPpm(pars$a_SFO1[s], pars$a_SW_h[s], length(sp_re),
            dref = 4.79, ref = FALSE)  # 4.79: distance water to TSP

        # phasing
        if (verbose == 3) {
            message("\tPhasing")
        }

        if (abs(min(sp_re[0:(length(sp_re)/3)])) > max(sp_re[0:(length(sp_re)/3)])) {
            sp_re <- sp_re * (-1)
        }
        sp_re <- .phaseTsp(sp_re, sp_im, ppm, seq(0, pi, by = 0.01), 0, idx_tsp = get.idx(c(-0.15,
            0.15), ppm) - 1)[, 1]
        if (abs(min(sp_re[0:(length(sp_re)/3)])) > max(sp_re[0:(length(sp_re)/3)])) {
            sp_re <- sp_re * (-1)
        }


        # calibration
        if (verbose == 3) {
            message("Calibrate w TSP")
        }
        ppm <- .calibTsp(sp_re, ppm)

        switch(return, absorption = {
            sp_out <- sp_re
        }, dispersion = {
            sp_out <- sp_im
        }, magnitude = {
            sp_out <- sp_mag
        })

        if (verbose == 3) {
            message("Approx. spec to common ppm scale")
        }
        fspec <- approxfun(ppm, sp_out)
        spec_out <- fspec(pref)
        if (verbose == 3) {
            message("###")
        }
        return(spec_out)
    }, FUN.VALUE = ppm_ref)


    out <- t(out)
    colnames(out) <- ppm_ref

    if (verbose == 3) {
        message("Prep rownames for X and ppm")
    }
    fnam <- strsplit(f_list[[1]], .Platform$file.sep)
    idx_keep <- which((apply(do.call(rbind, fnam), 2, function(x) length(unique(x)))) >
        1)
    fnam <- vapply(fnam, function(x, st = idx_keep) {
        paste(x[idx_keep], collapse = .Platform$file.sep)
    }, FUN.VALUE = "")

    rownames(out) <- fnam
    rownames(pars) <- fnam
    if (verbose > 0) {
        message("Adding objects X, ppm and meta to the global workspace.")
    }
    assign("X", out, envir = .GlobalEnv)
    assign("ppm", ppm_ref, envir = .GlobalEnv)
    assign("meta", pars, envir = .GlobalEnv)


}


