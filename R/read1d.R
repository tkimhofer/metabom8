#' @title Import 1D NMR spectra (TopSpin processed)
#'
#' @description
#' Imports Bruker TopSpin-processed 1D NMR spectra. Spectrometer and processing
#' parameters are extracted from `acqus` and `procs` files. Optionally filters
#' experiments using acquisition metadata.
#'
#' @param path Character. Directory path containing Bruker NMR experiments.
#' @param exp_type Named list. Filter experiments by acquisition parameters (e.g., list(pulprog = 'noesygppr1d')).
#' @param n_max Integer. Maximum number of spectra to import. Default: 1000.
#' @param filter Logical. Filter out experiments with incomplete file systems.
#' @param recursive Logical. Search path recursively. Default: TRUE.
#' @param verbose Logical or numeric. Verbosity level.
#'
#' @return
#' The following objects are assigned to the global environment:
#' \itemize{
#'   \item \code{X}: Numeric matrix of spectra (samples x variables)
#'   \item \code{ppm}: Numeric vector of chemical shift values
#'   \item \code{meta}: Data frame of acquisition/processing metadata
#' }
#' Also returns a (named) list invisibly for testing.
#'
#' @examples
#' path <- system.file("extdata", package = "metabom8")
#' read1d(path, exp_type = list(pulprog = "noesygppr1d"), n_max = 2)
#'
#' @family NMR
#' @export
read1d <- function(path,
                   exp_type = list(pulprog = "noesygppr1d"),
                   n_max = 1000,
                   filter = TRUE,
                   recursive = TRUE,
                   verbose = TRUE) {
  path <- path.expand(path)

  if (as.character(match.call()[[1]]) == "read1d") {
    warning("`read1d` will be deprecated; please use `read1d_proc()`.", call. = FALSE)
  }

  f_list <- .detect1d_procs(path, n_max = n_max, filter = filter, recursive = recursive, verbose = verbose)
  if (is.null(f_list)) stop("No valid experiments found.")

  pars <- .extract_pars1d(f_list)
  exp_filt <- .filterExp_files(pars, exp_type, f_list, n_max = n_max)
  f_list <- exp_filt[[1]]
  pars <- exp_filt[[2]]

  ppm_ref <- .chemShift(swidth = pars$a_SW[1], offset = pars$p_OFFSET[1], si = pars$p_SI[1])

  out <- vapply(seq_along(f_list[[1]]), function(s) {
    cs_ppm <- .chemShift(swidth = pars$a_SW[s], offset = pars$p_OFFSET[s], si = pars$p_SI[s])
    byteord <- c(little = 0, big = 1)
    spec <- readBin(
      f_list$f_1r[s], what = "int", n = pars$p_FTSIZE[s], size = 4,
      signed = TRUE, endian = names(byteord)[match(pars$p_BYTORDP[s], byteord)]
    )
    spec <- spec * (2 ^ pars$p_NC_proc[s])
    f_spec <- approxfun(x = cs_ppm, y = spec)
    f_spec(ppm_ref)
  }, FUN.VALUE = numeric(length(ppm_ref)))

  out <- t(out)
  colnames(out) <- ppm_ref

  fnam <- vapply(strsplit(f_list$f_1r, .Platform$file.sep), function(x) {
    paste(x, collapse = .Platform$file.sep)
  }, character(1))

  rownames(out) <- fnam
  rownames(pars) <- fnam
  out[is.na(out)] <- 0

  assign("X", out, envir = .GlobalEnv)
  assign("ppm", ppm_ref, envir = .GlobalEnv)
  assign("meta", pars, envir = .GlobalEnv)

  if (verbose) message("Imported ", nrow(out), " spectra.")

  invisible(list(X = out, ppm = ppm_ref, meta = pars))
}

#' @rdname read1d
#' @export
read1d_proc <- read1d


#' @title Read Bruker NMR Parameter Files
#'
#' @description
#' Helper function used by `read1d()` to extract acquisition (`acqus`) and processing (`procs`) parameters
#' from Bruker-formatted 1D NMR experiments.
#'
#' @param f_list List of file paths. Output from `.detect1d_procs()`.
#'
#' @return A data frame of extracted acquisition and processing metadata. Row names correspond to spectrum filenames.
#'
#' @keywords internal
.extract_pars1d <- function(f_list) {
  out <- lapply(seq_along(f_list[[1]]), function(i) {
    # Read procs file
    f_procs_lines <- readLines(f_list$f_procs[i], warn = FALSE)
    idx_procs <- grep("..", f_procs_lines, fixed = TRUE)
    f_procs_lines[idx_procs] <- vapply(idx_procs, function(j) {
      gsub(" .*", f_procs_lines[j + 1], f_procs_lines[j])
    }, FUN.VALUE = "")
    procs_vals <- strsplit(gsub("^##\\$", "", grep("^##\\$", f_procs_lines, value = TRUE), fixed = FALSE), "=")
    d_procs_val <- gsub("^ ", "", vapply(procs_vals, "[[", 2, FUN.VALUE = ""))
    names(d_procs_val) <- paste0("p_", vapply(procs_vals, "[[", 1, FUN.VALUE = ""))

    # Read acqus file
    f_acqus_lines <- readLines(f_list$f_acqus[i], warn = FALSE)
    idx_acqus <- grep("..", f_acqus_lines, fixed = TRUE)
    f_acqus_lines[idx_acqus] <- vapply(idx_acqus, function(j) {
      gsub(" .*", f_acqus_lines[j + 1], f_acqus_lines[j])
    }, FUN.VALUE = "")
    acqus_vals <- strsplit(gsub("^##\\$", "", grep("^##\\$", f_acqus_lines, value = TRUE), fixed = FALSE), "=")
    d_acqu_val <- gsub("^ ", "", vapply(acqus_vals, "[[", 2, FUN.VALUE = ""))
    names(d_acqu_val) <- paste0("a_", vapply(acqus_vals, "[[", 1, FUN.VALUE = ""))

    # Convert datetime
    idx_date <- grep("date", names(d_acqu_val), ignore.case = TRUE)
    d_acqu_val[idx_date] <- as.character(
      as.POSIXct("1970-01-01 00:00:00", tz = "UTC") + as.numeric(d_acqu_val[idx_date])
    )

    c(d_acqu_val, d_procs_val)
  })

  # Handle unequal metadata entries
  out_len <- vapply(out, length, FUN.VALUE = 1L)
  if (length(unique(out_len)) > 1) {
    all_names <- unique(unlist(lapply(out, names)))
    out <- lapply(out, function(x) {
      out_row <- setNames(rep(NA, length(all_names)), all_names)
      out_row[names(x)] <- x
      out_row
    })
  }

  out_df <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)

  is_num <- vapply(out_df, .is_numeric_trycatch, logical(1))
  # is_num <- vapply(out_df, function(col) all(!is.na(as.numeric(col))), logical(1))
  out_df[is_num] <- lapply(out_df[is_num], as.numeric)

  rownames(out_df) <- f_list[[2]]
  return(out_df)
}


#' @title Check for Intact File Systems - Helper Function for read1d
#' @description
#' Scans the specified directory recursively to find intact sets of Bruker NMR files:
#' 'procs', 'acqus', and '1r' files. Optionally filters incomplete experiments.
#'
#' @param datapath character. Path to directory containing spectra.
#' @param n_max integer. Maximum number of spectra to read (default 10).
#' @param filter logical. Whether to filter out incomplete file sets (default TRUE).
#' @param recursive logical. Whether to search directories recursively.
#' @param verbose integer. Verbosity level for messaging.
#'
#' @return
#' A list with elements:
#' \describe{
#'   \item{path}{Character vector of experiment folder paths (intact).}
#'   \item{exp_no}{Character vector of experiment folder IDs.}
#'   \item{f_procs}{Character vector of paths to 'procs' files.}
#'   \item{f_acqus}{Character vector of paths to 'acqus' files.}
#'   \item{f_1r}{Character vector of paths to '1r' files.}
#' }
#'
#' @keywords internal
.detect1d_procs <- function(datapath,
                            n_max = 10,
                            filter = TRUE,
                            recursive = TRUE,
                            verbose = 1) {

  stopifnot(is.character(datapath), length(datapath) == 1)
  stopifnot(is.numeric(n_max), length(n_max) == 1, n_max > 0)
  stopifnot(is.logical(filter), length(filter) == 1)
  stopifnot(is.logical(recursive), length(recursive) == 1)
  stopifnot(is.numeric(verbose), length(verbose) == 1)

  datapath <- gsub(paste0(.Platform$file.sep, "$"), "", datapath)

  # List files matching patterns recursively, if requested
  f_procs <- list.files(path = datapath, pattern = "^procs$", full.names = TRUE,
                        recursive = TRUE, ignore.case = TRUE)
  f_acqus <- list.files(path = datapath, pattern = "^acqus$", full.names = TRUE,
                        recursive = TRUE, ignore.case = TRUE)
  f_1r <- list.files(path = datapath, pattern = "^1r$", full.names = TRUE,
                     recursive = TRUE, ignore.case = TRUE)

  # Extract experiment IDs by trimming known suffixes and folder patterns
  pattern_acqus <- paste0("^", datapath, .Platform$file.sep, "|", .Platform$file.sep, "acqus$")
  pattern_procs <- paste0("^", datapath, .Platform$file.sep, "|", .Platform$file.sep, "pdata", .Platform$file.sep, ".*")
  pattern_1r <- paste0("^", datapath, .Platform$file.sep, "|", .Platform$file.sep, "pdata", .Platform$file.sep, ".*")

  id_a <- gsub(pattern_acqus, "", f_acqus)
  id_p <- gsub(pattern_procs, "", f_procs)
  id_f1r <- gsub(pattern_1r, "", f_1r)

  # Find matching indices (keep only matching files)
  idx_a <- which(id_a %in% id_f1r)
  idx_p <- which(id_p %in% id_f1r)
  idx_f1r <- which(id_f1r %in% id_a)

  # Filter incomplete sets if requested
  if (length(idx_a) != length(id_a) || length(idx_f1r) != length(id_f1r) || length(idx_p) != length(id_p)) {
    if (filter) {
      if (verbose > 1) {
        message("Reading experiments with matching procs, acqus and 1r files.")
      }
      f_acqus <- f_acqus[idx_a]
      id_a <- id_a[idx_a]
      f_procs <- f_procs[idx_p]
      id_p <- id_p[idx_p]
      f_1r <- f_1r[idx_f1r]
      id_f1r <- id_f1r[idx_f1r]
    } else {
      message("File system seems to be corrupt for some experiments. Consider setting `filter=TRUE`.")
      return(NULL)
    }
  }

  # Ensure matching between acqus, procs, and 1r files
  idm_f1r <- match(id_a, id_f1r)
  idm_p <- match(id_a, id_p)
  if (any(is.na(idm_f1r)) || any(is.na(idm_p))) {
    stop("Mismatch in matching file indices between acqus, procs, and 1r files.")
  }

  f_procs <- f_procs[idm_p]
  f_1r <- f_1r[idm_f1r]
  id_f1r <- id_f1r[idm_f1r]

  # Sanity check lengths are consistent
  if (length(unique(c(length(f_acqus), length(f_1r), length(f_procs)))) != 1) {
    stop("Lengths of filtered files (acqus, 1r, procs) do not match after filtering.")
  }

  # Limit to n_max spectra
  if (n_max < length(f_1r)) {
    idx_nmax <- seq_len(n_max)
    f_acqus <- f_acqus[idx_nmax]
    f_procs <- f_procs[idx_nmax]
    f_1r <- f_1r[idx_nmax]
    id_f1r <- id_f1r[idx_nmax]
    message("Reached n_max - not all spectra read-in.")
  }

  # Remove trailing "/acqus" from path for experiment folder
  p_intact <- gsub(paste0(.Platform$file.sep, "acqus$"), "", f_acqus)

  # Return list of intact experiment data
  return(list(
    path = p_intact,
    exp_no = id_f1r,
    f_procs = f_procs,
    f_acqus = f_acqus,
    f_1r = f_1r
  ))
}


#' @title Calculate Chemical Shift Axis
#'
#' @description
#' Computes a 1D NMR chemical shift axis based on sweep width, offset, and size.
#'
#' @param swidth Numeric. Sweep width (Hz).
#' @param offset Numeric. Frequency offset (ppm).
#' @param si Integer. Number of data points.
#'
#' @return Numeric vector of chemical shift values (ppm).
#'
#' @keywords internal
.chemShift <- function(swidth, offset, si) {
  dppm <- swidth / (si - 1)
  seq(offset, offset - swidth, by = -dppm)
}

