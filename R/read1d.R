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

  f_list <- .detect1d_procs(path, n_max = 1e6, filter = filter, recursive = recursive, verbose = verbose)
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

    # Convert date from Bruker format
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

  # Convert numeric columns
  is_num <- vapply(out_df, function(col) all(!is.na(as.numeric(col))), logical(1))
  out_df[is_num] <- lapply(out_df[is_num], as.numeric)

  rownames(out_df) <- f_list[[2]]
  return(out_df)
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

