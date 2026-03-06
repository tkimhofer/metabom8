#' @title Import 1D NMR spectra (TopSpin processed)
#'
#' @description
#' Imports TopSpin-processed 1D NMR spectra together with spectrometer acquisition
#' and TopSpin processing parameters (`acqus` and `procs`, respectively).
#'
#' @param path Character. Directory path containing Bruker NMR experiments.
#'
#' @param exp_type Named list. Optional filtering specification based on
#' acquisition or processing metadata. Each list element must correspond to a
#' metadata field (e.g. \code{pulprog}, \code{ns}, \code{rg}). Filtering supports:
#' \itemize{
#'   \item Exact match: \code{list(pulprog = "noesygppr1d")}
#'   \item Membership: \code{list(pulprog = c("zg30", "noesygppr1d"))}
#'   \item Numeric range: \code{list(ns = list(range = c(16, 128)))}
#'   \item Generic comparison: \code{list(ns = list(op = ">=", value = 32))}
#' }
#' Multiple fields are combined using logical AND.
#'
#' @param n_max Integer. Maximum number of spectra to import. Default: 1000.
#' @param filter Logical. Filter out experiments with incomplete file systems.
#' @param recursive Logical. Search \code{path} recursively. Default: TRUE.
#' @param verbose Logical or numeric. Verbosity level.
#' @param to_global Logical. If \code{TRUE}, the returned objects are additionally
#' assigned to the global environment.
#'
#' @return
#' A named list with three elements:
#' \describe{
#'   \item{X}{A numeric matrix of spectra (rows = samples, columns = ppm values).}
#'   \item{ppm}{A numeric vector of chemical shift values (ppm).}
#'   \item{meta}{A data frame of acquisition and processing metadata,
#'   row-aligned with \code{X}.}
#' }
#'
#' If \code{to_global = TRUE}, objects with the same names in the global
#' environment will be overwritten.
#'
#' @examples
#' path <- system.file("extdata", package = "metabom8")
#'
#' read1d_proc(path, exp_type = list(pulprog = "noesygppr1d"), n_max = 2)
#'
#' @importFrom fs path_common path_rel
#' @importFrom stats approxfun
#' @export
read1d <- function(path,
                   exp_type = list(pulprog = "noesygppr1d"),
                   n_max = 1000,
                   filter = TRUE,
                   recursive = TRUE,
                   verbose = 1,
                   to_global = FALSE) {
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
    f_spec <- approxfun(x = cs_ppm, y = spec, rule = 2)
    f_spec(ppm_ref)
  }, FUN.VALUE = numeric(length(ppm_ref)))

  out <- t(out)
  colnames(out) <- ppm_ref


  pars$exp_dir <- dirname(f_list$f_1r)

  if (length(f_list$f_1r) ==1){
    x <- strsplit(f_list$f_1r, .Platform$file.sep)[[1]]
    xr <- x[(length(x) - 3):(length(x)-1)]
    fnam <- paste(xr, collapse = .Platform$file.sep)
  } else{
    cpath <- path_common(f_list$f_1r)
    fnam <- path_rel(pars$exp_dir, cpath)
  }

  rownames(out) <- fnam
  rownames(pars) <- fnam
  out[is.na(out)] <- 0


  if (verbose) {
    if (nrow(out)==1) message("Imported ", nrow(out), " spectrum.") else{
      message("Imported ", nrow(out), " spectra.")
    }}

  prov <- list(
    source = list(
      vendor = "bruker",
      input = "topspin_proc",
      processed_upstream = TRUE
    ),
    path = list(
      root = path,
      experiments = f_list$path,
      recursive = recursive,
      filter_incomplete = filter,
      n_max = n_max
    ),
    time = list(
      import_time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    pipeline = list(
      fct = "metabom8::read1d_proc",
      mode = "absorption",
      steps = c("read_1r", "ppm_from_acqus_procs", "interpolation"),
      params = list( exp_type = exp_type)
    ),
    axis = list(
      ppm_source = "topspin_procs",
      common_grid = TRUE,
      interpolation = list(
        applied = TRUE,
        method = "linear",
        rule = 2,
        from = "cs_ppm",
        to = "ppm_ref",
        notes = "Each spectrum interpolated onto common ppm_ref grid."
      )
    ),
    session = list(
      metabom8 = as.character(utils::packageVersion("metabom8")),
      R = paste0(R.version$major, ".", R.version$minor)
    )
  )

  attr(out, "m8_axis") <- list(ppm = ppm_ref)
  attr(out, "m8_provenance") <- prov

  out <- .m8_stamp(
    out,
    step = "import ft, processed spectra (1r)",
    params = list(source_dir = path, exp_filter = .list_to_string(exp_type), n_spectra = nrow(out)),
    notes = "TopSpin-processed spectra imported and interpolated onto a common ppm grid."
  )

  collect_out <- list(X = out, ppm = ppm_ref, meta = pars)

  if (to_global) {
    if (verbose > 0) message("Assigning objects X, ppm, meta to global environment.")
    list2env(collect_out, envir = .GlobalEnv)
    # assign("X", out, envir = .GlobalEnv)
    # assign("ppm", ppm_ref, envir = .GlobalEnv)
    # assign("meta", pars, envir = .GlobalEnv)
  }

  invisible(collect_out)
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
#' @importFrom stats setNames
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
.detect1d_procs <- function(datapath, n_max = 10, filter = TRUE, recursive = TRUE, verbose = 1) {

  datapath <- gsub(paste0(.Platform$file.sep, "$"), "", path.expand(datapath))

  f_procs <- list.files(datapath, pattern="^procs$", full.names=TRUE, recursive=recursive, ignore.case=TRUE)
  f_acqus <- list.files(datapath, pattern="^acqus$", full.names=TRUE, recursive=recursive, ignore.case=TRUE)
  f_1r    <- list.files(datapath, pattern="^1r$",    full.names=TRUE, recursive=recursive, ignore.case=TRUE)

  if (length(f_acqus) == 0 || length(f_procs) == 0 || length(f_1r) == 0) return(NULL)

  # experiment directory = folder containing acqus
  expdir_a <- dirname(f_acqus)

  # for processed files, experiment directory = strip "/pdata/..." tail
  expdir_p <- sub(paste0(.Platform$file.sep, "pdata", .Platform$file.sep, ".*$"), "", f_procs)
  expdir_1 <- sub(paste0(.Platform$file.sep, "pdata", .Platform$file.sep, ".*$"), "", f_1r)

  # keep only experiments that have all three
  common <- Reduce(intersect, list(expdir_a, expdir_p, expdir_1))

  if (length(common) == 0) {
    if (filter) return(NULL)
    stop("No experiments with matching acqus, procs and 1r found.")
  }

  # map to common exp dirs
  ia <- match(common, expdir_a)
  ip <- match(common, expdir_p)
  i1 <- match(common, expdir_1)

  f_acqus <- f_acqus[ia]
  f_procs <- f_procs[ip]
  f_1r    <- f_1r[i1]

  # limit
  if (length(f_1r) > n_max) {
    keep <- seq_len(n_max)
    f_acqus <- f_acqus[keep]; f_procs <- f_procs[keep]; f_1r <- f_1r[keep]
    common <- common[keep]
    if (verbose) message("Reached n_max - not all spectra read-in.")
  }

  return(list(
    path   = common,
    exp_no = basename(common),   # "162" for your example
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

