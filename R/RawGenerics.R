#' @keywords internal
.em <- function(n, lb) {
  idx <- seq_len(n)
  out <- exp(-(idx * lb * pi) / n)
  minmax(out)
}

#' @keywords internal
.cosine <- function(n) {
  idx <- seq_len(n)
  out <- cos(pi * idx / n)
  minmax(out)
}

#' @keywords internal
.sine <- function(n) {
  idx <- seq_len(n)
  out <- sin(pi * idx / n)
  minmax(out)
}

#' @keywords internal
.sem <- function(n, lb = 1.5) {
  idx <- seq_len(n)
  out <- sin(pi * idx / n) * exp(-(idx * lb * pi) / n)
  minmax(out)
}

#' @keywords internal
.expGaus_resyG <- function(n, lb = -10, gb = 0.12, aq_t) {
  idx <- seq_len(n)
  b <- -lb / (2 * gb * aq_t)
  out <- exp((-lb * idx) - (b * idx^2))
  minmax(out)
}

#' @keywords internal
.gauss <- function(n, lb = -1.2, gb = 0.3, para) {
  sw_hz <- para$a_SW_h
  AQ <- para$a_TD / (2 * sw_hz)
  idx <- seq_len(n) - 1
  e <- (-pi * idx * lb) / sw_hz
  g <- ((lb * pi) / (2 * gb * AQ)) * (idx / sw_hz)
  out <- exp(e - g^2)
  minmax(out)
}

#' @title Apodisation function dispatcher
#' @description Applies one of several apodisation functions to FID data.
#' @param n Integer, number of data points.
#' @param pars A list containing the apodisation parameters: must include `fun` and any required shape parameters (e.g. `lb`, `gb`, `para`).
#' @return A numeric vector of apodisation weights.
#' @keywords internal
.fidApodisationFct <- function(n, pars) {
  if (is.null(pars$fun) || !pars$fun %in% c("uniform", "exponential", "cosine", "sine", "sem", "expGaus_resyG", "gauss")) {
    stop("Invalid apodisation function specified. Choose one of: 'uniform', 'exponential', 'cosine', 'sine', 'sem', 'expGaus_resyG', or 'gauss'.")
  }

  afun <- switch(pars$fun,
                 "uniform" = rep(1, n),
                 "exponential" = {
                   if (!"lb" %in% names(pars)) stop("Exponential apodisation requires 'lb' parameter.")
                   .em(n, pars$lb)
                 },
                 "cosine" = .cosine(n),
                 "sine" = .sine(n),
                 "sem" = {
                   if (!"lb" %in% names(pars)) stop("SEM apodisation requires 'lb' parameter.")
                   .sem(n, pars$lb)
                 },
                 "expGaus_resyG" = {
                   if (!all(c("lb", "gb", "aq_t") %in% names(pars))) {
                     stop("expGaus_resyG apodisation requires 'lb', 'gb', and 'aq_t' parameters.")
                   }
                   .expGaus_resyG(n, pars$lb, pars$gb, pars$aq_t)
                 },
                 "gauss" = {
                   if (!all(c("lb", "gb", "para") %in% names(pars))) {
                     stop("Gaussian apodisation requires 'lb', 'gb', and 'para' parameters.")
                   }
                   .gauss(n, pars$lb, pars$gb, pars$para)
                 }
  )

  if (!is.null(pars$plot) && isTRUE(pars$plot)) {
    plot(afun, type = "l", main = paste("Apodisation function:", pars$fun))
  }

  return(afun)
}


#' @title Filter Bruker NMR Experiments
#' @description Filters experiments based on acquisition parameters and limits the number returned.
#' @param pars Data frame. Parsed parameter data from acquisition files.
#' @param exp_type Named list. Filtering conditions for acquisition parameters.
#' @param f_list List. File paths of spectra and metadata (e.g., f_fid, f_1r).
#' @param n_max Integer. Maximum number of experiments to retain (mainly for debugging).
#' @return A list with filtered and sorted `f_list` and `pars`.
#' @keywords internal
#' @seealso \code{\link{.extract_acq_pars1d}}, \code{\link{.check1d_files_fid}}
#' @importFrom plyr ddply .
.filterExp_files <- function(pars, exp_type, f_list, n_max) {
  idx <- match(toupper(names(exp_type)), toupper(gsub("[ap]_", "", colnames(pars))))
  if (length(idx) == 0) {
    stop("No parameter(s) found that match the specification. Check 'exp_type' and parameter choices in 'acqus' and 'procs'.")
  }

  idx_na <- which(is.na(idx))
  if (length(idx_na) > 0) {
    if (length(idx_na) == length(idx)) {
      stop("No matching parameter names found. Check input argument 'exp_type'.")
    } else {
      message(sprintf("Experiment filter %s not in NMR acquisition list. Using remaining arguments to filter: %s",
                      paste(names(exp_type)[idx_na], collapse = ", "),
                      paste(names(exp_type)[-idx_na], collapse = ", ")))
      idx <- idx[-idx_na]
    }
  }

  fmat <- vapply(seq_along(idx), function(i) {
    vars <- gsub("^<|>$", "", pars[, idx[i]])
    vars %in% exp_type[[i]]
  }, logical(nrow(pars)))

  idx_filt <- apply(fmat == 1, 1, all)
  if (!any(idx_filt)) {
    stop("No files found that match the specified parameter specification levels.")
  }

  f_list <- lapply(f_list, function(x) x[idx_filt])
  idx_filt <- which(idx_filt)

  if (length(idx_filt) > n_max) {
    idx_filt <- idx_filt[seq_len(n_max)]
    f_list <- lapply(f_list, function(x) x[seq_len(n_max)])
  }

  pars <- pars[idx_filt, ]
  fnam <- strsplit(ifelse(!is.na(f_list$f_1r), f_list$f_1r, f_list$f_fid), .Platform$file.sep)
  fmat_comb <- do.call(rbind, fnam)
  uniq_cols <- vapply(seq_len(ncol(fmat_comb)), function(i) length(unique(fmat_comb[, i])) > 1, logical(1))
  idx_keep <- which(uniq_cols)

  if (length(idx_keep) > 1) {
    rack_ <- t(vapply(seq_along(fnam), function(i) {
      c(fnam[[i]][idx_keep[length(idx_keep) - 1]], pars$a_DATE[i])
    }, FUN.VALUE = c("", "")))

    colnames(rack_) <- c("a", "b")
    rack_order_ <- ddply(as.data.frame(rack_), .(a), function(x) mean(as.POSIXct(x$b)))
    rord_fac <- order(rack_order_$V1) * 1e5

    fnam1 <- vapply(fnam, function(x) x[idx_keep[length(idx_keep) - 1]], FUN.VALUE = "")
    rord_fac <- rord_fac[match(fnam1, rack_order_$a)]

    exp_ <- vapply(fnam, function(x) x[idx_keep[length(idx_keep)]], FUN.VALUE = "")
    if (any(is.na(as.numeric(exp_)))) {
      exp_ <- factor(exp_)
    }
    rr <- order(as.numeric(exp_) + rord_fac)
  } else {
    exp_ <- vapply(fnam, function(x) x[idx_keep], FUN.VALUE = "")
    if (any(is.na(as.numeric(exp_)))) {
      exp_ <- factor(exp_)
    }
    rr <- order(as.numeric(exp_))
  }

  pars <- pars[rr, ]
  f_list <- lapply(f_list, function(x) x[rr])
  pars$a_DATE <- as.POSIXct(pars$a_DATE)

  return(list(f_list, pars))
}



#' @title Extract Bruker NMR Acquisition Parameters
#' @description Parses `acqus` files to extract acquisition parameters for each experiment.
#' @param f_list List. Contains paths to NMR experiment folders (e.g., as returned by \code{.check1d_files_fid()}).
#' @return A data frame containing extracted acquisition parameters for each experiment.
#' @keywords internal
#' @seealso \code{\link{.filterExp_files}}, \code{\link{.check1d_files_fid}}
.extract_acq_pars1d <- function(f_list) {
  out <- lapply(f_list[[1]], function(fil) {
    f_acqu_path <- file.path(fil, "acqus")
    fhand <- file(f_acqu_path, open = "r")
    f_acqu <- readLines(fhand, warn = FALSE)
    close(fhand)

    idx <- grep("..", f_acqu, fixed = TRUE)
    f_acqu[idx] <- vapply(idx, function(i) {
      gsub(" .*", f_acqu[i + 1], f_acqu[i])
    }, FUN.VALUE = "")

    parsed <- strsplit(gsub("^##\\$", "", grep("^##\\$", f_acqu, value = TRUE)), "=")
    d_acqu_val <- gsub("^ ", "", vapply(parsed, `[[`, 2, FUN.VALUE = ""))
    names(d_acqu_val) <- paste0("a_", vapply(parsed, `[[`, 1, FUN.VALUE = ""))

    idx_date <- grep("date", names(d_acqu_val), ignore.case = TRUE)
    d_acqu_val[idx_date] <- as.character(
      as.POSIXct("1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S") +
        as.numeric(d_acqu_val[idx_date])
    )

    return(d_acqu_val)
  })

  if (is.list(out)) {
    le <- unique(vapply(out, length, FUN.VALUE = integer(1)))
    if (length(le) == 1) {
      out <- do.call(rbind, out)
    } else {
      unam <- unique(names(unlist(out)))
      out <- do.call(rbind, lapply(out, `[`, unam))
      colnames(out) <- unam
    }
  }

  if (nrow(out) != length(f_list[[1]])) {
    out <- t(out)
  }

  dtype_num <- vapply(out, function(x) !any(is.na(as.numeric(x))), logical(1))
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out[, dtype_num] <- lapply(out[, dtype_num, drop = FALSE], as.numeric)
  rownames(out) <- f_list[[2]]

  return(out)
}


#' @title Check for Intact Bruker NMR File Structures
#' @description Verifies that each NMR experiment has a matching `acqus` and `fid` file.
#' @param datapath Character. Directory containing NMR spectra.
#' @param procs_exp Integer. Processing experiment number (unused, reserved for future use).
#' @param n_max Integer. Maximum number of experiments to process.
#' @param filter Logical. Whether to filter out incomplete experiments.
#' @param recursive Logical. Whether to search directories recursively.
#' @param verbose Integer. Controls verbosity level (0 = silent, 1 = basic, 2 = detailed).
#' @return A list with paths to intact experiments and matching metadata.
#' @keywords internal
#' @seealso \code{\link{.extract_acq_pars1d}}, \code{\link{.filterExp_files}}
.check1d_files_fid <- function(datapath, n_max = 10, filter = TRUE, recursive, verbose) {
  datapath <- gsub(paste0(.Platform$file.sep, "$"), "", datapath)

  f_acqus <- list.files(datapath, pattern = "^acqus$", full.names = TRUE, recursive = recursive, ignore.case = TRUE)
  f_fid <- list.files(datapath, pattern = "^fid$", full.names = TRUE, recursive = recursive, ignore.case = TRUE)

  id_a <- gsub(paste("^", datapath, "/|/acqus$", sep = ""), "", f_acqus)
  id_fid <- gsub(paste("^", datapath, "/|/fid$", sep = ""), "", f_fid)

  idx_a <- which(id_a %in% id_fid)
  idx_fid <- which(id_fid %in% id_a)

  if (length(idx_a) != length(id_a) || length(idx_fid) != length(id_fid)) {
    if (filter) {
      if (verbose > 1) {
        message("Incomplete files detected - reading only experiments with intact folder structure.")
      }
      f_acqus <- f_acqus[idx_a]
      id_a <- id_a[idx_a]
      f_fid <- f_fid[idx_fid]
      id_fid <- id_fid[idx_fid]
    } else {
      message("File system appears to be corrupt for some experiments. Consider setting `filter = TRUE`.")
      return(NULL)
    }
  }

  idm_fid <- match(id_a, id_fid)
  if (any(is.na(idm_fid)) || any(diff(idm_fid) > 1)) {
    stop("Check matching of experiment folders.")
  }

  if (length(unique(c(length(f_acqus), length(f_fid)))) != 1) {
    stop("Mismatch in number of acquisitions and FID files.")
  }

  if (n_max < length(f_fid)) {
    idx_nmax <- seq_len(n_max)
    f_acqus <- f_acqus[idx_nmax]
    f_fid <- f_fid[idx_nmax]
    id_fid <- id_fid[idx_nmax]
    message("Reached n_max - not all spectra read-in.")
  }

  p_intact <- gsub("/acqus$", "", f_acqus)
  exp_no <- id_fid

  return(list(path = p_intact, exp_no = exp_no, f_acqus = f_acqus, f_fid = f_fid))
}

