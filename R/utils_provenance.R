# internal helper: append preprocessing step to a matrix attribute

.m8_stamp <- function(X, step, params = list(), notes = NULL) {
  rec <- list(
    step = step,
    params = params,
    notes = notes,
    time = as.character(Sys.time()),
    pkg = "metabom8",
    pkg_version = as.character(utils::packageVersion("metabom8"))
  )

  log <- attr(X, "m8_prep", exact = TRUE)
  if (is.null(log)) log <- list()
  log[[length(log) + 1L]] <- rec
  attr(X, "m8_prep") <- log
  X
}

# internal helper: safely carry attributes from old X to new X
.m8_copy_attrs <- function(from, to, keys = c("m8_prep", "m8_axis")) {
  for (k in keys) {
    val <- attr(from, k, exact = TRUE)
    if (!is.null(val)) attr(to, k) <- val
  }
  to
}

.format_params <- function(x, indent = 8, detail = FALSE, max_items = 8) {

  if (is.null(x) || length(x) == 0)
    return(character())

  space <- strrep(" ", indent)

  if (!is.list(x)) {
    return(paste0(space, .format_value(x)))
  }

  keys <- names(x)

  if (is.null(keys)) {
    return(vapply(seq_along(x), function(i) {
      paste0(space, "- ", .format_value(x[[i]]))
    }, character(1)))
  }

  max_width <- max(nchar(keys))

  out <- character()

  for (k in keys) {

    value <- x[[k]]

    if (!detail && .should_collapse(value, max_items = max_items)) {
      preview <- .preview_list(value, n = 3)
      out <- c(out, paste0(
        space,
        format(k, width = max_width, justify = "left"),
        " : ",
        sprintf("<%d items>", length(value)),
        if (nzchar(preview)) paste0(" [", preview, "]") else ""
      ))
      next
    }

    if (is.list(value)) {
      out <- c(out, paste0(
        space,
        format(k, width = max_width, justify = "left"),
        " :"
      ))
      out <- c(out, .format_params(value,
                                   indent = indent + 4,
                                   detail = detail,
                                   max_items = max_items
      ))
    } else {
      out <- c(out, paste0(
        space,
        format(k, width = max_width, justify = "left"),
        " : ",
        .format_value(value)
      ))
    }
  }

  out
}

.format_value <- function(v) {

  if (is.numeric(v) && length(v) > 1) {
    return(sprintf("[%s]",
                   paste(signif(v, 6), collapse = ", ")))
  }

  if (length(v) > 1) {
    return(sprintf("[%s]", paste(v, collapse = ", ")))
  }

  as.character(v)
}

#' @importFrom utils capture.output
.should_collapse <- function(x,
                             max_items = 8,
                             max_chars = 1200) {

  if (!is.list(x))
    return(FALSE)

  if (length(x) > max_items)
    return(TRUE)

  est_size <- nchar(paste(utils::capture.output(utils::str(x, max.level = 2)), collapse = ""))
  if (est_size > max_chars)
    return(TRUE)

  FALSE
}

.preview_list <- function(x, n = 3) {

  if (!is.list(x) || length(x) == 0)
    return("")

  x_sub <- x[seq_len(min(n, length(x)))]
  items <- vapply(x_sub, .format_value, character(1))
  preview <- paste(items, collapse = ", ")

  if (length(x) > n)
    preview <- paste0(preview, ", ...")

  preview
}


get_attr <- function(X, field) {
  attributes(X)[[field]]
}


#' Retrieve metabom8 provenance metadata
#'
#' Extracts preprocessing provenance stored in the \code{"m8_prep"} attribute.
#' Allows access to the full processing log, a specific step, or a specific
#' parameter within a step.
#'
#' @param x A metabom8 object (named list with element \code{X}) or a numeric
#'   matrix containing metabom8 provenance metadata.
#' @param step Optional. Either:
#'   \itemize{
#'     \item Numeric index of the preprocessing step
#'     \item Character string matching the recorded step name
#'   }
#'   If \code{NULL}, the full provenance log is returned.
#' @param param Optional character string specifying a parameter name within
#'   the selected step. If provided, only this parameter value is returned.
#'
#' @return
#' Depending on the arguments:
#' \itemize{
#'   \item Full provenance list (if \code{step = NULL})
#'   \item A single preprocessing step (if \code{step} specified)
#'   \item A single parameter value (if \code{param} specified)
#' }
#'
#' @details
#' Provenance metadata are recorded automatically by metabom8 preprocessing
#' functions and stored as structured attributes on the spectral matrix.
#' This function provides programmatic access to these records.
#'
#' @examples
#' data(hiit_raw)
#'
#' hiit_proc <- hiit_raw |>
#'   calibrate(type = "tsp") |>
#'   excise()
#'
#' # Retrieve full log
#' log <- get_provenance(hiit_proc)
#'
#' # Retrieve specific step
#' get_provenance(hiit_proc, step = 2)
#'
#' # Retrieve parameter from a named step
#' get_provenance(hiit_proc, step = "calibrate", param = "target")
#'
#' @family provenance
#' @export
get_provenance <- function(x, step = NULL, param = NULL) {

  obj  <- x$X %||% x
  prep <- attr(obj, "m8_prep")

  if (is.null(prep))
    stop("No metabom8 provenance found.", call. = FALSE)

  if (is.null(step))
    return(prep)

  if (is.numeric(step)) {
    if (step > length(prep))
      stop("Step index out of range.", call. = FALSE)

    sel <- prep[[step]]
  }

  else {
    idx <- which(vapply(prep, `[[`, character(1), "step") == step)

    if (length(idx) == 0)
      stop("Step not found.", call. = FALSE)

    sel <- prep[[idx]]
  }

  if (is.null(param))
    return(sel)

  if (!param %in% names(sel$params))
    stop("Parameter not found in selected step.", call. = FALSE)

  sel$params[[param]]
}



`%||%` <- function(a, b) if (!is.null(a)) a else b



.format_provenance <- function(x, detail = FALSE, max_items = 8) {
  prep <- attr(x$X %||% x, "m8_prep", exact = TRUE)

  if (is.null(prep) || length(prep) == 0) {
    return("No metabom8 preprocessing metadata found.")
  }

  out <- c(
    "metabom8 processing pipeline:",
    "=============================="
  )

  for (i in seq_along(prep)) {
    step <- prep[[i]]
    out <- c(out, sprintf("[%d] %s", i, step$step))

    if (!is.null(step$params)) {
      out <- c(out, "     params:")
      out <- c(out, .format_params(step$params, indent = 8, detail = detail, max_items = max_items))
    }

    if (!is.null(step$notes)) {
      out <- c(out, paste0("     notes: ", step$notes))
    }

    out <- c(out, "")
  }

  out
}


#' Print metabom8 preprocessing pipeline
#'
#' Displays the recorded preprocessing history attached to a metabom8
#' spectral matrix. The function reads the \code{"m8_prep"} attribute
#' and prints each processing step in chronological order, including
#' parameters and notes.
#'
#' The input can be either:
#' \itemize{
#'   \item A metabom8-style list containing \code{$X}, or
#'   \item A matrix with metabom8 provenance attributes attached.
#' }
#'
#' If no preprocessing metadata is found, a message is printed.
#'
#' @param x A metabom8 object or spectral matrix with attached
#'   \code{"m8_prep"} provenance metadata.
#' @param detail Prints full information trail, incl. timestamp / versioning
#' @param max_items Limits individual list entries (e.g., parameter) to specified
#' number
#'
#' @return Invisibly returns \code{NULL}. This function is called for
#'   its side effect of printing pipeline information.
#'
#' @details
#' metabom8 records preprocessing steps as an ordered list of
#' transformation descriptors stored in the \code{"m8_prep"} attribute.
#' Each step typically contains:
#' \itemize{
#'   \item \code{step}: Name of the preprocessing operation
#'   \item \code{params}: Parameters used
#'   \item \code{notes}: Optional description
#'   \item \code{time}: Timestamp
#'   \item \code{pkg}: Package name and version
#' }
#'
#' This function provides a compact audit trail of the processing
#' workflow, facilitating reproducibility and provenance inspection.
#'
#' @examples
#' params <- list(
#'   runtime  = "docker",
#'   image    = Sys.getenv("IMAGE", "docker-image-dummy"),
#'   workflow = Sys.getenv("M8_WORKFLOW", "std_prof-urine"),
#'   agent    = paste0("snakemake/", Sys.getenv("SNAKEMAKE_VERSION", "v?")),
#'   run_id   = Sys.getenv("M8_RUN_ID", "m8-2605-001")
#' )
#'
#' data(hiit_raw)
#' print_provenance(hiit_raw)
#'
#' hiit_proc <- hiit_raw |>
#'   calibrate(type = "tsp") |>
#'   excise() |>
#'   add_note('🔎 dilution-adaptive acquisition mode -> verify snr after normalisation',
#'     params)
#'
#' print_provenance(hiit_proc)
#' @family provenance
#' @export
print_provenance <- function(x, detail = FALSE, max_items = 8) {
  txt <- .format_provenance(x, detail = detail, max_items = max_items)
  writeLines(txt)
  invisible(txt)
}

#' Add user note to metabom8 provenance
#' Appends a user annotation to the \code{"m8_prep"} attribute.
#' The step title is formatted as \code{"note {username}"}.
#' The timestamp is stored in \code{params}, and the user message
#' is stored in \code{notes}.
#' @param x A metabom8 object or matrix with \code{"m8_prep"} metadata.
#' @param note Character string describing the annotation.
#' @param params Named list providing paramter key-value pairs.
#' @return The input object with updated provenance metadata.
#' @family provenance
#' @examples
#' params <- list(
#'   runtime  = "docker",
#'   image    = Sys.getenv("IMAGE", "docker-image-dummy"),
#'   workflow = Sys.getenv("M8_WORKFLOW", "std_prof-urine"),
#'   agent    = paste0("snakemake/", Sys.getenv("SNAKEMAKE_VERSION", "v?")),
#'   run_id   = Sys.getenv("M8_RUN_ID", "m8-2605-001")
#' )
#'
#' data(hiit_raw)
#' print_provenance(hiit_raw)
#'
#' hiit_proc <- hiit_raw |>
#'   calibrate(type = "tsp") |>
#'   excise() |>
#'   add_note('🔎 dilution-adaptive acquisition mode -> verify snr after normalisation',
#'     params)
#'
#' print_provenance(hiit_proc)
#'
#' @export
add_note <- function(x, note, params=NULL) {

  stopifnot(
    is.character(note), length(note) == 1
  )

  target <- if (is.list(x) && !is.null(x$X)) x$X else x

  target <- .m8_stamp(
    target,
    step   = "user note",
    params = params,
    notes  = note
  )

  if (is.list(x) && !is.null(x$X)) {
    x$X <- target
    return(x)
  }

  target
}


#' List available preprocessing steps
#'
#' Returns a named character vector describing the preprocessing
#' operations implemented in \pkg{metabom8}. For more information
#' on individual functionalities please refer to the function help
#' pages.
#'
#' @return A named character vector with preprocessing function
#'   identifiers as names and short descriptions as values.
#' @examples
#' list_preprocessing()
#' @export
list_preprocessing <- function() {
  c(
    calibrate        = "Chemical shift calibration to reference signal",
    excise           = "Remove spectral regions (e.g., residual water, urea)",
    baseline_correct = "Baseline correction",
    norm_erectic     = "Normalise spectra based on ERETIC signal",
    pqn              = "Apply Probabilistic Quantile Normalisation (PQN)",
    binning          = "Bin spectra",
    correct_lw       = "Apply line width correction (experimental)",
    align_segment    = "Spectral alignment of a single chemical shift segment",
    align_spectra    = "Automated segment-wise full spectrum alignment"
  )
}

#' List available preprocessing functions
#' Returns the preprocessing utilities provided by \pkg{metabom8}.
#' @return A named character vector describing preprocessing functions.
#' @examples
#' list_preprocessing()
#' @family preprocessing
#' @export
print_preprocessing <- function() {
  fns <- list_preprocessing()
  out <- c("Available preprocessing functions:", "")
  out <- c(out, paste0(" - ", format(names(fns), width = 18), fns))
  writeLines(out)
  invisible(fns)
}
