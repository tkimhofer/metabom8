#########################
### MODEL VALIDATION ###
#########################

setValidity("m8_model", function(object) {
  msg <- character()

  if (length(object@engine) != 1L || is.na(object@engine) || object@engine == "") {
    msg <- c(msg, "slot 'engine' must be a single non-empty character string.")
  } else if (!object@engine %in% c("pca", "pls", "opls")) {
    msg <- c(msg, "slot 'engine' must be one of: 'pca', 'pls', 'opls'.")
  }

  if (!is.list(object@dims) || is.null(object@dims$n) || is.null(object@dims$p)) {
    msg <- c(msg, "slot 'dims' must be a list containing elements 'n' and 'p'.")
  } else {
    n <- object@dims$n
    p <- object@dims$p
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0) {
      msg <- c(msg, "dims$n must be a single positive numeric value.")
    }
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p <= 0) {
      msg <- c(msg, "dims$p must be a single positive numeric value.")
    }
  }

  if (!is.list(object@ctrl)) {
    msg <- c(msg, "slot 'ctrl' must be a list.")
  } else {
    ctrl <- object@ctrl
    if (!is.null(ctrl$type)) {
      if (!is.character(ctrl$type) || length(ctrl$type) != 1L || is.na(ctrl$type)) {
        msg <- c(msg, "ctrl$type must be a single character string.")
      } else if (!ctrl$type %in% c("R", "DA", "DA-mY", "unsupervised")) {
        msg <- c(msg, "ctrl$type must be one of: 'R', 'DA', 'unsupervised'.")
      }
    } else {
      msg <- c(msg, "ctrl$type is missing for engine.")
    }

    n_sel <- ctrl$ncomp_selected
    n_tst <- ctrl$ncomp_tested

    if (!is.null(n_sel) && (!is.numeric(n_sel) || length(n_sel) != 1L || is.na(n_sel) || n_sel < 0)) {
      msg <- c(msg, "ctrl$ncomp_selected must be a single non-negative numeric value.")
    }
    if (!is.null(n_tst) && (!is.numeric(n_tst) || length(n_tst) != 1L || is.na(n_tst) || n_tst < 0)) {
      msg <- c(msg, "ctrl$ncomp_tested must be a single non-negative numeric value.")
    }

    if (!is.null(n_sel) && !is.null(n_tst) &&
        is.numeric(n_sel) && is.numeric(n_tst) &&
        !is.na(n_sel) && !is.na(n_tst)) {
      if (n_sel > n_tst) {
        msg <- c(msg, "ctrl$ncomp_selected cannot exceed ctrl$ncomp_tested.")
      }
      if (object@engine %in% c("pls", "opls") && n_tst < 1) {
        msg <- c(msg, "For PLS/OPLS models, ctrl$ncomp_tested should be >= 1.")
      }
    }

    check_len <- function(x, name, n_req, offset=0) {
      if (!is.null(x) && !is.null(n_req) && is.numeric(n_req) && !is.na(n_req) && n_req > 0) {
        if (!is.numeric(x) && !is.logical(x))
          return(paste0(name, " must be numeric (or logical NA)."))
        expected <- n_req - offset
        if (length(x) < expected)
          return(paste0(name, " has length ", length(x),
                        " but expected at least ", expected,
                        " based on ncomp_tested = ", n_req, "."))
      }
      NULL
    }

    if (!is.null(n_tst) && is.numeric(n_tst) && !is.na(n_tst) && n_tst > 0) {
      offset <- if (object@engine == "opls") 1 else 0
      m <- c(
        check_len(ctrl$q2, "ctrl$q2", n_tst, offset),
        check_len(ctrl$r2, "ctrl$r2", n_tst, offset),
        check_len(ctrl$aucs_te, "ctrl$aucs_te", n_tst, offset),
        check_len(ctrl$aucs_tr, "ctrl$aucs_tr", n_tst, offset)
      )
      m <- m[!vapply(m, is.null, logical(1))]
      if (length(m)) msg <- c(msg, m)
    }

    if (identical(object@engine, "opls")) {
      if (!is.null(ctrl$nc_pred) && (!is.numeric(ctrl$nc_pred) || length(ctrl$nc_pred) != 1L || is.na(ctrl$nc_pred) || ctrl$nc_pred > 1)) {
        msg <- c(msg, "For OPLS, ctrl$nc_pred must be a single numeric value >= 1.")
      }
      if (!is.null(ctrl$nc_orth) && (!is.numeric(ctrl$nc_orth) || length(ctrl$nc_orth) != 1L || is.na(ctrl$nc_orth) || ctrl$nc_orth < 0)) {
        msg <- c(msg, "For OPLS, ctrl$nc_orth must be a single non-negative numeric value.")
      }
      if (!is.null(ctrl$nc_pred) && !is.null(ctrl$nc_orth) && !is.null(ctrl$ncomp_selected) &&
          is.numeric(ctrl$nc_pred) && is.numeric(ctrl$nc_orth) && is.numeric(ctrl$ncomp_selected) &&
          !is.na(ctrl$nc_pred) && !is.na(ctrl$nc_orth) && !is.na(ctrl$ncomp_selected)) {
        if ((ctrl$nc_pred + ctrl$nc_orth) != ctrl$ncomp_selected) {
          msg <- c(msg, "For OPLS, ctrl$ncomp_selected must equal ctrl$nc_pred + ctrl$nc_orth.")
        }
      }
    }
  }

  if (!is.list(object@fit)) msg <- c(msg, "slot 'fit' must be a list.")
  if (!is.list(object@provenance)) msg <- c(msg, "slot 'provenance' must be a list.")
  if (!is.list(object@session)) msg <- c(msg, "slot 'session' must be a list.")

  if (!is.list(object@cv) && !is.null(object@cv)) {
    msg <- c(msg, "slot 'cv' must be a list'.")
  }
  if (length(msg)) msg else TRUE
})
