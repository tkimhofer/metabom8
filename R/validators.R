# collection of validators

.arg_check_scaling <- function(x){
  stopifnot(exprs = {
    is.list(x)
    !is.null(names(x))
    length(x) == 2
    all(names(x) %in% c('center', 'scale'))
    is.logical(x$center) && length(x$center) == 1
    is.character(x$scale) && length(x$scale) == 1
    tolower(x$scale) %in% c('none', 'uv', 'pareto')
  })
}

.check_preprocess <- function(x){
  stopifnot(exprs = {
    is.list(x)
    !is.null(names(x))
    length(x) == 4
    all(names(x) %in% c('center', 'scale', 'X_mean', 'X_sd'))
    length(X_mean) == length(X_sd)
    all(!is.infinite(x$X_mean) && !is.null(x$X_mean))
    all(!is.infinite(x$X_sd) && !is.null(x$X_sd))
  })
}

#' @importFrom stats runif
.check_probs <- function(cv_pars) {
  if (!"probs" %in% names(cv_pars)) {
    stop("quantile probs argument missing", call. = FALSE)
  }

  probs <- cv_pars$probs

  if (!is.numeric(probs)) {
    stop("probs must be numeric", call. = FALSE)
  }

  if (length(probs) < 2L) {
    stop("at least two quantile probs must be provided", call. = FALSE)
  }

  if (any(probs < 0 | probs > 1)) {
    stop("quantile probs need to be in [0,1]", call. = FALSE)
  }

  if (anyDuplicated(probs)) {
    stop("probs must be strictly increasing (no duplicates)", call. = FALSE)
  }

  if (!isTRUE(all.equal(sort(probs), probs))) {
    stop("probs must be sorted in increasing order", call. = FALSE)
  }

  if (probs[1] != 0 || tail(probs, 1) != 1) {
    stop("probs must start at 0 and end at 1", call. = FALSE)
  }
}

.check_split <- function(cv_pars) {
  if (!"split" %in% names(cv_pars)) {
    stop("split argument missing", call. = FALSE)
  }

  split <- cv_pars$split

  if (!is.numeric(split) || length(split) != 1L || is.na(split)) {
    stop("split must be a single numeric value", call. = FALSE)
  }

  if (!(split > 0 && split < 1)) {
    stop("split must be strictly between 0 and 1", call. = FALSE)
  }

  if (split < 0.05 || split > 0.95) {
    warning("Very small or large split may lead to unstable estimates")
  }
}

.arg_check_cv <- function(cv_pars, model_type=c('DA', 'R'), n, Y_prepped){

  if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)

  seed <- .Random.seed

  type <- match.arg(model_type, c("DA", "R"))

  k <- cv_pars$k

  stopifnot(
    is.list(cv_pars),
    !is.null(names(cv_pars)),
    all(c('method', 'k') %in% names(cv_pars)),
    length(k) == 1L && !is.na(k) && is.numeric(k) && k %% 1 == 0 && k >= 2L
  )

  out <- switch(
    cv_pars$method,

    "kfold" = {
      stopifnot("k must be <= number of samples" = k <= n)

      idc <- .cvSetsMethod(Y=Y_prepped, type, method = "k-fold", k = k)
      list(
        train    = idc,
        strategy = "KFold",
        n        = n,
        seed     = seed
      )
    },

    "StratifiedKFold" = {
      if (type == "R") {
        .check_probs(cv_pars)
        probs <- cv_pars$probs
      } else{ probs <- NULL}

      idc <- .cvSetsMethod(Y=Y_prepped, type, method = "k-fold_stratified", k = k, probs=probs)
      list(
        train    = idc,
        strategy = "StratifiedKFold",
        n        = n,
        seed     = seed
      )
    },

    "MonteCarlo" = {
      .check_split(cv_pars)
      idc <- .cvSetsMethod(Y_prepped, type, method = "MC", k = k, split = cv_pars$split)
      list(
        train    = idc,
        strategy = "MonteCarlo",
        n        = n,
        seed     = seed
      )
    },

    "BalancedMonteCarlo" = {
      if (type == "R") {
        .check_probs(cv_pars)
        probs <- cv_pars$probs
      } else{ probs <- NULL}
      .check_split(cv_pars)

      idc <- .cvSetsMethod(Y=Y_prepped, type=type, method = "MC_balanced", k = k, split = cv_pars$split, probs=probs)
      list(
        train    = idc,
        strategy = "BalancedMonteCarlo",
        n        = n,
        seed     = seed
      )
    },

    "BalancedBoot" = {
      if (type == "R") {
        .check_probs(cv_pars)
        probs <- cv_pars$probs
      } else{ probs <- NULL}
      .check_split(cv_pars)

      idc <- .cvSetsMethod(Y=Y_prepped, type=type, method = "Boot_balanced", k = k, split = cv_pars$split, probs=probs)
      list(
        train    = idc,
        strategy = "Boot_balanced",
        n        = n,
        seed     = seed
      )
    }

  )

  out
}


