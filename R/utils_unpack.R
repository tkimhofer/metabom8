.m8_unpack_dat <- function(x, ppm = NULL, meta = NULL) {
  if (is.list(x) && !is.null(x$X)) {
    if (is.null(ppm))  ppm  <- x$ppm
    if (is.null(meta)) meta <- x$meta
    x <- x$X
    if (!.check_X_ppm(x, ppm))
      stop("Dimensions of X and ppm do not match.", call. = FALSE)
  }


  list(X = x, ppm = ppm, meta = meta)
}
