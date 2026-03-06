#### ppm

test_that("chemical shift axis is computed correctly", {
  ppm <- .chemShift(swidth = 1000, offset = 10, si = 5)
  expected <- seq(10, 10 - 1000, by = -1000 / (5 - 1))
  expect_equal(ppm, expected, tolerance = 1e-12)
  expect_length(ppm, 5)
})


#### parsing metadata
test_that(".extract_pars1d correctly parses metadata", {
  temp_dir <- tempfile("10_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  acqus_path <- file.path(temp_dir, "acqus")
  procs_path <- file.path(temp_dir, "pdata", "1", "procs")
  dir.create(dirname(procs_path), recursive = TRUE, showWarnings = FALSE)

  f_1r_path <- file.path(temp_dir, "pdata", "1", "1r")
  dir.create(dirname(f_1r_path), recursive = TRUE, showWarnings = FALSE)
  writeLines("dummy data", f_1r_path)

  writeLines(c("##$DATE= 1640995200", "##$EXP= test"), acqus_path)
  writeLines(c("##$SI= 1024", "##$SW= 5000", "##$OFFSET= 4.7"), procs_path)

  expect_true(dir.exists(dirname(procs_path)))
  expect_true(dir.exists(dirname(f_1r_path)))
  expect_true(file.exists(acqus_path))
  expect_true(file.exists(procs_path))
  expect_true(file.exists(f_1r_path))

  f_list <- list(
    f_acqus = acqus_path,
    f_procs = procs_path,
    f_1r = f_1r_path,
    path = temp_dir,
    exp_no = "10"
  )

  out <- .extract_pars1d(f_list)

  expect_s3_class(out, "data.frame")
  expect_true("a_DATE" %in% names(out))
  expect_true("p_SI" %in% names(out))
  expect_equal(out$p_SI, 1024)
})


### read1d_procs
test_that("read1d_proc returns named list with X, ppm meta", {
  skip_if_not_installed("testthat")
  skip_if_not(file.exists(system.file("extdata", package = "metabom8")),
              "Example data not installed.")

  path <- system.file("extdata", package = "metabom8")

  dat <- read1d_proc(
    path,
    exp_type = list(pulprog = 'noesygppr1d'),
    n_max = 1,
    verbose = 1
    )

  expect_true(is.list(dat))
  expect_named(dat, c("X","ppm","meta"), ignore.order = TRUE)
  expect_s3_class(dat$meta, "data.frame")
  expect_type(dat$ppm, "double")
  expect_true(is.matrix(dat$X))
  expect_equal(ncol(dat$X), length(dat$ppm))
})



test_that("read1d_proc returns expected global objects", {
  skip_if_not(dir.exists(system.file("extdata", package = "metabom8")),
              "extdata not found / installed.")
  on.exit(rm(list=c("X","ppm","meta"), envir=.GlobalEnv), add=TRUE)

  path <- system.file("extdata", package = "metabom8")

  dat <- read1d_proc(
    path,
    exp_type = list(pulprog = 'noesygppr1d'),
    n_max = 1,
    verbose = 1,
    to_global = TRUE
    )

  expect_true(exists("X", envir = .GlobalEnv))
  expect_true(exists("ppm", envir = .GlobalEnv))
  expect_true(exists("meta", envir = .GlobalEnv))
  expect_s3_class(get("meta", envir = .GlobalEnv), "data.frame")
  expect_type(get("ppm", envir = .GlobalEnv), "double")
  expect_type(get("X", envir = .GlobalEnv), "double")

})


### read1d_raw



test_that("read1d_raw returns named list with X, ppm meta", {
  skip_if_not(dir.exists(system.file("extdata", package = "metabom8")),
              "extdata not found / installed.")

  path <- system.file("extdata", package = "metabom8")

  dat <- read1d_raw(
    path = path,
    exp_type = list(exp = "PROF_PLASMA_NOESY"),
    apodisation = list(fun = "exponential", lb = 0.2),
    zerofill = 1L,
    mode = "absorption",
    verbose = 0,
    n_max = 2,
    recursive = TRUE,
    filter = TRUE
  )

  expect_true(is.list(dat))
  expect_named(dat, c("X","ppm","meta"), ignore.order = TRUE)
  expect_s3_class(dat$meta, "data.frame")
  expect_type(dat$ppm, "double")
  expect_true(is.matrix(dat$X))

  expect_equal(ncol(dat$X), length(dat$ppm))
  expect_equal(nrow(dat$X), nrow(dat$meta))
})



test_that("read1d_raw returns expected global objects", {
  skip_if_not(dir.exists(system.file("extdata", package = "metabom8")),
              "extdata not found / installed.")
  on.exit(rm(list=c("X","ppm","meta"), envir=.GlobalEnv), add=TRUE)

  path <- system.file("extdata", package = "metabom8")

  read1d_raw(
    path = path,
    exp_type = list(exp = "PROF_PLASMA_NOESY"),
    apodisation = list(fun = "exponential", lb = 0.2),
    zerofill = 1L,
    mode = "absorption",
    verbose = 0,
    n_max = 2,
    recursive = TRUE,
    filter = TRUE,
    to_global = TRUE
  )


  expect_true(exists("X", envir = .GlobalEnv))
  expect_true(exists("ppm", envir = .GlobalEnv))
  expect_true(exists("meta", envir = .GlobalEnv))

  expect_s3_class(get("meta", envir = .GlobalEnv), "data.frame")
  expect_type(get("ppm", envir = .GlobalEnv), "double")
  expect_true(is.matrix(get("X", envir = .GlobalEnv)))

  expect_equal(ncol(get("X", envir = .GlobalEnv)),
               length(get("ppm", envir = .GlobalEnv)))
  expect_equal(nrow(get("X", envir = .GlobalEnv)),
               nrow(get("meta", envir = .GlobalEnv)))

})
