test_that(".chemShift returns correct ppm axis", {
  ppm <- .chemShift(swidth = 1000, offset = 10, si = 5)
  expected <- seq(10, 10 - 1000, by = -1000 / (5 - 1))
  expect_equal(ppm, expected)
  expect_length(ppm, 5)
})

test_that(".extract_pars1d correctly parses metadata", {
  # Mock f_list input
  temp_dir <- tempdir()
  acqus_path <- file.path(temp_dir, "acqus")
  procs_path <- file.path(temp_dir, "pdata", "1", "procs")
  dir.create(dirname(procs_path), recursive = TRUE)

  # Simulate minimal acqus and procs content
  writeLines(c("##$DATE= 1640995200", "##$EXP= test"), acqus_path)
  writeLines(c("##$SI= 1024", "##$SW= 5000", "##$OFFSET= 4.7"), procs_path)

  f_list <- list(
    f_acqus = c(acqus_path),
    f_procs = c(procs_path),
    f_1r = c(file.path(temp_dir, "pdata", "1", "1r")),
    path = temp_dir,
    exp_no = "exp1"
  )

  out <- .extract_pars1d(f_list)
  expect_s3_class(out, "data.frame")
  expect_equal(rownames(out), "exp1")
  expect_true(any(grepl("^a_DATE$", names(out))))
  expect_true(any(grepl("^p_SI$", names(out))))
  expect_type(out$p_SI, "double")
})

test_that("read1d returns expected global objects", {
  skip_if_not_installed("testthat")
  skip_if_not(file.exists(system.file("extdata", package = "metabom8")),
              "Example data not installed.")

  path <- system.file("extdata", package = "metabom8")
  suppressWarnings(read1d(path, exp_type = list(), n_max = 1, verbose = FALSE))

  expect_true(exists("X", envir = .GlobalEnv))
  expect_true(exists("ppm", envir = .GlobalEnv))
  expect_true(exists("meta", envir = .GlobalEnv))
  expect_s3_class(get("meta", envir = .GlobalEnv), "data.frame")
  expect_type(get("ppm", envir = .GlobalEnv), "double")
  expect_type(get("X", envir = .GlobalEnv), "double")
})
