test_that(".chemShift returns correct ppm axis", {
  ppm <- .chemShift(swidth = 1000, offset = 10, si = 5)
  expected <- seq(10, 10 - 1000, by = -1000 / (5 - 1))
  expect_equal(ppm, expected)
  expect_length(ppm, 5)
})

test_that(".extract_pars1d correctly parses metadata", {
  # Create unique temp directory
  temp_dir <- tempfile("10_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  acqus_path <- file.path(temp_dir, "acqus")
  procs_path <- file.path(temp_dir, "pdata", "1", "procs")
  dir.create(dirname(procs_path), recursive = TRUE, showWarnings = FALSE)

  # Create 1r dummy file directory & file without warnings
  f_1r_path <- file.path(temp_dir, "pdata", "1", "1r")
  dir.create(dirname(f_1r_path), recursive = TRUE, showWarnings = FALSE)
  writeLines("dummy data", f_1r_path)

  # Write minimal acqus and procs content
  writeLines(c("##$DATE= 1640995200", "##$EXP= test"), acqus_path)
  writeLines(c("##$SI= 1024", "##$SW= 5000", "##$OFFSET= 4.7"), procs_path)

  # Sanity check directories and files exist before function call
  expect_true(dir.exists(dirname(procs_path)))
  expect_true(dir.exists(dirname(f_1r_path)))
  expect_true(file.exists(acqus_path))
  expect_true(file.exists(procs_path))
  expect_true(file.exists(f_1r_path))

  # Construct file list input
  f_list <- list(
    f_acqus = acqus_path,
    f_procs = procs_path,
    f_1r = f_1r_path,
    path = temp_dir,
    exp_no = "10"
  )

  out <- .extract_pars1d(f_list)

  expect_s3_class(out, "data.frame")
  expect_true(any(grepl("^a_DATE$", names(out))))
  expect_true(any(grepl("^p_SI$", names(out))))
  expect_type(out$p_SI, "double")
})


test_that("read1d returns expected global objects", {
  skip_if_not_installed("testthat")
  skip_if_not(file.exists(system.file("extdata", package = "metabom8")),
              "Example data not installed.")

  path <- system.file("extdata", package = "metabom8")
  suppressWarnings(read1d(path, exp_type = list(pulprog = 'noesygppr1d'), n_max = 1, verbose = 1))

  expect_true(exists("X", envir = .GlobalEnv))
  expect_true(exists("ppm", envir = .GlobalEnv))
  expect_true(exists("meta", envir = .GlobalEnv))
  expect_s3_class(get("meta", envir = .GlobalEnv), "data.frame")
  expect_type(get("ppm", envir = .GlobalEnv), "double")
  expect_type(get("X", envir = .GlobalEnv), "double")
})
