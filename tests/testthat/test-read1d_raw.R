test_that("read1d_raw correctly reads and processes FID data", {
  # Path to example data included in the package
  test_path <- system.file("extdata", package = "metabom8")

  # Run the function on minimal subset
  read1d_raw(
    path = test_path,
    exp_type = list(exp = "PROF_PLASMA_NOESY"),
    apodisation = list(fun = "exponential", lb = 0.2),
    zerofil = 1,
    return = "absorption",
    verbose = 0,
    n_max = 2,
    recursive = TRUE,
    filter = TRUE
  )

  # Check if objects exist
  expect_true(exists("X", envir = .GlobalEnv))
  expect_true(exists("ppm", envir = .GlobalEnv))
  expect_true(exists("meta", envir = .GlobalEnv))

  # Check object types
  expect_true(is.matrix(X))
  expect_true(is.numeric(ppm))
  expect_true(is.data.frame(meta))

  # Check dimensions match
  expect_equal(nrow(X), nrow(meta))
  expect_equal(ncol(X), length(ppm))

  # Clean up
  rm(X, ppm, meta, envir = .GlobalEnv)
})
