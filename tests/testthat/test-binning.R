test_that("binning with width returns matrix of expected dimensions", {
  set.seed(1)
  X <- matrix(rnorm(1000), nrow = 10)
  ppm <- seq(10, 0, length.out = ncol(X))
  Xb <- binning(X, ppm, width = 0.1)

  expect_true(is.matrix(Xb))
  expect_equal(nrow(Xb), nrow(X))
  expect_lt(ncol(Xb), ncol(X))
  expect_false(any(is.na(Xb)))
})

test_that("binning with npoints returns matrix of expected shape", {
  set.seed(2)
  X <- matrix(rnorm(1000), nrow = 10)
  ppm <- seq(10, 0, length.out = ncol(X))
  Xb <- binning(X, ppm, npoints = 50)

  expect_true(is.matrix(Xb))
  expect_equal(nrow(Xb), nrow(X))
  expect_equal(ncol(Xb), 50)
  expect_false(any(is.na(Xb)))
})

test_that("binning throws error if both width and npoints are specified", {
  X <- matrix(rnorm(100), nrow = 10)
  ppm <- seq(10, 0, length.out = ncol(X))

  expect_error(binning(X, ppm, width = 0.1, npoints = 20),
               "Please specify one or the other")
})

test_that("binning throws error if neither width nor npoints specified", {
  X <- matrix(rnorm(100), nrow = 10)
  ppm <- seq(10, 0, length.out = ncol(X))

  expect_error(binning(X, ppm),
               "Define bin width in ppm or desired number of bins")
})

test_that("binning throws error if width too small", {
  X <- matrix(rnorm(100), nrow = 10)
  ppm <- seq(10, 0, length.out = ncol(X))

  expect_error(binning(X, ppm, width = 0.00001),
               "Bin width equals or is smaller")
})

test_that("binning throws error if npoints is too high", {
  X <- matrix(rnorm(100), nrow = 10)
  ppm <- seq(10, 0, length.out = ncol(X))

  expect_error(binning(X, ppm, npoints = 200),
               "npoints cannot be larger or equal")
})
