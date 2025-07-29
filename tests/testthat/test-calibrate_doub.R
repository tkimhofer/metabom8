test_that(".calibrate_doub returns correctly shaped matrix", {
  set.seed(123)
  ppm <- seq(0.5, 6, length.out = 1024)
  X <- matrix(rnorm(1024 * 10), nrow = 10)

  # add doublet to row 1 for glucose region
  X[1, which.min(abs(ppm - 5.226))] <- 1
  X[1, which.min(abs(ppm - 5.240))] <- 1

  result <- metabom8:::.calibrate_doub(X, ppm, type = "glu")
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(X))
})

test_that(".calibrate_doub works with 'ala' type", {
  set.seed(42)
  ppm <- seq(0.5, 2, length.out = 1024)
  X <- matrix(rnorm(1024 * 5), nrow = 5)
  X[2, which.min(abs(ppm - 1.475))] <- 1
  X[2, which.min(abs(ppm - 1.49))] <- 1

  result <- metabom8:::.calibrate_doub(X, ppm, type = "ala")
  expect_true(is.matrix(result))
})

test_that(".calibrate_doub works with numeric range", {
  ppm <- seq(0.5, 6, length.out = 512)
  X <- matrix(rnorm(512 * 3), nrow = 3)
  X[1, which.min(abs(ppm - 4.99))] <- 1
  X[1, which.min(abs(ppm - 5.005))] <- 1

  result <- metabom8:::.calibrate_doub(X, ppm, type = c(4.95, 5.05), j_const = c(0.01, 0.02))
  expect_equal(dim(result), dim(X))
})

test_that(".calibrate_doub warns when no peak is found", {
  ppm <- seq(0.5, 6, length.out = 1024)
  X <- matrix(0, nrow = 2, ncol = 1024)

  expect_warning({
    result <- metabom8:::.calibrate_doub(X, ppm, type = "glu")
  }, regexp = "Could not calibrate spectrum")
})
