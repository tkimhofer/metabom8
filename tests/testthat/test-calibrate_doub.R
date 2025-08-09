test_that(".calibrate_doub returns correctly shaped matrix", {
  set.seed(123)
  n_spec <- 5
  ppm <- seq(5.21, 5.25, length.out = 1024)
  cent_loc <- 5.233
  j_const <- c(0.0065, 0.007)
  lw <- 0.0002

  gauss <- function(x, mu, sd) exp(-((x - mu)^2) / (2 * sd^2))

  X <- t(sapply(1:n_spec, function(i) {
    jc <- 0.0065
    peak1 <- gauss(ppm, cent_loc - jc/2, lw)
    peak2 <- gauss(ppm, cent_loc + jc/2, lw)
    peak1 + peak2
  }))

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

  set.seed(123)
  n_spec <- 2
  ppm <- seq(4.8, 5.1, length.out = 1024)
  cent_loc <- 4.9
  j_const <- 0.008
  lw <- 0.0002

  gauss <- function(x, mu, sd) exp(-((x - mu)^2) / (2 * sd^2))

  X <- t(sapply(1:n_spec, function(i) {
    peak1 <- gauss(ppm, cent_loc - j_const/2, lw)
    peak2 <- gauss(ppm, cent_loc + j_const/2, lw)
    peak1 + peak2
  }))

  # matspec(X, ppm)

  result <- metabom8:::.calibrate_doub(X, ppm, type = c(4.8, 5.1), j_const = c(0.005, 0.01))
  expect_equal(dim(result), dim(X))
})

test_that(".calibrate_doub warns when no peak is found", {
  ppm <- seq(0.5, 6, length.out = 1024)
  X <- matrix(0, nrow = 2, ncol = 1024)

  expect_warning({
    result <- metabom8:::.calibrate_doub(X, ppm, type = "glu")
  }, regexp = "Could not calibrate")
})
