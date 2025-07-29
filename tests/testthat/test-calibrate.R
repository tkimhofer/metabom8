test_that("calibrate() aligns spectrum correctly using TSP", {
  set.seed(123)
  ppm <- seq(10, 0, length.out = 1000)
  X <- matrix(rnorm(10 * length(ppm), sd = 1e-2), nrow = 10)

  # Add TSP peak at 0.1 ppm (simulate shifted peak)
  tsp_idx <- which.min(abs(ppm - 0.1))
  X[, tsp_idx] <- X[, tsp_idx] + 10

  X_cal <- calibrate(X, ppm, type = "tsp")

  new_peak_pos <- which.max(X_cal[1, ])
  expect_true(abs(ppm[new_peak_pos] - 0) < 0.01)
})

test_that("calibrate() with numeric range recenters peak to mean of interval", {
  set.seed(1)
  ppm <- seq(10, 0, length.out = 1000)
  X <- matrix(rnorm(5 * length(ppm), sd = 0.001), nrow = 5)

  # Simulate strong peak at 2.1 ppm
  peak_idx <- which.min(abs(ppm - 2.1))
  X[, peak_idx] <- X[, peak_idx] + 10

  X_cal <- calibrate(X, ppm, type = c(2.0, 2.2))
  new_peak_pos <- which.max(X_cal[1, ])
  expect_true(abs(ppm[new_peak_pos] - mean(c(2.0, 2.2))) < 0.01)
})

test_that("calibrate() fails for incorrect input", {
  ppm <- seq(10, 0, length.out = 1000)
  X <- matrix(rnorm(1000), nrow = 1)

  expect_error(calibrate(X, ppm = NULL, type = "tsp"))
  expect_error(calibrate(X, ppm, type = c("glucose", "tsp")))
  expect_error(calibrate(X, ppm, type = c(2.0)))
  expect_error(calibrate("not a matrix", ppm, type = "tsp"))
})

test_that("calibrate() works with glucose and alanine types", {
  skip_if_not(exists(".calibrate_doub"), "requires internal .calibrate_doub()")

  ppm <- seq(10, 0, length.out = 1000)
  X <- matrix(rnorm(2 * length(ppm), sd = 1e-2), nrow = 2)

  X_glu <- calibrate(X, ppm, type = "glucose")
  expect_true(is.matrix(X_glu))
  expect_equal(dim(X_glu), dim(X))

  X_ala <- calibrate(X, ppm, type = "alanine")
  expect_true(is.matrix(X_ala))
})
