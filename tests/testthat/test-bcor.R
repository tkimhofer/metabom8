# tests/testthat/test-bcor.R

test_that("bcor returns matrix of same dimensions as input", {
  set.seed(1)
  mat <- matrix(rnorm(100), nrow = 10)
  result <- bcor(mat)
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(mat))
})

test_that("bcor handles a single spectrum (row vector)", {
  vec <- rnorm(100)
  result <- bcor(vec)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), length(vec))
})

test_that("bcor replaces NA with zero", {
  mat <- matrix(rnorm(100), nrow = 10)
  mat[1, 1] <- NA
  expect_message(res <- bcor(mat), "replacing with zeros")
  expect_false(any(is.na(res)))
})

# test_that("bcor and bline return identical output", {
#   mat <- matrix(rnorm(100), nrow = 5)
#   expect_equal(bline(mat), bcor(mat))
# })

test_that("bcor returns different values than uncorrected spectra", {
  set.seed(42)
  n <- 1000
  ppm <- seq(0, 10, length.out = n)
  genSpec <- function(h, nz) {
    baseline <- 0.01*(ppm - 5)^2 + 0.05*sin(2 * pi * ppm / 2) + 0.2
    peak <- h * exp(-((ppm - 5)^2) / (2 * 0.1^2))
    baseline + peak + rnorm(n, 0, nz)
  }
  X <- rbind(genSpec(40, 0.02), genSpec(40, 0.015))
  # matspec(X, ppm)
  X_bcor <- bcor(X, lambda=1e5)
  # matspec(X_bcor, ppm)
  expect_false(isTRUE(all.equal(X, X_bcor)))
})
