# tests/testthat/test-storm.R

test_that("storm selects a valid subset of spectra", {
  # Simulate test data
  set.seed(123)
  X <- matrix(rnorm(1000), nrow = 20)  # 20 spectra, 50 variables
  ppm <- seq(0.5, 2, length.out = ncol(X))

  # Inject a fake peak at 1.2 ppm for spectra 1:10
  peak_idx <- which.min(abs(ppm - 1.2))
  X[1:10, peak_idx + (-2:2)] <- X[1:10, peak_idx + (-2:2)] + 3

  # Mock scale_rcpp and get_idx if not already defined
  if (!exists("scale_rcpp")) {
    scale_rcpp <- function(x, center = TRUE, scale = FALSE) {
      if (center) x <- sweep(x, 2, colMeans(x), "-")
      if (scale) x <- sweep(x, 2, apply(x, 2, sd), "/")
      return(x)
    }
  }

  if (!exists("get_idx")) {
    get_idx <- function(range, ppm) {
      which(ppm >= min(range) & ppm <= max(range))
    }
  }

  subset_idx <- storm(X, ppm, b = 3, q = 0.05, idx.refSpec = 1, shift = c(1.15, 1.25))

  matspec(X, ppm)
  spec(X[1,], ppm)
  dim(X)
  # Expectations
  expect_type(subset_idx, "integer")
  expect_true(length(subset_idx) > 0)
  expect_true(all(subset_idx %in% seq_len(nrow(X))))
  expect_true(1 %in% subset_idx)  # reference spectrum should be included
})
