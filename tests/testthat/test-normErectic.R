test_that("normErectic returns correct normalization or integral", {
  set.seed(123)
  ppm <- seq(0.2, 16, length.out = 16000)
  X <- matrix(rnorm(10 * length(ppm), mean = 1e3, sd = 10), nrow = 10)
  colnames(X) <- round(ppm, 4)

  # Add synthetic ERETIC peak near 12 ppm
  eretic_idx <- which(ppm >= 11.95 & ppm <= 12.05)
  X[, eretic_idx] <- X[, eretic_idx] + 500

  # Test integral
  integrals <- normErectic(X, integr = TRUE)
  expect_true(is.numeric(integrals))
  expect_equal(length(integrals), nrow(X))

  # Test normalized output
  Xn <- normErectic(X, integr = FALSE)
  expect_true(is.matrix(Xn))
  expect_equal(dim(Xn), dim(X))
  expect_true(all(is.finite(Xn)))
})
