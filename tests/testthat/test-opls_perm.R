test_that("opls_perm returns valid output", {
  set.seed(123)
  X <- matrix(rnorm(200), nrow = 20)
  colnames(X) <- sprintf("ppm%.4f", seq(0, 1, length.out = 10))
  Y <- rep(c("A", "B"), each = 10)
  model <- opls(X, Y = Y, ncomp = 1, ortho = 0)

  res <- opls_perm(model, n = 3, plot = FALSE)

  expect_s3_class(res, "data.frame")
  expect_true("Non-permuted" %in% res$model)
  expect_equal(nrow(res), 4)  # 3 permutations + 1 original
})
