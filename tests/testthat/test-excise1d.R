test_that("excise1d removes expected ppm regions and returns correct format", {
  ppm <- seq(0, 10, length.out = 1000)
  X <- matrix(rnorm(100 * 1000), nrow = 100)

  result <- excise1d(X, ppm)

  expect_type(result, "list")
  expect_named(result, c("Xc", "ppc"))
  expect_true(is.matrix(result$Xc))
  expect_true(is.numeric(result$ppc))
  expect_equal(nrow(result$Xc), nrow(X))
  expect_equal(ncol(result$Xc), length(result$ppc))
  expect_lt(ncol(result$Xc), ncol(X))  # should have fewer columns
})
