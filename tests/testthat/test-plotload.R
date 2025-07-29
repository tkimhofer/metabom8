test_that("plotload returns ggplot object for PCA", {
  data(covid)
  model <- pca(X)
  g <- plotload(model, pc = 1)
  expect_s3_class(g, "ggplot")
})

test_that("plotload errors with invalid model", {
  expect_error(plotload(matrix(1:10), shift = c(0, 1)), "Input must be")
})

test_that("plotload errors with insufficient shift", {
  data(covid)
  model <- pca(X)
  expect_error(plotload(model, shift = c(20, 21)), "Insufficient")
})
