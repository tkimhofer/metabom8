test_that("specload returns ggplot object with PCA input", {
  data(covid)
  X <- covid$X
  an <- covid$an
  model <- pca(X)
  g <- specload(model, shift = c(1, 2), an = list(an$type), pc = 1)
  expect_s3_class(g, "ggplot")
})

test_that("specload stops for wrong model input", {
  mat <- matrix(rnorm(100), nrow = 10)
  expect_error(specload(mod = mat, shift = c(0, 1), an = list("A")), "Requires")
})

test_that("specload handles invalid ppm range", {
  data(covid)
  X <- covid$X
  an <- covid$an
  model <- pca(X)
  expect_error(specload(model, shift = c(20, 21), an = list(an$type)), "Insufficient")
})
