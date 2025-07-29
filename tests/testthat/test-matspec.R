test_that("matspec validates input", {
  X <- matrix(rnorm(100), nrow = 10)
  ppm <- seq(0, 10, length.out = 10)
  expect_error(matspec(X, ppm[-1]), "X and ppm dimensions do not match")
})

test_that("matspec runs silently in base mode", {
  X <- matrix(rnorm(100), nrow = 10)
  ppm <- seq(0, 10, length.out = 10)
  expect_silent(matspec(X, ppm, interactive = FALSE))
})

test_that("matspec returns a plotly object in interactive mode", {
  X <- matrix(rnorm(100), nrow = 10)
  ppm <- seq(0, 10, length.out = 10)
  plt <- matspec(X, ppm, interactive = TRUE)
  expect_s3_class(plt, "plotly")
})

test_that("matspec handles shift range correctly", {
  X <- matrix(rnorm(100), nrow = 10)
  ppm <- seq(0, 10, length.out = 10)
  expect_silent(matspec(X, ppm, shift = c(2, 5), interactive = FALSE))
})
