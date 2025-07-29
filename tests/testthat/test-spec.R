test_that("spec handles input validation correctly", {
  expect_error(spec(matrix(1:10, nrow = 2), 1:10), "x must be a numeric vector")
  expect_error(spec(1:10, 1:9), "x and ppm must have equal length")
})

test_that("spec returns plotly object in interactive mode", {
  ppm <- seq(0, 10, length.out = 100)
  x <- sin(ppm)
  plt <- spec(x, ppm, interactive = TRUE)
  expect_s3_class(plt, "plotly")
})

test_that("spec produces silent base plot and overlay", {
  ppm <- seq(0, 10, length.out = 100)
  x1 <- sin(ppm)
  x2 <- cos(ppm)

  expect_silent(spec(x1, ppm, interactive = FALSE))
  expect_silent(spec(x2, ppm, interactive = FALSE, add = TRUE, col = "red"))
})

test_that("spec respects shift range", {
  ppm <- seq(0, 10, length.out = 100)
  x <- sin(ppm)
  idx <- get_idx(c(2, 4), ppm)
  expect_silent(spec(x, ppm, shift = c(2, 4), interactive = FALSE))
  expect_equal(range(ppm[idx]), c(4, 2)) # reversed by default
})
