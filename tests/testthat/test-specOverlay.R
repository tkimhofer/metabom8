test_that("specOverlay returns ggplot object", {
  mat <- matrix(rnorm(500), nrow = 10)
  ppm <- seq(0.1, 5, length.out = ncol(mat))
  group <- rep(c("A", "B"), length.out = nrow(mat))

  p <- specOverlay(X = mat, ppm = ppm, shift = c(0.2, 0.4), an = list(Group = group))
  expect_s3_class(p, "ggplot")
})

test_that("specOverlay fails gracefully with wrong dimensions", {
  expect_error(specOverlay(X = matrix(1:10, nrow = 5), ppm = 1:3, an = list("A")),
               "Non-matching dimensions")
})

test_that("specOverlay works with color and linetype", {
  mat <- matrix(rnorm(500), nrow = 10)
  ppm <- seq(0.1, 5, length.out = ncol(mat))
  group <- rep(c("A", "B"), length.out = 10)
  run_order <- seq_len(10)

  p <- specOverlay(X = mat, ppm = ppm, shift = c(0.2, 0.4),
                   an = list(Facet = group, Color = run_order, Line = group))
  expect_s3_class(p, "ggplot")
})
