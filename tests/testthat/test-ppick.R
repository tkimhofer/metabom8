test_that("ppick returns expected output structure", {
  set.seed(1)
  X <- matrix(rnorm(5 * 100), nrow = 5)
  ppm <- seq(0.5, 9.5, length.out = 100)

  peaks_max <- ppick(X, ppm, type = "max")
  peaks_min <- ppick(X, ppm, type = "min")
  peaks_both <- ppick(X, ppm, type = "both")

  expect_length(peaks_max, 5)
  expect_type(peaks_max, "list")
  expect_true(all(sapply(peaks_max, is.data.frame)))
  expect_true(all(c("idc", "ppm", "Int", "Etype") %in% colnames(peaks_both[[1]])))

  expect_true(all(peaks_max[[1]]$Etype == -1))
  expect_true(all(peaks_min[[1]]$Etype == 1))
  expect_true(all(peaks_both[[1]]$Etype %in% c(-1, 1)))
})
