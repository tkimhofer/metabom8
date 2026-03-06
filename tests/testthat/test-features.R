set.seed(1)

## ppick

test_that("ppick returns expected output structure", {

  X <- matrix(rnorm(5 * 100), nrow = 5)
  ppm <- seq(0.5, 9.5, length.out = 100)

  peaks_max  <- ppick(X, ppm, type = "max")
  peaks_min  <- ppick(X, ppm, type = "min")
  peaks_both <- ppick(X, ppm, type = "both")

  expect_length(peaks_max, 5)
  expect_true(all(vapply(peaks_max, is.data.frame, logical(1))))

  expect_named(peaks_both[[1]], c("idc","ppm","Int","Etype"))

  expect_true(all(peaks_max[[1]]$Etype == -1))
  expect_true(all(peaks_min[[1]]$Etype == 1))
  expect_true(all(peaks_both[[1]]$Etype %in% c(-1,1)))

  expect_equal(peaks_both[[1]]$ppm,
               ppm[peaks_both[[1]]$idc])

})

## ppick2

test_that("ppick2 validates input parameters", {

  X <- matrix(rnorm(100), nrow = 1)
  ppm <- seq_len(100)

  expect_error(ppick2(X, ppm, fil_n = 4))
  expect_error(ppick2(X, ppm, fil_n = 3))
  expect_error(ppick2(X, ppm, fil_p = 5, fil_n = 5))

})

test_that("ppick2 returns correct structure", {

  ppm <- seq(0, 10, length.out = 200)

  peak1 <- exp(-(ppm - 3)^2 / 0.01)
  peak2 <- exp(-(ppm - 7)^2 / 0.02)

  X <- matrix(peak1 + peak2, nrow = 1)

  peaks <- ppick2(X, ppm, type = "max", min_snr = 1)

  expect_type(peaks, "list")
  expect_length(peaks, 1)

  expect_true(is.data.frame(peaks[[1]]))

  expect_named(peaks[[1]],
               c("idc","ppm","Int","Etype",
                 "height","snr","curvature","prominence"),
               ignore.order = TRUE)

})

test_that("ppick2 detects known peaks", {

  ppm <- seq(0, 10, length.out = 500)

  peak1 <- exp(-(ppm - 3)^2 / 0.02)
  peak2 <- exp(-(ppm - 7)^2 / 0.02)

  X <- matrix(peak1 + peak2, nrow = 1)

  peaks <- ppick2(X, ppm, type = "max", min_snr = 1)

  detected <- peaks[[1]]$ppm

  expect_true(any(abs(detected - 3) < 0.05))
  expect_true(any(abs(detected - 7) < 0.05))

})

test_that("ppick2 infers ppm from colnames", {

  ppm <- seq(0, 10, length.out = 200)

  X <- matrix(exp(-(ppm - 4)^2 / 0.02), nrow = 1)

  colnames(X) <- as.character(ppm)

  peaks <- ppick2(X, ppm = NULL, min_snr = 1)

  expect_true(is.list(peaks))
})
