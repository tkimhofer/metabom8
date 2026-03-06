expect_valid_spectra <- function(x) {

  expect_true(is.list(x))

  expect_named(x, c("X","ppm","meta"), ignore.order = TRUE)

  expect_true(is.matrix(x$X))
  expect_type(x$ppm, "double")
  expect_s3_class(x$meta, "data.frame")

  expect_equal(ncol(x$X), length(x$ppm))

  expect_true(all(diff(x$ppm) != 0))

  expect_false(anyNA(x$X))
}

data("hiit_raw")
hiit_raw = list(X = hiit_raw$X[2:3,], ppm = hiit_raw$ppm, meta = hiit_raw$meta[2:3,])

test_that("normalize_pqn preserves spectral structure", {

  res <- pqn(hiit_raw)

  expect_valid_spectra(res)

})



test_that("preprocessing pipeline preserves invariants", {

  hiit <-
    hiit_raw |>
    calibrate(type = "tsp") |>
    excise() |>
    correct_baseline() |>
    align_spectra() |>
    pqn()

  expect_valid_spectra(hiit)

})

preproc_funs <- list(
  calibrate = function(x) calibrate(x),
  correct_lw = function(x) correct_lw(x),
  excise = function(x) excise(x),
  correct_baseline_lin = function(x) correct_baseline(x, method='linear'),
  correct_baseline_als = function(x) correct_baseline(x, method='asls'),
  binning = function(x) binning(x, width=0.01),
  norm_eretic     = function(x) norm_eretic(x),
  norm_pqn     = function(x) pqn(x),
  align_segment = function(x) align_segment(x, shift = c(3, 3.1)),
  align_spectra = function(x) align_spectra(x)
)

for (nm in names(preproc_funs)) {

  test_that(paste(nm, "preserves spectral invariants"), {
    res <- preproc_funs[[nm]](hiit_raw)
    expect_valid_spectra(res)
  })

}

