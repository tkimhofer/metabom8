# tests/testthat/test-storm.R

test_that("storm selects a valid subset of spectra", {
  set.seed(153)

  # Parameters
  n_spectra <- 10    # number of spectra (rows)
  n_vars <- 500      # number of variables (columns)
  ppm <- seq(0.5, 2, length.out = n_vars)  # ppm scale

  # Baseline noise matrix
  X <- matrix(rnorm(n_spectra * n_vars, 0, 0.3), nrow = n_spectra, ncol = n_vars)

  # Gaussian peak parameters
  peak_center <- 1.2
  peak_width <- 0.05

  # Peak amplitudes vary per spectrum (positive values only)
  peak_amplitude <- rnorm(n_spectra, mean = 500, sd = 100)
  peak_amplitude[peak_amplitude < 0] <- 0  # no negative amplitudes

  # Gaussian peak shape (length n_vars)
  gauss_shape <- exp(-0.5 * ((ppm - peak_center) / peak_width)^2)

  # Inject peaks with different amplitudes per spectrum
  for (i in 1:5) {
    X[i, ] <- X[i, ] + peak_amplitude[i] * gauss_shape
  }

  matspec(X, ppm)

  subset_idx <- storm(X, ppm, b = 30, q = 0.05, idx.refSpec = 1, shift = c(1.15, 1.25))

  # Expectations
  expect_type(subset_idx, "integer")
  expect_true(all(subset_idx %in% seq_len(nrow(X))))
  expect_true(1 %in% subset_idx)  # ref should be included
  expect_equal(sort(subset_idx), 1:5) # acc to synth ds
})
