test_that("es_cdelta returns expected range and handles input correctly", {
  ref <- rnorm(100, mean = 0)
  comp <- rnorm(100, mean = 1)

  result <- es_cdelta(ref, comp)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= -1 && result <= 1)

  # Symmetry check (reverse comparison changes sign)
  expect_equal(es_cdelta(comp, ref), -result, tolerance = 1e-6)

  # Error for too few elements
  expect_error(es_cdelta(1:3, 1:3), "at least 5")

  # Error for non-numeric
  expect_error(es_cdelta(ref, as.character(comp)), "must be numeric")
})
