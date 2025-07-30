test_that("alignSegment returns matrix of correct dimensions", {
  mat <- matrix(rnorm(100), nrow = 10)
  result <- alignSegment(mat)
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(mat))
})

test_that("alignSegment keeps low-correlation rows unshifted", {
  set.seed(42)
  mat <- matrix(rnorm(100), nrow = 10)
  mat[5, ] <- rnorm(10)  # deliberately uncorrelated row
  result <- alignSegment(mat, clim = 0.99)
  expect_equal(result[5, ], mat[5, ])  # Should be unchanged
})

test_that("median reference is used correctly", {
  mat <- matrix(rnorm(100), nrow = 10)
  result <- alignSegment(mat, med = TRUE)
  expect_equal(nrow(result), 10)
})

test_that("normalization does not crash", {
  mat <- matrix(rnorm(100), nrow = 10)
  expect_silent(alignSegment(mat, norm = TRUE))
})

test_that("alignment improves similarity for correlated rows", {
  mat <- matrix(rnorm(100), nrow = 10)
  mat[6, ] <- mat[1, ]  # perfectly correlated with row 1
  aligned <- alignSegment(mat, idx_ref = 1, clim = 0.5)
  expect_true(cor(aligned[6, ], aligned[1, ]) > 0.9)
})

# edge cases

test_that("out-of-bounds idx_ref throws error", {
  mat <- matrix(rnorm(100), nrow = 10)
  expect_error(alignSegment(mat, idx_ref = 20, med=FALSE))
})

test_that("non-matrix input throws error", {
  expect_error(alignSegment(list(1, 2, 3)))
})

test_that("too few rows returns error", {
  mat <- matrix(rnorm(10), nrow = 1)
  expect_error(alignSegment(mat))
})

test_that("all rows shift when clim = 0", {
  mat <- matrix(rnorm(100), nrow = 10)
  result <- alignSegment(mat, clim = 0)
  expect_false(all(result == mat))  # Some rows should have shifted
})
