test_that("pca() runs and returns valid PCA_metabom8 object", {
  # Sample input matrix
  set.seed(123)
  X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  colnames(X) <- paste0("V", 1:10)
  rownames(X) <- paste0("S", 1:100)

  # Run PCA with default NIPALS
  model <- pca(X, pc = 3, center = TRUE, scale = "UV", method = "nipals")

  # Basic structural checks
  expect_s4_class(model, "PCA_metabom8")
  expect_true(is.matrix(model@t))
  expect_true(is.matrix(model@p))
  expect_equal(ncol(model@t), 3)
  expect_equal(nrow(model@p), 3)
  expect_equal(nrow(model@X), 100)
  expect_named(model@Parameters, c("center", "scale", "nPC", "R2", "tssx"))
})

test_that("pca() handles too many components gracefully", {
  X <- matrix(rnorm(30), nrow = 5)
  expect_warning(model <- pca(X, pc = 10), regexp = "Too many components")
  expect_s4_class(model, "PCA_metabom8")
})

test_that("pca() errors on invalid input types", {
  X <- matrix(rnorm(100), nrow = 10)

  expect_error(pca(X, scale = "invalid"), regexp = "Check scale parameter")
  expect_error(pca(X, center = "yes"), regexp = "must be logical")
  expect_error(pca(X, method = 42), regexp = "must be a character")
})
