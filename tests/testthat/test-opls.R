test_that("opls returns a valid OPLS_metabom8 object", {
  skip_if_not_installed("metabom8")
  skip_on_cran()

  set.seed(123)

  # Simulate a test dataset
  n <- 30
  p <- 10
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("V", seq_len(p))
  Y <- sample(c("A", "B"), size = n, replace = TRUE)

  # Run model
  model <- opls(X, Y, plotting = FALSE, center = TRUE, scale = "UV",
                cv = list(method = "k-fold_stratified", k = 5, split = 2/3),
                maxPCo = 3)

  # Check class
  expect_s4_class(model, "OPLS_metabom8")

  # Check slots
  expect_true(!is.null(model@t_pred))
  expect_true(!is.null(model@p_pred))
  expect_true(!is.null(model@X_scaled))
  expect_true(!is.null(model@summary))
  expect_true(model@nPC >= 0)

  # Check dimensions
  expect_equal(nrow(model@X_scaled), nrow(X))
  expect_equal(ncol(model@X_scaled), ncol(X))

  # Check residual matrix
  expect_equal(dim(model@X_res), dim(model@X_scaled))
})
