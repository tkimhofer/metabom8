test_that("predict_opls returns correct output structure", {
  # Simulate simple model
  X <- matrix(rnorm(100 * 10), 100, 10)
  Y <- rep(c("A", "B"), each = 50)
  model <- opls(X, Y = factor(Y), plotting = FALSE)

  # Predict
  pred <- predict_opls(model, X)

  expect_type(pred, "list")
  expect_named(pred, c("Y_predicted", "t_pred", "t_orth", "t_orth_pca"), ignore.order = TRUE)
  expect_equal(length(pred$Y_predicted), nrow(X))
  expect_equal(nrow(pred$t_pred), nrow(X))
  expect_equal(nrow(pred$t_orth), nrow(X))
})
