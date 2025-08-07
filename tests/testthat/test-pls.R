# test-pls.R
test_that("PLS model runs correctly with numeric input", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  Y <- iris$Species

  model <- pls(X, Y, center = TRUE, scale = "UV",
               cv = list(method = "k-fold", k = 3, split = 2/3),
               maxPCo = 5, plotting = FALSE)

  expect_s4_class(model, "PLS_metabom8")
  expect_true(all(c("t_pred", "p_pred", "w_pred") %in% slotNames(model)))
  expect_equal(nrow(model@t_pred), nrow(X))
  expect_equal(nrow(model@Y_res), nrow(X))
})

test_that("PLS returns valid summary metrics", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  Y <- iris$Species
  model <- pls(X, Y, plotting = FALSE)

  expect_type(model@summary, "list")
  expect_named(model@summary, c("PC_pred","PC_orth", "R2X", "AUROC", "AUROC_CV"), ignore.order = TRUE)
})

test_that("PLS handles factor and numeric Y properly", {
  X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  Y <- factor(rep(c("A", "B"), length.out = 100))
  model <- pls(X, Y, plotting = FALSE)

  expect_s4_class(model, "PLS_metabom8")
  expect_equal(length(model@Y$ori), 100)
})
