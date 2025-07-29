# tests/testthat/test_plotscores.R

test_that("plotscores() returns ggplot for PCA model", {
  skip_if_not_installed("ggplot2")

  # Simulate simple PCA model
  data(iris)
  X <- as.matrix(iris[, 1:4])
  model <- pca(X, pc = 2)

  # Test with default components and no annotations
  p <- plotscores(model)
  expect_s3_class(p, "ggplot")

  # Test with 1-level annotation
  an1 <- list(Species = iris$Species)
  p1 <- plotscores(model, an = an1)
  expect_s3_class(p1, "ggplot")

  # Test with 2-level annotation
  an2 <- list(Species = iris$Species, Shape = iris$Species)
  p2 <- plotscores(model, an = an2)
  expect_s3_class(p2, "ggplot")

  # Test with 3-level annotation (labels)
  an3 <- list(Species = iris$Species, Shape = iris$Species, Label = rownames(iris))
  p3 <- plotscores(model, an = an3)
  expect_s3_class(p3, "ggplot")

  # Test with QC points
  qc_idx <- c(1, 10, 20)
  p_qc <- plotscores(model, qc = qc_idx)
  expect_s3_class(p_qc, "ggplot")
})

test_that("plotscores() works with OPLS model and CV", {
  skip_if_not_installed("ggplot2")

  # Simulate simple OPLS model
  data(iris)
  X <- as.matrix(iris[, 1:4])
  Y <- factor(iris$Species)
  mod_opls <- opls(X, Y, center = TRUE, scale = "UV", cv = list(method = "k-fold", k = 5, split = 2/3))

  # Use only two groups for clearer model
  idx <- Y %in% c("versicolor", "virginica")
  mod_opls2 <- opls(X[idx, ], Y[idx], center = TRUE, scale = "UV", cv = list(method = "k-fold", k = 3, split = 2/3))

  # Test with CV scores
  p_cv <- plotscores(mod_opls2, cv = TRUE)
  expect_s3_class(p_cv, "ggplot")
})
