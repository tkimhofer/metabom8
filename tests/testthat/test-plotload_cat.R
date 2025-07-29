test_that("plotload_cat returns ggplot object for PCA_metabom8", {
  # Example data
  data(covid)

  # Fit PCA model
  model <- pca(X = X, pc = 2)

  # Call function with 1 annotation (color)
  g <- plotload_cat(model, pc = c(1, 2), an = list(Group = rep("A", ncol(X))))

  expect_s3_class(g, "ggplot")
  expect_true("ggplot" %in% class(g))
})


test_that("plotload_cat returns ggplot object for OPLS_metabom8", {
  data(covid)

  model <- opls(X = X, Y = an$type, plotting = FALSE)

  # Use first orthogonal component
  g <- plotload_cat(model, pc = c("1", "o1"), an = list(Group = rep("B", ncol(X))))

  expect_s3_class(g, "ggplot")
  expect_true("ggplot" %in% class(g))
})


test_that("plotload_cat fails gracefully with missing pc and an", {
  data(covid)
  model <- pca(X = X, pc = 2)

  # Test default behaviour with no pc or an
  expect_s3_class(plotload_cat(model), "ggplot")
})
