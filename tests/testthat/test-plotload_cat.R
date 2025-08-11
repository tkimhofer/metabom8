# test_that("plotload_cat returns ggplot object for PCA_metabom8", {
#   # Example data
#   data(covid)
#   X <- covid$X
#   an <- covid$an
#   # Fit PCA model
#   model <- pca(X = X, pc = 2)
#
#   # Call function with 1 annotation (color)
#   g <- plotload_cat(model, pc = c(1, 2), an = list(Group = rep("A", ncol(X))))
#
#   expect_s3_class(g, "ggplot")
#   expect_true("ggplot" %in% class(g))
# })
#
# test_that("plotload_cat works with default pc and an when features < 150", {
#   X_small <- iris[1:149, 1:4]
#   model_small <- pca(X = X_small, pc = 2)
#
#   p <- plotload_cat(model_small)
#
#   expect_s3_class(p, "ggplot")
# })
#
#
# test_that("plotload_cat fails gracefully with missing pc and an", {
#   data(covid)
#   X <- covid$X
#   an <- covid$an
#
#   model <- pca(X = X, pc = 2)
#
#   # Expect error due to missing 'an' with large feature set
#   expect_error(
#     plotload_cat(model),
#     "Annotation list 'an' is missing and too many features to set a default."
#   )
#
# })
