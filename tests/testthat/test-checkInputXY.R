#
# test_that(".checkXclassNas passes for valid numeric matrix", {
#   X <- matrix(rnorm(10), nrow = 5)
#   expect_silent(.checkXclassNas(X))
# })
#
# test_that(".checkXclassNas fails for non-matrix input", {
#   x_df <- data.frame(a = 1:5, b = 6:10)
#   expect_error(.checkXclassNas(x_df), "Input X must be a numeric matrix")
# })
#
# test_that(".checkXclassNas fails with NA/NaN/Inf values", {
#   X <- matrix(rnorm(9), nrow = 3)
#   X[1, 1] <- NA
#   expect_error(.checkXclassNas(X), "contains missing")
#
#   X[1, 1] <- NaN
#   expect_error(.checkXclassNas(X), "contains missing")
#
#   X[1, 1] <- Inf
#   expect_error(.checkXclassNas(X), "contains missing")
# })
#
# # ----------- .checkDimXY Tests -----------
#
# test_that(".checkDimXY passes for valid X and Y", {
#   X <- matrix(rnorm(30), nrow = 10)
#   Y <- matrix(rnorm(20), nrow = 10)
#   expect_silent(.checkDimXY(X, Y))
# })
#
# test_that(".checkDimXY fails if row counts differ", {
#   X <- matrix(rnorm(20), nrow = 5)
#   Y <- matrix(rnorm(24), nrow = 6)
#   expect_error(.checkDimXY(X, Y), "Dimensions of input X and Y do not match")
# })
#
# test_that(".checkDimXY fails if X has fewer columns than Y", {
#   X <- matrix(rnorm(30), nrow = 10, ncol = 2)
#   Y <- matrix(rnorm(30), nrow = 10, ncol = 3)
#   expect_error(.checkDimXY(X, Y), "Number of variables .* should be higher than")
# })
#
# # ----------- .checkYclassNas Tests -----------
#
# test_that(".checkYclassNas works for numeric Y (regression)", {
#   Y <- rnorm(10)
#   result <- .checkYclassNas(Y)
#   expect_true(is.matrix(result[[1]]))
#   expect_equal(ncol(result[[1]]), 1)
#   expect_equal(result[[3]], "R")
# })
#
# test_that(".checkYclassNas works for binary factor Y (classification)", {
#   Y <- factor(rep(c("A", "B"), each = 5))
#   result <- .checkYclassNas(Y)
#   expect_true(is.matrix(result[[1]]))
#   expect_equal(ncol(result[[1]]), 1)
#   expect_equal(result[[3]], "DA")
# })
#
# test_that(".checkYclassNas works for multi-class factor Y", {
#   Y <- factor(c("A", "B", "C", "A", "B", "C"))
#   result <- .checkYclassNas(Y)
#   expect_true(is.matrix(result[[1]]))
#   expect_gt(ncol(result[[1]]), 1)
#   expect_equal(result[[3]], "DA-mY")
# })
#
# test_that(".checkYclassNas errors for constant class", {
#   Y <- rep("A", 5)
#   expect_error(.checkYclassNas(Y), "only a single level")
# })
#
# test_that(".checkYclassNas errors on NA/NaN/Inf", {
#   Y <- c(1, 2, NA, 4)
#   expect_error(.checkYclassNas(Y), "contains NA")
#
#   Y <- c(1, 2, NaN, 4)
#   expect_error(.checkYclassNas(Y), "contains NA")
#
#   Y <- c(1, 2, Inf, 4)
#   expect_error(.checkYclassNas(Y), "contains NA")
# })
#
#
# # ----------- .prepareInputs Tests -----------
#
# test_that("Regression task: numeric Y vector and matrix", {
#   set.seed(1)
#   X <- matrix(rnorm(50), nrow = 10)
#   Y <- rnorm(10)
#
#   result <- .prepareInputs(X, Y, center = TRUE, scale = "None")
#   expect_true(is.matrix(result$X))
#   expect_equal(result$type, "R")
#   expect_true(is.numeric(result$Y))
#   expect_true(is.matrix(result$Y))
#   expect_true(ncol(result2$Y) == 1)
#   # expect_true(is.numeric(result$Y))
#   # expect_null(result$Ydummy)
#
#   # Y as single-column matrix
#   result2 <- .prepareInputs(X, matrix(Y, ncol = 1), center = TRUE, scale = "UV")
#   expect_equal(result2$type, "R")
#   expect_true(is.numeric(result$Y))
#   expect_true(is.matrix(result2$Y))
#   expect_true(ncol(result2$Y) == 1)
# })
#
# test_that("Binary classification: factor Y vector and matrix", {
#   set.seed(2)
#   X <- matrix(rnorm(50), nrow = 10)
#   Y <- factor(sample(c("A", "B"), 10, replace = TRUE))
#
#   result <- .prepareInputs(X, Y, center = FALSE, scale = "Pareto")
#   expect_equal(result$type, "DA")
#   expect_true(is.numeric(result$Y))
#   expect_equal(ncol(result$Ydummy), 2)
#
#   # Y as single-column matrix
#   Y_mat <- matrix(as.character(Y), ncol = 1)
#   result2 <- .prepareInputs(X, Y_mat, center = FALSE, scale = "None")
#   expect_equal(result2$type, "DA")
#   expect_equal(ncol(result2$Ydummy), 2)
# })
#
# test_that("Multiclass classification: 3-level factor", {
#   set.seed(3)
#   X <- matrix(rnorm(60), nrow = 10)
#   Y <- factor(sample(c("D1", "B1", "D2"), 10, replace = TRUE))
#
#   result <- .prepareInputs(X, Y, center = TRUE, scale = "UV")
#   expect_equal(result$type, "DA-mY")
#   expect_equal(ncol(result$Ydummy), 2)
#   expect_equal(nrow(result$Ydummy), 3)
# })
#
# test_that("Invalid center or scale throws errors", {
#   X <- matrix(rnorm(20), nrow = 5)
#   Y <- rnorm(5)
#
#   expect_error(.prepareInputs(X, Y, center = "yes", scale = "UV"), "Check center")
#   expect_error(.prepareInputs(X, Y, center = TRUE, scale = "BadScale"), "Check scale")
# })
#
#
# test_that(".prepareY handles numeric vector (regression case)", {
#   y <- c(1.2, 3.4, 5.6)
#   result <- .prepareY(y)
#
#   expect_type(result, "list")
#   expect_true(is.matrix(result[[1]]))
#   expect_equal(dim(result[[1]]), c(3, 1))
#   expect_equal(result[[1]][, 1], y)
#   expect_equal(ncol(result[[1]]), 1)
#   expect_equal(result[[2]], data.frame())  # mapping should be empty
# })
#
# test_that(".prepareY handles binary classification (factor)", {
#   y <- factor(c("A", "B", "A", "B"))
#   result <- .prepareY(y)
#
#   expect_type(result, "list")
#   expect_true(is.matrix(result[[1]]))
#   expect_equal(dim(result[[1]]), c(4, 1))
#   expect_equal(unique(as.vector(result[[1]])), c(1, 2))
#   expect_true(all(result[[2]]$Level %in% c("A", "B")))
#   expect_equal(sort(result[[2]]$Code), c(1, 2))
# })
#
# test_that(".prepareY handles binary classification (character)", {
#   y <- c("yes", "no", "yes", "no")
#   result <- .prepareY(y)
#
#   expect_type(result, "list")
#   expect_true(is.matrix(result[[1]]))
#   expect_equal(dim(result[[1]]), c(4, 1))
#   expect_equal(sort(unique(as.vector(result[[1]]))), c(1, 2))
#   expect_equal(nrow(result[[2]]), 2)
# })
#
# test_that(".prepareY handles multiclass classification", {
#   y <- factor(c("cat", "dog", "cat", "mouse", "dog"))
#   result <- .prepareY(y)
#
#   y_mat <- result[[1]]
#   mapping <- result[[2]]
#
#   expect_type(result, "list")
#   expect_true(is.matrix(y_mat))
#   expect_equal(dim(y_mat), c(5, 3))
#   expect_true(all(y_mat %in% c(1, -1)))
#   expect_equal(mapping$Level, c("cat", "dog", "mouse"))
#   expect_equal(mapping$Column, 1:3)
# })
#
# test_that(".prepareY returns correct colnames for multiclass", {
#   y <- c("low", "medium", "high", "medium", "high")
#   result <- .prepareY(y)
#   expect_equal(colnames(result[[1]]), c("high", "low", "medium"))
# })
#
#
# test_that(".prepareY errors on multi-column numeric input", {
#   Y <- matrix(1:6, nrow = 3, ncol = 2)
#   expect_error(.prepareY(Y), "Only single-column numeric input is supported for Y at this time")
# })
