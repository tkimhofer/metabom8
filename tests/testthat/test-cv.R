test_that(".kFold works with valid input", {
  Y <- matrix(rnorm(100), nrow = 20)
  folds <- .kFold(k = 5, Y)
  expect_type(folds, "list")
  expect_length(folds, 5)
  expect_true(all(sapply(folds, is.integer)))
})

test_that(".kFold defaults to LOO-CV with invalid k", {
  Y <- matrix(rnorm(10), nrow = 10)
  folds <- .kFold(k = 100, Y)
  expect_length(folds, 10)
})

test_that(".kFoldStratified works with DA", {
  Y <- matrix(rep(c("A", "B", "C"), length.out = 30), ncol = 1)
  strat <- list("DA", Y, NULL)
  folds <- .kFoldStratified(3, strat)
  expect_type(folds, "list")
  expect_length(folds, 3)
})

test_that(".kFoldStratified returns NULL on imbalanced data", {
  Y <- matrix(c(rep("A", 28), "B", "B"), ncol = 1)
  strat <- list("DA", Y, NULL)
  expect_message(folds <- .kFoldStratified(3, strat))
  expect_null(folds)
})

test_that(".mc works with valid input", {
  Y <- matrix(rnorm(100), ncol = 1)
  sets <- .mc(k = 5., Y = Y, split = 0.7)
  expect_type(sets, "list")
  expect_length(sets, 5)
  expect_true(all(sapply(sets, function(x) length(x) == floor(nrow(Y) * 0.7))))
})

test_that(".mcBalanced works with DA input", {
  Y <- matrix(rep(c("A", "B", "C"), length.out = 90), ncol = 1)
  stratified <- list("DA", Y, NULL)
  sets <- .mcBalanced(k = 5, split = 0.5, stratified)
  expect_type(sets, "list")
  expect_length(sets, 5)
})

test_that(".cvSetsMethod works for each method", {
  Y <- matrix(rep(1:100), ncol = 1)

  # k-fold
  sets1 <- .cvSetsMethod(Y, type = "R", method = "k-fold", k = 5)
  expect_type(sets1, "list")
  expect_length(sets1, 5)

  # k-fold_stratified
  sets2 <- .cvSetsMethod(Y, type = "R", method = "k-fold_stratified", k = 3)
  expect_type(sets2, "list")
  expect_length(sets2, 3)

  # MC
  sets3 <- .cvSetsMethod(Y, type = "R", method = "MC", k = 4, split = 0.7)
  expect_type(sets3, "list")
  expect_length(sets3, 4)

  # MC_balanced with categorical Y
  Yc <- matrix(rep(c("A", "B", "C"), length.out = 90), ncol = 1)
  sets4 <- .cvSetsMethod(Yc, type = "DA", method = "MC_balanced", k = 5, split = 0.6)
  expect_type(sets4, "list")
  expect_length(sets4, 5)
})

test_that(".cvSetsMethod fails with invalid method", {
  Y <- matrix(1:10, ncol = 1)
  expect_error(.cvSetsMethod(Y, type = "R", method = "invalid", k = 3))
})
