set.seed(15)

###constructors

test_that("kfold constructor returns valid object", {

  cv <- kfold(k=5)
  expect_type(cv, "list")
  expect_equal(cv$k, 5)

})

test_that("kfold rejects invalid k", {

  expect_error(kfold(1))
  expect_error(kfold(NA))
})

test_that("stratified_kfold constructor works", {

  cv <- stratified_kfold(k = 5, type = "R", probs = c(0,0.5,1))

  expect_type(cv, "list")

  expect_equal(cv$k, 5)
  expect_equal(cv$type, "R")
  expect_equal(cv$probs, c(0,0.5,1))
})

test_that("mc constructor works", {

  cv <- mc(k = 10, split = 0.7)

  expect_type(cv, "list")

  expect_equal(cv$k, 10)
  expect_equal(cv$split, 0.7)

})

test_that("mc validates split", {

  expect_error(mc(10, -0.2))
  expect_error(mc(10, 1))
  expect_error(mc(10, 0))

})

test_that("balanced_mc constructor works", {

  cv <- balanced_mc(k = 5, split = 0.7, type = "DA")

  expect_type(cv, "list")

  expect_equal(cv$k, 5)
  expect_equal(cv$split, 0.7)
  expect_equal(cv$type, "DA")

})


### instantiation

test_that("kfold instantiates correct folds", {

  n = 100
  Y <- matrix(rnorm(n), ncol = 1)

  cv <- kfold(5)
  inst <- .arg_check_cv(cv, model_type='R', n=n, Y_prepped = Y)

  folds <- inst$train

  expect_length(folds, 5)
  expect_true(all(sapply(folds, is.integer)))
  expect_true(all(sapply(folds, function(i) all(i <= nrow(Y)))))
  expect_true(all(sapply(folds, length) > 0))
  expect_true(all(unlist(folds) >= 1))
  expect_true(all(unlist(folds) <= nrow(Y)))
  expect_true(length(unique(unlist(folds))) == nrow(Y))
})


test_that("stratified_kfold preserves class balance", {

  n = 100
  Y <- matrix(sample(c(0,1), n, replace = TRUE, prob = c(0.8,0.2)), ncol=1)

  cv <- stratified_kfold(5, type="DA")

  inst <- .arg_check_cv(cv, model_type='DA', n=n, Y_prepped = Y)
  folds <- inst$train

  props <- sapply(folds, function(i) mean(Y[i] == 1))

  expect_true(sd(props) < 0.1)
  expect_true(all(unlist(folds) >= 1))
  expect_true(all(unlist(folds) <= nrow(Y)))
  expect_true(all(props > 0))

})


test_that("mc generates correct training sizes", {

  n = 100
  Y <- matrix(rnorm(n), ncol = 1)
  cv <- mc(10, split = 0.7)

  inst <- .arg_check_cv(cv, model_type='R', n=n, Y_prepped = Y)
  folds <- inst$train

  sizes <- sapply(folds, length)

  expect_true(all(abs(sizes - 70) <= 2))
  expect_true(all(unlist(folds) >= 1))
  expect_true(all(sizes == floor(n * 0.7)))
})


test_that("balanced_mc returns correct number of splits", {

  Y <- c(rep(0,80), rep(1,20))
  Y <- matrix(Y, ncol = 1)

  strat <- list("DA", Y, NULL)

  sets <- .mcBalanced(k = 10, split = 0.7, stratified = strat)

  expect_length(sets, 10)
  expect_true(all(unlist(sets) >= 1))
  expect_true(all(unlist(sets) <= nrow(Y)))

})

test_that("balanced_mc produces balanced class samples", {

  Y <- c(rep(0,80), rep(1,20))
  Y <- matrix(Y, ncol = 1)

  strat <- list("DA", Y, NULL)

  sets <- .mcBalanced(10, 0.7, strat)

  props <- sapply(sets, function(i) mean(Y[i] == 1))

  expect_true(sd(props) < 0.05)

})

test_that("balanced_mc samples without replacement", {

  Y <- c(rep(0,80), rep(1,20))
  Y <- matrix(Y, ncol = 1)

  strat <- list("DA", Y, NULL)

  sets <- .mcBalanced(5, 0.7, strat)

  expect_false(any(sapply(sets, function(i) any(duplicated(i)))))
})

test_that("balanced_mc detects extreme imbalance", {

  Y <- c(rep(0,95), rep(1,5))
  Y <- matrix(Y, ncol = 1)

  strat <- list("DA", Y, NULL)

  expect_error(.mcBalanced(10, 0.7, strat))

})

test_that("bootstrap sampling contains duplicates", {

  n = 100
  Y <- matrix(rnorm(n), ncol = 1)

  cv <- balanced_boot(5, split = 0.7, type = 'R', probs=c(0, 0.5, 1))

  inst <- .arg_check_cv(cv, model_type='R', n=n, Y_prepped = Y)
  folds <- inst$train

  expect_true(any(sapply(folds, function(i) any(duplicated(i)))))
  expect_true(all(unlist(folds) >= 1))
  expect_true(all(unlist(folds) <= nrow(Y)))
})
