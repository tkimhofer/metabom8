set.seed(1)

test_that(".plsComponentCv returns valid structure for first component", {

  n <- 30
  p <- 5
  q <- 1

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)

  cv.set <- list(1:15, 16:30)
  Ycs_fold <- list(Y, Y)

  acc <- list(
    sum_test = matrix(0, n, q),
    n_test   = rep(0, n),
    sum_train = matrix(0, n, q),
    n_train   = rep(0, n)
  )

  res <- .plsComponentCv(
    X = X,
    cv.set = cv.set,
    Ycs_fold = Ycs_fold,
    nc = 1,
    mod.cv = NULL,
    acc = acc
  )

  expect_true(is.list(res))
  expect_true(all(c("mod.cv","acc") %in% names(res)))

  expect_length(res$mod.cv, length(cv.set))

  expect_true(all(vapply(res$mod.cv, function(m) {
    all(c("t_xp","y_pred_train","y_pred_test","x_res") %in% names(m))
  }, logical(1))))

})


test_that(".plsComponentCv accumulates predictions", {

  n <- 20
  p <- 4

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)

  cv.set <- list(1:10, 11:20)
  Ycs_fold <- list(Y, Y)

  acc <- list(
    sum_test = matrix(0, n, 1),
    n_test = rep(0, n),
    sum_train = matrix(0, n, 1),
    n_train = rep(0, n)
  )

  res <- .plsComponentCv(
    X,
    cv.set,
    Ycs_fold,
    nc = 1,
    mod.cv = NULL,
    acc = acc
  )

  expect_true(any(res$acc$n_test > 0))
  expect_true(any(res$acc$n_train > 0))

})


test_that(".plsComponentCv appends components when nc > 1", {

  n <- 20
  p <- 4

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)

  cv.set <- list(1:10, 11:20)
  Ycs_fold <- list(Y, Y)

  acc <- list(
    sum_test = matrix(0, n, 1),
    n_test = rep(0, n),
    sum_train = matrix(0, n, 1),
    n_train = rep(0, n)
  )

  res1 <- .plsComponentCv(X, cv.set, Ycs_fold, 1, NULL, acc)

  res2 <- .plsComponentCv(
    X,
    cv.set,
    Ycs_fold,
    nc = 2,
    mod.cv = res1$mod.cv,
    acc = res1$acc
  )

  expect_true(ncol(res2$mod.cv[[1]]$t_xp) == 2)

})

test_that(".plsComponentCv produces matrices with correct dimensions", {

  n <- 30
  p <- 6

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)

  cv.set <- list(1:15, 16:30)
  Ycs_fold <- list(Y, Y)

  acc <- list(
    sum_test = matrix(0, n, 1),
    n_test = rep(0, n),
    sum_train = matrix(0, n, 1),
    n_train = rep(0, n)
  )

  res <- .plsComponentCv(X, cv.set, Ycs_fold, 1, NULL, acc)

  m <- res$mod.cv[[1]]

  expect_equal(nrow(m$t_xp), n)
  expect_equal(nrow(m$x_res), n)
  expect_equal(ncol(m$x_res), p)

})


test_that(".oplsComponentCv returns valid structure for first component", {

  n <- 30
  p <- 6

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(sample(c(0,1), n, replace = TRUE), n, 1)

  cv.set <- list(1:15, 16:30)
  Ycs_fold <- list(Y, Y)

  acc <- list(
    sum_test  = matrix(0, n, 1),
    n_test    = rep(0, n),
    sum_train = matrix(0, n, 1),
    n_train   = rep(0, n)
  )

  res <- .oplsComponentCv(
    X = X,
    Ycs_fold = Ycs_fold,
    cv.set = cv.set,
    nc = 1,
    mod.cv = NULL,
    acc = acc
  )

  expect_true(is.list(res))
  expect_true(all(c("mod.cv","acc") %in% names(res)))

  expect_length(res$mod.cv, length(cv.set))

  expect_true(all(vapply(res$mod.cv, function(m) {
    all(c("t_xo","t_xp","y_pred_train","y_pred_test","x_res") %in% names(m))
  }, logical(1))))

})


test_that(".oplsComponentCv produces matrices with correct dimensions", {

  n <- 40
  p <- 8

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)

  cv.set <- list(1:20, 21:40)
  Ycs_fold <- list(Y, Y)

  acc <- list(
    sum_test  = matrix(0, n, 1),
    n_test    = rep(0, n),
    sum_train = matrix(0, n, 1),
    n_train   = rep(0, n)
  )

  res <- .oplsComponentCv(X, Ycs_fold, cv.set, 1, NULL, acc)

  m <- res$mod.cv[[1]]

  expect_equal(nrow(m$t_xo), n)
  expect_equal(nrow(m$t_xp), n)
  expect_equal(nrow(m$x_res), n)

  expect_equal(ncol(m$x_res), p)

})


test_that(".oplsComponentCv updates acc counters", {

  n <- 20
  p <- 4

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)

  cv.set <- list(1:10, 11:20)
  Ycs_fold <- list(Y, Y)

  acc <- list(
    sum_test  = matrix(0, n, 1),
    n_test    = rep(0, n),
    sum_train = matrix(0, n, 1),
    n_train   = rep(0, n)
  )

  res <- .oplsComponentCv(X, Ycs_fold, cv.set, 1, NULL, acc)

  expect_true(any(res$acc$n_test > 0))
  expect_true(any(res$acc$n_train > 0))

})


test_that(".oplsComponentCv appends components when nc > 1", {

  n <- 30
  p <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)

  cv.set <- list(1:15, 16:30)
  Ycs_fold <- list(Y, Y)

  acc <- list(
    sum_test  = matrix(0, n, 1),
    n_test    = rep(0, n),
    sum_train = matrix(0, n, 1),
    n_train   = rep(0, n)
  )

  res1 <- .oplsComponentCv(X, Ycs_fold, cv.set, 1, NULL, acc)

  res2 <- .oplsComponentCv(
    X,
    Ycs_fold,
    cv.set,
    nc = 2,
    mod.cv = res1$mod.cv,
    acc = res1$acc
  )

  expect_equal(ncol(res2$mod.cv[[1]]$t_xo), 2)
  expect_equal(ncol(res2$mod.cv[[1]]$t_xp), 2)

})


test_that(".oplsComponentCv stores predictions in correct rows", {

  n <- 20
  p <- 4

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)

  cv.set <- list(1:10, 11:20)
  Ycs_fold <- list(Y, Y)

  acc <- list(
    sum_test  = matrix(0, n, 1),
    n_test    = rep(0, n),
    sum_train = matrix(0, n, 1),
    n_train   = rep(0, n)
  )

  res <- .oplsComponentCv(X, Ycs_fold, cv.set, 1, NULL, acc)

  m <- res$mod.cv[[1]]

  expect_true(any(m$y_pred_test[-cv.set[[1]], ] != 0))
  expect_true(any(m$y_pred_train[cv.set[[1]], ] != 0))

})
