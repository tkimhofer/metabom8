### .scaleMatRcpp

test_that("mean is computed on training rows only", {

  X <- matrix(1:20, nrow = 5)
  idc <- 0:2   # rows 1,2,3 in R

  res <- .scaleMatRcpp(X, idc, center = TRUE, scale_type = 0)

  m_train <- colMeans(X[idc+1,,drop=FALSE])

  expect_equal(res$mean, m_train)

})

test_that("UV scaling matches base R", {

  set.seed(1)
  X <- matrix(rnorm(100), nrow = 10)
  idc <- 0:6

  res <- .scaleMatRcpp(X, idc, center = TRUE, scale_type = 1)

  m <- colMeans(X[idc+1,,drop=FALSE])
  s <- apply(X[idc+1,,drop=FALSE],2,sd)
  s[s < 1e-12] <- 1

  X_ref <- sweep(X,2,m,"-")
  X_ref <- sweep(X_ref,2,s,"/")

  expect_equal(res$X_prep, X_ref, tolerance = 1e-10)

})

test_that("Pareto scaling matches definition", {

  set.seed(1)
  X <- matrix(rnorm(100), nrow = 10)
  idc <- 0:6

  res <- .scaleMatRcpp(X, idc, center = TRUE, scale_type = 2)

  m <- colMeans(X[idc+1,,drop=FALSE])
  s <- apply(X[idc+1,,drop=FALSE],2,sd)

  denom <- sqrt(s)
  denom[denom < 1e-6] <- 1

  X_ref <- sweep(X,2,m,"-")
  X_ref <- sweep(X_ref,2,denom,"/")

  expect_equal(res$X_prep, X_ref, tolerance = 1e-10)

})

test_that("constant column becomes zero after UV scaling", {

  X <- cbind(rnorm(10), rep(5,10))
  idc <- 0:7

  res <- .scaleMatRcpp(X, idc, center = TRUE, scale_type = 1)

  expect_true(all(abs(res$X_prep[,2]) < 1e-12))

})


test_that("n_train = 1 is safe", {

  X <- matrix(rnorm(50), nrow = 5)
  idc <- 0L

  res <- .scaleMatRcpp(X, idc, center = TRUE, scale_type = 1)

  expect_false(any(is.nan(res$X_prep)))
  expect_false(any(is.infinite(res$X_prep)))

})

test_that("scaled matrix has same dimensions as input", {

  X <- matrix(rnorm(100), nrow = 10)
  idc <- 0:5

  res <- .scaleMatRcpp(X, idc, TRUE, 1)

  expect_equal(dim(res$X_prep), dim(X))

})

test_that("mean and sd have correct length", {

  X <- matrix(rnorm(100), nrow = 10)
  idc <- 0:5

  res <- .scaleMatRcpp(X, idc, TRUE, 1)

  expect_equal(length(res$mean), ncol(X))
  expect_equal(length(res$sd),   ncol(X))

})

test_that("single column remains matrix", {

  X <- matrix(rnorm(10), ncol = 1)
  idc <- 0:5

  res <- .scaleMatRcpp(X, idc, TRUE, 1)

  expect_true(is.matrix(res$X_prep))
  expect_equal(dim(res$X_prep), dim(X))

})

test_that("constant column does not drop dimensions", {

  X <- cbind(rnorm(10), rep(5,10))
  idc <- 0:5

  res <- .scaleMatRcpp(X, idc, TRUE, 1)

  expect_true(is.matrix(res$X_prep))
  expect_equal(ncol(res$X_prep), 2)

})

test_that("single row remains matrix", {

  X <- matrix(rnorm(5), nrow = 1)
  idc <- 0L

  res <- .scaleMatRcpp(X, idc, TRUE, 1)

  expect_true(is.matrix(res$X_prep))
  expect_equal(dim(res$X_prep), dim(X))

})

test_that("test rows use training estimates", {

  X <- matrix(rnorm(50), nrow = 5)
  idc <- 0:2

  res <- .scaleMatRcpp(X, idc, center = TRUE, scale_type = 1)

  m <- colMeans(X[idc+1,,drop=FALSE])
  s <- apply(X[idc+1,,drop=FALSE],2,sd)
  s[s < 1e-12] <- 1

  i_test <- 5

  expect_equal(
    res$X_prep[i_test,],
    (X[i_test,] - m)/s,
    tolerance = 1e-10
  )

})


### .sdRcpp

test_that("sd_rcpp agrees with base R mean/sd", {

  set.seed(1)
  X <- matrix(rnorm(1000), nrow = 100)

  res <- .sdRcpp(X)

  expect_equal(res$mean, colMeans(X), tolerance = 1e-12)
  expect_equal(res$sd,   apply(X, 2, sd), tolerance = 1e-12)

})

test_that("sd_rcpp handles constant columns", {

  X <- matrix(5, nrow = 10, ncol = 3)

  res <- .sdRcpp(X)

  expect_equal(res$mean, rep(5, 3))
  expect_equal(res$sd,   rep(0, 3))

})

test_that("sd_rcpp handles n = 1 correctly", {

  X <- matrix(c(1,2,3), nrow = 1)

  res <- .sdRcpp(X)

  expect_equal(res$mean, c(1,2,3))
  expect_true(all(is.na(res$sd)))

})

test_that("sd_rcpp is numerically stable for large values", {

  set.seed(1)
  X <- matrix(1e6 + rnorm(1000, sd = 1e-3), nrow = 100)

  res <- .sdRcpp(X)

  expect_equal(res$mean, colMeans(X), tolerance = 1e-8)
  expect_equal(res$sd,   apply(X, 2, sd), tolerance = 1e-8)

})


test_that("sd_rcpp returns finite output when possible", {

  X <- matrix(rnorm(1000), nrow = 100)

  res <- .sdRcpp(X)

  expect_false(any(is.nan(res$mean)))
  expect_false(any(is.infinite(res$mean)))
  expect_false(any(is.nan(res$sd)))
  expect_false(any(is.infinite(res$sd)))

})

test_that("scaling without centering works", {

  set.seed(1)
  X <- matrix(rnorm(100), nrow = 10)
  idc <- 0:6

  res <- .scaleMatRcpp(X, idc, center = FALSE, scale_type = 1)

  s <- apply(X[idc+1,,drop=FALSE],2,sd)
  s[s < 1e-12] <- 1

  X_ref <- sweep(X,2,s,"/")

  expect_equal(res$X_prep, X_ref, tolerance = 1e-10)

})

test_that("center only scaling works", {

  set.seed(1)
  X <- matrix(rnorm(100), nrow = 10)
  idc <- 0:6

  res <- .scaleMatRcpp(X, idc, center = TRUE, scale_type = 0)

  m <- colMeans(X[idc+1,,drop=FALSE])

  X_ref <- sweep(X,2,m,"-")

  expect_equal(res$X_prep, X_ref, tolerance = 1e-10)

})
