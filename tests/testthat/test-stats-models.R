set.seed(15)
### pca

test_that("pca() runs and returns valid PCA_metabom8 object", {

  X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  colnames(X) <- paste0("V", 1:10)
  rownames(X) <- paste0("S", 1:100)

  sc <- UVScaling(center = TRUE)
  model <- pca(X, ncomp = 2, scaling = sc, method = "nipals")

  Tx <- scores(model)
  Px <- loadings(model)

  expect_s4_class(model, "m8_model")

  expect_true(is.matrix(Tx))
  expect_true(is.matrix(Px))
  expect_equal(nrow(Tx), nrow(X))
  expect_equal(ncol(Tx), 2)
  expect_equal(ncol(Px), ncol(X))

  expect_true(all(c("method","type","ncomp_selected","ncomp_tested") %in%
                    names(model@ctrl)))

  R <- cor(Tx)
  expect_true(all(abs(R[lower.tri(R)]) < 1e-3))

  expect_true(all(model@ctrl$r2x_comp >= 0))
  expect_true(all(model@ctrl$r2x_comp <= 1))
})


test_that("pca() handles too many components gracefully", {
  X <- matrix(rnorm(30), nrow = 5)

  sc <- UVScaling(center = TRUE)
  expect_warning(
    model <- pca(X, ncomp = 10, scaling = sc, method = "nipals"),
    regexp = "Too many"
    )

  expect_s4_class(model, "m8_model")
  expect_lte(model@ctrl$ncomp_selected, min(nrow(X)-1, ncol(X)))
})


# pls
test_that("PLS-DA model runs correctly with numeric input", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  Y <- iris$Species

  sc <- UVScaling(center = TRUE)
  kf <- kfold(3)
  model <- pls(X, Y, scaling = sc, validation_strategy = kf)

  Tx <- scores(model)
  Px <- loadings(model)
  Wx <- weights(model)

  expect_s4_class(model, "m8_model")

  expect_true(is.matrix(Tx))
  expect_true(is.matrix(Px))
  expect_true(is.matrix(Wx))

  expect_equal(nrow(Tx), nrow(X))
  expect_equal(ncol(Px), ncol(X))
  expect_equal(ncol(Wx), ncol(X))

  expect_equal(ncol(Tx), model@ctrl$ncomp_selected)
  expect_equal(nrow(Px), model@ctrl$ncomp_selected)
  expect_equal(nrow(Wx), model@ctrl$ncomp_selected)

  expect_true(all(c("type","ncomp_selected","ncomp_tested", "stop_reason") %in%
                    names(model@ctrl)))

  expect_true(model@ctrl$type %in% c("DA","DA-mY"))

})


test_that("PLS regression model runs correctly", {

  n <- 80
  p <- 10

  X <- matrix(rnorm(n*p), n, p)

  beta <- c(1.2, -0.8, 0.5, rep(0, p-3))
  Y <- X %*% beta + rnorm(n, sd = 0.2)

  sc <- UVScaling(center = TRUE)
  kf <- kfold(3)

  model <- pls(X, Y, scaling = sc, validation_strategy = kf)

  Tx <- scores(model)
  Px <- loadings(model)
  Wx <- weights(model)

  expect_s4_class(model, "m8_model")

  expect_true(is.matrix(Tx))
  expect_true(is.matrix(Px))
  expect_true(is.matrix(Wx))

  expect_equal(nrow(Tx), nrow(X))
  expect_equal(ncol(Px), ncol(X))
  expect_equal(ncol(Wx), ncol(X))

  expect_equal(ncol(Tx), model@ctrl$ncomp_selected)
  expect_true(model@ctrl$q2[1] > 0)
})


test_that("OPLS-DA model runs correctly", {

  data(iris)

  X <- as.matrix(iris[,1:4])
  Y <- iris$Species

  sc <- UVScaling(center = TRUE)
  kf <- kfold(3)

  model <- opls(X, Y, scaling = sc, validation_strategy = kf)

  Tp <- scores(model)
  To <- scores(model, orth = TRUE)

  Pp  <- loadings(model)
  Po  <- loadings(model, orth = TRUE)
  Wx  <- weights(model)

  expect_s4_class(model, "m8_model")

  expect_true(is.matrix(Tp))
  expect_true(is.matrix(Pp))
  expect_true(is.matrix(Wx))

  expect_equal(nrow(Tp), nrow(X))
  expect_equal(ncol(Pp), ncol(X))
  expect_equal(ncol(Wx), ncol(X))

  expect_equal(ncol(Tp), model@ctrl$ncomp_selected-1)

  expect_true(all(c("type","ncomp_selected","ncomp_tested","stop_reason") %in%
                    names(model@ctrl)))

})


test_that("OPLS pred component correlates with Y, orthogonal comp is independent of Y", {

  data(iris)

  X <- as.matrix(iris[,1:4])
  Y <- iris$Species

  sc <- UVScaling(center = TRUE)
  kf <- kfold(3)

  model <- opls(X, Y, scaling = sc, validation_strategy = kf)
  Tp <- scores(model)

  Ynum <- as.numeric(Y)
  expect_true(abs(cor(Tp[,1], Ynum)) > 0.2)

  Txo <- scores(model, orth = TRUE)

  if (!is.null(Txo)) {
    r <- cor(Txo[,1], Ynum)
    expect_true(abs(r) < 0.3)
  }

})
