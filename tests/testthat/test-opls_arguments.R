test_that("opls cv-method kfold", {
  data(iris)

  n=nrow(iris)
  m=ncol(iris)-1

  smod=opls(X=iris[,1:4], Y=iris[,5], cv=list(method='k-fold',  k=7, split=2/3), plotting = FALSE)

  expect_equal(class(smod)[1], 'OPLS_metabom8')
  expect_equal(nrow(smod@t_pred), n)
  expect_equal(nrow(smod@t_pred_cv), n)

  expect_equal(nrow(smod@t_orth), n)
  expect_equal(nrow(smod@t_orth_cv), n)

  expect_equal(ncol(smod@p_pred), m)
  expect_equal(nrow(smod@w_pred), m) # this need change in later versions

  expect_equal(ncol(smod@p_orth), m)


  expect_false(any(is.na(smod@t_pred)))
  expect_false(any(is.na(smod@t_orth)))

  expect_false(any(is.na(smod@t_orth_cv)))
  expect_false(any(is.na(smod@t_pred_cv)))

  expect_false(any(is.na(smod@p_pred)))
  expect_false(any(is.na(smod@p_orth)))

  expect_false(any(is.na(smod@Y$dummy)))

})



test_that("opls cv-method loocv", {
  data(iris)

  n=nrow(iris)
  m=ncol(iris)-1

  smod=opls(X=iris[,1:4], Y=iris[,5], cv=list(method='k-fold', k=n, split=2/3), plotting = FALSE)

  expect_equal(class(smod)[1], 'OPLS_metabom8')
  expect_equal(nrow(smod@t_pred), n)
  expect_equal(nrow(smod@t_pred_cv), n)

  expect_equal(nrow(smod@t_orth), n)
  expect_equal(nrow(smod@t_orth_cv), n)

  expect_equal(ncol(smod@p_pred), m)
  expect_equal(nrow(smod@w_pred), m) # this need change in later versions

  expect_equal(ncol(smod@p_orth), m)


  expect_false(any(is.na(smod@t_pred)))
  expect_false(any(is.na(smod@t_orth)))

  expect_false(any(is.na(smod@t_orth_cv)))
  expect_false(any(is.na(smod@t_pred_cv)))

  expect_false(any(is.na(smod@p_pred)))
  expect_false(any(is.na(smod@p_orth)))

  expect_false(any(is.na(smod@Y$dummy)))

})



test_that("opls cv-method k-fold_stratified", {
  data(iris)

  n=nrow(iris)
  m=ncol(iris)-1

  smod=opls(X=iris[,1:4], Y=iris[,5], cv=list(method='k-fold_stratified', k=7, split=2/3), plotting = FALSE)

  expect_equal(class(smod)[1], 'OPLS_metabom8')
  expect_equal(nrow(smod@t_pred), n)
  expect_equal(nrow(smod@t_pred_cv), n)

  expect_equal(nrow(smod@t_orth), n)
  expect_equal(nrow(smod@t_orth_cv), n)

  expect_equal(ncol(smod@p_pred), m)
  expect_equal(nrow(smod@w_pred), m) # this need change in later versions

  expect_equal(ncol(smod@p_orth), m)


  expect_false(any(is.na(smod@t_pred)))
  expect_false(any(is.na(smod@t_orth)))

  expect_false(any(is.na(smod@t_orth_cv)))
  expect_false(any(is.na(smod@t_pred_cv)))

  expect_false(any(is.na(smod@p_pred)))
  expect_false(any(is.na(smod@p_orth)))

  expect_false(any(is.na(smod@Y$dummy)))

})



test_that("opls cv-method k-fold_stratified, k too high/low", {
  data(iris)

  n=nrow(iris)
  m=ncol(iris)-1

  k_high=nrow(iris)-10
  smod=opls(X=iris[,1:4], Y=iris[,5], cv=list(method='k-fold_stratified', k=k_high, split=2/3), plotting = FALSE)
  expect_true(smod@Parameters$cv$k<(k_high))

  k_high=nrow(iris)+10
  smod=opls(X=iris[,1:4], Y=iris[,5], cv=list(method='k-fold_stratified', k=k_high, split=2/3), plotting = FALSE)
  expect_true(smod@Parameters$cv$k<(k_high))

  k_low=0
  smod=opls(X=iris[,1:4], Y=iris[,5], cv=list(method='k-fold_stratified', k=k_low, split=2/3), plotting = FALSE)
  expect_true(smod@Parameters$cv$k>(k_low))



})


test_that("opls cv-method k-fold, k too high/low", {
  data(iris)

  n=nrow(iris)
  m=ncol(iris)-1

  k_high=n+300
  smod=opls(X=iris[,1:4], Y=iris[,5], cv=list(method='k-fold', k=k_high, split=2/3), plotting = FALSE)
  expect_true(smod@Parameters$cv$k<(k_high))

  k_low=0
  smod=opls(X=iris[,1:4], Y=iris[,5], cv=list(method='k-fold', k=k_high, split=2/3), plotting = FALSE)
  expect_true(smod@Parameters$cv$k==n)

})






test_that("opls cv-method MC, split too low/high", {
  data(iris)
  n=nrow(iris)
  m=ncol(iris)-1

  expect_error(opls(X=iris[,1:4], Y=iris[,5], cv=list(method='MC', k=10, split=-1.5), plotting = FALSE))
  expect_error(opls(X=iris[,1:4], Y=iris[,5], cv=list(method='MC', k=10, split=1.5), plotting = FALSE))
})


test_that("opls cv-method MC", {
  data(iris)

  n=nrow(iris)
  m=ncol(iris)-1

  smod=opls(X=iris[,1:4], Y=iris[,5], cv=list(method='MC', k=10, split=2/3), plotting = FALSE)

  expect_equal(class(smod)[1], 'OPLS_metabom8')
  expect_equal(nrow(smod@t_pred), n)
  expect_equal(nrow(smod@t_pred_cv), n)

  expect_equal(nrow(smod@t_orth), n)
  expect_equal(nrow(smod@t_orth_cv), n)

  expect_equal(ncol(smod@p_pred), m)
  expect_equal(nrow(smod@w_pred), m) # this need change in later versions

  expect_equal(ncol(smod@p_orth), m)


  expect_false(any(is.na(smod@t_pred)))
  expect_false(any(is.na(smod@t_orth)))

  #expect_false(any(is.na(smod@t_orth_cv)))
  #expect_false(any(is.na(smod@t_pred_cv)))

  expect_false(any(is.na(smod@p_pred)))
  expect_false(any(is.na(smod@p_orth)))

  expect_false(any(is.na(smod@Y$dummy)))

})




test_that("opls cv-method MC_balanced", {
  data(iris)

  n=nrow(iris)
  m=ncol(iris)-1

  smod=opls(X=iris[,1:4], Y=iris[,5], cv=list(method='MC_balanced', k=10, split=2/3), plotting = FALSE)

  expect_equal(class(smod)[1], 'OPLS_metabom8')
  expect_equal(nrow(smod@t_pred), n)
  expect_equal(nrow(smod@t_pred_cv), n)

  expect_equal(nrow(smod@t_orth), n)
  expect_equal(nrow(smod@t_orth_cv), n)

  expect_equal(ncol(smod@p_pred), m)
  expect_equal(nrow(smod@w_pred), m) # this need change in later versions

  expect_equal(ncol(smod@p_orth), m)


  expect_false(any(is.na(smod@t_pred)))
  expect_false(any(is.na(smod@t_orth)))

  #expect_false(any(is.na(smod@t_orth_cv)))
  #expect_false(any(is.na(smod@t_pred_cv)))

  expect_false(any(is.na(smod@p_pred)))
  expect_false(any(is.na(smod@p_orth)))

  expect_false(any(is.na(smod@Y$dummy)))

})





