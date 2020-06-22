test_that("cv kfold", {
  Y=cbind(seq_len(12))
  cv_idc=.kFold(3, Y)
  cv_test=lapply(cv_idc, function(idc, k=3, y_le=12){
    which(!seq_len(y_le) %in% idc)
  })
  cv_test=unique(unlist(cv_test))
  expect_equal(length(cv_test), 12)
  expect_true(all(cv_test %in% 1:12))
})


test_that("cv kfoldStratified (num Y)", {

  strat=list(type='R', Y=cbind(seq_len(12)), probs=c(0, 0.33, 0.66, 1))

  cv_idc=.kFoldStratified(k=3, stratified = strat)
  #browser()
  #print(cv_idc)
  cv_test=lapply(cv_idc, function(idc, k=3, y_le=12){
    which(!seq_len(y_le) %in% idc)
  })

  cv_test=unique(unlist(cv_test))
  #print((cv_test))
  expect_equal(length(cv_test), 12)
  expect_true(all(cv_test %in% 1:12))
})


#
# test_that("cv kfold_stratified (cat Y)", {
#   Y=cbind(factor(c(rep('A',6), rep('B', 6))))
#   cv_idc=.kFold(3, Y)
#   cv_test=lapply(cv_idc, function(idc, k=3, y_le=12){
#     which(!seq_len(y_le) %in% idc)
#   })
#   cv_test=unique(unlist(cv_test))
#   expect_equal(length(cv_test), 12)
#   expect_true(all(cv_test %in% 1:12))
# })



test_that("cv kfold_stratified (cat Y)", {
  Y=cbind(factor(rep(c(1,2), each=20)))
  n=nrow(Y)
  k=2
  lev_min=floor((min(table(Y))/k))

  if(lev_min>k) lev_min = (lev_min-(k/2))


  expect_true(lev_min>0)

  cv_idc=.kFoldStratified(2, list('DA', Y, c(1,1)))
  expect_length(cv_idc, k)

  cv_test=lapply(cv_idc, function(idc, k=2, y_le=n){
    which(!seq_len(y_le) %in% idc)
  })

  cv_ele=unique(unlist(cv_test))
  expect_equal(length(cv_ele), n) # check le
  expect_true(all(cv_ele %in% seq(n))) # check indices

  #check if each cv set has a minimum % of each group in it
  min_ele=sapply(cv_idc, function(idc, y=Y[,1], lmin=lev_min){
    y_cv=y[idc]
    all(unique(y) %in% y_cv) &  all(table(y_cv)>=lmin)
  })
  expect_true(all(min_ele))
})



test_that("cv mc", {
  Y=cbind(factor(c(rep('A',6), rep('B', 6))))
  cv_idc=.mc(k=10L, Y, split=2/3)
  cv_test=lapply(cv_idc, function(idc, k=3, y_le=12){
    length(idc)
  })
  expect_equal(nrow(Y)*2/3, unique(unlist(cv_test)))
})

test_that("cv mcBalanced", {
  Y=cbind(factor(c(rep('A',6), rep('B', 6))))
  cv_idc=.mc(k=10L, Y, split=2/3)
  cv_test=lapply(cv_idc, function(idc, k=3, y_le=12){
    length(idc)
  })
  expect_equal(nrow(Y)*2/3, unique(unlist(cv_test)))
})

