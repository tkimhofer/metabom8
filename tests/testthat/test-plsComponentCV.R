# test-plsComponentCv.R

test_that(".plsComponentCv produces valid outputs", {
  set.seed(123)
  X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  Y <- model.matrix(~ factor(rep(c("A", "B"), each = 50)))  # dummy matrix

  cv.set <- split(1:100, rep(1:5, each = 20))  # 5-fold CV

  out <- .plsComponentCv(X, Y, cv.set = cv.set, nc = 1, mod.cv = NULL)

  expect_type(out, "list")
  expect_equal(length(out), length(cv.set))
  expect_true(all(sapply(out, function(x) is.list(x))))
  expect_true(all(sapply(out, function(x) "t_xp" %in% names(x))))
  expect_true(all(sapply(out, function(x) nrow(x$t_xp) == 100)))
})
