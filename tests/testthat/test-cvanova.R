test_that("cvanova returns expected structure and valid p-value", {
  skip_if_not_installed("metabom8")

  # print(ls())
  # message("class of X: ", class(X))
  # message("typeof X: ", typeof(X))
  # print(class(X))   # should be "matrix" or "data.frame"
  # print(typeof(X))  # should be "double" or "integer"
  data(covid, package = "metabom8")
  # print(ls())
  # print(class(X))   # should be "matrix" or "data.frame"
  # print(typeof(X))  # should be "double" or "integer"
  X_data <- as.matrix(X)

  # print(class(X))
  # print(dim(X))
  # print(sapply(X, class))  # if data.frame
  # X <- as.matrix(X)
  # print(class(X))          # should be "matrix"
  # print(mode(X))

  model <- opls(X_data, Y = factor(an$type), plotting=FALSE)

  result <- cvanova(model)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
  expect_equal(colnames(result), c("SS", "DF", "MS", "F_value", "p_value"))
  expect_true(is.na(result$p_value[1]))
  expect_false(is.na(result$p_value[4]))
})
