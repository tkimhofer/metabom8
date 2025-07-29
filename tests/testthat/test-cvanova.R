test_that("cvanova returns expected structure and valid p-value", {
  skip_if_not_installed("metabom8")

  data(covid, package = "metabom8")
  model <- opls(X, Y = factor(an$type))

  result <- cvanova(model)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
  expect_equal(colnames(result), c("SS", "DF", "MS", "F_value", "p_value"))
  expect_true(is.na(result$p_value[1]))
  expect_false(is.na(result$p_value[4]))
})
