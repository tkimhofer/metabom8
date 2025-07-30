test_that("eruption returns valid plot and data frame", {
  skip_if_not_installed("metabom8")

  data(covid, package = "metabom8")
  model <- opls(X, Y = factor(an$type), plotting=FALSE)

  result <- eruption(model)

  expect_type(result, "list")
  expect_named(result, c("data", "plot"))

  expect_s3_class(result$data, "data.frame")
  expect_s3_class(result$plot, "ggplot")

  expect_true(all(c("Cd", "pval") %in% colnames(result$data)))
})
