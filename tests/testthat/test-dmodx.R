test_that("dmodx returns expected structure", {
  data(covid)
  model <- opls(X, Y = an$type, plotting = FALSE)

  result <- dmodx(model, plot = FALSE)

  expect_s3_class(result, "data.frame")
  expect_named(result, c("ID", "DmodX", "passedT.test"))
  expect_equal(nrow(result), nrow(X))
  expect_type(result$DmodX, "double")
  expect_type(result$passedT.test, "logical")
})

test_that("dmodx handles incorrect input gracefully", {
  expect_error(dmodx(list()), "OPLS_metabom8")
})

test_that("dmodx returns the same values with plot = TRUE", {
  data(covid)
  model <- opls(X, Y = an$type, plotting = FALSE)

  expect_silent({
    result <- dmodx(model, plot = TRUE)
  })

  expect_s3_class(result, "data.frame")
  expect_named(result, c("ID", "DmodX", "passedT.test", "Group"))
})

test_that("dmodx DmodX values are non-negative", {
  data(covid)
  model <- opls(X, Y = an$type, plotting = FALSE)

  result <- dmodx(model, plot = FALSE)
  expect_true(all(result$DmodX >= 0))
})
