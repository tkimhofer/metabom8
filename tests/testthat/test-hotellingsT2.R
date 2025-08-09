test_that(".hotellingsT2 returns correct ellipse coordinates", {
  set.seed(123)
  x <- rnorm(100)
  y <- 0.5 * x + rnorm(100, sd = 0.5)

  df <- .hotellingsT2(x, y, alpha = 0.95)

  expect_true(is.data.frame(df))
  expect_named(df, c("V1", "V2"))
  expect_true(nrow(df) > 0)
  expect_true(all(is.finite(df$V1)))
  expect_true(all(is.finite(df$V2)))
})
