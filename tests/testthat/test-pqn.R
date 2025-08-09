test_that("pqn normalises data as expected", {
  set.seed(42)
  X <- matrix(runif(100), nrow = 10)
  rownames(X) <- paste0("S", 1:10)

  # Basic PQN
  Xn <- pqn(X, add_DilF = NULL)
  expect_equal(dim(Xn), dim(X))
  expect_true(all(!is.na(Xn)))

  # Check that dilution factors scale input
  df <- NULL
  Xn2 <- pqn(X, add_DilF = "df")
  expect_true(exists("df", envir = .GlobalEnv))
  expect_equal(length(get("df", envir = .GlobalEnv)), nrow(X))

  # Clean up
  rm(df, envir = .GlobalEnv)
})
