test_that("stocsy output structure", {
  data(covid)
  X <- covid$X
  # an <- covid$an
  ppm <- covid$ppm
  mod <- stocsy(X, ppm, driver = 5.233, plotting = FALSE)
  expect_s4_class(mod, "stocsy1d_metabom8")
  expect_true(all(c("r", "cov", "ppm", "driver") %in% slotNames(mod)))
})

test_that("plotStocsy returns ggplot", {
  data(covid)
  X <- covid$X
  # an <- covid$an
  ppm <- covid$ppm
  mod <- stocsy(X, ppm, driver = 5.233, plotting = FALSE)
  p <- plotStocsy(mod, shift = c(5.1, 5.3))
  expect_s3_class(p, "gg")
})
