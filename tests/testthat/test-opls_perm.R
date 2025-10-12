# test_that("opls_perm returns valid output", {
#   set.seed(123)
#   n <- 10; p <- 20
#   t_pred <- rnorm(n,0,2)
#   p_pred <- c(runif(5,1,2), rep(0,p-5))
#   t_orth <- rnorm(n)
#   p_orth <- c(rep(0,5), runif(5,0.5,1), rep(0,p-10))
#   X <- outer(t_pred,p_pred) + outer(t_orth,p_orth) + matrix(rnorm(n*p,0,0.1), n, p)
#   y <- t_pred + rnorm(n,0,0.1)
#
#   model <- opls(X, y, maxPCo = 1, plotting = FALSE)
#   res <- opls_perm(model, n = 3, plot = FALSE)
#   # res <- opls_perm(model, n = 300)
#   expect_s3_class(res, "data.frame")
#   expect_true("Non-permuted" %in% res$model)
#   expect_equal(nrow(res), 4)  # 3 permutations + 1 original
# })
