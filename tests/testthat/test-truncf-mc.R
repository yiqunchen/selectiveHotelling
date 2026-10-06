test_that("truncated-F Monte Carlo methods return valid p-values", {
  set.seed(11)
  X <- matrix(rnorm(120), 60, 2)
  km <- kmeans_with_trajectory(X, k = 3); cl <- km$cluster; init <- km$init_centers_idx
  cf <- function(x) kmeans_with_fixed_init(x, init)$cluster
  b <- selective_truncated_F(X, 1, 2, clustering = "custom", method = "bc", cl_fun = cf,
                             cl_obs = cl, M = 50L)
  i <- selective_truncated_F(X, 1, 2, clustering = "custom", method = "is", cl_fun = cf,
                             cl_obs = cl, n_is = 200L)
  expect_true(b$pval >= 0 && b$pval <= 1)
  expect_true(i$pval >= 0 && i$pval <= 1)
  expect_equal(length(b$draws), 51L)
  expect_true(all(b$draws > 0 & b$draws < pi / 2))
})
