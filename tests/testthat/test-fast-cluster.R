# =============================================================================
# Fast clustering routines (R/fast_cluster.R): the same clusterings as the
# package routines, and the Monte Carlo entry points use them.
# =============================================================================

test_that("fast routines reproduce kmeans_with_fixed_init() and hclust_fun()", {
  for (s in 1:40) {
    set.seed(s); n <- sample(c(30, 90), 1); p <- sample(c(3, 200), 1)
    X <- matrix(stats::rnorm(n * p), n, p) + outer(rep(c(0, 2, 4), length.out = n), c(1, rep(0, p - 1)))
    init <- sample.int(n, 3)
    expect_identical(make_fast_kmeans_cl_fun(init)(X), as.integer(kmeans_with_fixed_init(X, init)$cluster))
    expect_identical(make_fast_hclust_cl_fun(3)(X), hclust_fun(X, 3, square = TRUE)$cluster)
  }
})

test_that("fast routines agree on separated, tied-free and high-dimensional data", {
  set.seed(77)
  for (p in c(2, 50, 600)) {
    X <- matrix(stats::rnorm(120 * p), 120, p) + outer(rep(c(0, 4, 8), each = 40), c(1, rep(0, p - 1)))
    init <- c(1, 41, 81)
    expect_identical(make_fast_kmeans_cl_fun(init)(X), as.integer(kmeans_with_fixed_init(X, init)$cluster))
    expect_identical(make_fast_hclust_cl_fun(3)(X), hclust_fun(X, 3, square = TRUE)$cluster)
    expect_identical(make_hclust_cl_fun(3, "average", TRUE)(X), hclust_fun(X, 3, square = TRUE)$cluster)
  }
})

test_that("selective_hotelling_test gives the p-value of the original k-means route", {
  set.seed(5); X <- matrix(stats::rnorm(90 * 3), 90, 3) + outer(rep(c(0, 2.5, 5), each = 30), c(1, 0, 0))
  set.seed(105); a <- selective_hotelling_test(X, 1, 2, K = 3, M = 100)
  set.seed(105); km <- kmeans_with_trajectory(X, k = 3)
  b <- selective_hotelling_bc(X, function(x) kmeans_with_fixed_init(x, km$init_centers_idx)$cluster,
                              1, 2, M = 100, cl_obs = as.integer(km$cluster))
  expect_identical(a$pval, b$pval)
})
