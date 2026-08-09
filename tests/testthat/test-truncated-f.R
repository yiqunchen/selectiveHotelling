# =============================================================================
# selective_truncated_F: closed-form one-dimensional truncated-F selective test.
# Exact k-means (Lloyd trajectory) + exact hierarchical (pair-preserved) engines,
# and the algorithm-agnostic grid engine.
# =============================================================================

hclust_cl_fun <- function(K = 3L, method = "average", square = TRUE) {
  force(K); force(method); force(square)
  function(Xmat) {
    d <- stats::dist(Xmat); if (isTRUE(square)) d <- d^2
    as.integer(stats::cutree(stats::hclust(d, method = method), k = K))
  }
}
pick_big_pair <- function(cl) {
  sizes <- table(cl)
  big <- as.integer(names(sizes)[sizes >= 3L])
  if (length(big) < 2L) return(NULL)
  list(k1 = big[1], k2 = big[2])
}

test_that("exact hclust route: valid p-value, S3 object, exact-vs-grid agree", {
  set.seed(101)
  X <- matrix(stats::rnorm(30 * 3), 30, 3)
  cl <- hclust_cl_fun(3L)(X)
  pk <- pick_big_pair(cl); skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))

  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "hclust", K = 3L)
  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval >= 0 && res$pval <= 1)
  expect_identical(res$condition_on, "pair-preserved")
  expect_identical(res$sampler, "truncated-F")
  expect_true(is.na(res$se))                 # exact: no Monte Carlo error
  expect_identical(res$df1, ncol(X))
  expect_true(res$n_intervals >= 1L)

  # exact closed-form arc and grid engine condition on the SAME event => agree
  res_g <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "hclust", K = 3L,
                               method = "grid", n_grid = 600L)
  expect_equal(res$pval, res_g$pval, tolerance = 5e-3)

  expect_output(print(res), "truncated-F")
  expect_output(summary(res), "Valid draws")
})

test_that("exact k-means route: valid p-value and T2 = (r-2) cot^2(theta) identity", {
  set.seed(303)
  X <- matrix(stats::rnorm(45 * 3), 45, 3)
  init_idx <- c(1L, 16L, 31L)
  ref <- kmeans_with_fixed_init(X, init_centers_idx = init_idx)
  cl <- ref$cluster
  pk <- pick_big_pair(cl); skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))

  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "kmeans", K = 3L,
                             init_centers_idx = init_idx, ridge = 0)
  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval >= 0 && res$pval <= 1)
  expect_identical(res$condition_on, "kmeans-trajectory")

  r <- res$df2 + res$df1 + 1L                # = r (pair size)
  T_from_angle <- (r - 2) * (cos(res$theta_obs) / sin(res$theta_obs))^2
  expect_equal(res$T_obs, T_from_angle, tolerance = 1e-6)
})

test_that("trivial (always-feasible) selection reduces to the classical F p-value", {
  set.seed(202)
  X <- matrix(stats::rnorm(30 * 3), 30, 3)
  cl <- hclust_cl_fun(3L)(X)
  pk <- pick_big_pair(cl); skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))

  const_cl <- function(Xmat) cl          # ignores its input => pair always preserved
  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "custom",
                             cl_fun = const_cl, cl_obs = cl)
  expect_equal(res$n_intervals, 1L)
  expect_equal(res$pval, res$p_uncond, tolerance = 1e-6)
})

test_that("custom cl_fun (grid) works and is algorithm-agnostic", {
  set.seed(404)
  X <- matrix(stats::rnorm(45 * 3), 45, 3)
  init_idx <- c(1L, 16L, 31L)
  km_cl_fun <- function(Xmat) kmeans_with_fixed_init(Xmat, init_centers_idx = init_idx)$cluster
  cl <- km_cl_fun(X)
  pk <- pick_big_pair(cl); skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))
  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "custom",
                             cl_fun = km_cl_fun, cl_obs = cl)
  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval >= 0 && res$pval <= 1)
})

test_that("errors when the pair is too small relative to dimension (r-1 <= p)", {
  set.seed(11)
  X <- matrix(stats::rnorm(12 * 6), 12, 6)
  cl <- c(rep(1L, 3), rep(2L, 3), rep(3L, 6))
  const_cl <- function(Xmat) cl
  expect_error(selective_truncated_F(X, 1L, 2L, clustering = "custom",
                                   cl_fun = const_cl, cl_obs = cl), "r-1 > p")
})

test_that("small exactness sanity: exact hclust null p-values ~ uniform", {
  set.seed(7)
  reps <- 120L
  clf <- hclust_cl_fun(3L)
  pvals <- rep(NA_real_, reps)
  for (i in seq_len(reps)) {
    X <- matrix(stats::rnorm(30 * 3), 30, 3)
    cl <- clf(X)
    pk <- pick_big_pair(cl)
    if (is.null(pk)) next
    if (!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2))) next
    res <- tryCatch(selective_truncated_F(X, pk$k1, pk$k2, clustering = "hclust",
                                        K = 3L, cl_obs = cl),
                    error = function(e) NULL)
    if (!is.null(res)) pvals[i] <- res$pval
  }
  ok <- pvals[is.finite(pvals)]
  expect_gt(length(ok), 30L)
  expect_gt(mean(ok), 0.34)
  expect_lt(mean(ok), 0.66)
  expect_lt(mean(ok <= 0.05), 0.16)
})
