# =============================================================================
# selective_campbellF: closed-form one-dimensional truncated-F selective test.
# =============================================================================

hclust_cl_fun <- function(K = 3L, method = "average") {
  force(K); force(method)
  function(Xmat) {
    hc <- stats::hclust(stats::dist(Xmat), method = method)
    as.integer(stats::cutree(hc, k = K))
  }
}

pick_big_pair <- function(cl) {
  sizes <- table(cl)
  big <- as.integer(names(sizes)[sizes >= 3L])
  if (length(big) < 2L) return(NULL)
  list(k1 = big[1], k2 = big[2])
}

test_that("shape/validity: hclust cl_fun gives an exact p-value and S3 object", {
  set.seed(101)
  X <- matrix(stats::rnorm(30 * 3), 30, 3)   # p = 3 so r-1 > p for typical pairs
  clf <- hclust_cl_fun(K = 3L)
  cl <- clf(X)
  pk <- pick_big_pair(cl)
  skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))

  res <- selective_campbellF(X, clf, pk$k1, pk$k2, cl_obs = cl)

  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval >= 0 && res$pval <= 1)
  expect_identical(res$condition_on, "pair-preserved")
  expect_identical(res$sampler, "campbell-truncated-F")
  expect_true(is.finite(res$T_obs))
  expect_true(is.na(res$se))                 # exact: no Monte Carlo error
  expect_identical(res$df1, ncol(X))
  expect_true(res$n_intervals >= 1L)
  expect_identical(res$K, length(unique(cl)))

  # print()/summary() run without error
  expect_output(print(res), "Selective Hotelling")
  expect_output(print(res), "campbell-truncated-F")
  expect_output(summary(res), "Valid draws")
})

test_that("exact angle/statistic identity: T_obs = (r-2) cot^2(theta_obs)", {
  set.seed(7)
  X <- matrix(stats::rnorm(30 * 3), 30, 3)
  clf <- hclust_cl_fun(K = 3L)
  cl <- clf(X)
  pk <- pick_big_pair(cl)
  skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))

  res <- selective_campbellF(X, clf, pk$k1, pk$k2, cl_obs = cl, ridge = 0)
  r  <- res$df2 + res$df1 + 1L                # r - p - 1 + p + 1 = r
  T_from_angle <- (r - 2) * (cos(res$theta_obs) / sin(res$theta_obs))^2
  expect_equal(res$T_obs, T_from_angle, tolerance = 1e-6)
})

test_that("trivial (always-feasible) selection reduces to the classical F p-value", {
  set.seed(202)
  X <- matrix(stats::rnorm(30 * 3), 30, 3)
  clf <- hclust_cl_fun(K = 3L)
  cl <- clf(X)
  pk <- pick_big_pair(cl)
  skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))

  # A cl_fun that IGNORES its input and always returns the observed labels: the
  # pair is preserved for EVERY theta, so I_E = [0, pi/2] and the truncated-F
  # p-value must equal the unconditional (classical) F p-value.
  const_cl <- function(Xmat) cl
  res <- selective_campbellF(X, const_cl, pk$k1, pk$k2, cl_obs = cl)
  expect_equal(res$n_intervals, 1L)
  expect_equal(res$pval, res$p_uncond, tolerance = 1e-6)
})

test_that("algorithm-agnostic: SAME function works with fixed-init k-means cl_fun", {
  set.seed(303)
  X <- matrix(stats::rnorm(45 * 3), 45, 3)
  init_idx <- c(1L, 16L, 31L)
  km_cl_fun <- function(Xmat) kmeans_with_fixed_init(Xmat, init_centers_idx = init_idx)$cluster
  cl <- km_cl_fun(X)
  pk <- pick_big_pair(cl)
  skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))

  res <- selective_campbellF(X, km_cl_fun, pk$k1, pk$k2, cl_obs = cl)
  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval >= 0 && res$pval <= 1)
})

test_that("errors when the pair is too small relative to dimension (r-1 <= p)", {
  set.seed(11)
  X <- matrix(stats::rnorm(12 * 6), 12, 6)   # p = 6, tiny clusters
  # Force a labelling with two 3-member clusters => r = 6, r-1 = 5 <= p = 6.
  cl <- c(rep(1L, 3), rep(2L, 3), rep(3L, 6))
  const_cl <- function(Xmat) cl
  expect_error(selective_campbellF(X, const_cl, 1L, 2L, cl_obs = cl),
               "r-1 > p")
})

test_that("small exactness sanity: null p-values ~ uniform, rejection not absurd", {
  set.seed(7)
  reps <- 120L
  clf <- hclust_cl_fun(K = 3L)
  pvals <- rep(NA_real_, reps)
  for (i in seq_len(reps)) {
    X <- matrix(stats::rnorm(30 * 3), 30, 3)
    cl <- clf(X)
    pk <- pick_big_pair(cl)
    if (is.null(pk)) next
    if (!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2))) next
    res <- tryCatch(selective_campbellF(X, clf, pk$k1, pk$k2, cl_obs = cl, n_grid = 200L),
                    error = function(e) NULL)
    if (!is.null(res)) pvals[i] <- res$pval
  }
  ok <- pvals[is.finite(pvals)]
  expect_gt(length(ok), 30L)
  expect_gt(mean(ok), 0.34)          # exact => mean near 0.5
  expect_lt(mean(ok), 0.66)
  expect_lt(mean(ok <= 0.05), 0.16)  # selective Type-I near 0.05, not inflated
})
