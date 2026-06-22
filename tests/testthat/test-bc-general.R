# =============================================================================
# selective_hotelling_bc: algorithm-agnostic Besag-Clifford selective test.
# =============================================================================

# A deterministic, data-pure hierarchical-clustering routine (average linkage).
hclust_cl_fun <- function(K = 3L, method = "average") {
  force(K); force(method)
  function(Xmat) {
    hc <- stats::hclust(stats::dist(Xmat), method = method)
    as.integer(stats::cutree(hc, k = K))
  }
}

# Pick a pair of clusters, each with >= 3 members.
pick_big_pair <- function(cl) {
  sizes <- table(cl)
  big <- as.integer(names(sizes)[sizes >= 3L])
  if (length(big) < 2L) return(NULL)
  list(k1 = big[1], k2 = big[2])
}

test_that("validity/shape smoke: hclust cl_fun gives a valid p-value, guards pass", {
  set.seed(101)
  X <- matrix(stats::rnorm(30 * 2), 30, 2)
  clf <- hclust_cl_fun(K = 3L)
  cl <- clf(X)
  pk <- pick_big_pair(cl)
  expect_false(is.null(pk))

  res <- selective_hotelling_bc(X, clf, pk$k1, pk$k2, M = 80L, L = 3L,
                                cl_obs = cl)

  # p-value in (0, 1]
  expect_true(res$pval > 0 && res$pval <= 1)
  # GUARDS
  expect_true(isTRUE(res$feas_ok))
  expect_lt(res$stiefel_max, 1e-8)
  # S3 object + primary fields
  expect_s3_class(res, "selectiveHotelling")
  expect_identical(res$condition_on, "pair-preserved")
  expect_identical(res$sampler, "besag-clifford")
  expect_true(is.finite(res$T_obs))
  expect_equal(length(res$draws), 81L)
  expect_true(is.finite(res$se))
  expect_identical(res$K, length(unique(cl)))

  # print()/summary() run without error
  expect_output(print(res), "Selective Hotelling")
  expect_output(summary(res), "Valid draws")
})

test_that("algorithm-agnostic: SAME function works with fixed-init k-means cl_fun", {
  set.seed(202)
  X <- matrix(stats::rnorm(36 * 2), 36, 2)

  # Pure (data-only-dependent) k-means via a FIXED initialization: no RNG.
  init_idx <- c(1L, 13L, 25L)
  km_cl_fun <- function(Xmat) {
    kmeans_with_fixed_init(Xmat, init_centers_idx = init_idx)$cluster
  }
  cl <- km_cl_fun(X)
  pk <- pick_big_pair(cl)
  skip_if(is.null(pk), "no fixed-init k-means pair with >= 3 members each")

  res <- selective_hotelling_bc(X, km_cl_fun, pk$k1, pk$k2, M = 80L, L = 3L,
                                cl_obs = cl)

  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval > 0 && res$pval <= 1)
  expect_true(isTRUE(res$feas_ok))
  expect_lt(res$stiefel_max, 1e-8)
})

test_that("oracle_Sigma path returns a valid extra p-value", {
  set.seed(7)
  X <- matrix(stats::rnorm(30 * 2), 30, 2)
  clf <- hclust_cl_fun(K = 3L)
  cl <- clf(X)
  pk <- pick_big_pair(cl)
  skip_if(is.null(pk))

  res <- selective_hotelling_bc(X, clf, pk$k1, pk$k2, M = 60L, L = 3L,
                                cl_obs = cl, oracle_Sigma = diag(2))
  expect_true(is.finite(res$p_oracle))
  expect_true(res$p_oracle > 0 && res$p_oracle <= 1)
})

test_that("small exactness sanity: null mean p-value ~ uniform, rejection not absurd", {
  set.seed(7)
  reps <- 36L
  M <- 60L
  pvals <- rep(NA_real_, reps)
  clf <- hclust_cl_fun(K = 3L)

  for (i in seq_len(reps)) {
    X <- matrix(stats::rnorm(30 * 2), 30, 2)
    cl <- clf(X)
    pk <- pick_big_pair(cl)
    if (is.null(pk)) next
    # require finite Hotelling T^2 for the chosen pair
    if (!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2))) next
    res <- tryCatch(
      selective_hotelling_bc(X, clf, pk$k1, pk$k2, M = M, L = 3L, cl_obs = cl),
      error = function(e) NULL)
    if (!is.null(res)) pvals[i] <- res$pval
  }

  ok <- pvals[is.finite(pvals)]
  expect_gt(length(ok), 10L)            # enough usable reps
  # Loose bounds: BC is exact, so the mean should be near 0.5 and rejection near 0.05.
  expect_gt(mean(ok), 0.30)
  expect_lt(mean(ok), 0.70)
  expect_lt(mean(ok <= 0.05), 0.25)
})
