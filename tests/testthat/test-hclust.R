# =============================================================================
# Hierarchical-clustering support: hclust_fun() helper and the hclust route of
# selective_hotelling_test() (which dispatches to the exact Besag-Clifford
# engine). Kept FAST (small M / few reps) so the suite stays quick.
# =============================================================================

# Pick a pair of hclust clusters, each with >= 3 members.
pick_big_pair_h <- function(cl) {
  sizes <- table(cl)
  big <- as.integer(names(sizes)[sizes >= 3L])
  if (length(big) < 2L) return(NULL)
  list(k1 = big[1], k2 = big[2])
}

test_that("hclust_fun returns the documented structure (no init_centers_idx)", {
  set.seed(11)
  X <- matrix(stats::rnorm(30 * 2), 30, 2)

  out <- hclust_fun(X, k = 3)

  # $cluster: integer vector of the right length, labels in 1:3
  expect_true(is.integer(out$cluster))
  expect_length(out$cluster, nrow(X))
  expect_true(all(out$cluster %in% 1:3))

  # $trajectory: a length-1 list whose single element is the final assignment
  expect_true(is.list(out$trajectory))
  expect_length(out$trajectory, 1L)
  expect_identical(out$trajectory[[1]], out$cluster)

  # CRUCIALLY: no init_centers_idx field (pure black box, never re-seeded)
  expect_null(out$init_centers_idx)
  expect_false("init_centers_idx" %in% names(out))
})

test_that("hclust_fun square=TRUE (default) can differ from square=FALSE", {
  # Squared vs raw Euclidean distances can give different average-linkage merges
  # (averaging is not preserved under squaring). Use a dataset/k where they are
  # known to differ on this platform.
  set.seed(2)
  X <- matrix(stats::rnorm(40 * 3), 40, 3)

  cl_sq  <- hclust_fun(X, k = 3, square = TRUE)$cluster
  cl_raw <- hclust_fun(X, k = 3, square = FALSE)$cluster

  # Both are valid labelings of the right length.
  expect_length(cl_sq, nrow(X))
  expect_length(cl_raw, nrow(X))

  # A permutation-equivalent partition would have exactly one nonzero per row
  # and column in the cross-tabulation; here the partitions genuinely differ.
  tab <- table(cl_sq, cl_raw)
  is_perm <- (nrow(tab) == ncol(tab)) &&
    all(rowSums(tab > 0) == 1L) && all(colSums(tab > 0) == 1L)
  expect_false(is_perm)
})

test_that("selective_hotelling_test(clustering='hclust') works end-to-end", {
  set.seed(101)
  X <- matrix(stats::rnorm(30 * 2), 30, 2)

  # Pick a pair from the ACTUAL hclust clustering (squared-distance default),
  # each with >= 3 members.
  cl <- hclust_fun(X, k = 3)$cluster
  pk <- pick_big_pair_h(cl)
  expect_false(is.null(pk))

  res <- selective_hotelling_test(
    X, pk$k1, pk$k2, K = 3, clustering = "hclust", M = 80L
  )

  expect_s3_class(res, "selectiveHotelling")
  expect_true(is.finite(res$pval))
  expect_gt(res$pval, 0)
  expect_lte(res$pval, 1)
  expect_identical(res$sampler, "besag-clifford")
  expect_identical(res$condition_on, "pair-preserved")
  expect_true(is.finite(res$T_obs))

  # print()/summary() run without error.
  expect_output(print(res), "Selective Hotelling")
  expect_output(summary(res), "Valid draws")
})

test_that("small exactness sanity: hclust route mean p-value ~ uniform", {
  skip_on_cran()
  set.seed(7)
  reps <- 36L
  M <- 60L
  pvals <- rep(NA_real_, reps)

  for (i in seq_len(reps)) {
    X <- matrix(stats::rnorm(30 * 2), 30, 2)
    cl <- hclust_fun(X, k = 3)$cluster
    pk <- pick_big_pair_h(cl)
    if (is.null(pk)) next
    if (!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2))) next
    res <- tryCatch(
      selective_hotelling_test(X, pk$k1, pk$k2, K = 3,
                               clustering = "hclust", M = M),
      error = function(e) NULL)
    if (!is.null(res)) pvals[i] <- res$pval
  }

  ok <- pvals[is.finite(pvals)]
  expect_gt(length(ok), 10L)
  # Loose bounds: BC is exact, so the mean is near 0.5 and rejection near 0.05.
  expect_gt(mean(ok), 0.30)
  expect_lt(mean(ok), 0.70)
  expect_lt(mean(ok <= 0.05), 0.25)
})

test_that("backward-compat: clustering='kmeans' defaults still use the MCMC sampler", {
  set.seed(2024)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat

  res <- selective_hotelling_test(
    X, 1, 2, K = 3, clustering = "kmeans", sampler = "fast",
    burn = 100, B = 200, thin = 2,
    verbose = FALSE, show_progress = FALSE
  )

  expect_s3_class(res, "selectiveHotelling")
  expect_identical(res$sampler, "fast")
  expect_identical(res$condition_on, "final")
  expect_identical(res$K, 3)
  expect_true(is.finite(res$pval))
  expect_gt(res$pval, 0)
  expect_lte(res$pval, 1)
})
