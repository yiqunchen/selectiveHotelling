# =============================================================================
# kmeans_with_fixed_init determinism.
#
# With a fixed initialization and no RNG use, repeated calls on the same data
# must produce byte-identical results.
# =============================================================================

test_that("kmeans_with_fixed_init is deterministic across repeated calls", {
  set.seed(101)
  X <- generate_exchangeable_null_data(40, 3, rho = 0.15)$dat
  idx <- c(1L, 15L, 30L)

  r1 <- kmeans_with_fixed_init(X, idx)
  r2 <- kmeans_with_fixed_init(X, idx)

  expect_identical(r1$cluster, r2$cluster)
  expect_identical(r1$trajectory, r2$trajectory)
  expect_identical(r1$n_iter, r2$n_iter)
  expect_identical(r1$init_centers_idx, r2$init_centers_idx)
})

test_that("kmeans_with_fixed_init determinism is unaffected by intervening RNG use", {
  # Even if the global RNG advances between calls, the fixed-init clustering must
  # be identical (it consumes no randomness).
  set.seed(202)
  X <- generate_exchangeable_null_data(40, 3, rho = 0.15)$dat
  idx <- c(2L, 16L, 31L)

  c1 <- kmeans_with_fixed_init(X, idx)$cluster
  invisible(stats::rnorm(100))
  c2 <- kmeans_with_fixed_init(X, idx)$cluster

  expect_identical(c1, c2)
})

test_that("kmeans_with_fixed_init infers k from the length of init_centers_idx", {
  set.seed(303)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat

  res2 <- kmeans_with_fixed_init(X, c(1L, 16L))
  res4 <- kmeans_with_fixed_init(X, c(1L, 8L, 16L, 24L))

  expect_setequal(unique(res2$cluster), 1:2)
  expect_true(all(res4$cluster %in% 1:4))
})
