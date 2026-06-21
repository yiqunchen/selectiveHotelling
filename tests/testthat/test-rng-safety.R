# =============================================================================
# RNG SAFETY -- the headline regression suite.
#
# kmeans_with_fixed_init() is the loop-safe re-clustering routine used inside the
# MCMC samplers for every per-step feasibility check. It must NOT disturb the
# global RNG stream that drives the Cayley/HAR proposals. The old footgun (a
# per-step set.seed() / per-step random init) pinned every proposal to one fixed
# rotation, freezing the chain. These tests guard against any regression of that
# kind.
# =============================================================================

test_that("kmeans_with_fixed_init does not consume or reset the global RNG", {
  set.seed(123)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat
  idx <- c(1L, 11L, 21L)

  # rnorm sequence WITH an invisible kmeans_with_fixed_init() in the middle.
  set.seed(7)
  before_a <- stats::rnorm(5)
  invisible(kmeans_with_fixed_init(X, idx))
  after_a <- stats::rnorm(5)
  seq_with <- c(before_a, after_a)

  # The identical rnorm sequence WITHOUT the call.
  set.seed(7)
  before_b <- stats::rnorm(5)
  after_b <- stats::rnorm(5)
  seq_without <- c(before_b, after_b)

  # The function must neither consume draws nor reset .Random.seed.
  expect_equal(seq_with, seq_without)
  expect_lt(max(abs(seq_with - seq_without)), 1e-12)
})

test_that("kmeans_with_fixed_init leaves .Random.seed byte-identical", {
  set.seed(321)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat
  idx <- c(2L, 12L, 22L)

  # Ensure .Random.seed exists in the global env.
  set.seed(11)
  invisible(stats::rnorm(1))
  seed_before <- get(".Random.seed", envir = .GlobalEnv)

  invisible(kmeans_with_fixed_init(X, idx))

  seed_after <- get(".Random.seed", envir = .GlobalEnv)
  expect_identical(seed_before, seed_after)
})

test_that("two cayley_left_step() calls separated by kmeans_with_fixed_init() DIFFER", {
  # This is the direct regression test for the old per-step set.seed() footgun:
  # if the intervening re-clustering reset the RNG, the two proposals from the
  # same Q would be identical and the chain would never move.
  cayley_left_step <- sh_internal("cayley_left_step")

  set.seed(99)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat
  idx <- c(1L, 11L, 21L)

  # A fixed Stiefel point to step away from.
  set.seed(2024)
  Q <- qr.Q(qr(matrix(stats::rnorm(6 * 2), 6, 2)))

  set.seed(555)
  step1 <- cayley_left_step(Q, 0.1)
  invisible(kmeans_with_fixed_init(X, idx))
  step2 <- cayley_left_step(Q, 0.1)

  # The two proposals must be genuinely different draws.
  expect_gt(max(abs(step1 - step2)), 1e-8)
})

test_that("kmeans_with_trajectory(seed=) restores the global RNG on exit", {
  # The documented RNG policy: when `seed` is supplied it is applied LOCALLY and
  # the global RNG stream is restored on exit, so calling this inside an MCMC
  # loop can never freeze the chain.
  kmeans_with_trajectory <- sh_internal("kmeans_with_trajectory")

  set.seed(202)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat

  set.seed(13)
  invisible(stats::rnorm(1))
  seed_before <- get(".Random.seed", envir = .GlobalEnv)

  invisible(kmeans_with_trajectory(X, k = 3, seed = 42))

  seed_after <- get(".Random.seed", envir = .GlobalEnv)
  expect_identical(seed_before, seed_after)
})

test_that("kmeans_with_trajectory(init_centers_idx=) does not touch the global RNG", {
  # The loop-safe path actually used by the samplers: an explicit fixed
  # initialization means no RNG is consumed at all.
  kmeans_with_trajectory <- sh_internal("kmeans_with_trajectory")

  set.seed(404)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat

  set.seed(17)
  invisible(stats::rnorm(1))
  seed_before <- get(".Random.seed", envir = .GlobalEnv)

  invisible(kmeans_with_trajectory(X, k = 3, init_centers_idx = c(1L, 11L, 21L)))

  seed_after <- get(".Random.seed", envir = .GlobalEnv)
  expect_identical(seed_before, seed_after)
})

test_that("kmeans_with_trajectory() with default (random) init advances the RNG", {
  # Documents the real contract: the DEFAULT initialization is data-independent
  # but RANDOM (sample.int), so it deliberately consumes the caller's RNG. This
  # is by design -- reproducibility is the caller's responsibility (seed once at
  # the call site). Crucially the loop-safe routine used per MCMC step is
  # kmeans_with_fixed_init(), tested above, NOT this default path.
  kmeans_with_trajectory <- sh_internal("kmeans_with_trajectory")

  set.seed(606)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat

  set.seed(19)
  invisible(stats::rnorm(1))
  seed_before <- get(".Random.seed", envir = .GlobalEnv)

  invisible(kmeans_with_trajectory(X, k = 3))

  seed_after <- get(".Random.seed", envir = .GlobalEnv)
  expect_false(identical(seed_before, seed_after))
})
