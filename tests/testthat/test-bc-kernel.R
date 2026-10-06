# =============================================================================
# selective_hotelling_bc(kernel = ...): the optional Givens proposal.
# =============================================================================

.km_fixed_cl <- function(init_idx) function(Xmat)
  kmeans_with_fixed_init(Xmat, init_centers_idx = init_idx)$cluster

test_that("both kernels give a p-value in [0, 1] and a feasible Stiefel path", {
  set.seed(303)
  X <- matrix(stats::rnorm(36 * 2), 36, 2)
  clf <- .km_fixed_cl(c(1L, 13L, 25L))
  cl <- clf(X)
  sz <- table(cl); big <- as.integer(names(sz)[sz >= 4L])
  skip_if(length(big) < 2L, "need two clusters with >= 4 members")
  for (kern in c("cayley", "givens")) {
    set.seed(7)
    res <- selective_hotelling_bc(X, clf, big[1], big[2], M = 40L, L = 2L,
                                  cl_obs = cl, kernel = kern,
                                  guard_feas_checks = 41L)
    expect_identical(res$kernel, kern)
    expect_true(res$pval >= 0 && res$pval <= 1)
    expect_true(isTRUE(res$feas_ok))           # every path point re-verified
    expect_lt(res$stiefel_max, 1e-8)
    expect_equal(length(res$path_u), 41L)
    expect_true(res$acc >= 0 && res$acc <= 1)
  }
})

test_that("the default kernel is givens and matches an explicit kernel = \"givens\" call", {
  set.seed(304)
  X <- matrix(stats::rnorm(30 * 2), 30, 2)
  clf <- .km_fixed_cl(c(1L, 11L, 21L))
  cl <- clf(X)
  sz <- table(cl); big <- as.integer(names(sz)[sz >= 4L])
  skip_if(length(big) < 2L, "need two clusters with >= 4 members")
  set.seed(9); a <- selective_hotelling_bc(X, clf, big[1], big[2], M = 30L, cl_obs = cl)
  set.seed(9); b <- selective_hotelling_bc(X, clf, big[1], big[2], M = 30L, cl_obs = cl,
                                           kernel = "givens")
  expect_identical(a$kernel, "givens")
  expect_identical(a$path_u, b$path_u)
  expect_identical(a$pval, b$pval)
})

test_that("the Givens proposal stays on the Stiefel manifold and moves two rows", {
  set.seed(305)
  Q <- qr.Q(qr(matrix(stats::rnorm(6 * 2), 6, 2)))
  Qp <- .bc_givens_step(Q, eps = 0.3)
  expect_lt(max(abs(crossprod(Qp) - diag(2))), 1e-12)
  expect_lte(sum(rowSums(abs(Qp - Q)) > 1e-14), 2L)
})

test_that("the shrinking-bracket Givens kernel gives a valid path and counts cl_fun calls", {
  set.seed(306)
  X <- matrix(stats::rnorm(36 * 2), 36, 2)
  clf <- .km_fixed_cl(c(1L, 13L, 25L))
  cl <- clf(X)
  sz <- table(cl); big <- as.integer(names(sz)[sz >= 4L])
  skip_if(length(big) < 2L, "need two clusters with >= 4 members")
  set.seed(8)
  res <- selective_hotelling_bc(X, clf, big[1], big[2], M = 40L, L = 2L,
                                cl_obs = cl, kernel = "givens_slice",
                                guard_feas_checks = 41L)
  expect_identical(res$kernel, "givens_slice")
  expect_true(is.na(res$eps))                  # no step size
  expect_true(res$pval >= 0 && res$pval <= 1)
  expect_true(isTRUE(res$feas_ok))
  expect_lt(res$stiefel_max, 1e-8)
  expect_gt(res$move_frac, 0.9)                # the bracket move nearly always moves
  expect_gte(res$n_path_calls, 40L * 2L)       # >= one call per step
  # all calls = path + guard (c) + guard (b) (cl_obs supplied, eps not tuned)
  expect_identical(res$n_cluster_calls, res$n_path_calls + 1L + 41L)
  # eps is ignored: a different eps gives the identical path
  set.seed(8)
  res2 <- selective_hotelling_bc(X, clf, big[1], big[2], M = 40L, L = 2L,
                                 cl_obs = cl, kernel = "givens_slice", eps = 0.01,
                                 guard_feas_checks = 41L)
  expect_identical(res$path_u, res2$path_u)
})

test_that("Metropolis kernels use exactly one cl_fun call per step", {
  set.seed(307)
  X <- matrix(stats::rnorm(30 * 2), 30, 2)
  clf <- .km_fixed_cl(c(1L, 11L, 21L))
  cl <- clf(X)
  sz <- table(cl); big <- as.integer(names(sz)[sz >= 4L])
  skip_if(length(big) < 2L, "need two clusters with >= 4 members")
  for (kern in c("cayley", "givens")) {
    set.seed(10)
    res <- selective_hotelling_bc(X, clf, big[1], big[2], M = 25L, L = 3L,
                                  kernel = kern, guard_feas_checks = 5L)
    expect_identical(res$n_path_calls, 25L * 3L)
    # cl_obs computed internally (+1), guard (c) (+1), guard (b) (+5)
    expect_identical(res$n_cluster_calls, 25L * 3L + 7L)
  }
})

test_that("the shrinking-bracket step stays in a known feasible set on the Stiefel", {
  set.seed(308)
  feas <- function(Q) sum(Q[1:3, ]^2) <= 0.25
  repeat { Q <- qr.Q(qr(matrix(stats::rnorm(20 * 2), 20, 2))); if (feas(Q)) break }
  moved <- 0L
  for (s in 1:200) {
    st <- .bc_givens_slice_step(Q, feas)
    expect_true(feas(st$Q))
    expect_gte(st$n_calls, 1L)
    moved <- moved + st$moved
    Q <- st$Q
  }
  expect_lt(max(abs(crossprod(Q) - diag(2))), 1e-10)
  expect_gt(moved, 190L)
})
