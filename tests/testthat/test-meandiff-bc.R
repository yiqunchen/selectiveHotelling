# =============================================================================
# selective_meandiff_bc: Besag-Clifford test on the O(r-1) orbit, any p.
# =============================================================================

km_cl_fun <- function(init) { force(init); function(x) kmeans_with_fixed_init(x, init)$cluster }

test_that("the Hotelling statistic reproduces selective_hotelling_bc exactly (p < r)", {
  set.seed(3); X <- matrix(stats::rnorm(45 * 4), 45, 4)
  clf <- km_cl_fun(c(1, 2, 3))
  set.seed(10); a <- selective_hotelling_bc(X, clf, 1, 2, M = 20L)
  set.seed(10); b <- selective_meandiff_bc(X, clf, 1, 2, M = 20L, stats = "hotelling", adapt = NULL)
  expect_equal(unname(b$pvals["hotelling"]), a$pval, tolerance = 1e-12)
  n_eff <- b$n1 * b$n2 / (b$n1 + b$n2)
  expect_equal(unname(b$path[, "hotelling"]) / n_eff, a$path_u, tolerance = 1e-8)
})

test_that("runs for p > n1 + n2: guards pass, observed statistics match the data", {
  set.seed(4); X <- matrix(stats::rnorm(45 * 200), 45, 200)
  clf <- km_cl_fun(c(1, 2, 3)); cl <- clf(X)
  expect_lt(sum(cl == 1) + sum(cl == 2), 200)
  res <- selective_meandiff_bc(X, clf, 1, 2, M = 60L, stats = c("l2", "diag", "max", "ridge", "hotelling"))
  expect_true(res$feas_ok); expect_lt(res$scatter_max, 1e-8)
  expect_true(all(res$pvals[c("l2", "diag", "max", "ridge", "adapt")] > 0 &
                  res$pvals[c("l2", "diag", "max", "ridge", "adapt")] <= 1))
  expect_true(is.na(res$pvals["hotelling"]))           # constant when p >= r - 1
  d <- colMeans(X[cl == 1, ]) - colMeans(X[cl == 2, ])
  n_eff <- res$n1 * res$n2 / (res$n1 + res$n2)
  expect_equal(unname(res$obs["l2"]), n_eff * sum(d^2), tolerance = 1e-8)
  s0 <- apply(X[cl %in% 1:2, ], 2, var)
  expect_equal(unname(res$obs["max"]), n_eff * max(d^2 / s0), tolerance = 1e-8)
})

test_that("the reconstruction map returns the data at Z_obs", {
  set.seed(5); X <- matrix(stats::rnorm(30 * 50), 30, 50)
  clf <- km_cl_fun(c(1, 2, 3)); cl <- clf(X)
  st <- sh_internal(".hd_setup")(X, cl, 1, 2, clf)
  expect_equal(st$X_of_Z(st$Z_obs), X, tolerance = 1e-10)
  expect_true(st$feasible(st$Z_obs))
})

test_that("null p-values for p > r are roughly uniform", {
  pv <- vapply(1:120, function(s) {
    set.seed(1000 + s); X <- matrix(stats::rnorm(24 * 60), 24, 60)
    clf <- km_cl_fun(c(1, 9, 17)); cl <- clf(X)
    if (sum(cl == 1) < 2 || sum(cl == 2) < 2) return(NA_real_)
    selective_meandiff_bc(X, clf, 1, 2, M = 30L, stats = c("l2", "max"),
                          adapt = c("l2", "max"), cl_obs = cl)$pvals["adapt"]
  }, numeric(1))
  pv <- pv[!is.na(pv)]
  expect_gt(stats::ks.test(pv, "punif")$p.value, 0.001)
})

test_that("cq, diag_w and proj match their formulas and leave the path unchanged", {
  set.seed(6); X <- matrix(stats::rnorm(45 * 120), 45, 120)
  clf <- km_cl_fun(c(1, 2, 3)); cl <- clf(X)
  set.seed(1); a <- selective_meandiff_bc(X, clf, 1, 2, M = 40L, stats = c("l2", "max"), adapt = NULL)
  set.seed(1); b <- selective_meandiff_bc(X, clf, 1, 2, M = 40L, adapt = NULL,
                                          stats = c("l2", "max", "cq", "diag_w", "proj"))
  expect_equal(a$path[, "l2"], b$path[, "l2"])
  X1 <- X[cl == 1, ]; X2 <- X[cl == 2, ]; n1 <- nrow(X1); n2 <- nrow(X2)
  d <- colMeans(X1) - colMeans(X2); n_eff <- n1 * n2 / (n1 + n2)
  S1 <- stats::cov(X1); S2 <- stats::cov(X2)
  expect_equal(unname(b$obs["cq"]),
               n_eff * (sum(d^2) - sum(diag(S1)) / n1 - sum(diag(S2)) / n2), tolerance = 1e-8)
  Sw <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  expect_equal(unname(b$obs["diag_w"]), n_eff * sum(d^2 / diag(Sw)), tolerance = 1e-8)
  expect_true(is.finite(b$obs["proj"]) && b$pvals["proj"] > 0)
})

test_that("t2 is the regular two-sample T^2 and ranks like hotelling when p < r", {
  set.seed(8); X <- matrix(stats::rnorm(60 * 8), 60, 8)
  clf <- km_cl_fun(c(1, 2, 3)); cl <- clf(X)
  b <- selective_meandiff_bc(X, clf, 1, 2, M = 50L, stats = c("hotelling", "t2"), adapt = NULL)
  expect_equal(unname(b$obs["t2"]), hotelling_T2(X, cl, 1, 2), tolerance = 1e-8)
  expect_identical(rank(b$path[, "t2"]), rank(b$path[, "hotelling"]))
  set.seed(9); X <- matrix(stats::rnorm(40 * 100), 40, 100)   # p > r: not constant
  cl <- clf(X)
  b <- selective_meandiff_bc(X, clf, 1, 2, M = 50L, stats = "t2", adapt = NULL)
  expect_equal(unname(b$obs["t2"]), hotelling_T2(X, cl, 1, 2), tolerance = 1e-6)
  expect_gt(length(unique(signif(b$path[, "t2"], 8))), 1L)
})

test_that("several min-p combinations are reported by name", {
  set.seed(10); X <- matrix(stats::rnorm(45 * 80), 45, 80)
  b <- selective_meandiff_bc(X, make_fast_kmeans_cl_fun(c(1, 2, 3)), 1, 2, M = 30L,
                             stats = c("l2", "diag_w", "max"),
                             adapt = list(adapt_l2_max = c("l2", "max"), adapt_dw_max = c("diag_w", "max")))
  expect_true(all(c("adapt_l2_max", "adapt_dw_max") %in% names(b$pvals)))
})
