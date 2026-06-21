# =============================================================================
# mc_tail_pvalue: Monte Carlo upper-tail p-value with add-one correction.
#   pval = (1 + #{draws >= T_obs}) / (n_ok + 1)
# =============================================================================

mc_tail_pvalue <- sh_internal("mc_tail_pvalue")

test_that("mc_tail_pvalue returns a probability in (0, 1]", {
  set.seed(1)
  draws <- stats::rnorm(200)
  res <- mc_tail_pvalue(draws, T_obs = 0)
  expect_true(is.finite(res$pval))
  expect_gt(res$pval, 0)
  expect_lte(res$pval, 1)
  expect_identical(res$n, 200L)
})

test_that("mc_tail_pvalue is monotone: larger T_obs gives a smaller-or-equal p", {
  set.seed(2)
  draws <- stats::rnorm(300)
  grid <- c(-2, -1, 0, 1, 2, 3)
  pvals <- vapply(grid, function(t) mc_tail_pvalue(draws, t)$pval, numeric(1))
  # Non-increasing in T_obs.
  expect_true(all(diff(pvals) <= 1e-12))
  # Strictly smaller for a clearly larger T_obs.
  expect_lt(mc_tail_pvalue(draws, 2)$pval, mc_tail_pvalue(draws, 0)$pval)
})

test_that("mc_tail_pvalue uses the add-one form (minimum p = 1/(n+1))", {
  set.seed(3)
  draws <- stats::rnorm(200)
  # T_obs above every draw: no exceedances => p = 1/(n+1).
  res <- mc_tail_pvalue(draws, T_obs = max(draws) + 100)
  expect_equal(res$pval, 1 / (200 + 1), tolerance = 1e-12)
})

test_that("mc_tail_pvalue equals 1 when T_obs is below every draw", {
  set.seed(4)
  draws <- stats::rnorm(150)
  # Every draw exceeds T_obs => (1 + n)/(n + 1) = 1.
  res <- mc_tail_pvalue(draws, T_obs = min(draws) - 100)
  expect_equal(res$pval, 1, tolerance = 1e-12)
})

test_that("mc_tail_pvalue ignores non-finite draws and uses only the finite count", {
  set.seed(5)
  draws <- c(stats::rnorm(100), NA_real_, NA_real_, Inf, -Inf, NaN)
  res <- mc_tail_pvalue(draws, T_obs = 0)
  expect_identical(res$n, 100L)
  expect_gt(res$pval, 0)
  expect_lte(res$pval, 1)
})

test_that("mc_tail_pvalue returns NA when there are no finite draws or T_obs is non-finite", {
  expect_true(is.na(mc_tail_pvalue(c(NA_real_, Inf, NaN), 1)$pval))
  expect_identical(mc_tail_pvalue(c(NA_real_, Inf), 1)$n, 0L)
  set.seed(6)
  expect_true(is.na(mc_tail_pvalue(stats::rnorm(50), NA_real_)$pval))
})
