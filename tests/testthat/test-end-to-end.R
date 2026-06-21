# =============================================================================
# selective_hotelling_test: end-to-end pipeline on small null data.
#
# Kept intentionally small (modest burn/B) so the suite runs quickly while still
# exercising the full path: cluster -> polar reduction -> MCMC -> p-value, plus
# the S3 print/summary methods.
# =============================================================================

test_that("selective_hotelling_test returns a valid selectiveHotelling object on null data", {
  set.seed(2024)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat

  # Sanity: clusters 1 and 2 must each have >= 2 members for a finite T^2.
  cl <- kmeans_fun(X, k = 3)$cluster
  expect_gte(sum(cl == 1), 2L)
  expect_gte(sum(cl == 2), 2L)

  res <- selective_hotelling_test(
    X, 1, 2, K = 3, sampler = "fast",
    burn = 150, B = 300, thin = 2,
    verbose = FALSE, show_progress = FALSE
  )

  expect_s3_class(res, "selectiveHotelling")
  expect_true(is.list(res))

  # p-value is finite and a valid probability in (0, 1].
  expect_true(is.finite(res$pval))
  expect_gt(res$pval, 0)
  expect_lte(res$pval, 1)

  # Bookkeeping fields are populated as documented.
  expect_true(is.finite(res$T_obs))
  expect_gte(res$T_obs, 0)
  expect_identical(res$k1, 1)
  expect_identical(res$k2, 2)
  expect_identical(res$K, 3)
  expect_identical(res$condition_on, "final")
  expect_identical(res$sampler, "fast")
})

test_that("print and summary methods for selectiveHotelling run without error", {
  set.seed(2024)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat

  res <- selective_hotelling_test(
    X, 1, 2, K = 3, sampler = "fast",
    burn = 100, B = 200, thin = 2,
    verbose = FALSE, show_progress = FALSE
  )

  # Methods print to the console; capture output and confirm they run, return
  # the object invisibly, and emit the expected headers.
  expect_invisible(print(res))
  out_print <- capture.output(print(res))
  expect_true(any(grepl("Selective Hotelling", out_print)))
  expect_true(any(grepl("Selective p-val", out_print)))

  expect_invisible(summary(res))
  out_summary <- capture.output(summary(res))
  expect_true(any(grepl("Valid draws", out_summary)))
  expect_true(any(grepl("Acceptance rate", out_summary)))
})

test_that("selective_hotelling_mcmc_fast (exported low-level sampler) returns a finite p-value", {
  set.seed(2024)
  X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat
  cl_fun <- function(x) kmeans_fun(x, k = 3)

  res <- selective_hotelling_mcmc_fast(
    X, cl_fun, 1, 2,
    burn = 100, B = 200, thin = 2,
    condition_on = "final",
    verbose = FALSE, show_progress = FALSE
  )

  expect_true(is.finite(res$pval))
  expect_gt(res$pval, 0)
  expect_lte(res$pval, 1)
  expect_true(is.finite(res$T_obs))
  expect_length(res$draws, 200L)
})
