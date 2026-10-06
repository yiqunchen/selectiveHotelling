# =============================================================================
# selective_hotelling_test: end-to-end pipeline (cluster -> polar reduction ->
# Besag-Clifford path with the random Givens kernel -> p-value), plus the S3
# print/summary methods. Kept small (short paths) so the suite runs quickly.
# =============================================================================

three_blobs <- function(seed) {
  set.seed(seed)
  rbind(matrix(stats::rnorm(20), 10, 2),
        matrix(stats::rnorm(20) + 6, 10, 2),
        matrix(stats::rnorm(20) - 6, 10, 2))
}

test_that("selective_hotelling_test (k-means default) returns a valid object", {
  X <- three_blobs(2)
  set.seed(11)
  res <- selective_hotelling_test(X, 1, 2, K = 3, M = 100L)
  expect_s3_class(res, "selectiveHotelling")
  expect_true(is.finite(res$pval))
  expect_gt(res$pval, 0)
  expect_lte(res$pval, 1)
  expect_true(is.finite(res$T_obs))
  expect_identical(res$sampler, "besag-clifford")
  expect_identical(res$condition_on, "pair-preserved")
  expect_identical(as.integer(res$K), 3L)
})

test_that("print and summary methods for selectiveHotelling run without error", {
  X <- three_blobs(2)
  set.seed(11)
  res <- selective_hotelling_test(X, 1, 2, K = 3, M = 50L)
  expect_invisible(print(res))
  out_print <- capture.output(print(res))
  expect_true(any(grepl("Selective Hotelling", out_print)))
  expect_true(any(grepl("Selective p-val", out_print)))
  expect_invisible(summary(res))
  out_summary <- capture.output(summary(res))
  expect_true(any(grepl("Valid draws", out_summary)))
})
