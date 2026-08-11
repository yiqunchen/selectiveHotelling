# =============================================================================
# selective_truncated_F: closed-form one-dimensional truncated-F selective test.
# Exact k-means (Lloyd trajectory) + exact hierarchical (pair-preserved) engines,
# and the algorithm-agnostic grid engine.
# =============================================================================

hclust_cl_fun <- function(K = 3L, method = "average", square = TRUE) {
  force(K); force(method); force(square)
  function(Xmat) {
    d <- stats::dist(Xmat); if (isTRUE(square)) d <- d^2
    as.integer(stats::cutree(stats::hclust(d, method = method), k = K))
  }
}
pick_big_pair <- function(cl) {
  sizes <- table(cl)
  big <- as.integer(names(sizes)[sizes >= 3L])
  if (length(big) < 2L) return(NULL)
  list(k1 = big[1], k2 = big[2])
}

test_that("exact hclust route: valid p-value, S3 object, exact-vs-grid agree", {
  set.seed(101)
  X <- matrix(stats::rnorm(30 * 3), 30, 3)
  cl <- hclust_cl_fun(3L)(X)
  pk <- pick_big_pair(cl); skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))

  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "hclust", K = 3L)
  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval >= 0 && res$pval <= 1)
  expect_identical(res$condition_on, "pair-preserved")
  expect_identical(res$sampler, "truncated-F")
  expect_true(is.na(res$se))                 # exact: no Monte Carlo error
  expect_identical(res$df1, ncol(X))
  expect_true(res$n_intervals >= 1L)

  # exact closed-form arc and grid engine condition on the SAME event => agree
  res_g <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "hclust", K = 3L,
                               method = "grid", n_grid = 600L)
  expect_equal(res$pval, res_g$pval, tolerance = 5e-3)

  expect_output(print(res), "truncated-F")
  expect_output(summary(res), "Valid draws")
})

test_that("exact k-means route: valid p-value and T2 = (r-2) cot^2(theta) identity", {
  set.seed(303)
  X <- matrix(stats::rnorm(45 * 3), 45, 3)
  init_idx <- c(1L, 16L, 31L)
  ref <- kmeans_with_fixed_init(X, init_centers_idx = init_idx)
  cl <- ref$cluster
  pk <- pick_big_pair(cl); skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))

  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "kmeans", K = 3L,
                             init_centers_idx = init_idx, ridge = 0)
  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval >= 0 && res$pval <= 1)
  expect_identical(res$condition_on, "kmeans-trajectory")

  r <- res$df2 + res$df1 + 1L                # = r (pair size)
  T_from_angle <- (r - 2) * (cos(res$theta_obs) / sin(res$theta_obs))^2
  expect_equal(res$T_obs, T_from_angle, tolerance = 1e-6)
})

test_that("trivial (always-feasible) selection reduces to the classical F p-value", {
  set.seed(202)
  X <- matrix(stats::rnorm(30 * 3), 30, 3)
  cl <- hclust_cl_fun(3L)(X)
  pk <- pick_big_pair(cl); skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))

  const_cl <- function(Xmat) cl          # ignores its input => pair always preserved
  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "custom",
                             cl_fun = const_cl, cl_obs = cl)
  expect_equal(res$n_intervals, 1L)
  expect_equal(res$pval, res$p_uncond, tolerance = 1e-6)
})

test_that("custom cl_fun (grid) works and is algorithm-agnostic", {
  set.seed(404)
  X <- matrix(stats::rnorm(45 * 3), 45, 3)
  init_idx <- c(1L, 16L, 31L)
  km_cl_fun <- function(Xmat) kmeans_with_fixed_init(Xmat, init_centers_idx = init_idx)$cluster
  cl <- km_cl_fun(X)
  pk <- pick_big_pair(cl); skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))
  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "custom",
                             cl_fun = km_cl_fun, cl_obs = cl)
  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval >= 0 && res$pval <= 1)
})

test_that("errors when the pair is too small relative to dimension (r-1 <= p)", {
  set.seed(11)
  X <- matrix(stats::rnorm(12 * 6), 12, 6)
  cl <- c(rep(1L, 3), rep(2L, 3), rep(3L, 6))
  const_cl <- function(Xmat) cl
  expect_error(selective_truncated_F(X, 1L, 2L, clustering = "custom",
                                   cl_fun = const_cl, cl_obs = cl), "r-1 > p")
})

test_that("small exactness sanity: exact hclust null p-values ~ uniform", {
  # OPT-IN: this block is far too slow for a routine suite run. A single
  # selective_truncated_F() call on the "hclust" exact route takes MINUTES even at
  # n=30, p=3 -- the cost is the pair-preserved arc solver's event-driven cell
  # traversal, not the (shared, fast) truncated-F p-value: the same p-value code
  # runs 12 replicates of the k-means route at n=400, p=10 in seconds. At 120
  # replicates this block alone runs for hours and pushes the whole suite past any
  # sane time limit, so it is gated rather than silently making `test_dir` unusable.
  #
  #   SELHOT_SLOW_TESTS=1 Rscript -e 'testthat::test_file(...)'
  #     (or sbatch slurm/run_tests_slow.sbatch)
  #
  # If the hclust arc solver is ever optimised, drop the gate.
  skip_if(!identical(Sys.getenv("SELHOT_SLOW_TESTS"), "1"),
          "slow: hclust exact-route uniformity check (set SELHOT_SLOW_TESTS=1)")
  set.seed(7)
  reps <- 120L
  clf <- hclust_cl_fun(3L)
  pvals <- rep(NA_real_, reps)
  for (i in seq_len(reps)) {
    X <- matrix(stats::rnorm(30 * 3), 30, 3)
    cl <- clf(X)
    pk <- pick_big_pair(cl)
    if (is.null(pk)) next
    if (!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2))) next
    res <- tryCatch(selective_truncated_F(X, pk$k1, pk$k2, clustering = "hclust",
                                        K = 3L, cl_obs = cl),
                    error = function(e) NULL)
    if (!is.null(res)) pvals[i] <- res$pval
  }
  ok <- pvals[is.finite(pvals)]
  expect_gt(length(ok), 30L)
  expect_gt(mean(ok), 0.34)
  expect_lt(mean(ok), 0.66)
  expect_lt(mean(ok <= 0.05), 0.16)
})

# ---- hclust_event = "trace": the lifetime/endpoint characterisation -----------
# Guards the theorem itself, not just the plumbing: the O(n^2) lifetime constraints
# must hold at t EXACTLY when the first n-K merges are unchanged at t. Without this
# a refactor of the Lance-Williams recursion breaks exactness silently.
test_that("trace lifetime constraints characterise the merge trace exactly", {
  e <- asNamespace("selectiveHotelling")
  # canonical, label-free signature of the first n-K merges
  trace_sig <- function(X, K) {
    n <- nrow(X); hc <- stats::hclust(stats::dist(X)^2, method = "average")
    mg <- hc$merge; node <- vector("list", n - 1L); sig <- character(0)
    for (s in seq_len(n - 1L)) {
      a <- mg[s, 1]; b <- mg[s, 2]
      Gi <- if (a < 0) -a else node[[a]]; Hi <- if (b < 0) -b else node[[b]]
      node[[s]] <- c(Gi, Hi)
      if (s <= n - K)
        sig <- c(sig, paste0(paste(sort(Gi), collapse = "."), "|",
                             paste(sort(Hi), collapse = ".")))
    }
    paste(sig, collapse = " ")
  }
  set.seed(7)
  for (K in c(3L, 4L)) {
    n <- 22L; p <- 3L
    X  <- matrix(stats::rnorm(n * p), n, p)
    cl <- as.integer(stats::cutree(stats::hclust(stats::dist(X)^2, method = "average"), k = K))
    tb <- sort(table(cl), decreasing = TRUE)
    if (length(tb) < 2L || tb[2] < 3L) next
    k1 <- as.integer(names(tb)[1]); k2 <- as.integer(names(tb)[2])
    red <- e$.truncF_reduce(X, cl, k1, k2)
    geo <- red$geo
    PUW <- e$.geodesic_PUW(red$X_of_Q, red$prep$X_fixed,
                           geo$M_const, geo$M_cos, geo$M_sin)
    hc  <- stats::hclust(stats::dist(X)^2, method = "average")
    co  <- e$.hclust_trace_constraints(
             e$.hclust_pointpair_harmonics(PUW$P, PUW$U, PUW$W), hc, K)
    expect_gt(length(co$A), 0L)
    ref <- trace_sig(X, K)
    ts  <- seq(1e-6, pi / 2 - 1e-6, length.out = 60L)
    pred <- vapply(ts, function(t)
      max(co$A * cos(2*t) + co$B * sin(2*t) + co$C * cos(t) +
          co$D * sin(t) + co$E) <= 1e-9, logical(1))
    act <- vapply(ts, function(t)
      trace_sig(red$X_of_Q(geo$Qtheta(t)), K) == ref, logical(1))
    expect_identical(pred, act)
    # theta_obs must always be feasible: the observed data is in its own event
    expect_true(max(co$A * cos(2*geo$theta_obs) + co$B * sin(2*geo$theta_obs) +
                    co$C * cos(geo$theta_obs) + co$D * sin(geo$theta_obs) +
                    co$E) <= 1e-8)
  }
})

test_that("trace conditioning refuses non-average linkage and labels its event", {
  # Average linkage chains off singletons, so a fixed seed often has no pair with
  # >= 3 members. Search rather than skip: skipping would make this test vacuous.
  X <- cl <- k1 <- k2 <- NULL
  for (s in 1:60) {
    set.seed(s)
    Xc <- matrix(stats::rnorm(24L * 3L), 24L, 3L)
    clc <- as.integer(stats::cutree(stats::hclust(stats::dist(Xc)^2, method = "average"), k = 3L))
    tb  <- sort(table(clc), decreasing = TRUE)
    if (length(tb) >= 2L && tb[2] >= 3L) {
      X <- Xc; cl <- clc
      k1 <- as.integer(names(tb)[1]); k2 <- as.integer(names(tb)[2])
      break
    }
  }
  expect_false(is.null(X))          # 60 seeds without a testable pair = real problem
  # the theorem is average-linkage only -- must error, not return a wrong arc set
  for (lk in c("centroid", "single", "complete"))
    expect_error(selective_truncated_F(X, k1 = k1, k2 = k2, clustering = "hclust",
                   method = "exact", K = 3L, cl_obs = cl, hclust_method = lk,
                   hclust_event = "trace"), "average linkage only")
  out <- selective_truncated_F(X, k1 = k1, k2 = k2, clustering = "hclust",
           method = "exact", K = 3L, cl_obs = cl, hclust_method = "average",
           hclust_event = "trace")
  # "trace" and "pair" are DIFFERENT events and must not share a label
  expect_identical(out$condition_on, "hclust-merge-trace")
  expect_true(is.na(out$pval) || (out$pval >= 0 && out$pval <= 1))
})
