# =============================================================================
# selective_truncated_F: the one-dimensional selective test. Besag-Clifford route
# on the angle (default, any clustering) and the closed forms for the k-means
# trajectory and the hclust merge order. No grid approximation.
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

test_that("hclust route (default method = 'bc'): valid p-value and S3 object", {
  set.seed(101)
  X <- matrix(stats::rnorm(30 * 3), 30, 3)
  cl <- hclust_cl_fun(3L)(X)
  pk <- pick_big_pair(cl); skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))
  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "hclust", K = 3L, M = 200L)
  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval > 0 && res$pval <= 1)
  expect_identical(res$mc_method, "bc")
  expect_identical(res$onedim_kernel, "givens")
  expect_identical(res$df1, ncol(X))
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

  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "kmeans", method = "exact", K = 3L,
                             init_centers_idx = init_idx, ridge = 0)
  expect_s3_class(res, "selectiveHotelling")
  expect_true(res$pval >= 0 && res$pval <= 1)
  expect_identical(res$condition_on, "kmeans-trajectory")

  r <- res$df2 + res$df1 + 1L                # = r (pair size)
  T_from_angle <- (r - 2) * (cos(res$theta_obs) / sin(res$theta_obs))^2
  expect_equal(res$T_obs, T_from_angle, tolerance = 1e-6)
})

test_that("trivial (always-feasible) selection reduces to the classical F p-value", {
  # custom cl_fun routes to the exact Besag-Clifford path on the angle (no grid).
  # With an event that always holds, the path targets the unconditional angle law,
  # so its p-value is a Monte Carlo estimate of the classical F p-value. Labels are
  # assigned at random so that the classical p-value is moderate and resolvable.
  set.seed(202)
  X <- matrix(stats::rnorm(30 * 3), 30, 3)
  cl <- rep(1:3, length.out = 30)
  const_cl <- function(Xmat) cl          # ignores its input => pair always preserved
  res <- selective_truncated_F(X, 1L, 2L, clustering = "custom",
                             cl_fun = const_cl, cl_obs = cl, M = 5000L)
  expect_identical(res$mc_method, "bc")
  expect_gt(res$p_uncond, 0.05)
  expect_lt(abs(res$pval - res$p_uncond), 0.06)
})

test_that("custom cl_fun works and is algorithm-agnostic", {
  set.seed(404)
  X <- matrix(stats::rnorm(45 * 3), 45, 3)
  init_idx <- c(1L, 16L, 31L)
  km_cl_fun <- function(Xmat) kmeans_with_fixed_init(Xmat, init_centers_idx = init_idx)$cluster
  cl <- km_cl_fun(X)
  pk <- pick_big_pair(cl); skip_if(is.null(pk))
  skip_if(!is.finite(hotelling_T2(X, cl, pk$k1, pk$k2)))
  res <- selective_truncated_F(X, pk$k1, pk$k2, clustering = "custom",
                             cl_fun = km_cl_fun, cl_obs = cl, M = 200L)
  expect_s3_class(res, "selectiveHotelling")
  expect_error(selective_truncated_F(X, pk$k1, pk$k2, clustering = "custom",
                                     cl_fun = km_cl_fun, method = "exact"), "exact")
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

test_that("merge-order conditioning refuses other linkages and labels its event", {
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
  # (single linkage has its own Kruskal-edge trace event, tested below)
  for (lk in c("centroid", "complete"))
    expect_error(selective_truncated_F(X, k1 = k1, k2 = k2, clustering = "hclust",
                   method = "exact", K = 3L, cl_obs = cl, hclust_method = lk,
                   hclust_event = "trace"), "average and single")
  out <- selective_truncated_F(X, k1 = k1, k2 = k2, clustering = "hclust",
           method = "exact", K = 3L, cl_obs = cl, hclust_method = "average",
           hclust_event = "trace")
  # the merge-order event is labelled as such
  expect_identical(out$condition_on, "hclust-merge-trace")
  expect_true(is.na(out$pval) || (out$pval >= 0 && out$pval <= 1))
})

# ---- single linkage, hclust_event = "trace": the Kruskal-edge event ----------
# The closed-form arcs must agree with a brute-force check that the first n-K
# Kruskal (MST) edges are the observed point pairs, in order, at every grid angle.
test_that("single-linkage Kruskal-edge arcs match brute force", {
  e <- asNamespace("selectiveHotelling")
  kedges <- function(X, q) {
    D2 <- as.matrix(stats::dist(X))^2; n <- nrow(X)
    ut <- which(upper.tri(D2)); o <- ut[order(D2[ut])]
    comp <- seq_len(n); out <- character(0)
    for (k in o) {
      i <- (k - 1L) %% n + 1L; j <- (k - 1L) %/% n + 1L
      if (comp[i] != comp[j]) {
        out <- c(out, paste(min(i, j), max(i, j)))
        comp[comp == comp[j]] <- comp[i]
        if (length(out) == q) break
      }
    }
    paste(out, collapse = ",")
  }
  n_checked <- 0L
  for (s in 1:40) {
    set.seed(s)
    n <- 24L; p <- 2L
    X <- matrix(stats::rnorm(n * p), n, p); X[, 1] <- X[, 1] + rep(c(-5, 0, 5), 8)
    cl <- hclust_fun(X, k = 3, method = "single")$cluster
    tb <- sort(table(cl), decreasing = TRUE)
    if (length(tb) < 2L || tb[2] < 3L) next
    k1 <- as.integer(names(tb)[1]); k2 <- as.integer(names(tb)[2])
    out <- selective_truncated_F(X, k1 = k1, k2 = k2, clustering = "hclust",
             method = "exact", K = 3L, cl_obs = cl, hclust_method = "single", hclust_event = "trace")
    expect_identical(out$condition_on, "hclust-single-kruskal")
    expect_true(out$pval >= 0 && out$pval <= 1)
    arcs <- out$intervals
    inside <- function(t) any(t >= arcs[, 1] & t <= arcs[, 2])
    expect_true(inside(out$theta_obs))
    red <- e$.truncF_reduce(X, cl, k1, k2); geo <- red$geo
    ref <- kedges(X, n - 3L)
    w <- arcs[, 2] - arcs[, 1]
    ts <- c(seq(1e-6, pi / 2 - 1e-6, length.out = 150L),
            seq(arcs[1, 1] - w[1], arcs[1, 2] + w[1], length.out = 61L))
    ts <- ts[ts > 0 & ts < pi / 2]
    # skip points within 1e-9 of a boundary, where the brute force is a coin flip
    ts <- ts[vapply(ts, function(t) min(abs(t - arcs)) > 1e-9, logical(1))]
    act  <- vapply(ts, function(t) kedges(red$X_of_Q(geo$Qtheta(t)), n - 3L) == ref, logical(1))
    pred <- vapply(ts, inside, logical(1))
    expect_identical(pred, act)
    n_checked <- n_checked + 1L
    if (n_checked >= 3L) break
  }
  expect_gte(n_checked, 3L)
})

test_that("sweep arc solver agrees with the scalar solver", {
  e <- asNamespace("selectiveHotelling")
  set.seed(11)
  for (rep in 1:5) {
    nC <- 40L
    co <- list(A = stats::rnorm(nC), B = stats::rnorm(nC), C = stats::rnorm(nC),
               D = stats::rnorm(nC), E = stats::rnorm(nC) - 3)
    a1 <- e$.harmonic_feasible_arcs(co, 0, pi / 2)
    a2 <- e$.harmonic_feasible_arcs_sweep(co, 0, pi / 2)
    expect_equal(unname(a1), unname(a2), tolerance = 1e-8)
  }
})
