# =============================================================================
# test-conditioning-set.R -- BRUTE-FORCE validation of the conditioning set.
#
# Every selective p-value here is a ratio of Beta masses over a feasible set S.
# If S is wrong the p-value is wrong, and nothing else in the pipeline will notice:
# the test suite passed for months while the hclust march was silently truncating S
# and returning a fabricated 1.0. These tests therefore check S itself, against the
# only ground truth available -- actually re-running the clustering at many theta
# and asking whether the conditioned-on event holds.
#
# One test per engine. Each asserts BOTH directions:
#   * no theta inside the computed arcs where the event FAILS  (S is not too big:
#     an inflated S inflates the denominator and makes the test ANTI-conservative)
#   * no theta outside them where the event HOLDS              (S is not too small:
#     a deflated S is what destroys power -- see the merge-trace negative result)
# =============================================================================

# reconstruct the geodesic pieces the arc engines are given
.cs_geo <- function(X, cl, k1, k2) {
  e   <- asNamespace("selectiveHotelling")
  red <- e$.truncF_reduce(X, cl, k1, k2)
  geo <- red$geo
  list(red = red, geo = geo,
       PUW = e$.geodesic_PUW(red$X_of_Q, red$prep$X_fixed,
                             geo$M_const, geo$M_cos, geo$M_sin))
}
.in_arcs <- function(arcs, t, tol = 1e-9) {
  if (is.null(arcs) || !nrow(arcs)) return(FALSE)
  any(t >= arcs[, 1] - tol & t <= arcs[, 2] + tol)
}
# a testable pair for average-linkage hclust (it chains, so search seeds)
.find_hclust_case <- function(n, p, K, seeds = 1:80) {
  for (s in seeds) {
    set.seed(s)
    X  <- matrix(stats::rnorm(n * p), n, p)
    cl <- as.integer(stats::cutree(stats::hclust(stats::dist(X)^2, method = "average"), k = K))
    tb <- sort(table(cl), decreasing = TRUE)
    if (length(tb) >= 2L && tb[2] >= 3L)
      return(list(X = X, cl = cl, k1 = as.integer(names(tb)[1]),
                  k2 = as.integer(names(tb)[2])))
  }
  NULL
}

test_that("hclust PAIR-PRESERVED arcs match brute-force re-clustering", {
  e <- asNamespace("selectiveHotelling")
  K <- 3L
  cse <- .find_hclust_case(20L, 3L, K)
  skip_if(is.null(cse), "no testable average-linkage pair in the search range")
  g <- .cs_geo(cse$X, cse$cl, cse$k1, cse$k2)
  arcs <- suppressWarnings(
    e$.hclust_pairpreserved_arcs(g$PUW$P, g$PUW$U, g$PUW$W, g$red$X_of_Q,
                                 g$geo$Qtheta, K, cse$cl, cse$k1, cse$k2,
                                 method = "average", square = TRUE))
  skip_if(isTRUE(attr(arcs, "capped")), "march capped; arc set is truncated by design")
  ts <- seq(1e-6, pi / 2 - 1e-6, length.out = 150L)
  held <- vapply(ts, function(t) {
    clt <- as.integer(stats::cutree(
      stats::hclust(stats::dist(g$red$X_of_Q(g$geo$Qtheta(t)))^2, method = "average"), k = K))
    e$.bc_pair_preserved(cse$cl, clt, cse$k1, cse$k2)
  }, logical(1))
  inside <- vapply(ts, function(t) .in_arcs(arcs, t), logical(1))
  expect_identical(inside, held)
  expect_true(.in_arcs(arcs, g$geo$theta_obs))   # observed data is always feasible
})

test_that("hclust MERGE-TRACE arcs match brute-force re-clustering", {
  e <- asNamespace("selectiveHotelling")
  K <- 3L
  cse <- .find_hclust_case(20L, 3L, K)
  skip_if(is.null(cse), "no testable average-linkage pair in the search range")
  g  <- .cs_geo(cse$X, cse$cl, cse$k1, cse$k2)
  hc <- stats::hclust(stats::dist(cse$X)^2, method = "average")
  arcs <- e$.hclust_trace_arcs(g$PUW$P, g$PUW$U, g$PUW$W, K, hc)
  # ground truth: the first n-K merges, as label-free member sets
  sig <- function(X) {
    n <- nrow(X); mg <- stats::hclust(stats::dist(X)^2, method = "average")$merge
    node <- vector("list", n - 1L); out <- character(0)
    for (s in seq_len(n - 1L)) {
      a <- mg[s, 1]; b <- mg[s, 2]
      Gi <- if (a < 0) -a else node[[a]]; Hi <- if (b < 0) -b else node[[b]]
      node[[s]] <- c(Gi, Hi)
      if (s <= n - K) out <- c(out, paste0(paste(sort(Gi), collapse = "."), "|",
                                           paste(sort(Hi), collapse = ".")))
    }
    paste(out, collapse = " ")
  }
  ref <- sig(cse$X)
  # the trace arc is narrow, so sample AROUND it rather than across [0, pi/2]
  w  <- if (nrow(arcs)) sum(arcs[, 2] - arcs[, 1]) else 1e-3
  lo <- max(g$geo$theta_obs - 20 * w, 1e-9)
  hi <- min(g$geo$theta_obs + 20 * w, pi / 2 - 1e-9)
  ts <- seq(lo, hi, length.out = 200L)
  held   <- vapply(ts, function(t) sig(g$red$X_of_Q(g$geo$Qtheta(t))) == ref, logical(1))
  inside <- vapply(ts, function(t) .in_arcs(arcs, t), logical(1))
  expect_identical(inside, held)
  expect_true(.in_arcs(arcs, g$geo$theta_obs))
})

test_that("k-means Lloyd-trajectory arcs match brute-force re-running Lloyd", {
  e <- asNamespace("selectiveHotelling")
  K <- 3L
  set.seed(4)
  X  <- matrix(stats::rnorm(30L * 3L), 30L, 3L)
  km <- kmeans_fun(X, k = K)
  idx <- as.integer(km$init_centers_idx)
  ref <- kmeans_with_fixed_init(X, idx)
  tb  <- sort(table(ref$cluster), decreasing = TRUE)
  skip_if(length(tb) < 2L || tb[2] < 3L, "no testable k-means pair for this seed")
  k1 <- as.integer(names(tb)[1]); k2 <- as.integer(names(tb)[2])
  g   <- .cs_geo(X, ref$cluster, k1, k2)
  co  <- e$.kmeans_traj_coeffs(g$PUW$P, g$PUW$U, g$PUW$W, idx, ref$trajectory, K)
  arcs <- e$.harmonic_feasible_arcs(co, 0, pi / 2)
  # ground truth: re-run Lloyd from the SAME init and demand the identical trajectory
  same_traj <- function(Xt) {
    r <- tryCatch(kmeans_with_fixed_init(Xt, idx), error = function(e) NULL)
    if (is.null(r)) return(FALSE)
    length(r$trajectory) == length(ref$trajectory) &&
      all(vapply(seq_along(ref$trajectory),
                 function(i) identical(as.integer(r$trajectory[[i]]),
                                       as.integer(ref$trajectory[[i]])), logical(1)))
  }
  ts <- seq(1e-6, pi / 2 - 1e-6, length.out = 150L)
  held   <- vapply(ts, function(t) same_traj(g$red$X_of_Q(g$geo$Qtheta(t))), logical(1))
  inside <- vapply(ts, function(t) .in_arcs(arcs, t), logical(1))
  expect_identical(inside, held)
  expect_true(.in_arcs(arcs, g$geo$theta_obs))
})

test_that("a capped hclust march is reported, never silently truncated", {
  # The march has a cell budget. If it runs out, the arc set is a strict subset of
  # the truth and the p-value must NOT be reported as a number -- historically this
  # returned exactly 1.0, which reads as a confident non-rejection.
  e <- asNamespace("selectiveHotelling")
  K <- 3L
  cse <- .find_hclust_case(20L, 3L, K)
  skip_if(is.null(cse), "no testable average-linkage pair in the search range")
  g <- .cs_geo(cse$X, cse$cl, cse$k1, cse$k2)
  arcs <- suppressWarnings(
    e$.hclust_pairpreserved_arcs(g$PUW$P, g$PUW$U, g$PUW$W, g$red$X_of_Q,
                                 g$geo$Qtheta, K, cse$cl, cse$k1, cse$k2,
                                 method = "average", square = TRUE,
                                 max_cells = 3L))          # force the cap
  expect_true(isTRUE(attr(arcs, "capped")))
  res <- e$.truncF_result(g$red, arcs, cse$X, cse$cl, cse$k1, cse$k2,
                          "pair-preserved", 1e-8)
  expect_true(is.na(res$pval))
  expect_true(isTRUE(res$capped))
})
