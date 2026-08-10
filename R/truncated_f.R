# =============================================================================
# truncated_f.R -- CLOSED-FORM one-dimensional truncated-F selective Hotelling test.
#
# The strong-conditioning, Monte-Carlo-free companion of selective_hotelling_bc.
# Where the Besag-Clifford test leaves the polar factor Q UNIFORM on a
# high-dimensional feasible subset of the Stiefel manifold (needing simulation),
# this test conditions ADDITIONALLY on every angular nuisance coordinate of Q
# except the single Hotelling angle theta. After that, theta is the ONLY random
# coordinate, its law is a 1-D truncated F, and the selective p-value is a
# closed-form ratio of Beta masses -- exact, with NO sampling and NO mixing.
#
# THREE feasible-set engines (all share the truncated-Beta p-value):
#   * clustering = "kmeans", method = "exact"  -> EXACT closed-form arc on the
#       fixed-init Lloyd TRAJECTORY (harmonic_arc.R::.kmeans_traj_coeffs), the
#       same E_tau the k-means high-dim-Q test conditions on.
#   * clustering = "hclust", method = "exact"  -> EXACT closed-form arc on the
#       PAIR-PRESERVED (cutree-K partition) event via event-driven cell traversal
#       (harmonic_arc.R::.hclust_pairpreserved_arcs). This is the event
#       selective_hotelling_bc and Gao et al. (clusterpval) condition on; the
#       stronger agglomeration-trajectory event over-conditions and loses power.
#   * method = "grid" (or clustering = "custom")  -> algorithm-agnostic feasible
#       set of an arbitrary cl_fun via dense grid + bisection (pair-preserved).
# =============================================================================

# -----------------------------------------------------------------------------
# .truncF_geodesic(Q_obs, v, tol): the Hotelling-angle 1-D slice.
# With a = v/||v|| the unit contrast direction, cos^2 theta = ||Q'a||^2 and the
# geodesic Q(theta) = M_const + cos theta*M_cos + sin theta*M_sin holds every
# angular nuisance fixed while sweeping theta, satisfying Q(theta_obs) = Q_obs,
# Q(theta) in St(r-1,p), and ||Q(theta)'a||^2 = cos^2 theta.
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.truncF_geodesic <- function(Q_obs, v, tol = 1e-9) {
  cvar  <- sum(v^2)                                   # = 1/n1 + 1/n2
  a_vec <- v / sqrt(cvar)
  w_a   <- as.numeric(crossprod(Q_obs, a_vec))        # ||w_a|| = cos theta_obs
  cos0  <- min(sqrt(sum(w_a^2)), 1)
  sin0  <- sqrt(max(0, 1 - cos0^2))
  theta_obs <- acos(cos0)
  degenerate <- (cos0 < tol) || (sin0 < tol)

  gamma <- if (cos0 > tol) w_a / cos0 else { g <- w_a; g[1] <- g[1] + 1; g / sqrt(sum(g^2)) }
  u_col <- as.numeric(Q_obs %*% gamma)
  b_vec <- if (sin0 > tol) (cos0 * a_vec - u_col) / sin0 else a_vec

  M_cos   <- outer(a_vec,  gamma)
  M_sin   <- outer(-b_vec, gamma)
  M_const <- Q_obs - cos0 * M_cos - sin0 * M_sin

  list(a_vec = a_vec, cos0 = cos0, sin0 = sin0, theta_obs = theta_obs,
       degenerate = degenerate, cvar = cvar,
       M_const = M_const, M_cos = M_cos, M_sin = M_sin,
       Qtheta = function(th) M_const + cos(th) * M_cos + sin(th) * M_sin)
}

# ---- closed-form upper-tail selective p-value from truncated Beta(H=cos^2 th) --
#' @keywords internal
#' @noRd
.truncF_pvalue <- function(intervals, theta_obs, p, r) {
  alpha <- p / 2; beta <- (r - 1 - p) / 2

  # LOG-SPACE Beta mass. The natural form,
  #     pbeta(cos(lo)^2) - pbeta(cos(hi)^2),
  # is a difference of two nearly-equal CDF values. The feasible theta-arc
  # shrinks as n grows (the Lloyd trajectory contributes ~n(K-1) constraints per
  # iteration), so both terms converge, the difference underflows, and the
  # p-value degenerates to 0/0 = NA. Measured on null iso p=10 data:
  #     n=50   fine
  #     n=100  total feasible mass already 1e-13 .. 1e-15  (5/6 replicates survive)
  #     n=200  total feasible mass EXACTLY 0                (0/6 survive)
  # which wiped out all six n=200 cells of the unified sweep.
  #
  # log(e^la - e^lb) = la + log1p(-exp(lb - la)) is algebraically identical and
  # stays accurate at log-masses far below where the raw difference underflows,
  # because pbeta(log.p = TRUE) is computed in log space throughout. Only the
  # RATIO is needed, so the common scale cancels exactly.
  # TAIL CHOICE MATTERS. theta near 0.5 rad puts cos^2(theta) ~ 0.72, while a
  # Beta(p/2, (r-1-p)/2) with large r has mean ~0.08 -- so BOTH endpoint CDFs sit
  # at essentially 1, and a lower-tail log difference collapses (log(1)-log(1)).
  # Whichever tail holds the SMALL probability is the one that carries the
  # resolution, so pick it by where the interval sits. For an arc so narrow that
  # even the correct tail cancels, fall back to density x width, which is exact
  # to first order and cannot cancel at all. Only the RATIO is used and both
  # numerator and denominator go through this same routine, so any residual
  # first-order error largely cancels.
  lbmass <- function(lo, hi) {
    if (hi <= lo) return(-Inf)
    xl <- cos(hi)^2                      # smaller x  (theta larger)
    xu <- cos(lo)^2                      # larger  x  (theta smaller)
    if (!(xu > xl)) return(-Inf)
    mid <- (xl + xu) / 2
    upper <- stats::pbeta(mid, alpha, beta) > 0.5
    if (upper) {                          # mass = Q(xl) - Q(xu),  Q = upper tail
      la <- stats::pbeta(xl, alpha, beta, lower.tail = FALSE, log.p = TRUE)
      lb <- stats::pbeta(xu, alpha, beta, lower.tail = FALSE, log.p = TRUE)
    } else {                              # mass = F(xu) - F(xl)
      la <- stats::pbeta(xu, alpha, beta, log.p = TRUE)
      lb <- stats::pbeta(xl, alpha, beta, log.p = TRUE)
    }
    if (is.finite(la) && lb < la) {
      out <- la + log1p(-exp(lb - la))
      if (is.finite(out)) return(out)
    }
    # narrow-arc fallback: log( f(mid) * (xu - xl) )
    stats::dbeta(mid, alpha, beta, log = TRUE) + log(xu - xl)
  }
  logsumexp <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(-Inf)
    m <- max(v); m + log(sum(exp(v - m)))
  }

  ld <- numeric(0); ln <- numeric(0)
  if (!is.null(intervals) && nrow(intervals)) {
    for (k in seq_len(nrow(intervals))) {
      lo <- intervals[k, 1]; hi <- intervals[k, 2]
      ld <- c(ld, lbmass(lo, hi))
      if (lo < theta_obs) ln <- c(ln, lbmass(lo, min(hi, theta_obs)))
    }
  }
  lden <- logsumexp(ld); lnum <- logsumexp(ln)

  # The numerator integrates a SUBSET of the denominator's arcs, so the ratio is
  # in [0,1] mathematically; min(0, .) in log space enforces that against
  # rounding, and the [0,1] clamp is a belt-and-braces guard on the exponential.
  pval <- if (is.finite(lden)) exp(min(0, lnum - lden)) else NA_real_
  if (is.finite(pval)) pval <- max(0, min(1, pval))
  # NOTE: a genuinely uncomputable case stays NA rather than being coerced to 0
  # or 1. Coercing would fabricate a p-value and silently bias Type-I and power;
  # an NA is counted and reported by the aggregators instead.

  list(pval = pval, log_denom = lden, log_numer = lnum,
       denom_mass = exp(lden), numer_mass = exp(lnum),
       alpha = alpha, beta = beta)
}

# ---- feasible theta-set of a black-box predicate: grid + bisection ------------
#' @keywords internal
#' @noRd
.truncF_bisect_boundary <- function(f, a, b, tol) {
  fa <- isTRUE(f(a))
  while (b - a > tol) {
    mid <- 0.5 * (a + b)
    if (isTRUE(f(mid)) == fa) a <- mid else b <- mid
  }
  0.5 * (a + b)
}

#' @keywords internal
#' @noRd
.truncF_intervals_grid <- function(feasible_fn, theta_obs,
                                     n_grid = 400L, tol = 1e-6,
                                     lo = 0, hi = pi / 2) {
  grid <- sort(unique(c(seq(lo, hi, length.out = n_grid), theta_obs)))
  fe   <- vapply(grid, function(t) isTRUE(feasible_fn(t)), logical(1))
  n <- length(grid); s <- 1L; out <- list()
  while (s <= n) {
    if (!fe[s]) { s <- s + 1L; next }
    e <- s
    while (e + 1L <= n && fe[e + 1L]) e <- e + 1L
    left  <- if (s == 1L) grid[1] else .truncF_bisect_boundary(feasible_fn, grid[s - 1L], grid[s], tol)
    right <- if (e == n)  grid[n] else .truncF_bisect_boundary(feasible_fn, grid[e], grid[e + 1L], tol)
    out[[length(out) + 1L]] <- c(left, right)
    s <- e + 1L
  }
  if (!length(out)) return(matrix(numeric(0), ncol = 2))
  m <- do.call(rbind, out); colnames(m) <- c("t_lo", "t_hi"); m
}

# ---- shared reduction: prep, X_of_Q, contrast, truncated-F geodesic -----------
#' @keywords internal
#' @noRd
.truncF_reduce <- function(X, cl_obs, k1, k2) {
  prep <- prepare_pair_objects(X, cl_obs, k1, k2)
  m <- prep$m; r <- prep$r; p <- ncol(X)
  if (r - 1 <= p)
    stop(sprintf("truncated-F test needs pair size r-1 > p (got r=%d, p=%d).", r, p))
  X_of_Q <- function(Q) {
    AQS <- prep$A %*% Q %*% prep$S_sqrt
    XQ <- prep$X_fixed
    XQ[prep$idx1, ] <- XQ[prep$idx1, ] + AQS[1:m, , drop = FALSE]
    XQ[prep$idx2, ] <- XQ[prep$idx2, ] + AQS[(m + 1):r, , drop = FALSE]
    XQ
  }
  contrast <- c(rep(1 / m, m), rep(-1 / (r - m), r - m))
  v   <- as.numeric(crossprod(prep$A, contrast))
  geo <- .truncF_geodesic(prep$Q_obs, v)
  list(prep = prep, X_of_Q = X_of_Q, geo = geo, p = p, r = r)
}

# ---- assemble the S3 result from feasible arcs --------------------------------
#' @keywords internal
#' @noRd
.truncF_result <- function(red, arcs, X, cl_obs, k1, k2, condition_on,
                             ridge, extra = list()) {
  geo <- red$geo; p <- red$p; r <- red$r
  pv <- .truncF_pvalue(arcs, geo$theta_obs, p, r)
  T_obs    <- hotelling_T2(X, cl_obs, k1, k2, ridge)
  F_obs    <- ((r - p - 1) / p) * (cos(geo$theta_obs) / sin(geo$theta_obs))^2
  p_uncond <- 1 - stats::pf(F_obs, p, r - p - 1)
  out <- c(list(
    pval = pv$pval, se = NA_real_, T_obs = T_obs,
    draws = NULL, acc = NA_real_, n_unique_draws = NA_integer_,
    k1 = k1, k2 = k2, K = length(unique(cl_obs)),
    condition_on = condition_on, sampler = "truncated-F",
    p_uncond = p_uncond, theta_obs = geo$theta_obs, F_obs = F_obs,
    df1 = p, df2 = r - p - 1, intervals = arcs,
    denom_mass = pv$denom_mass, n_intervals = nrow(arcs),
    degenerate = geo$degenerate), extra)
  class(out) <- "selectiveHotelling"
  out
}

# ---- EXACT k-means (Lloyd trajectory E_tau) -----------------------------------
#' @keywords internal
#' @noRd
.truncF_exact_kmeans <- function(X, k1, k2, K, init_centers_idx, ridge, verbose) {
  if (is.null(init_centers_idx)) init_centers_idx <- kmeans_fun(X, k = K)$init_centers_idx
  init_centers_idx <- as.integer(init_centers_idx)
  ref  <- kmeans_with_fixed_init(X, init_centers_idx)
  cl_obs <- ref$cluster; traj <- ref$trajectory
  if (sum(cl_obs == k1) < 3L || sum(cl_obs == k2) < 3L)
    stop("Tested clusters must each have >= 3 members.")
  red <- .truncF_reduce(X, cl_obs, k1, k2); geo <- red$geo
  PUW <- .geodesic_PUW(red$X_of_Q, red$prep$X_fixed, geo$M_const, geo$M_cos, geo$M_sin)
  co  <- .kmeans_traj_coeffs(PUW$P, PUW$U, PUW$W, init_centers_idx, traj, K)
  arcs <- .harmonic_feasible_arcs(co, 0, pi / 2)
  if (verbose) cat(sprintf("[truncated-Fkmeans-exact] r=%d p=%d theta=%.4f #arcs=%d\n",
                           red$r, red$p, geo$theta_obs, nrow(arcs)))
  .truncF_result(red, arcs, X, cl_obs, k1, k2, "kmeans-trajectory", ridge)
}

# ---- EXACT hclust (pair-preserved via event-driven cell traversal) ------------
#' @keywords internal
#' @noRd
.truncF_exact_hclust <- function(X, k1, k2, K, hclust_method, hclust_square,
                                   cl_obs, ridge, verbose) {
  if (is.null(cl_obs))
    cl_obs <- hclust_fun(X, k = K, method = hclust_method, square = hclust_square)$cluster
  cl_obs <- as.integer(cl_obs)
  if (sum(cl_obs == k1) < 3L || sum(cl_obs == k2) < 3L)
    stop("Tested clusters must each have >= 3 members.")
  # SAFETY: closed-form harmonic breakpoints (.hclust_traj_coeffs) exist for the
  # single-harmonic linkages (average, centroid) and the min/max-envelope linkages
  # (single, complete). For any OTHER linkage (median/McQuitty/Ward etc.), fall
  # back to the grid engine, which is linkage-correct because it re-runs the true
  # clustering rule at each probe. (Never silently substitute a wrong linkage.)
  if (!(hclust_method %in% c("average", "centroid", "single", "complete"))) {
    if (isTRUE(verbose))
      cat(sprintf("[truncated-Fhclust-exact] linkage '%s' has no closed-form arc; using grid engine.\n",
                  hclust_method))
    cf <- function(x) hclust_fun(x, k = K, method = hclust_method, square = hclust_square)$cluster
    return(.truncF_grid(X, cf, k1, k2, cl_obs,
                          n_grid = 400L, boundary_tol = 1e-6, ridge, verbose))
  }
  red <- .truncF_reduce(X, cl_obs, k1, k2); geo <- red$geo
  PUW <- .geodesic_PUW(red$X_of_Q, red$prep$X_fixed, geo$M_const, geo$M_cos, geo$M_sin)
  arcs <- .hclust_pairpreserved_arcs(PUW$P, PUW$U, PUW$W, red$X_of_Q, geo$Qtheta,
                                     K, cl_obs, k1, k2, hclust_method, hclust_square)
  if (verbose) cat(sprintf("[truncated-Fhclust-exact] r=%d p=%d theta=%.4f #arcs=%d\n",
                           red$r, red$p, geo$theta_obs, nrow(arcs)))
  .truncF_result(red, arcs, X, cl_obs, k1, k2, "pair-preserved", ridge)
}

# ---- GRID (algorithm-agnostic, pair-preserved) --------------------------------
#' @keywords internal
#' @noRd
.truncF_grid <- function(X, cl_fun, k1, k2, cl_obs, n_grid, boundary_tol,
                           ridge, verbose) {
  .as_labels <- function(z) if (is.list(z) && !is.null(z$cluster)) z$cluster else z
  if (is.null(cl_obs)) cl_obs <- .as_labels(cl_fun(X))
  cl_obs <- as.integer(cl_obs)
  if (length(cl_obs) != nrow(X)) stop("cl_fun(X) length != nrow(X).")
  if (sum(cl_obs == k1) < 3L || sum(cl_obs == k2) < 3L)
    stop("Tested clusters must each have >= 3 members.")
  red <- .truncF_reduce(X, cl_obs, k1, k2); geo <- red$geo
  feasible_theta <- function(th) {
    cl_new <- tryCatch(.as_labels(cl_fun(red$X_of_Q(geo$Qtheta(th)))), error = function(e) NULL)
    if (is.null(cl_new) || length(cl_new) != length(cl_obs)) return(FALSE)
    .bc_pair_preserved(cl_obs, cl_new, k1, k2)
  }
  if (!isTRUE(feasible_theta(geo$theta_obs)))
    stop("GUARD FAILED: observed theta not feasible (is cl_fun a pure function of its data?).")
  arcs <- .truncF_intervals_grid(feasible_theta, geo$theta_obs,
                                   n_grid = as.integer(n_grid), tol = boundary_tol)
  if (verbose) cat(sprintf("[truncated-Fgrid] r=%d p=%d theta=%.4f #arcs=%d\n",
                           red$r, red$p, geo$theta_obs, nrow(arcs)))
  .truncF_result(red, arcs, X, cl_obs, k1, k2, "pair-preserved", ridge)
}

#' Closed-form one-dimensional truncated-F selective Hotelling test (unknown covariance)
#'
#' @description
#' An **exact, Monte-Carlo-free** selective test of equality of the means of a
#' discovered cluster pair `(k1, k2)`, with the covariance \eqn{\Sigma} treated as
#' **unknown**. It is the strong-conditioning, closed-form companion of
#' [selective_hotelling_bc()]: instead of a Besag-Clifford Stiefel path it
#' conditions additionally on every angular nuisance coordinate of the polar
#' factor \eqn{Q} except the single Hotelling angle \eqn{\theta}, whose law is
#' then a 1-D truncated F.
#'
#' @details
#' Reduce the tested pair via the polar decomposition
#' \eqn{X_{sub} = X_{fixed} + A Q S^{1/2}} (see [prepare_pair_objects()]). With
#' \eqn{a} the unit two-sample-contrast direction, \eqn{\cos^2\theta = \|Q'a\|^2},
#' the two-sample Hotelling statistic is \eqn{T^2 = (r-2)\cot^2\theta} and
#' \eqn{F = \tfrac{r-p-1}{p}\cot^2\theta \sim F_{p,\,r-p-1}} unconditionally,
#' while under the null \eqn{\cos^2\theta \sim \mathrm{Beta}(p/2, (r-1-p)/2)}
#' independently of the other angular coordinates of \eqn{Q} (\eqn{r = n_1 + n_2}
#' is the pair size). Holding those coordinates fixed, \eqn{Q(\theta)} traces a
#' geodesic on which the reconstructed data \eqn{X(\theta)} is a harmonic curve,
#' so the clustering-preservation event restricts \eqn{\theta} to a feasible set
#' \eqn{I_E \subseteq [0,\pi/2]} that is computed **in closed form** for k-means
#' and hierarchical clustering (and by grid+bisection for an arbitrary `cl_fun`).
#' Since \eqn{T^2 \ge T^2_{obs} \iff \theta \le \theta_{obs}}, the selective
#' p-value is the truncated-Beta ratio
#' \deqn{p = \frac{\mathrm{mass}\{\theta \in I_E,\ \theta \le \theta_{obs}\}}{\mathrm{mass}\{\theta \in I_E\}}.}
#' No sampling, no burn-in, no mixing assumption.
#'
#' **Conditioning events.** `clustering = "kmeans"` conditions on the fixed-init
#' Lloyd **trajectory** \eqn{E_\tau} (the k-means selection event); `"hclust"`
#' conditions on the **cutree-K partition preserving the pair** -- the same event
#' [selective_hotelling_bc()] and Gao et al. use (the stronger
#' agglomeration-trajectory event over-conditions and loses power, so it is not
#' used). Both closed-form feasible sets are verified against brute-force
#' re-clustering.
#'
#' **What conditioning buys and costs.** Fixing the contrast *direction* makes the
#' p-value a pivot free of \eqn{\Sigma}, so on this fiber knowing \eqn{\Sigma}
#' cannot improve the test; the price is power (the high-dimensional
#' [selective_hotelling_bc()] keeps directional freedom and can be lifted by a
#' known \eqn{\Sigma}). In exchange this test is exact, closed-form, and needs no
#' path-length or step-size tuning.
#'
#' @param X A numeric \eqn{n \times p} data matrix with observations in rows.
#' @param k1,k2 The two tested cluster labels. Each cluster must have at least
#'   three members, and the pair size must exceed the dimension
#'   (\eqn{r - 1 > p}) for a truncated F to exist.
#' @param clustering One of `"kmeans"`, `"hclust"`, or `"custom"` (use `cl_fun`).
#' @param method `"exact"` (closed-form arc; default) or `"grid"` (dense grid +
#'   bisection on the black-box pair-preserved predicate). `"custom"` always uses
#'   the grid engine.
#' @param K Number of clusters (for the built-in `"kmeans"`/`"hclust"` routines).
#' @param init_centers_idx Fixed initial-centre row indices for the `"kmeans"`
#'   route (conditions on that exact fixed-init run); drawn via [kmeans_fun()] if
#'   `NULL`.
#' @param hclust_method,hclust_square Agglomeration method and squared-distance
#'   flag for the `"hclust"` route (defaults match [hclust_fun()]). The exact
#'   closed-form arc engine is derived for **average** linkage on squared
#'   distances; for any other linkage the `"exact"` route automatically falls
#'   back to the linkage-correct `"grid"` engine (single/complete linkage make
#'   each merge criterion a min/max envelope of harmonics, and
#'   Ward/centroid/median/McQuitty use different Lance-Williams weights, so the
#'   single-harmonic reduction does not apply as-is).
#' @param cl_fun A black-box clustering routine `function(Xmat) -> labels` for
#'   `clustering = "custom"` (must be a pure function of its data). Also used as
#'   the predicate for `method = "grid"`.
#' @param cl_obs Optional precomputed observed clustering (saves one call).
#' @param n_grid,boundary_tol Grid size and bisection tolerance for the grid
#'   engine.
#' @param ridge Small non-negative ridge for the reported observed \eqn{T^2}.
#' @param verbose If `TRUE`, print a one-line diagnostic.
#'
#' @return An S3 object of class `"selectiveHotelling"` with primary fields
#'   `pval` (closed-form truncated-F selective p-value), `T_obs`, `k1`, `k2`,
#'   `K`, `condition_on` (`"kmeans-trajectory"` or `"pair-preserved"`), and
#'   `sampler` (`"truncated-F"`). Additional fields: `p_uncond` (the
#'   classical, selection-ignoring F p-value), `theta_obs`, `F_obs`, `df1`
#'   (\eqn{= p}), `df2` (\eqn{= r-p-1}), `intervals` (feasible \eqn{\theta}-arcs),
#'   `denom_mass`, `se` (`NA`; exact, not Monte Carlo), and `n_intervals`.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(30 * 3), 30, 3)
#' # exact hierarchical-clustering route:
#' cl <- hclust_fun(X, k = 3)$cluster
#' big <- as.integer(names(table(cl))[table(cl) >= 3])
#' res <- selective_truncated_F(X, big[1], big[2], clustering = "hclust", K = 3)
#' res$pval
#' print(res)
#'
#' @seealso [selective_hotelling_bc()], [selective_hotelling_test()],
#'   [hotelling_T2()]
#' @export
selective_truncated_F <- function(X, k1, k2,
                                clustering = c("kmeans", "hclust", "custom"),
                                method = c("exact", "grid"),
                                K = 3L, init_centers_idx = NULL,
                                hclust_method = "average", hclust_square = TRUE,
                                cl_fun = NULL, cl_obs = NULL,
                                n_grid = 400L, boundary_tol = 1e-6,
                                ridge = 1e-8, verbose = FALSE) {
  stopifnot(is.matrix(X))
  clustering <- match.arg(clustering)
  method     <- match.arg(method)

  if (clustering == "custom") {
    if (is.null(cl_fun)) stop("clustering = 'custom' requires a `cl_fun`.")
    return(.truncF_grid(X, cl_fun, k1, k2, cl_obs, n_grid, boundary_tol, ridge, verbose))
  }
  if (clustering == "kmeans") {
    if (method == "exact")
      return(.truncF_exact_kmeans(X, k1, k2, K, init_centers_idx, ridge, verbose))
    cf <- if (!is.null(cl_fun)) cl_fun else {
      ii <- if (is.null(init_centers_idx)) kmeans_fun(X, k = K)$init_centers_idx else init_centers_idx
      function(x) kmeans_with_fixed_init(x, init_centers_idx = ii)$cluster
    }
    return(.truncF_grid(X, cf, k1, k2, cl_obs, n_grid, boundary_tol, ridge, verbose))
  }
  # clustering == "hclust"
  if (method == "exact")
    return(.truncF_exact_hclust(X, k1, k2, K, hclust_method, hclust_square,
                                  cl_obs, ridge, verbose))
  cf <- if (!is.null(cl_fun)) cl_fun else
    function(x) hclust_fun(x, k = K, method = hclust_method, square = hclust_square)$cluster
  .truncF_grid(X, cf, k1, k2, cl_obs, n_grid, boundary_tol, ridge, verbose)
}
