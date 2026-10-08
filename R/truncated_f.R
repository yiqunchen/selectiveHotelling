# =============================================================================
# truncated_f.R -- the one-dimensional selective Hotelling test p_polar,1d.
#
# Conditions on every coordinate of the polar factor Q except the Hotelling angle
# theta; the selective law of theta is the null Beta law truncated to the angles at
# which the clustering returns the tested pair, and the p-value is the ratio of
# Beta masses below and above the observed angle. Two ways to compute it:
#   * method = "bc" (DEFAULT, any clustering; Algorithm 2 of the paper):
#       Besag-Clifford path on the angle (truncated_f_mc.R), random-step kernel
#       by default. Exact for any M. Conditions on the PAIR event.
#   * method = "exact" (closed form, no Monte Carlo), for two stronger events:
#       - clustering = "kmeans": the whole fixed-init Lloyd trajectory
#         (harmonic_arc.R::.kmeans_traj_coeffs);
#       - clustering = "hclust": the merge order up to the K-cut. Average linkage:
#         the observed first n-K merges (.hclust_trace_arcs). Single linkage: the
#         observed first n-K Kruskal (MST) edges (.hclust_single_kruskal_arcs).
#       Conditioning on these events costs power (aug-task 0v.16, 0v.30).
#   * method = "is": importance sampling on the angle (consistent, not exact).
# There is no grid approximation anywhere (author rule, 2026-10-05).
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
.truncF_lbmass <- function(alpha, beta) {

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
  function(lo, hi) {
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
}

.truncF_pvalue <- function(intervals, theta_obs, p, r) {
  alpha <- p / 2; beta <- (r - 1 - p) / 2
  lbmass <- .truncF_lbmass(alpha, beta)
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

# ---- EXACT hclust (merge order up to the K-cut) ----------------------------------
#' @keywords internal
#' @noRd
.truncF_exact_hclust <- function(X, k1, k2, K, hclust_method, hclust_square,
                                   cl_obs, ridge, verbose) {
  if (is.null(cl_obs))
    cl_obs <- hclust_fun(X, k = K, method = hclust_method, square = hclust_square)$cluster
  cl_obs <- as.integer(cl_obs)
  if (sum(cl_obs == k1) < 3L || sum(cl_obs == k2) < 3L)
    stop("Tested clusters must each have >= 3 members.")
  # The merge-order constraints are derived for AVERAGE linkage (lifetime/endpoint
  # theorem, .hclust_trace_constraints) and for SINGLE linkage via its Kruskal edges
  # (.hclust_single_kruskal_constraints). Other linkages have no closed form here:
  # use method = "bc", which reruns the clustering.
  if (!(hclust_method %in% c("average", "single")))
    stop(sprintf(paste0("method = 'exact' for hclust is available for average and single ",
                        "linkage only (got '%s'); use method = 'bc'."), hclust_method))
  red <- .truncF_reduce(X, cl_obs, k1, k2); geo <- red$geo
  PUW <- .geodesic_PUW(red$X_of_Q, red$prep$X_fixed, geo$M_const, geo$M_cos, geo$M_sin)
  single <- identical(hclust_method, "single")
  arcs <- if (single) {
    # single linkage is invariant to squaring, so the edges are found on squared
    # distances whatever hclust_square says
    D2 <- as.matrix(stats::dist(X))^2
    hc_obs <- stats::hclust(stats::as.dist(D2), method = "single")
    .hclust_single_kruskal_arcs(PUW$P, PUW$U, PUW$W, K, hc_obs, D2)
  } else {
    hc_obs <- stats::hclust(if (isTRUE(hclust_square)) stats::dist(X)^2 else stats::dist(X),
                            method = hclust_method)
    .hclust_trace_arcs(PUW$P, PUW$U, PUW$W, K, hc_obs)
  }
  if (verbose) cat(sprintf("[truncated-Fhclust-exact] r=%d p=%d theta=%.4f #arcs=%d\n",
                           red$r, red$p, geo$theta_obs, nrow(arcs)))
  .truncF_result(red, arcs, X, cl_obs, k1, k2,
                 if (single) "hclust-single-kruskal" else "hclust-merge-trace", ridge)
}

#' One-dimensional selective Hotelling test (unknown covariance)
#'
#' @description
#' The one-dimensional test \eqn{p_{\mathrm{polar,1d}}} of the paper: it conditions on
#' every coordinate of the polar factor \eqn{Q} except the Hotelling angle
#' \eqn{\theta}, so that the selective null law of \eqn{\theta} is a truncated Beta
#' (equivalently, a truncated F for the scaled \eqn{T^2}).
#'
#' @details
#' With \eqn{\cos^2\theta = u(Q)}, \eqn{T^2 = (r-2)\cot^2\theta}, and under the null
#' \eqn{\cos^2\theta \sim \mathrm{Beta}(p/2, (r-1-p)/2)} independently of the other
#' coordinates of \eqn{Q}. Given those coordinates, the clustering restricts
#' \eqn{\theta} to a set \eqn{S_\theta}, and the p-value is the Beta mass of
#' \eqn{S_\theta \cap (0, \theta_{obs}]} divided by the Beta mass of \eqn{S_\theta}.
#'
#' `method = "bc"` (default) computes it for any clustering, on the event that the
#' clustering returns the tested pair, with the Besag--Clifford construction on the
#' angle (Algorithm 2 of the paper); it is exact for any `M`. `method = "exact"`
#' gives a closed form, with no Monte Carlo, for two stronger events: the whole
#' Lloyd trajectory (`clustering = "kmeans"`) and the merge order up to the cut
#' (`clustering = "hclust"`, average or single linkage). These events condition on
#' more and cost power. `method = "is"` uses importance sampling (consistent, not
#' exact).
#'
#' @param X A numeric \eqn{n \times p} data matrix with observations in rows.
#' @param k1,k2 The two tested cluster labels. Each cluster must have at least
#'   three members, and \eqn{r - 1 > p}.
#' @param clustering `"kmeans"`, `"hclust"`, or `"custom"` (use `cl_fun`).
#' @param method `"bc"` (default), `"exact"` (closed form; `"kmeans"` and `"hclust"`
#'   only), or `"is"`.
#' @param K Number of clusters (built-in `"kmeans"`/`"hclust"`).
#' @param init_centers_idx Fixed initial-centre row indices for `"kmeans"`; drawn via
#'   [kmeans_fun()] if `NULL`.
#' @param hclust_method,hclust_square Linkage and squared-distance flag for `"hclust"`
#'   (defaults match [hclust_fun()]).
#' @param hclust_event Kept for compatibility; the only value is `"trace"`, the merge
#'   order used by `method = "exact"`.
#' @param cl_fun A clustering routine `function(Xmat) -> labels` (a pure function of
#'   its data) for `clustering = "custom"`; optional otherwise.
#' @param cl_obs Optional precomputed observed clustering (saves one call).
#' @param M,L Path length and number of kernel steps between path points for
#'   `method = "bc"`.
#' @param onedim_kernel For `method = "bc"`: `"givens"` (default; random step on the
#'   angle, step size drawn log-uniformly on `[1e-3, 1]` per block of `L` steps) or
#'   `"slice"` (slice sampling with shrinkage).
#' @param n_is,is_mix,is_scale Number of draws, weight of the null law in the
#'   defensive mixture, and the standard deviation of its normal part, for
#'   `method = "is"`.
#' @param ridge Small non-negative ridge for the reported observed \eqn{T^2}.
#' @param verbose If `TRUE`, print a one-line diagnostic.
#'
#' @return An S3 object of class `"selectiveHotelling"` with `pval`, `T_obs`, `k1`,
#'   `k2`, `K`, `condition_on`, `sampler`, `p_uncond` (the classical F p-value),
#'   `theta_obs`, `F_obs`, `df1`, `df2`, and, for `method = "exact"`, the feasible
#'   arcs `intervals`.
#'
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(40), 20, 2), matrix(rnorm(40) + 5, 20, 2),
#'            matrix(rnorm(40) - 5, 20, 2))
#' cl_fun <- function(x) hclust_fun(x, k = 3)$cluster
#' res <- selective_truncated_F(X, 1, 2, clustering = "custom", cl_fun = cl_fun, M = 200)
#' res$pval
#'
#' @seealso [selective_hotelling_bc()], [selective_hotelling_test()], [hotelling_T2()]
#' @export
selective_truncated_F <- function(X, k1, k2,
                                clustering = c("kmeans", "hclust", "custom"),
                                method = c("bc", "exact", "is"),
                                K = 3L, init_centers_idx = NULL,
                                hclust_method = "average", hclust_square = TRUE,
                                hclust_event = "trace",
                                cl_fun = NULL, cl_obs = NULL,
                                M = 5000L, L = 2L, onedim_kernel = c("givens", "slice"),
                                n_is = 2000L, is_mix = 0.5, is_scale = NULL,
                                ridge = 1e-8, verbose = FALSE) {
  stopifnot(is.matrix(X))
  clustering    <- match.arg(clustering)
  method        <- match.arg(method)
  onedim_kernel <- match.arg(onedim_kernel)
  if (!identical(hclust_event, "trace"))
    stop("hclust_event must be 'trace' (the pair-event exact march was removed; use method = 'bc').")
  mc <- function(cf) .truncF_mc(X, cf, k1, k2, cl_obs, method, as.integer(M), as.integer(L),
                                as.integer(n_is), is_mix, is_scale, ridge, verbose,
                                onedim_kernel = onedim_kernel)
  if (clustering == "custom") {
    if (is.null(cl_fun)) stop("clustering = 'custom' requires a `cl_fun`.")
    if (method == "exact") stop("method = 'exact' needs clustering = 'kmeans' or 'hclust'.")
    return(mc(cl_fun))
  }
  if (clustering == "kmeans") {
    if (method == "exact")
      return(.truncF_exact_kmeans(X, k1, k2, K, init_centers_idx, ridge, verbose))
    cf <- if (!is.null(cl_fun)) cl_fun else {
      ii <- if (is.null(init_centers_idx)) kmeans_fun(X, k = K)$init_centers_idx else init_centers_idx
      make_fast_kmeans_cl_fun(ii)                    # kmeans_with_fixed_init(), BLAS
    }
    return(mc(cf))
  }
  # clustering == "hclust"
  if (method == "exact")
    return(.truncF_exact_hclust(X, k1, k2, K, hclust_method, hclust_square,
                                cl_obs, ridge, verbose))
  cf <- if (!is.null(cl_fun)) cl_fun else
    make_hclust_cl_fun(K, hclust_method, hclust_square)   # hclust_fun(), fast when square
  mc(cf)
}
