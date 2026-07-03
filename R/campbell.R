# =============================================================================
# campbell.R -- CLOSED-FORM one-dimensional truncated-F selective Hotelling test.
#
# The strong-conditioning, Monte-Carlo-free companion of selective_hotelling_bc.
# Where the Besag-Clifford test leaves the polar factor Q UNIFORM on a
# high-dimensional feasible subset of the Stiefel manifold (needing simulation),
# this test conditions ADDITIONALLY on every angular nuisance coordinate of Q
# except the single Hotelling angle theta. After that, theta is the ONLY random
# coordinate, its law is a 1-D truncated F, and the selective p-value is a
# closed-form ratio of Beta masses -- exact, with NO sampling and NO mixing.
# =============================================================================

# -----------------------------------------------------------------------------
# .campbell_geodesic(Q_obs, v, tol): the Hotelling-angle 1-D slice.
# With a = v/||v|| the unit contrast direction, cos^2 theta = ||Q'a||^2 and the
# geodesic Q(theta) = M_const + cos theta*M_cos + sin theta*M_sin holds every
# angular nuisance fixed while sweeping theta, satisfying Q(theta_obs) = Q_obs,
# Q(theta) in St(r-1,p), and ||Q(theta)'a||^2 = cos^2 theta.
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.campbell_geodesic <- function(Q_obs, v, tol = 1e-9) {
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
       Qtheta = function(th) M_const + cos(th) * M_cos + sin(th) * M_sin)
}

# -----------------------------------------------------------------------------
# .campbell_truncated_pvalue(intervals, theta_obs, p, r): closed-form upper-tail
# (large-T^2) selective p-value from the truncated Beta law of H = cos^2 theta ~
# Beta(p/2, (r-1-p)/2). T^2 >= T^2_obs  <=>  theta <= theta_obs.
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.campbell_truncated_pvalue <- function(intervals, theta_obs, p, r) {
  alpha <- p / 2; beta <- (r - 1 - p) / 2
  bmass <- function(lo, hi) {
    if (hi <= lo) return(0)
    stats::pbeta(cos(lo)^2, alpha, beta) - stats::pbeta(cos(hi)^2, alpha, beta)
  }
  denom <- 0; numer <- 0
  if (!is.null(intervals) && nrow(intervals)) {
    for (k in seq_len(nrow(intervals))) {
      lo <- intervals[k, 1]; hi <- intervals[k, 2]
      denom <- denom + bmass(lo, hi)
      if (lo < theta_obs) numer <- numer + bmass(lo, min(hi, theta_obs))
    }
  }
  list(pval = if (denom > 0) as.numeric(numer / denom) else NA_real_,
       denom_mass = as.numeric(denom), numer_mass = as.numeric(numer),
       alpha = alpha, beta = beta)
}

# -----------------------------------------------------------------------------
# Feasible theta-set of a black-box predicate on [lo,hi]: dense grid + bisection.
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.campbell_bisect_boundary <- function(f, a, b, tol) {
  fa <- isTRUE(f(a))
  while (b - a > tol) {
    mid <- 0.5 * (a + b)
    if (isTRUE(f(mid)) == fa) a <- mid else b <- mid
  }
  0.5 * (a + b)
}

#' @keywords internal
#' @noRd
.campbell_intervals_grid <- function(feasible_fn, theta_obs,
                                     n_grid = 400L, tol = 1e-6,
                                     lo = 0, hi = pi / 2) {
  grid <- sort(unique(c(seq(lo, hi, length.out = n_grid), theta_obs)))
  fe   <- vapply(grid, function(t) isTRUE(feasible_fn(t)), logical(1))
  n <- length(grid); s <- 1L; out <- list()
  while (s <= n) {
    if (!fe[s]) { s <- s + 1L; next }
    e <- s
    while (e + 1L <= n && fe[e + 1L]) e <- e + 1L
    left  <- if (s == 1L) grid[1] else .campbell_bisect_boundary(feasible_fn, grid[s - 1L], grid[s], tol)
    right <- if (e == n)  grid[n] else .campbell_bisect_boundary(feasible_fn, grid[e], grid[e + 1L], tol)
    out[[length(out) + 1L]] <- c(left, right)
    s <- e + 1L
  }
  if (!length(out)) return(matrix(numeric(0), ncol = 2))
  m <- do.call(rbind, out); colnames(m) <- c("t_lo", "t_hi"); m
}

#' Closed-form one-dimensional truncated-F selective Hotelling test (unknown covariance)
#'
#' @description
#' An **exact, Monte-Carlo-free** selective test of equality of the means of a
#' discovered cluster pair `(k1, k2)`, valid after clustering `X` with an
#' **arbitrary** black-box clustering routine `cl_fun`, treating the covariance
#' \eqn{\Sigma} as unknown. It is the strong-conditioning, closed-form companion
#' of [selective_hotelling_bc()]: it conditions on exactly the same
#' `"pair-preserved"` selection event, but instead of a Besag-Clifford Stiefel
#' path it conditions additionally on every angular nuisance coordinate of the
#' polar factor \eqn{Q} except the single Hotelling angle \eqn{\theta}, whose law
#' is then a 1-D truncated F.
#'
#' @details
#' Reduce the tested pair via the polar decomposition
#' \eqn{X_{sub} = X_{fixed} + A Q S^{1/2}} (see [prepare_pair_objects()]). Let
#' \eqn{v = A'(1/n_1,\dots,-1/n_2,\dots)} be the mean-difference contrast,
#' \eqn{c = \|v\|^2 = 1/n_1 + 1/n_2}, and \eqn{a = v/\sqrt{c}} the unit contrast
#' direction. Define the Hotelling angle by \eqn{\cos^2\theta = \|Q'a\|^2}. Then
#' the two-sample Hotelling statistic is \eqn{T^2 = (r-2)\cot^2\theta} and
#' \eqn{F = \tfrac{r-p-1}{p}\cot^2\theta \sim F_{p,\,r-p-1}} unconditionally,
#' while under the null \eqn{\cos^2\theta \sim \mathrm{Beta}(p/2, (r-1-p)/2)}
#' independently of every other angular coordinate of \eqn{Q} (\eqn{r = n_1+n_2}
#' is the pair size). Holding those coordinates fixed, \eqn{Q(\theta)} traces a
#' geodesic on which the reconstructed data \eqn{X(\theta)} is a simple harmonic
#' curve; the feasible set
#' \eqn{I_E = \{\theta \in [0,\pi/2] : cl\_fun(X(\theta))\ \mathrm{preserves}\ (k_1,k_2)\}}
#' is found by a dense grid with bisection boundary refinement. Since
#' \eqn{T^2 \ge T^2_{obs} \iff \theta \le \theta_{obs}}, the selective p-value is
#' the closed-form truncated-Beta ratio
#' \deqn{p = \frac{\mathrm{mass}\{\theta \in I_E,\ \theta \le \theta_{obs}\}}{\mathrm{mass}\{\theta \in I_E\}},}
#' with each interval's Beta mass evaluated via [stats::pbeta()] on
#' \eqn{H = \cos^2\theta}. No sampling, no burn-in, no mixing assumption.
#'
#' **What conditioning buys and costs.** Because the extra conditioning fixes the
#' contrast *direction*, the p-value is a pivot that does not involve
#' \eqn{\Sigma} at all -- so, remarkably, on this fiber knowing \eqn{\Sigma}
#' cannot improve the test. The price is paid as power: the high-dimensional
#' [selective_hotelling_bc()] retains directional freedom (and can therefore be
#' lifted by a known \eqn{\Sigma}), whereas this test trades some power for an
#' exact, closed-form, \eqn{\Sigma}-free p-value that needs no path-length tuning.
#'
#' @param X A numeric \eqn{n \times p} data matrix with observations in rows.
#' @param cl_fun A black-box clustering routine
#'   `function(Xmat) -> integer label vector` (or a list carrying `$cluster`). It
#'   **must be a pure function of its data argument** (only data-independent
#'   randomness, e.g. a fixed initialization or a deterministic algorithm such as
#'   hierarchical clustering). Called once on `X` and many times on reconstructed
#'   data \eqn{X(\theta)} along the fiber.
#' @param k1,k2 The two tested cluster labels, in the labelling of `cl_fun(X)`.
#'   Each cluster must have at least three members, and the pair size must exceed
#'   the dimension (\eqn{r - 1 > p}) for a truncated F to exist.
#' @param cl_obs Optional precomputed `cl_fun(X)` (saves one call).
#' @param n_grid Number of grid points on \eqn{[0,\pi/2]} used to locate the
#'   feasible interval(s); boundaries are then refined by bisection. Larger is
#'   more accurate at feasibility boundaries.
#' @param boundary_tol Bisection tolerance (radians) for interval endpoints.
#' @param ridge Small non-negative ridge for the reported observed \eqn{T^2} via
#'   [hotelling_T2()].
#' @param verbose If `TRUE`, print a one-line diagnostic.
#'
#' @return An S3 object of class `"selectiveHotelling"`: a list whose primary
#'   fields are `pval` (the closed-form truncated-F selective p-value), `T_obs`
#'   (observed two-sample Hotelling \eqn{T^2}), `k1`, `k2`, `K` (number of
#'   distinct labels in `cl_fun(X)`), `condition_on` (`"pair-preserved"`), and
#'   `sampler` (`"campbell-truncated-F"`). Additional fields: `p_uncond` (the
#'   classical, selection-ignoring F p-value), `theta_obs`, `F_obs`, `df1` (=
#'   \eqn{p}), `df2` (= \eqn{r-p-1}), `intervals` (feasible \eqn{\theta}-arcs),
#'   `denom_mass`, `se` (`NA`; the test is exact, not Monte Carlo), and
#'   `n_intervals`.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(30 * 2), 30, 2)
#' hclust_cl <- function(Xmat) {
#'   hc <- stats::hclust(stats::dist(Xmat), method = "average")
#'   as.integer(stats::cutree(hc, k = 3))
#' }
#' cl <- hclust_cl(X)
#' sizes <- table(cl); big <- as.integer(names(sizes)[sizes >= 3])
#' res <- selective_campbellF(X, hclust_cl, big[1], big[2])
#' res$pval
#' print(res)
#'
#' @seealso [selective_hotelling_bc()], [selective_hotelling_test()],
#'   [hotelling_T2()]
#' @export
selective_campbellF <- function(X, cl_fun, k1, k2, cl_obs = NULL,
                                n_grid = 400L, boundary_tol = 1e-6,
                                ridge = 1e-8, verbose = FALSE) {
  stopifnot(is.matrix(X), n_grid >= 8L)
  .as_labels <- function(z) if (is.list(z) && !is.null(z$cluster)) z$cluster else z
  if (is.null(cl_obs)) cl_obs <- .as_labels(cl_fun(X))
  cl_obs <- as.integer(cl_obs)
  if (length(cl_obs) != nrow(X)) stop("cl_fun(X) length != nrow(X).")
  if (sum(cl_obs == k1) < 3L || sum(cl_obs == k2) < 3L)
    stop("Tested clusters must each have >= 3 members for a Hotelling T^2.")

  prep <- prepare_pair_objects(X, cl_obs, k1, k2)
  m <- prep$m; r <- prep$r; p <- ncol(X)
  if (r - 1 <= p)
    stop(sprintf("Campbell test needs pair size r-1 > p (got r=%d, p=%d).", r, p))
  Q_obs <- prep$Q_obs

  X_of_Q <- function(Q) {
    AQS <- prep$A %*% Q %*% prep$S_sqrt
    XQ <- prep$X_fixed
    XQ[prep$idx1, ] <- XQ[prep$idx1, ] + AQS[1:m, , drop = FALSE]
    XQ[prep$idx2, ] <- XQ[prep$idx2, ] + AQS[(m + 1):r, , drop = FALSE]
    XQ
  }
  contrast <- c(rep(1 / m, m), rep(-1 / (r - m), r - m))
  v   <- as.numeric(crossprod(prep$A, contrast))
  geo <- .campbell_geodesic(Q_obs, v)
  theta_obs <- geo$theta_obs

  feasible_theta <- function(th) {
    cl_new <- tryCatch(.as_labels(cl_fun(X_of_Q(geo$Qtheta(th)))), error = function(e) NULL)
    if (is.null(cl_new) || length(cl_new) != length(cl_obs)) return(FALSE)
    .bc_pair_preserved(cl_obs, cl_new, k1, k2)
  }
  # GUARD: the observed point sits at theta_obs and must be feasible (X(theta_obs)
  # == X up to rounding, so cl_fun must reproduce cl_obs).
  if (!isTRUE(feasible_theta(theta_obs)))
    stop("GUARD FAILED: observed theta not feasible (is cl_fun a pure function of its data?).")

  intervals <- .campbell_intervals_grid(feasible_theta, theta_obs,
                                        n_grid = as.integer(n_grid), tol = boundary_tol)
  pv <- .campbell_truncated_pvalue(intervals, theta_obs, p, r)

  T_obs    <- hotelling_T2(X, cl_obs, k1, k2, ridge)
  F_obs    <- ((r - p - 1) / p) * (cos(theta_obs) / sin(theta_obs))^2
  p_uncond <- 1 - stats::pf(F_obs, p, r - p - 1)

  if (verbose)
    cat(sprintf("[campbellF] r=%d p=%d theta_obs=%.4f F_obs=%.4g df=(%d,%d) #arcs=%d p=%.4f (uncond %.4f)\n",
                r, p, theta_obs, F_obs, p, r - p - 1, nrow(intervals), pv$pval, p_uncond))

  out <- list(
    pval = pv$pval, se = NA_real_, T_obs = T_obs,
    draws = NULL, acc = NA_real_, n_unique_draws = NA_integer_,
    k1 = k1, k2 = k2, K = length(unique(cl_obs)),
    condition_on = "pair-preserved", sampler = "campbell-truncated-F",
    p_uncond = p_uncond, theta_obs = theta_obs, F_obs = F_obs,
    df1 = p, df2 = r - p - 1, intervals = intervals,
    denom_mass = pv$denom_mass, n_intervals = nrow(intervals),
    degenerate = geo$degenerate)
  class(out) <- "selectiveHotelling"
  out
}
