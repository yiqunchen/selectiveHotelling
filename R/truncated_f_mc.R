# -----------------------------------------------------------------------------
# truncated_f_mc.R -- Monte Carlo versions of the one-dimensional test
# p_polar,1d for a black-box clustering rule (no closed-form feasible set).
#
# Along the one-dimensional path Q_N(phi) of the Hotelling angle, the null law of
# the angle is f(phi) = 2/B(p/2,(r-1-p)/2) cos^{p-1}(phi) sin^{r-p-2}(phi) on
# (0, pi/2), and the selective law restricts it to the angles at which the
# clustering still returns the tested pair. The p-value is
#   P(theta <= theta_obs | theta feasible),
# because a smaller angle means a larger Hotelling T^2. Two estimators:
#
#   method = "bc": EXACT for any M. A Besag-Clifford random-index path whose kernel
#     is a slice sampler with shrinkage (Neal 2003, Sec. 4.2) for the target
#     f(phi) 1{phi feasible}. The initial bracket is the whole interval (0, pi/2),
#     which is fixed, so the kernel is reversible with respect to the target and
#     every step ends in a move.
#   method = "is": CONSISTENT (not exact at finite N). Importance sampling as in
#     Gao et al. (2024, clusterpval::test_clusters_approx): draws from a defensive
#     mixture of the null law f and a normal centred at theta_obs, weights
#     f/g * 1{feasible}.
#
# No grid approximation is used. Each feasibility check reruns the clustering.
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.truncF_mc_setup <- function(X, cl_fun, k1, k2, cl_obs) {
  .as_labels <- function(z) if (is.list(z) && !is.null(z$cluster)) z$cluster else z
  if (is.null(cl_obs)) cl_obs <- .as_labels(cl_fun(X))
  cl_obs <- as.integer(cl_obs)
  if (length(cl_obs) != nrow(X)) stop("cl_fun(X) length != nrow(X).")
  red <- .truncF_reduce(X, cl_obs, k1, k2); geo <- red$geo
  n_calls <- 0L
  feasible <- function(th) {
    n_calls <<- n_calls + 1L
    cl_new <- tryCatch(.as_labels(cl_fun(red$X_of_Q(geo$Qtheta(th)))), error = function(e) NULL)
    if (is.null(cl_new) || length(cl_new) != length(cl_obs)) return(FALSE)
    .bc_pair_preserved(cl_obs, cl_new, k1, k2)
  }
  if (!isTRUE(feasible(geo$theta_obs)))
    stop("GUARD FAILED: observed angle not feasible (is cl_fun a pure function of its data?).")
  p <- red$p; r <- red$r
  a <- p / 2; b <- (r - 1 - p) / 2
  logf <- function(ph) log(2) - lbeta(a, b) + (p - 1) * log(cos(ph)) + (r - p - 2) * log(sin(ph))
  list(red = red, geo = geo, cl_obs = cl_obs, feasible = feasible, logf = logf,
       a = a, b = b, calls = function() n_calls)
}

# One slice-sampling step with shrinkage on (0, pi/2), target f * 1{feasible}.
#' @keywords internal
#' @noRd
.truncF_slice_step <- function(ph, st) {
  ly <- st$logf(ph) + log(stats::runif(1))
  lo <- 0; hi <- pi / 2
  repeat {
    prop <- stats::runif(1, lo, hi)
    if (st$logf(prop) > ly && isTRUE(st$feasible(prop))) return(prop)
    if (prop < ph) lo <- prop else hi <- prop
    if (hi - lo < 1e-12) return(ph)
  }
}

# One random-step (Givens-style) Metropolis step on (0, pi/2), target f * 1{feasible}:
# propose phi' = phi + N(0, (pi * eps)^2) (symmetric), stay if phi' leaves (0, pi/2),
# otherwise accept with probability min(1, f(phi')/f(phi)) when the clustering is
# reproduced (feasibility is checked only after the density test passes).
#' @keywords internal
#' @noRd
.truncF_rw_step <- function(ph, st, eps) {
  prop <- ph + stats::rnorm(1, 0, pi * eps)
  if (prop <= 0 || prop >= pi / 2) return(ph)
  if (log(stats::runif(1)) > st$logf(prop) - st$logf(ph)) return(ph)
  if (isTRUE(st$feasible(prop))) prop else ph
}

#' @keywords internal
#' @noRd
.truncF_bc <- function(st, M, L, kernel = "slice", eps_range = c(1e-3, 1)) {
  th0 <- st$geo$theta_obs
  J <- sample.int(M + 1L, 1L) - 1L
  path <- numeric(M + 1L); path[J + 1L] <- th0
  ph <- th0
  leps <- log(eps_range)
  move <- function(ph) {
    if (kernel == "slice") { for (s in seq_len(L)) ph <- .truncF_slice_step(ph, st); return(ph) }
    e <- exp(stats::runif(1, leps[1], leps[2]))      # step size per block of L steps
    for (s in seq_len(L)) ph <- .truncF_rw_step(ph, st, e)
    ph
  }
  if (J < M) for (j in (J + 1L):M) { ph <- move(ph); path[j + 1L] <- ph }
  ph <- th0
  if (J > 0L) for (j in (J - 1L):0L) { ph <- move(ph); path[j + 1L] <- ph }
  xi <- stats::runif(1)
  pval <- (sum(path < th0) + xi * sum(path == th0)) / (M + 1L)
  list(pval = pval, path = path, J = J)
}

#' @keywords internal
#' @noRd
.truncF_is <- function(st, N, mix, scale) {
  th0 <- st$geo$theta_obs
  s <- if (is.null(scale)) max(0.02, 0.25 * min(th0, pi / 2 - th0)) else scale
  # normal centred at theta_obs, truncated to (0, pi/2), drawn by inverse CDF
  plo <- stats::pnorm(0, th0, s); phi_ <- stats::pnorm(pi / 2, th0, s)
  rtn <- function(k) stats::qnorm(stats::runif(k, plo, phi_), th0, s)
  dtn <- function(x) stats::dnorm(x, th0, s) / (phi_ - plo)
  rnull <- function(k) acos(sqrt(stats::rbeta(k, st$a, st$b)))
  use_null <- stats::runif(N) < mix
  ph <- numeric(N); ph[use_null] <- rnull(sum(use_null)); ph[!use_null] <- rtn(sum(!use_null))
  g <- mix * exp(st$logf(ph)) + (1 - mix) * dtn(ph)
  fe <- vapply(ph, function(x) isTRUE(st$feasible(x)), logical(1))
  w <- ifelse(fe, exp(st$logf(ph)) / g, 0)
  if (!any(w > 0)) return(list(pval = NA_real_, se = NA_real_, ess = 0, n_feasible = 0L))
  num <- sum(w[ph <= th0]); den <- sum(w)
  pval <- num / den
  # delta-method standard error of the ratio estimator
  z <- w * ((ph <= th0) - pval)
  se <- sqrt(sum(z^2)) / den
  list(pval = pval, se = se, ess = den^2 / sum(w^2), n_feasible = sum(fe), scale = s)
}

#' @keywords internal
#' @noRd
.truncF_mc <- function(X, cl_fun, k1, k2, cl_obs, method, M, L, n_is, is_mix,
                       is_scale, ridge, verbose, onedim_kernel = "slice") {
  st <- .truncF_mc_setup(X, cl_fun, k1, k2, cl_obs)
  t0 <- proc.time()[["elapsed"]]
  res <- if (method == "bc") .truncF_bc(st, M, L, kernel = onedim_kernel) else .truncF_is(st, n_is, is_mix, is_scale)
  secs <- proc.time()[["elapsed"]] - t0
  out <- .truncF_result(st$red, matrix(c(0, pi / 2), 1L,
                                       dimnames = list(NULL, c("lo", "hi"))),
                        X, st$cl_obs, k1, k2, "pair-preserved", ridge)
  out$pval <- res$pval
  out$se <- if (!is.null(res$se)) res$se else sqrt(res$pval * (1 - res$pval) / (M + 1L))
  out$intervals <- NULL; out$denom_mass <- NA_real_; out$n_intervals <- NA_integer_
  out$sampler <- if (method == "bc") sprintf("truncated-F, Besag-Clifford %s path (exact)", onedim_kernel)
                 else "truncated-F, importance sampling (consistent)"
  out$mc_method <- method
  out$n_cluster_calls <- st$calls()
  out$seconds <- secs
  if (method == "bc") { out$M <- M; out$L <- L; out$draws <- res$path; out$onedim_kernel <- onedim_kernel }
  else { out$n_is <- n_is; out$ess <- res$ess; out$n_feasible <- res$n_feasible; out$is_scale <- res$scale }
  if (verbose) cat(sprintf("[truncated-F %s] r=%d p=%d theta_obs=%.4f p-value=%.4f calls=%d (%.1fs)\n",
                           method, st$red$r, st$red$p, st$geo$theta_obs, out$pval, out$n_cluster_calls, secs))
  out
}
