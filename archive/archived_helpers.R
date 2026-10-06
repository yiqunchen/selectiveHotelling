# Archived 2026-10-05: helpers used only by the add-one MCMC samplers (samplers.R).
# ---- from polar.R ----
#' Calibrate the Cayley step size for a minimum acceptance rate
#'
#' Repeatedly probes the feasibility of Cayley proposals from `Q` and shrinks
#' `eps` until the empirical pilot acceptance rate reaches `lower` (or `eps`
#' hits `min_eps`).
#'
#' @param Q Current Stiefel point.
#' @param eps Initial step size.
#' @param feasible_Q A predicate returning `TRUE` if a proposed `Q` is feasible.
#' @param n_probe Number of proposals per calibration round.
#' @param lower Target minimum acceptance rate.
#' @param shrink Multiplicative shrink factor for `eps`.
#' @param min_eps Lower bound on `eps`.
#' @param max_rounds Maximum number of calibration rounds.
#'
#' @return A list with elements `eps` (calibrated step) and `pilot_acceptance`.
#'
#' @keywords internal
calibrate_cayley_eps <- function(Q, eps, feasible_Q,
                                 n_probe = 80, lower = 0.05,
                                 shrink = 0.25, min_eps = 1e-8,
                                 max_rounds = 12) {
  eps <- max(eps, min_eps)
  last_rate <- NA_real_

  for (round in seq_len(max_rounds)) {
    accepted <- 0L
    for (probe in seq_len(n_probe)) {
      accepted <- accepted + isTRUE(feasible_Q(cayley_left_step(Q, eps)))
    }
    last_rate <- accepted / n_probe
    if (last_rate >= lower || eps <= min_eps) break
    eps <- max(eps * shrink, min_eps)
  }

  list(eps = eps, pilot_acceptance = last_rate)
}

#' Adapt the Cayley step size toward a target acceptance rate
#'
#' Multiplicatively rescales `eps` based on a windowed acceptance rate so that
#' the chain moves toward the `target` acceptance rate, clamped to
#' `[min_eps, max_eps]`.
#'
#' @param eps Current step size.
#' @param window_rate Observed acceptance rate over the recent window.
#' @param target Target acceptance rate.
#' @param min_eps,max_eps Bounds on the returned step size.
#'
#' @return The adapted step size.
#'
#' @keywords internal
adapt_cayley_eps <- function(eps, window_rate, target = 0.4,
                             min_eps = 1e-8, max_eps = 0.5) {
  if (!is.finite(window_rate)) return(eps)

  factor <- if (window_rate == 0) {
    0.20
  } else if (window_rate < 0.05) {
    0.35
  } else if (window_rate < target / 2) {
    0.70
  } else if (window_rate > 0.80) {
    1.25
  } else {
    exp(0.5 * (window_rate - target))
  }

  max(min(eps * factor, max_eps), min_eps)
}

# ---- from utils.R ----
#' Hotelling T-squared on the selection-mapped cluster pair
#'
#' Resolves the (possibly relabeled) clusters that correspond to the originally
#' selected pair `(k1, k2)` via [selected_pair_labels()] and then evaluates
#' [hotelling_T2()] on those labels.
#'
#' @param X A numeric data matrix with observations in rows.
#' @param cl0 The original cluster assignment vector.
#' @param cl1 The new cluster assignment vector (possibly relabeled).
#' @param k1,k2 The originally selected cluster labels (in `cl0`).
#' @param ridge Ridge passed to [hotelling_T2()].
#'
#' @return The \eqn{T^2} value for the mapped pair, or `NA_real_` if the pair is
#'   not preserved.
#'
#' @keywords internal
hotelling_T2_selected <- function(X, cl0, cl1, k1, k2, ridge = 1e-8) {
  labels <- selected_pair_labels(cl0, cl1, k1, k2)
  if (is.null(labels)) return(NA_real_)
  hotelling_T2(X, cl1, labels[["k1"]], labels[["k2"]], ridge)
}

#' Monte Carlo upper-tail p-value with batched standard error
#'
#' Computes the Monte Carlo upper-tail p-value \eqn{(1 + \#\{draws \ge
#' T_{obs}\})/(n + 1)} from a vector of sampled statistics, together with a
#' batched-means standard error estimate.
#'
#' @param draws Numeric vector of sampled \eqn{T^2} values (may contain `NA`).
#' @param T_obs The observed \eqn{T^2} value.
#'
#' @return A list with elements `pval`, `se`, and `n` (the number of finite
#'   draws used).
#'
#' @keywords internal
mc_tail_pvalue <- function(draws, T_obs) {
  ok <- is.finite(draws)
  n_ok <- sum(ok)
  if (n_ok == 0L || !is.finite(T_obs)) return(list(pval = NA_real_, se = NA_real_, n = n_ok))

  exceed <- draws[ok] >= T_obs
  pval <- (1 + sum(exceed)) / (n_ok + 1)

  nb <- max(10L, floor(n_ok / 50))
  bsz <- floor(n_ok / nb)
  se <- NA_real_
  if (bsz >= 1) {
    bm <- sapply(seq_len(nb), function(j) {
      st <- (j - 1) * bsz + 1
      en <- j * bsz
      (1 + sum(exceed[st:en])) / (bsz + 1)
    })
    se <- stats::sd(bm) / sqrt(nb)
  }

  list(pval = pval, se = se, n = n_ok)
}

# ---- from kmeans.R ----
#' Check preservation of the full k-means trajectory for a cluster pair
#'
#' Returns `TRUE` only if, at *every* k-means iteration, the observations that
#' were in clusters `k1` and `k2` of the original trajectory `traj0` form exactly
#' two distinct clusters in the candidate trajectory `traj1`. This is the
#' (experimental) trajectory-conditioning feasibility check.
#'
#' @param traj0,traj1 Lists of per-iteration cluster-assignment vectors.
#' @param k1,k2 The two original cluster labels being tested.
#'
#' @return Logical scalar.
#'
#' @keywords internal
preserve_trajectory <- function(traj0, traj1, k1, k2) {
  # Check that ALL intermediate cluster assignments are preserved for clusters k1, k2.
  # traj0, traj1: lists of cluster assignment vectors (one per k-means iteration)
  # Returns TRUE only if every iteration's assignments match for the tested clusters.

  # Must have same number of iterations
  if (length(traj0) != length(traj1)) return(FALSE)

  # Check each iteration
  for (t in seq_along(traj0)) {
    cl0_t <- traj0[[t]]
    cl1_t <- traj1[[t]]

    # Find which observations were in k1, k2 at iteration t in original
    i1 <- which(cl0_t == k1)
    i2 <- which(cl0_t == k2)

    # In the new trajectory, these same observations must form exactly two clusters
    L1 <- unique(cl1_t[i1])
    L2 <- unique(cl1_t[i2])

    # Each set must map to exactly one cluster label
    if (length(L1) != 1 || length(L2) != 1) return(FALSE)

    # The two clusters must be different
    if (L1 == L2) return(FALSE)

    # All obs in new cluster L1 must be exactly i1, and similarly for L2
    if (!all(cl1_t[i1] == L1)) return(FALSE)
    if (!all(cl1_t[i2] == L2)) return(FALSE)
    if (!all(which(cl1_t == L1) %in% i1)) return(FALSE)
    if (!all(which(cl1_t == L2) %in% i2)) return(FALSE)
  }
  TRUE
}

#' Map originally selected clusters to their current labels
#'
#' Feasibility checks permit label permutations, so downstream \eqn{T^2}
#' calculations must use the *mapped* labels rather than reusing `k1`/`k2`. This
#' returns the current labels of the originally selected clusters, or `NULL` if
#' the pair is not preserved.
#'
#' @param cl0,cl1 Original and candidate cluster-assignment vectors.
#' @param k1,k2 The two original cluster labels being tested.
#'
#' @return A named length-2 vector `c(k1 = L1, k2 = L2)`, or `NULL`.
#'
#' @keywords internal
selected_pair_labels <- function(cl0, cl1, k1, k2) {
  # Return the current labels corresponding to the originally selected clusters.
  # Feasibility checks allow label permutations, so downstream T^2 calculations
  # must use the mapped labels rather than blindly reusing k1/k2.
  i1 <- which(cl0 == k1); i2 <- which(cl0 == k2)
  if (length(i1) == 0 || length(i2) == 0) return(NULL)

  L1 <- unique(cl1[i1]); L2 <- unique(cl1[i2])
  if (length(L1) != 1 || length(L2) != 1 || L1 == L2) return(NULL)
  if (!all(cl1[i1] == L1) || !all(cl1[i2] == L2)) return(NULL)
  if (!all(which(cl1 == L1) %in% i1)) return(NULL)
  if (!all(which(cl1 == L2) %in% i2)) return(NULL)

  c(k1 = L1, k2 = L2)
}

# ---- from besag_clifford.R (calibrate = TRUE tuned eps from Q_obs, which can break exactness) ----
# -----------------------------------------------------------------------------
# .bc_calibrate_eps(Q_obs, feasible, ...): adapt the Cayley step eps toward a
# target acceptance rate, by short pilot bursts from Q_obs. Only affects POWER,
# never validity. Returns the chosen eps.
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.bc_calibrate_eps <- function(Q_obs, feasible, eps0 = 0.3,
                              target_lo = 0.20, target_hi = 0.60,
                              n_probe = 40L, rounds = 8L,
                              min_eps = 1e-4, max_eps = 2.0,
                              kernel = "cayley") {
  eps <- eps0
  for (rd in seq_len(rounds)) {
    acc <- mean(vapply(seq_len(n_probe), function(.)
      isTRUE(feasible(.bc_propose(Q_obs, eps, kernel))), logical(1)))
    if (acc < target_lo)      eps <- max(eps * 0.5, min_eps)
    else if (acc > target_hi) eps <- min(eps * 1.5, max_eps)
    else break
    if (eps <= min_eps || eps >= max_eps) break
  }
  eps
}

