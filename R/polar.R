# ---------- Core reduction & polar objects (no qr.Q; exact polar) ----------

#' Build polar-decomposition objects for a selected cluster pair
#'
#' Constructs the fixed reduction used by all samplers. For the rows of `X` in
#' clusters `k1` and `k2`, it forms the projector onto the orthogonal complement
#' of the all-ones vector, an orthonormal basis `A` of that complement, the
#' Gram matrix `S = X_sub' P_perp X_sub` with its symmetric square root, the
#' observed polar factor `Q_obs` on the Stiefel manifold, and the part of `X`
#' that is held fixed while only `Q` varies.
#'
#' @param X A numeric data matrix with observations in rows.
#' @param cl Cluster-assignment vector for the rows of `X`.
#' @param k1,k2 The two cluster labels defining the tested pair.
#'
#' @return A list with elements `idx1`, `idx2`, `m`, `n`, `r`, `A`, `P_perp`,
#'   `S_sqrt`, `X_fixed`, and `Q_obs`.
#'
#' @keywords internal
prepare_pair_objects <- function(X, cl, k1, k2) {
  idx1 <- which(cl == k1); idx2 <- which(cl == k2)
  X_sub <- rbind(X[idx1, , drop = FALSE], X[idx2, , drop = FALSE])
  m <- length(idx1); n <- length(idx2); r <- m + n

  # Projector onto v1^⊥
  v1 <- matrix(1, nrow = r, ncol = 1)
  P_perp <- diag(r) - v1 %*% t(v1) / as.numeric(t(v1) %*% v1)

  # A with AA' = P_perp and A'A = I_{r-1} (via eigendecomp)
  es <- eigen(P_perp, symmetric = TRUE)
  A <- es$vectors[, 1:(r-1), drop = FALSE] %*% diag(sqrt(es$values[1:(r-1)]), r-1)

  # S = X_sub^T P_perp X_sub; polar: M = A^T X_sub = Q S^{1/2}
  S_mat <- t(X_sub) %*% P_perp %*% X_sub
  Sfac  <- chol_psd_sqrt(S_mat)

  M     <- t(A) %*% X_sub
  Q_obs <- M %*% Sfac$invsqrt  # exact polar factor; already has orthonormal columns

  # Fixed part: remove P_perp X_sub from those rows
  AQS_obs <- P_perp %*% X_sub
  X_fixed <- X
  X_fixed[idx1, ] <- X_fixed[idx1, ] - AQS_obs[1:m, , drop = FALSE]
  X_fixed[idx2, ] <- X_fixed[idx2, ] - AQS_obs[(m+1):r, , drop = FALSE]

  list(idx1=idx1, idx2=idx2, m=m, n=n, r=r,
       A=A, P_perp=P_perp, S_sqrt=Sfac$sqrt,
       X_fixed=X_fixed, Q_obs=Q_obs)
}

# ---------- Cayley left-multiplicative step on Stiefel ---------------------

#' One left-multiplicative Cayley random-walk step on the Stiefel manifold
#'
#' Proposes a new Stiefel point by left-multiplying `Q` by an orthogonal matrix
#' generated from a random skew-symmetric matrix via the Cayley transform, with
#' step size `eps`.
#'
#' @param Q An `r` by `p` matrix with orthonormal columns.
#' @param eps Step size controlling proposal scale.
#'
#' @return A new `r` by `p` matrix with orthonormal columns.
#'
#' @keywords internal
cayley_left_step <- function(Q, eps) {
  r <- nrow(Q)
  G <- matrix(stats::rnorm(r * r), r, r)
  K <- G - t(G)  # skew-symmetric
  I <- diag(r)
  U <- solve(I - 0.5 * eps * K, I + 0.5 * eps * K)  # orthogonal via Cayley
  U %*% Q
}

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
