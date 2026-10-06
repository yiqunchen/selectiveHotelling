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

# ---- givens_left (moved from samplers.R on 2026-10-05; used by the Besag-Clifford kernels) ----
#' Givens rotation acting on two rows of Q
#'
#' Applies a Givens rotation of angle `theta` to rows `i` and `j` of `Q` (a left
#' action on those two rows), used by the geodesic hit-and-run / slice samplers.
#'
#' @param Q An `r` by `p` matrix with orthonormal columns.
#' @param i,j Row indices to rotate.
#' @param theta Rotation angle.
#'
#' @return The rotated matrix.
#'
#' @keywords internal
givens_left <- function(Q, i, j, theta) {
  Qi <- Q[i, , drop = FALSE]
  Qj <- Q[j, , drop = FALSE]
  cth <- cos(theta); sth <- sin(theta)
  Q[i, ] <-  cth * Qi + (-sth) * Qj
  Q[j, ] <-  sth * Qi +  cth  * Qj
  Q
}


# ---- .har_givens_PUW (moved from har.R on 2026-10-05; used by sim/path_bc_hclust.R) ----
# -----------------------------------------------------------------------------
# .har_givens_PUW(X_of_Q, X_fixed, Q, i, j): decompose the reconstructed data
# X(theta) = X_of_Q(givens_left(Q, i, j, theta)) as the harmonic curve
#     X(theta) = P + cos theta * U + sin theta * W .
# givens_left(Q,i,j,theta): row i = cos*Q[i] - sin*Q[j]; row j = sin*Q[i] + cos*Q[j].
# Hence Q(theta) = Qc + cos theta * Mcos + sin theta * Msin with
#   Mcos: rows i,j = Q[i,], Q[j,]  ;  Msin: row i = -Q[j,], row j = Q[i,]  ;
#   Qc  : Q with rows i,j zeroed.
# Feeds the package-native (global-free) .geodesic_PUW.
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.har_givens_PUW <- function(X_of_Q, X_fixed, Q, i, j) {
  rq <- nrow(Q); p <- ncol(Q)
  Mcos <- matrix(0, rq, p); Mcos[i, ] <- Q[i, ]; Mcos[j, ] <- Q[j, ]
  Msin <- matrix(0, rq, p); Msin[i, ] <- -Q[j, ]; Msin[j, ] <- Q[i, ]
  Qc <- Q; Qc[i, ] <- 0; Qc[j, ] <- 0
  .geodesic_PUW(X_of_Q, X_fixed, Qc, Mcos, Msin)
}

