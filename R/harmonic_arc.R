# =============================================================================
# harmonic_arc.R -- EXACT closed-form feasible arcs on a one-parameter geodesic.
#
# Along any geodesic Q(t) = Q_c + cos t * M_cos + sin t * M_sin of the polar
# factor, the affine reconstruction X(t) = X_of_Q(Q(t)) is a harmonic curve
#       X(t) = P + cos t * U + sin t * W        (each n x p),
# so every squared distance -- point-to-centroid (k-means) or point-to-point
# (hierarchical) -- is a trigonometric polynomial
#       g(t) = A cos 2t + B sin 2t + C cos t + D sin t + E,
# whose roots on [0, 2pi) are found in closed form (Weierstrass s = tan(t/2)).
# A clustering-preservation event is an intersection of {g_j(t) <= 0}, i.e. a
# union of feasible arcs whose endpoints are constraint roots. This module is the
# self-contained, GLOBAL-FREE package engine for those arcs (the sim twin lives
# in sim/analytic_arc.R and is verified against a brute-force grid). It powers
#   * selective_truncated_F(..., method = "exact")  -- the 1-D truncated-F test, and
#   * the exact-arc geodesic hit-and-run samplers,
# for BOTH k-means (Lloyd trajectory) and hierarchical (agglomeration trajectory)
# conditioning.
# =============================================================================

# ---- roots of a single harmonic on [0, 2pi) (verified twin: analytic_arc.R) ---
#' @keywords internal
#' @noRd
# ---- batched harmonic roots ---------------------------------------------------
# Roots of A cos2t + B sin2t + C cos t + D sin t + E = 0 on [0, 2pi) for MANY
# constraints at once. Same algebra as .harmonic_roots (Weierstrass s = tan(t/2)
# gives a quartic, then Newton polish), but the polish runs vectorised over all
# candidates from all constraints instead of scalar loops with per-call closures.
# Returns list(t, j) with j the constraint index of each accepted root.
#' @keywords internal
#' @noRd
.harmonic_roots_batch <- function(A, B, C, D, E, tol = 1e-9) {
  nC <- length(A)
  if (!nC) return(list(t = numeric(0), j = integer(0)))
  amp <- abs(A) + abs(B) + abs(C) + abs(D) + abs(E)

  p0 <- A + C + E; p1 <- 4 * B + 2 * D; p2 <- -6 * A + 2 * E
  p3 <- -4 * B + 2 * D; p4 <- A - C + E
  cand_t <- numeric(0); cand_j <- integer(0)
  for (j in seq_len(nC)) {
    if (amp[j] == 0) next
    coef <- c(p0[j], p1[j], p2[j], p3[j], p4[j])
    nz <- suppressWarnings(max(which(abs(coef) > tol * (sum(abs(coef)) + 1))))
    if (!is.finite(nz)) nz <- 1L
    if (nz < 2L) next
    rr <- tryCatch(polyroot(coef[seq_len(nz)]), error = function(e) complex(0))
    if (!length(rr)) next
    s <- Re(rr); s <- s[is.finite(s)]
    if (!length(s)) next
    cand_t <- c(cand_t, (2 * atan(s)) %% (2 * pi))
    cand_j <- c(cand_j, rep.int(j, length(s)))
  }

  if (length(cand_t)) {                       # vectorised Newton over ALL candidates
    Aj <- A[cand_j]; Bj <- B[cand_j]; Cj <- C[cand_j]; Dj <- D[cand_j]; Ej <- E[cand_j]
    t <- cand_t
    for (it in seq_len(8L)) {
      c2 <- cos(2 * t); s2 <- sin(2 * t); ct <- cos(t); st <- sin(t)
      gv  <- Aj * c2 + Bj * s2 + Cj * ct + Dj * st + Ej
      gpv <- -2 * Aj * s2 + 2 * Bj * c2 - Cj * st + Dj * ct
      ok <- abs(gpv) >= 1e-14
      if (!any(ok)) break
      step <- ifelse(ok, gv / gpv, 0)
      t <- t - step
      if (max(abs(step)) < 1e-13) break
    }
    t <- t %% (2 * pi)
    gv <- A[cand_j] * cos(2 * t) + B[cand_j] * sin(2 * t) +
          C[cand_j] * cos(t) + D[cand_j] * sin(t) + E[cand_j]
    keep <- abs(gv) <= 1e-7 * amp[cand_j] + 1e-12
    cand_t <- t[keep]; cand_j <- cand_j[keep]
  }

  # t = pi is a root whenever g(pi) = A - C + E is (numerically) zero
  gpi <- A - C + E
  jp  <- which(amp > 0 & abs(gpi) <= 1e-7 * amp + 1e-12)
  if (length(jp)) { cand_t <- c(cand_t, rep.int(pi, length(jp))); cand_j <- c(cand_j, jp) }
  list(t = cand_t, j = cand_j)
}

.harmonic_roots <- function(A, B, C, D, E, tol = 1e-9) {
  amp <- abs(A) + abs(B) + abs(C) + abs(D) + abs(E)
  if (amp == 0) return(numeric(0))
  g  <- function(t) A * cos(2 * t) + B * sin(2 * t) + C * cos(t) + D * sin(t) + E
  gp <- function(t) -2 * A * sin(2 * t) + 2 * B * cos(2 * t) - C * sin(t) + D * cos(t)
  # (1+s^2)^2 * g(t) is a real quartic in s = tan(t/2):
  p0 <- A + C + E; p1 <- 4 * B + 2 * D; p2 <- -6 * A + 2 * E
  p3 <- -4 * B + 2 * D; p4 <- A - C + E
  coef <- c(p0, p1, p2, p3, p4)
  nz <- max(which(abs(coef) > tol * (sum(abs(coef)) + 1)))
  if (!is.finite(nz)) nz <- 1L
  cand_s <- numeric(0)
  if (nz >= 2L) {
    rr <- tryCatch(polyroot(coef[seq_len(nz)]), error = function(e) complex(0))
    if (length(rr)) cand_s <- Re(rr)
  }
  roots_t <- numeric(0)
  for (s in cand_s) {
    if (!is.finite(s)) next
    t <- (2 * atan(s)) %% (2 * pi)
    for (it in seq_len(8L)) {
      gv <- g(t); gpv <- gp(t)
      if (abs(gpv) < 1e-14) break
      step <- gv / gpv; t <- t - step
      if (abs(step) < 1e-13) break
    }
    t <- t %% (2 * pi)
    if (abs(g(t)) <= 1e-7 * amp + 1e-12) roots_t <- c(roots_t, t)
  }
  if (abs(g(pi)) <= 1e-7 * amp + 1e-12) roots_t <- c(roots_t, pi)
  roots_t
}

# ---- feasible arcs on [lo, hi] of { max_j g_j(t) <= 0 } ------------------------
# coeffs: list with numeric vectors A,B,C,D,E (one entry per constraint).
# Returns a k x 2 matrix of [t_lo, t_hi] arcs where every constraint holds.
#' @keywords internal
#' @noRd
.harmonic_feasible_arcs <- function(coeffs, lo = 0, hi = 2 * pi,
                                    prune = TRUE, tol_mid = 1e-10) {
  A <- coeffs$A; B <- coeffs$B; C <- coeffs$C; D <- coeffs$D; E <- coeffs$E
  nC <- length(A)
  if (prune && nC > 0L) {
    # |g - E| <= M classifies most constraints without any root work:
    #   E + M < 0  -> g < 0 everywhere, constraint is inactive, DROP it
    #   E - M > 0  -> g > 0 everywhere, constraint can never hold, set is EMPTY
    M <- sqrt(A^2 + B^2) + sqrt(C^2 + D^2)
    if (any((E - M) > 1e-9)) {
      out <- matrix(numeric(0), ncol = 2); colnames(out) <- c("t_lo", "t_hi")
      return(out)
    }
    keep <- (E + M) >= -1e-9
    A <- A[keep]; B <- B[keep]; C <- C[keep]; D <- D[keep]; E <- E[keep]
    nC <- length(A)
  }
  g_all <- function(t) A * cos(2 * t) + B * sin(2 * t) + C * cos(t) + D * sin(t) + E
  maxg  <- function(t) if (nC == 0L) -Inf else max(g_all(t))

  bks <- c(lo, hi)
  if (nC > 0L) for (j in seq_len(nC)) {
    r <- .harmonic_roots(A[j], B[j], C[j], D[j], E[j])
    if (length(r)) bks <- c(bks, r[r >= lo - 1e-12 & r <= hi + 1e-12])
  }
  bks <- sort(unique(c(pmin(pmax(bks, lo), hi), lo, hi)))
  nseg <- length(bks) - 1L
  keep_lo <- numeric(0); keep_hi <- numeric(0)
  for (s in seq_len(nseg)) {
    a <- bks[s]; b <- bks[s + 1L]
    if (b - a <= 0) next
    if (maxg(0.5 * (a + b)) <= tol_mid) { keep_lo <- c(keep_lo, a); keep_hi <- c(keep_hi, b) }
  }
  arcs <- matrix(numeric(0), ncol = 2)
  if (length(keep_lo)) {
    cur_lo <- keep_lo[1]; cur_hi <- keep_hi[1]
    if (length(keep_lo) >= 2L) for (s in 2:length(keep_lo)) {
      if (abs(keep_lo[s] - cur_hi) <= 1e-12) cur_hi <- keep_hi[s]
      else { arcs <- rbind(arcs, c(cur_lo, cur_hi)); cur_lo <- keep_lo[s]; cur_hi <- keep_hi[s] }
    }
    arcs <- rbind(arcs, c(cur_lo, cur_hi))
  }
  colnames(arcs) <- c("t_lo", "t_hi")
  arcs
}

# ---- harmonic of || dP + cos t dU + sin t dW ||^2 (rows = many differences) ----
# Returns a list of A,B,C,D,E vectors (one per row of the difference matrices).
#' @keywords internal
#' @noRd
.sqdist_harmonic <- function(dP, dU, dW) {
  pp <- rowSums(dP * dP); uu <- rowSums(dU * dU); ww <- rowSums(dW * dW)
  pu <- rowSums(dP * dU); pw <- rowSums(dP * dW); uw <- rowSums(dU * dW)
  list(A = (uu - ww) / 2, B = uw, C = 2 * pu, D = 2 * pw, E = pp + (uu + ww) / 2)
}

# ---- geodesic data decomposition X(t) = P + cos t*U + sin t*W ------------------
# Q(t) = Qc + cos t*Mcos + sin t*Msin. X_of_Q affine => P/U/W below.
#' @keywords internal
#' @noRd
.geodesic_PUW <- function(X_of_Q, X_fixed, Qc, Mcos, Msin) {
  list(P = X_of_Q(Qc), U = X_of_Q(Mcos) - X_fixed, W = X_of_Q(Msin) - X_fixed)
}

# ---- center history along a fixed Lloyd trajectory (for k-means constraints) --
# Centers USED at each iteration t: t=1 -> init rows; t>=2 -> means of traj[[t-1]]
# clusters (empty clusters carry their previous center forward).
#' @keywords internal
#' @noRd
.center_history <- function(Xmat, init_idx, traj, K) {
  L <- length(traj)
  centers <- Xmat[init_idx, , drop = FALSE]
  CH <- vector("list", L)
  for (t in seq_len(L)) {
    CH[[t]] <- centers
    if (t < L) {
      lab <- traj[[t]]
      for (j in seq_len(K)) {
        mem <- which(lab == j)
        if (length(mem)) centers[j, ] <- colMeans(Xmat[mem, , drop = FALSE])
      }
    }
  }
  CH
}

# ---- k-means Lloyd-trajectory constraint harmonics ----------------------------
# For each (iteration t, point k, competitor label b != observed z):
#   || x_k(t) - c_z(t) ||^2 - || x_k(t) - c_b(t) ||^2 <= 0.
#' @keywords internal
#' @noRd
.kmeans_traj_coeffs <- function(P, U, W, init_idx, traj, K) {
  CHp <- .center_history(P, init_idx, traj, K)
  CHu <- .center_history(U, init_idx, traj, K)
  CHw <- .center_history(W, init_idx, traj, K)
  n <- nrow(P); L <- length(traj)
  A <- B <- C <- D <- E <- numeric(0)
  for (t in seq_len(L)) {
    Pc <- CHp[[t]]; Uc <- CHu[[t]]; Wc <- CHw[[t]]; lab <- traj[[t]]
    for (b in seq_len(K)) {
      act <- which(lab != b)                 # points whose observed label != b
      if (!length(act)) next
      z <- lab[act]
      # d_z = (P[k]-Pc[z]) + cos*(U[k]-Uc[z]) + sin*(W[k]-Wc[z]); d_b analogous
      ez0 <- P[act, , drop = FALSE] - Pc[z, , drop = FALSE]
      ez1 <- U[act, , drop = FALSE] - Uc[z, , drop = FALSE]
      ez2 <- W[act, , drop = FALSE] - Wc[z, , drop = FALSE]
      eb0 <- P[act, , drop = FALSE] - Pc[rep(b, length(act)), , drop = FALSE]
      eb1 <- U[act, , drop = FALSE] - Uc[rep(b, length(act)), , drop = FALSE]
      eb2 <- W[act, , drop = FALSE] - Wc[rep(b, length(act)), , drop = FALSE]
      hz <- .sqdist_harmonic(ez0, ez1, ez2)
      hb <- .sqdist_harmonic(eb0, eb1, eb2)
      A <- c(A, hz$A - hb$A); B <- c(B, hz$B - hb$B); C <- c(C, hz$C - hb$C)
      D <- c(D, hz$D - hb$D); E <- c(E, hz$E - hb$E)
    }
  }
  list(A = A, B = B, C = C, D = D, E = E)
}

# ---- point-pair squared-distance harmonics (n x n, CELL-INDEPENDENT) ----------
# These depend only on the geodesic (P, U, W), not on the cell being probed, so
# they are built ONCE per arc traversal rather than once per cell.
#' @keywords internal
#' @noRd
.hclust_pointpair_harmonics <- function(P, U, W) {
  n <- nrow(P)
  A <- B <- C <- D <- E <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) {
    js <- (i + 1L):n
    ii <- rep(i, length(js))
    h  <- .sqdist_harmonic(P[js, , drop = FALSE] - P[ii, , drop = FALSE],
                           U[js, , drop = FALSE] - U[ii, , drop = FALSE],
                           W[js, , drop = FALSE] - W[ii, , drop = FALSE])
    A[i, js] <- h$A; A[js, i] <- h$A;  B[i, js] <- h$B; B[js, i] <- h$B
    C[i, js] <- h$C; C[js, i] <- h$C;  D[i, js] <- h$D; D[js, i] <- h$D
    E[i, js] <- h$E; E[js, i] <- h$E
  }
  list(A = A, B = B, C = C, D = D, E = E)
}

# ---- ONE-SHOT arcs for a FIXED average-linkage merge trace --------------------
# Conditions on the OBSERVED merge sequence up to the K-cluster cut, rather than on
# the weaker "cutree-K preserves the tested pair". That is strictly more
# conditioning -- R_T is one component of the union S = U_T R_T over traces whose
# cut preserves the pair -- and it must be reported as such. What it buys is that
# the event becomes a FIXED conjunction, so there is no cell march and no cell
# budget: it cannot hit max_cells, hence cannot return the NA that a truncated
# march now yields (see .truncF_result).
#
# WHAT IT COSTS: EVERYTHING. Measured 2026-08-12 (job 2044756, 72 power cells x 300
# reps, n=200): this event has NO POWER. Its p-values are uniform under the
# alternative, so power equals the nominal level in every cell and at every delta up
# to 12 (0.03-0.13, against 1.000 for the BC test). Clinching cell -- cs, p=10,
# delta=12, cluster recovery 1.00: BC median p 0.0011 / max 0.0342 / power 1.000;
# trace median p 0.506 / power 0.063, KS-uniform p 0.187.
#
# WHY, structurally: the trace-valid region is a single tiny arc containing
# theta_obs, so the p-value is essentially the relative position of theta_obs inside
# its OWN arc -- uniform whether or not a signal is present. Conditioning on the
# exact merge sequence conditions away the quantity being tested. No algorithmic
# improvement recovers this; it is a property of the event.
#
# So this engine is EXACT and VALID (Type-I 0.044-0.071 across 12 cells, KS p >=
# 0.12, including p=50 at n=200 where the march cannot run) but is NOT an
# inferential method. Keep it as a negative result about over-conditioning; use
# selective_hotelling_bc(). See aug-task.md 1m.
#
# LINKAGE: derived and verified for AVERAGE linkage only. The `upd` recursion
# below is the unweighted Lance-Williams update, so passing a centroid/single/
# complete dendrogram produces the WRONG constraint set. Callers must guard.
#
# THE LIFETIME REDUCTION. The march re-compares every active cluster pair at every
# merge step, which is sum_s C(n_a, 2) ~ n^3/6 constraints (1.3M at n=200). But a
# pair that loses at several consecutive steps only needs its LAST losing step:
# with q = n - K, W_s the pair that merges at step s, h_s(t) = d_t(W_s) its
# time-varying height, and u(P) the last step at which P is active and loses,
#
#     trace T holds at t   <=>   d_P(t) > h_{u(P)}(t)  for every losing pair P.
#
# Necessity is immediate. Sufficiency: the endpoint constraints force the heights
# to be nondecreasing. For consecutive steps s, s+1 -- if W_{s+1} already existed
# at step s it lost there, so h_{s+1} = d_{W_{s+1}} > h_s; otherwise
# W_{s+1} = {A_s u B_s, C} and average linkage gives h_{s+1} as the |A|,|B|-weighted
# mean of d(A_s,C) and d(B_s,C), both of which die at step s and are therefore
# endpoint-constrained above h_s, so their weighted mean is too. Then any P losing
# at step s has u(P) >= s, so d_P(t) > h_{u(P)}(t) >= h_s(t) and W_s really is
# globally closest at every step. Note the heights are NOT assumed monotone a
# priori -- the constraints certify monotonicity exactly on the trace-valid region,
# which is what makes this work on a path where every row moves differently.
#
# Each pair identity is born once and dies once, so there are at most (n-1)^2
# constraints (39,601 at n=200, ~34x fewer than the stepwise count) and they are
# built in O(n^2) after the O(n^2 p) point-pair precomputation. d_P and h_u are
# both harmonics, so each constraint is a harmonic and stays quartic-solvable.
#' @keywords internal
#' @noRd
.hclust_trace_constraints <- function(pairH, hc_obs, K, method = "average") {
  # AVERAGE LINKAGE ONLY. The sufficiency step needs h_{s+1} to be the size-weighted
  # MEAN of two step-s endpoint-constrained distances. Centroid subtracts
  # |G||H|d(G,H)/(|G|+|H|)^2 from that mean and can fall below h_s; single/complete
  # are min/max envelopes, not a single harmonic. Guarded here as well as at the
  # caller so the recursion below cannot be reached with the wrong update rule.
  if (!identical(method, "average"))
    stop(sprintf("lifetime/endpoint trace constraints are derived for average linkage only (got '%s')",
                 method))
  n  <- nrow(pairH$A)
  q  <- n - K
  MA <- pairH$A; MB <- pairH$B; MC <- pairH$C; MD <- pairH$D; ME <- pairH$E
  sz   <- rep(1, n)
  born <- matrix(1L, n, n)          # step at which each pair identity became active
  cA <- cB <- cC <- cD <- cE <- vector("list", 2L * q + 2L)
  np <- 0L
  # constraint is d_P(t) > h_u(t), i.e. g = h_u - d_P <= 0
  push <- function(h, dA, dB, dC, dD, dE) {
    if (!length(dA)) return(invisible(NULL))
    np <<- np + 1L
    cA[[np]] <<- h[1] - dA; cB[[np]] <<- h[2] - dB; cC[[np]] <<- h[3] - dC
    cD[[np]] <<- h[4] - dD; cE[[np]] <<- h[5] - dE
  }
  mg <- hc_obs$merge
  node <- vector("list", n - 1L); lf <- rt <- vector("list", n - 1L)
  for (s in seq_len(n - 1L)) {
    a <- mg[s, 1]; b <- mg[s, 2]
    Gi <- if (a < 0) -a else node[[a]]
    Hi <- if (b < 0) -b else node[[b]]
    lf[[s]] <- Gi; rt[[s]] <- Hi; node[[s]] <- c(Gi, Hi)
  }
  akey <- seq_len(n); prevh <- NULL
  for (s in seq_len(q)) {
    Gi <- lf[[s]]; Hi <- rt[[s]]
    ix <- match(min(Gi), akey); iy <- match(min(Hi), akey)
    na <- length(akey)
    h  <- c(MA[ix, iy], MB[ix, iy], MC[ix, iy], MD[ix, iy], ME[ix, iy])
    oth <- seq_len(na)[-c(ix, iy)]
    # every {ix,C} and {iy,C} is active at step s, loses, and dies here: u = s
    for (w in c(ix, iy))
      push(h, MA[w, oth], MB[w, oth], MC[w, oth], MD[w, oth], ME[w, oth])
    # the winning pair itself last lost at s-1, if it was already active then
    if (s >= 2L && !is.null(prevh) && born[ix, iy] <= s - 1L)
      push(prevh, h[1], h[2], h[3], h[4], h[5])
    prevh <- h
    wG <- sz[ix]; wH <- sz[iy]; wn <- wG + wH
    upd <- function(M) { nr <- ((wG * M[ix, ] + wH * M[iy, ]) / wn)[oth]
                         rbind(cbind(M[oth, oth, drop = FALSE], nr), c(nr, 0)) }
    MA <- upd(MA); MB <- upd(MB); MC <- upd(MC); MD <- upd(MD); ME <- upd(ME)
    nb   <- rep(s + 1L, length(oth))
    born <- rbind(cbind(born[oth, oth, drop = FALSE], nb), c(nb, s + 1L))
    sz   <- c(sz[oth], wn)
    akey <- c(akey[oth], min(min(Gi), min(Hi)))
  }
  # pairs still active after the cut never lose again: u = q
  na <- length(akey)
  if (na >= 2L && !is.null(prevh)) {
    ut <- upper.tri(matrix(0, na, na))
    push(prevh, MA[ut], MB[ut], MC[ut], MD[ut], ME[ut])
  }
  length(cA) <- np; length(cB) <- np; length(cC) <- np
  length(cD) <- np; length(cE) <- np
  list(A = unlist(cA), B = unlist(cB), C = unlist(cC),
       D = unlist(cD), E = unlist(cE))
}

# One-shot trace-valid arcs: build the O(n^2) lifetime constraints, intersect.
#' @keywords internal
#' @noRd
.hclust_trace_arcs <- function(P, U, W, K, hc_obs, lo = 0, hi = pi / 2,
                               pairH = NULL, method = "average") {
  if (is.null(pairH)) pairH <- .hclust_pointpair_harmonics(P, U, W)
  co <- .hclust_trace_constraints(pairH, hc_obs, K, method = method)
  arcs <- .harmonic_feasible_arcs(co, lo, hi)
  attr(arcs, "n_constraints") <- length(co$A)
  attr(arcs, "capped") <- FALSE          # fixed conjunction: nothing to run out of
  arcs
}

# ---- feasible arcs by a SWEEP over sorted roots (large constraint sets) --------
# Same output as .harmonic_feasible_arcs, but it scales to tens of thousands of
# constraints. The scalar path solves one quartic per constraint and then tests
# EVERY constraint at the midpoint of EVERY segment, which is O(n_roots x n_C).
# Here the roots come from .harmonic_roots_batch, all constraints are evaluated
# once at the first midpoint, and at each later breakpoint only the constraints
# that have a root there are re-evaluated. A constraint can change sign only at
# one of its own roots, so the count of violated constraints stays exact.
#' @keywords internal
#' @noRd
.harmonic_feasible_arcs_sweep <- function(coeffs, lo = 0, hi = 2 * pi,
                                          tol_mid = 1e-10) {
  A <- coeffs$A; B <- coeffs$B; C <- coeffs$C; D <- coeffs$D; E <- coeffs$E
  empty <- matrix(numeric(0), ncol = 2); colnames(empty) <- c("t_lo", "t_hi")
  if (length(A)) {
    M <- sqrt(A^2 + B^2) + sqrt(C^2 + D^2)      # same amplitude prune as above
    if (any((E - M) > 1e-9)) return(empty)
    keep <- (E + M) >= -1e-9
    A <- A[keep]; B <- B[keep]; C <- C[keep]; D <- D[keep]; E <- E[keep]
  }
  nC <- length(A)
  if (!nC) { out <- matrix(c(lo, hi), 1L); colnames(out) <- c("t_lo", "t_hi"); return(out) }
  gj <- function(j, t) A[j] * cos(2 * t) + B[j] * sin(2 * t) + C[j] * cos(t) +
                       D[j] * sin(t) + E[j]
  rb <- .harmonic_roots_batch(A, B, C, D, E)
  in_rng <- rb$t > lo & rb$t < hi
  rt <- rb$t[in_rng]; rj <- rb$j[in_rng]
  bks <- sort(unique(c(lo, rt, hi)))
  owner <- split(rj, factor(match(rt, bks), levels = seq_along(bks)))
  nseg <- length(bks) - 1L
  viol <- gj(seq_len(nC), 0.5 * (bks[1] + bks[2])) > tol_mid
  nviol <- sum(viol)
  feas <- logical(nseg); feas[1] <- (nviol == 0L)
  if (nseg >= 2L) for (s in 2:nseg) {
    js <- owner[[s]]                           # constraints with a root at bks[s]
    if (length(js)) {
      js <- unique(js)
      nv <- gj(js, 0.5 * (bks[s] + bks[s + 1L])) > tol_mid
      nviol <- nviol + sum(nv) - sum(viol[js]); viol[js] <- nv
    }
    feas[s] <- (nviol == 0L)
  }
  if (!any(feas)) return(empty)
  # merge runs of feasible segments into maximal arcs
  r <- rle(feas); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  out <- cbind(bks[starts[r$values]], bks[ends[r$values] + 1L])
  colnames(out) <- c("t_lo", "t_hi")
  out
}

# ---- SINGLE linkage: the observed Kruskal-edge event ---------------------------
# Single linkage merges exactly along Kruskal's minimum spanning tree, so the
# first q = n - K merges are the first q Kruskal edges e_1, ..., e_q (point pairs;
# on squared distances the order is the same as on distances). We condition on
# those EDGES, as specific point pairs and in order. With t*(a, b) the first step
# at which a and b lie in the same forest component (capped at q), the event is
#   (i)  d(e_s) <= d(e_{s+1})               for s < q;
#   (ii) d(a, b) >= d(e_{t*(a, b)})          for every non-edge pair (a, b).
# Sufficiency: at step s every pair still crossing two components has
# t* >= s, so by (ii) and (i) d(a, b) >= d(e_{t*}) >= d(e_s), and e_s is the
# shortest available edge. Necessity: (a, b) crosses two components at every
# step up to t*, so it must lose to each of e_1, ..., e_{t*}.
#
# Each inequality is a difference of two point-pair harmonics, so this is a fixed
# conjunction of at most n(n-1)/2 - 1 harmonic constraints (no cell march, no
# budget). The event can be slightly SMALLER than "same first q cluster merges":
# a merge can keep the same two clusters while the closest point pair between
# them changes. It is contained in the merge-order event and so in the pair
# event, which keeps the pair weight constant on it.
#' @keywords internal
#' @noRd
.hclust_single_kruskal_constraints <- function(pairH, hc_obs, K, D2) {
  n <- nrow(pairH$A)
  q <- n - K
  mg <- hc_obs$merge
  node <- vector("list", q)
  ea <- eb <- integer(q)
  tstar <- matrix(q, n, n)                 # pairs never joined by step q: t* = q
  for (s in seq_len(q)) {
    a <- mg[s, 1]; b <- mg[s, 2]
    Gi <- if (a < 0) -a else node[[a]]
    Hi <- if (b < 0) -b else node[[b]]
    node[[s]] <- c(Gi, Hi)
    # the merge is realised by the closest observed point pair across G x H
    k <- which.min(D2[Gi, Hi, drop = FALSE])
    ea[s] <- Gi[(k - 1L) %% length(Gi) + 1L]
    eb[s] <- Hi[(k - 1L) %/% length(Gi) + 1L]
    tstar[Gi, Hi] <- s; tstar[Hi, Gi] <- s
  }
  ed <- cbind(ea, eb)
  hA <- pairH$A[ed]; hB <- pairH$B[ed]; hC <- pairH$C[ed]
  hD <- pairH$D[ed]; hE <- pairH$E[ed]
  # (ii) every non-edge pair stays at least as long as its joining edge
  ise <- matrix(FALSE, n, n); ise[ed] <- TRUE; ise[ed[, 2:1, drop = FALSE]] <- TRUE
  ut <- which(upper.tri(ise) & !ise)
  ts <- tstar[ut]
  A <- hA[ts] - pairH$A[ut]; B <- hB[ts] - pairH$B[ut]; C <- hC[ts] - pairH$C[ut]
  D <- hD[ts] - pairH$D[ut]; E <- hE[ts] - pairH$E[ut]
  # (i) the edge heights are nondecreasing
  if (q >= 2L) {
    s <- seq_len(q - 1L)
    A <- c(A, hA[s] - hA[s + 1L]); B <- c(B, hB[s] - hB[s + 1L])
    C <- c(C, hC[s] - hC[s + 1L]); D <- c(D, hD[s] - hD[s + 1L])
    E <- c(E, hE[s] - hE[s + 1L])
  }
  list(A = A, B = B, C = C, D = D, E = E, edges = ed)
}

# One-shot Kruskal-edge arcs for single linkage: build, then sweep.
#' @keywords internal
#' @noRd
.hclust_single_kruskal_arcs <- function(P, U, W, K, hc_obs, D2, lo = 0, hi = pi / 2,
                                        pairH = NULL) {
  if (is.null(pairH)) pairH <- .hclust_pointpair_harmonics(P, U, W)
  co <- .hclust_single_kruskal_constraints(pairH, hc_obs, K, D2)
  arcs <- .harmonic_feasible_arcs_sweep(co, lo, hi)
  attr(arcs, "n_constraints") <- length(co$A)
  attr(arcs, "capped") <- FALSE          # fixed conjunction: nothing to run out of
  arcs
}
