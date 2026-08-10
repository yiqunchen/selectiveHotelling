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
  if (prune && nC > 0L) {                       # drop constraints that never bind
    M <- sqrt(A^2 + B^2) + sqrt(C^2 + D^2)
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

# ---- observed agglomeration structure (for hierarchical constraints) ----------
# Runs hclust ONCE on the observed X; returns, for each of the first (n-K) merge
# steps, the merged pair (as original-index sets) and the list of competitor
# cluster pairs (all other active-cluster pairs). This structure is FIXED (we
# condition on the observed trajectory) and does not depend on theta.
#' @keywords internal
#' @noRd
.hclust_merge_structure <- function(X, K, method = "average", square = TRUE,
                                    hc = NULL) {
  n <- nrow(X)
  # `hc` lets the caller pass the hclust it already built for the cell labelling.
  if (is.null(hc)) {
    d <- stats::dist(X); if (isTRUE(square)) d <- d^2
    hc <- stats::hclust(d, method = method)
  }
  mg <- hc$merge
  node <- vector("list", n - 1L)             # index set formed at each merge
  left_of <- right_of <- vector("list", n - 1L)
  for (s in seq_len(n - 1L)) {
    a <- mg[s, 1]; b <- mg[s, 2]
    Gi <- if (a < 0) -a else node[[a]]
    Hi <- if (b < 0) -b else node[[b]]
    left_of[[s]] <- Gi; right_of[[s]] <- Hi
    node[[s]] <- c(Gi, Hi)
  }
  # active clusters just before each step; only the first (n-K) merges matter.
  # Clusters are identified by min(members): they are disjoint, so it is unique,
  # and it avoids a string key built O(n^2) times per step.
  active <- as.list(seq_len(n))              # start: singletons {1},...,{n}
  akey   <- vapply(active, min, integer(1))  # parallel vector of integer ids
  steps  <- vector("list", n - K)
  for (s in seq_len(n - K)) {
    Gi <- left_of[[s]]; Hi <- right_of[[s]]
    kG <- min(Gi); kH <- min(Hi)
    keep_other <- !(akey %in% c(kG, kH))
    others     <- active[keep_other]
    na    <- length(active)
    npair <- if (na >= 2L) na * (na - 1L) / 2L - 1L else 0L   # all pairs minus (G,H)
    comp  <- vector("list", max(npair, 0L))
    nc <- 0L
    if (na >= 2L) for (x in seq_len(na - 1L)) for (y in (x + 1L):na) {
      kx <- akey[x]; ky <- akey[y]
      if ((kx == kG && ky == kH) || (kx == kH && ky == kG)) next
      nc <- nc + 1L
      comp[[nc]] <- list(G = active[[x]], H = active[[y]])
    }
    length(comp) <- nc
    steps[[s]] <- list(mG = Gi, mH = Hi, comp = comp)
    active <- c(others, list(c(Gi, Hi)))     # perform the merge
    akey   <- c(akey[keep_other], min(kG, kH))
  }
  list(steps = steps, n = n, K = K)
}

# ---- hierarchical agglomeration-trajectory constraint harmonics ---------------
# For each merge step, "observed merged-pair linkage <= competitor linkage" is a
# harmonic <= 0. This holds whenever the linkage distance d(G,H) reduces to a
# SINGLE harmonic along the geodesic. Two such (SINGLE-HARMONIC) linkages:
#   * "average"  (UPGMA): d(G,H) = MEAN over G x H of the point-pair squared-
#                 distance harmonics.
#   * "centroid" (UPGMC, defined on SQUARED Euclidean distances): d(G,H) =
#                 ||cbar_G - cbar_H||^2 with cbar the member centroid -> the
#                 squared-distance harmonic of the centroid gap (member centroids
#                 are affine in the geodesic, so this is one harmonic).
#
# SAFETY: single/complete linkage make d(G,H) a min/max ENVELOPE of harmonics, and
# median/McQuitty/Ward use different Lance-Williams weights or a non-squared
# convention -- for those the reductions below are the WRONG constraint and would
# silently misplace the arc boundaries. We HARD-GUARD to the supported set; callers
# must route any other linkage to the (linkage-correct) grid engine.
#' @keywords internal
#' @noRd
.hclust_traj_coeffs <- function(P, U, W, mstruct, method = "average") {
  supported <- c("average", "centroid", "single", "complete")
  if (!(method %in% supported))
    stop(sprintf(paste0("closed-form hclust breakpoint harmonics are implemented for ",
                        "{%s} (got method='%s'); route the exact hclust test to the ",
                        "grid engine for this linkage."),
                 paste(supported, collapse = ", "), method))
  n <- nrow(P)

  # Point-pair harmonics (n x n) are needed by average and by the single/complete
  # envelope; centroid uses member centroids instead.
  if (method %in% c("average", "single", "complete")) {
    Amat <- Bmat <- Cmat <- Dmat <- Emat <- matrix(0, n, n)
    for (i in seq_len(n - 1L)) {
      js <- (i + 1L):n
      dP <- P[js, , drop = FALSE] - P[rep(i, length(js)), , drop = FALSE]
      dU <- U[js, , drop = FALSE] - U[rep(i, length(js)), , drop = FALSE]
      dW <- W[js, , drop = FALSE] - W[rep(i, length(js)), , drop = FALSE]
      h  <- .sqdist_harmonic(dP, dU, dW)
      Amat[i, js] <- h$A; Amat[js, i] <- h$A
      Bmat[i, js] <- h$B; Bmat[js, i] <- h$B
      Cmat[i, js] <- h$C; Cmat[js, i] <- h$C
      Dmat[i, js] <- h$D; Dmat[js, i] <- h$D
      Emat[i, js] <- h$E; Emat[js, i] <- h$E
    }
  }

  # ---- SINGLE-HARMONIC linkages (average, centroid): one harmonic per comparison
  if (method %in% c("average", "centroid")) {
    link <- if (identical(method, "average"))
      function(G, H) c(mean(Amat[G, H]), mean(Bmat[G, H]), mean(Cmat[G, H]),
                       mean(Dmat[G, H]), mean(Emat[G, H]))
    else
      function(G, H) {                     # d(G,H) = ||cbar_G - cbar_H||^2
        dP <- colMeans(P[G, , drop = FALSE]) - colMeans(P[H, , drop = FALSE])
        dU <- colMeans(U[G, , drop = FALSE]) - colMeans(U[H, , drop = FALSE])
        dW <- colMeans(W[G, , drop = FALSE]) - colMeans(W[H, , drop = FALSE])
        h  <- .sqdist_harmonic(matrix(dP, 1L), matrix(dU, 1L), matrix(dW, 1L))
        c(h$A, h$B, h$C, h$D, h$E)
      }
    A <- B <- C <- D <- E <- numeric(0)
    for (st in mstruct$steps) {
      lm <- link(st$mG, st$mH)
      for (cp in st$comp) {
        lc <- link(cp$G, cp$H)
        A <- c(A, lm[1] - lc[1]); B <- c(B, lm[2] - lc[2]); C <- c(C, lm[3] - lc[3])
        D <- c(D, lm[4] - lc[4]); E <- c(E, lm[5] - lc[5])
      }
    }
    return(list(A = A, B = B, C = C, D = D, E = E))
  }

  # ---- ENVELOPE linkages (single = min, complete = max over member pairs) -------
  # d(G,H) is a min/max envelope of the |G||H| member-pair harmonics, so it is NOT
  # a single harmonic. A merge decision (order change) flips exactly where two
  # linkage distances cross, i.e. where min(S)=min(T) [single] or max(S)=max(T)
  # [complete] for competing cluster-pairs S, T -- and BOTH occur at a root of some
  # cross-difference h_s - h_t (s in S, t in T). So the roots of all cross-
  # differences are a SUPERSET of the true order-change breakpoints; the cell
  # traversal re-clusters for each cell's label, so an over-complete breakpoint set
  # still yields the EXACT arcs. We prune differences that can never cross zero
  # (|E| > amplitude -> no root) to keep the count manageable. (Same breakpoint set
  # serves single and complete; only the per-cell label -- from the actual hclust
  # run with `method` -- distinguishes them.)
  pair_harm <- function(G, H) {                # |G||H| x 5 harmonic matrix
    idx <- cbind(rep(G, each = length(H)), rep(H, times = length(G)))
    cbind(Amat[idx], Bmat[idx], Cmat[idx], Dmat[idx], Emat[idx])
  }
  A <- B <- C <- D <- E <- numeric(0)
  for (st in mstruct$steps) {
    S <- pair_harm(st$mG, st$mH)
    for (cp in st$comp) {
      Tt <- pair_harm(cp$G, cp$H)
      nS <- nrow(S); nT <- nrow(Tt)
      si <- rep(seq_len(nS), each = nT); tj <- rep(seq_len(nT), times = nS)
      dA <- S[si, 1] - Tt[tj, 1]; dB <- S[si, 2] - Tt[tj, 2]; dC <- S[si, 3] - Tt[tj, 3]
      dD <- S[si, 4] - Tt[tj, 4]; dE <- S[si, 5] - Tt[tj, 5]
      amp  <- sqrt(dA^2 + dB^2) + sqrt(dC^2 + dD^2)   # max |harmonic - E|
      keep <- abs(dE) <= amp + 1e-9                   # can cross zero
      A <- c(A, dA[keep]); B <- c(B, dB[keep]); C <- c(C, dC[keep])
      D <- c(D, dD[keep]); E <- c(E, dE[keep])
    }
  }
  list(A = A, B = B, C = C, D = D, E = E)
}

# ---- EXACT hclust PAIR-PRESERVED arc via event-driven cell traversal ----------
# Conditioning on the full agglomeration trajectory over-conditions and destroys
# power (empirically ~5x lower); the scientifically correct hclust event -- the
# one selective_hotelling_bc and Gao et al. (clusterpval) use -- is that the
# cutree-K PARTITION preserves the tested pair. Along the geodesic that partition
# is piecewise-constant with harmonic-root boundaries. We enumerate its cells
# EXACTLY: at each cell probe the current agglomeration order (one hclust run),
# build its first-(n-K) merge constraints, and jump to the smallest harmonic root
# beyond the current theta (the next order change). Each cell is labelled by
# whether the pair is preserved; the feasible set is the union of preserved cells.
# Finite hclust runs (one per cell, typically a handful), exact boundaries.
#' @keywords internal
#' @noRd
.hclust_pairpreserved_arcs <- function(P, U, W, X_of_Q, Qtheta, K,
                                       cl_obs, k1, k2, method = "average",
                                       square = TRUE, lo = 0, hi = pi / 2,
                                       eps = 1e-7, max_cells = 500L) {
  labels_at <- function(th) {
    Xt <- X_of_Q(Qtheta(th)); d <- stats::dist(Xt); if (isTRUE(square)) d <- d^2
    as.integer(stats::cutree(stats::hclust(d, method = method), k = K))
  }
  cells_lo <- numeric(0); cells_hi <- numeric(0); cells_pp <- logical(0)
  th <- lo; guard <- 0L
  while (th < hi - 1e-12 && guard < max_cells) {
    guard <- guard + 1L
    probe <- min(th + eps, hi - 1e-12)
    Xp <- X_of_Q(Qtheta(probe))
    # one dist+hclust per cell, shared with .hclust_merge_structure below
    d_cell <- stats::dist(Xp); if (isTRUE(square)) d_cell <- d_cell^2
    hc_cell <- stats::hclust(d_cell, method = method)
    cl_cell <- as.integer(stats::cutree(hc_cell, k = K))
    ms_cell <- .hclust_merge_structure(Xp, K, method, square, hc = hc_cell)
    co <- .hclust_traj_coeffs(P, U, W, ms_cell, method = method)
    nb <- hi
    if (length(co$A)) {
      # |g_j| <= M + |E| with M = sqrt(A^2+B^2) + sqrt(C^2+D^2), so |E| > M means
      # g_j never crosses zero and contributes no breakpoint. Output-preserving:
      # nb is a min over roots and pruned constraints have none.
      Aj <- co$A; Bj <- co$B; Cj <- co$C; Dj <- co$D; Ej <- co$E
      Mj <- sqrt(Aj * Aj + Bj * Bj) + sqrt(Cj * Cj + Dj * Dj)
      cand <- which(abs(Ej) <= Mj + 1e-12)
      for (jj in cand) {
        rts <- .harmonic_roots(Aj[jj], Bj[jj], Cj[jj], Dj[jj], Ej[jj])
        rts <- rts[rts > th + 1e-9 & rts <= hi + 1e-9]
        if (length(rts)) nb <- min(nb, min(rts))
      }
    }
    if (nb <= th + 1e-12) nb <- min(th + 1e-4, hi)   # numeric floor (never stall)
    cells_lo <- c(cells_lo, th); cells_hi <- c(cells_hi, nb)
    cells_pp <- c(cells_pp, .bc_pair_preserved(cl_obs, cl_cell, k1, k2))
    th <- nb
  }
  # merge adjacent preserved cells into maximal arcs
  arcs <- matrix(numeric(0), ncol = 2)
  i <- 1L; nC <- length(cells_lo)
  while (i <= nC) {
    if (!cells_pp[i]) { i <- i + 1L; next }
    lo_i <- cells_lo[i]; hi_i <- cells_hi[i]
    while (i + 1L <= nC && cells_pp[i + 1L] && abs(cells_lo[i + 1L] - hi_i) <= 1e-9) {
      hi_i <- cells_hi[i + 1L]; i <- i + 1L
    }
    arcs <- rbind(arcs, c(lo_i, hi_i)); i <- i + 1L
  }
  colnames(arcs) <- c("t_lo", "t_hi")
  arcs
}
