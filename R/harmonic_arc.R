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

# ---- observed agglomeration structure (for hierarchical constraints) ----------
# Runs hclust ONCE on the observed X; returns, for each of the first (n-K) merge
# steps, the merged pair (as original-index sets) and the list of competitor
# cluster pairs (all other active-cluster pairs). This structure is FIXED (we
# condition on the observed trajectory) and does not depend on theta.
#' @keywords internal
#' @noRd
.hclust_merge_structure <- function(X, K, method = "average", square = TRUE,
                                    hc = NULL, want_comp = FALSE) {
  n <- nrow(X)
  # `hc` lets the caller pass the hclust it already built for the cell labelling;
  # `want_comp` materialises the O(n^3) competitor-pair list, which only the
  # envelope linkages (single/complete) still need.
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
    na <- length(active)
    # POSITIONS of the merging clusters in the active list. The Lance-Williams
    # walker replays the same remove-two-append-one bookkeeping, so positions are
    # enough and the O(n^2)-per-step `comp` list need not be materialised at all
    # (it is O(n^3) overall, and only the envelope linkages still need it).
    ix <- match(min(Gi), akey); iy <- match(min(Hi), akey)
    st <- list(mG = Gi, mH = Hi, ix = ix, iy = iy, na = na)
    if (isTRUE(want_comp)) {
      npair <- if (na >= 2L) na * (na - 1L) / 2L - 1L else 0L
      comp <- vector("list", max(npair, 0L)); nc <- 0L
      if (na >= 2L) for (x in seq_len(na - 1L)) for (y in (x + 1L):na) {
        if ((x == ix && y == iy) || (x == iy && y == ix)) next
        nc <- nc + 1L
        comp[[nc]] <- list(G = active[[x]], H = active[[y]])
      }
      length(comp) <- nc; st$comp <- comp
    }
    steps[[s]] <- st
    drop   <- c(ix, iy)
    active <- c(active[-drop], list(c(Gi, Hi)))     # perform the merge
    akey   <- c(akey[-drop], min(min(Gi), min(Hi)))
  }
  list(steps = steps, n = n, K = K)
}

# log(e^a + e^b), NA-safe in the -Inf ("zero mass") case.
#' @keywords internal
#' @noRd
.lse2 <- function(a, b) {
  if (!is.finite(a)) return(b)
  if (!is.finite(b)) return(a)
  m <- max(a, b); m + log(exp(a - m) + exp(b - m))
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

# ---- next merge-order breakpoint after `th` (average / centroid) --------------
# Replaces the "materialise every constraint, then solve" path for the two
# SINGLE-HARMONIC linkages. Cluster-pair linkage harmonics obey the same
# Lance-Williams recursion as the distances they encode (the identity holds at
# every t, and both sides are harmonics of the same form, so it holds coefficient
# by coefficient), so the na x na coefficient matrices are UPDATED across merges
# instead of recomputed from members. Per step the constraints are then a
# vectorised matrix subtraction, and only the ones that can cross zero are sent to
# the root solver. Returns the smallest root in (th, hi], or Inf if none.
#' @keywords internal
#' @noRd
.hclust_traj_next_root <- function(pairH, mstruct, th, hi, method = "average") {
  if (!(method %in% c("average", "centroid")))
    stop("Lance-Williams breakpoint path supports average/centroid only")
  MA <- pairH$A; MB <- pairH$B; MC <- pairH$C; MD <- pairH$D; ME <- pairH$E
  sz <- rep(1, mstruct$n)
  best <- Inf
  for (st in mstruct$steps) {
    ix <- st$ix; iy <- st$iy; na <- nrow(MA)
    if (na >= 2L) {
      ut <- upper.tri(MA)
      # merged-pair linkage minus every competitor's: g_j(t) <= 0 must hold.
      dA <- MA[ix, iy] - MA[ut]; dB <- MB[ix, iy] - MB[ut]
      dC <- MC[ix, iy] - MC[ut]; dD <- MD[ix, iy] - MD[ut]; dE <- ME[ix, iy] - ME[ut]
      # LIPSCHITZ PRUNE. Root-finding costs ~150us per constraint (a companion-
      # matrix solve) and there are O(n^3) of them, so the amplitude test |E| <= M
      # -- which only rejects ~15% -- is nowhere near selective enough. Instead use
      # |g_j'| <= L_j = 2 sqrt(A^2+B^2) + sqrt(C^2+D^2): from g_j(th) <= 0,
      #     g_j(t) <= g_j(th) + L_j (t - th)   for t >= th,
      # so a constraint whose bound stays negative across (th, ub] has NO root
      # there and can be dropped outright. Cells are narrow (ub - th is small), so
      # this rejects almost everything. It can never discard a real crossing.
      g0 <- dA * cos(2 * th) + dB * sin(2 * th) + dC * cos(th) + dD * sin(th) + dE
      Lj <- 2 * sqrt(dA * dA + dB * dB) + sqrt(dC * dC + dD * dD)
      ub <- min(best, hi)
      cand <- which(Lj > 1e-12 & g0 + Lj * (ub - th) >= -1e-12)
      # Solve in chunks, largest g_j(th) first: those cross soonest, so `best`
      # tightens early and the remaining chunks re-prune against it.
      if (length(cand)) {
        ord <- cand[order(g0[cand], decreasing = TRUE)]
        i <- 1L
        while (i <= length(ord)) {
          j  <- ord[i:min(i + 511L, length(ord))]
          ub <- min(best, hi)
          j  <- j[g0[j] + Lj[j] * (ub - th) >= -1e-12]
          if (length(j)) {
            rb  <- .harmonic_roots_batch(dA[j], dB[j], dC[j], dD[j], dE[j])
            rts <- rb$t[rb$t > th + 1e-9 & rb$t <= ub + 1e-9]
            if (length(rts)) best <- min(best, min(rts))
          }
          i <- i + 512L
        }
      }
    }
    wG <- sz[ix]; wH <- sz[iy]; wn <- wG + wH
    keep <- seq_len(na)[-c(ix, iy)]
    upd <- function(M) {
      nr <- (wG * M[ix, ] + wH * M[iy, ]) / wn
      # centroid (UPGMC) carries the extra -wG wH d(G,H)/wn^2 Lance-Williams term
      if (identical(method, "centroid")) nr <- nr - wG * wH * M[ix, iy] / (wn * wn)
      nr <- nr[keep]
      rbind(cbind(M[keep, keep, drop = FALSE], nr), c(nr, 0))
    }
    MA <- upd(MA); MB <- upd(MB); MC <- upd(MC); MD <- upd(MD); ME <- upd(ME)
    sz <- c(sz[keep], wn)
  }
  best
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
.hclust_traj_coeffs <- function(P, U, W, mstruct, method = "average",
                                pairH = NULL) {
  supported <- c("average", "centroid", "single", "complete")
  if (!(method %in% supported))
    stop(sprintf(paste0("closed-form hclust breakpoint harmonics are implemented for ",
                        "{%s} (got method='%s'); route the exact hclust test to the ",
                        "grid engine for this linkage."),
                 paste(supported, collapse = ", "), method))
  n <- nrow(P)

  # Point-pair harmonics (n x n) are needed by average and by the single/complete
  # envelope; centroid uses member centroids instead. `pairH` lets a caller that
  # traverses many cells on ONE geodesic build them once.
  if (method %in% c("average", "single", "complete")) {
    if (is.null(pairH)) pairH <- .hclust_pointpair_harmonics(P, U, W)
    Amat <- pairH$A; Bmat <- pairH$B; Cmat <- pairH$C
    Dmat <- pairH$D; Emat <- pairH$E
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

# ---- ONE-SHOT arcs for a FIXED average-linkage merge trace --------------------
# Conditions on the OBSERVED merge sequence up to the K-cluster cut, rather than on
# the weaker "cutree-K preserves the tested pair". That is strictly more
# conditioning -- R_T is one component of the union S = U_T R_T over traces whose
# cut preserves the pair -- and it must be reported as such. What it buys is that
# the event becomes a FIXED conjunction, so there is no cell march, no cell budget,
# and hence no max_cells truncation and no fabricated p = 1.
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
.hclust_trace_constraints <- function(pairH, hc_obs, K) {
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
                               pairH = NULL) {
  if (is.null(pairH)) pairH <- .hclust_pointpair_harmonics(P, U, W)
  co <- .hclust_trace_constraints(pairH, hc_obs, K)
  arcs <- .harmonic_feasible_arcs(co, lo, hi)
  attr(arcs, "n_constraints") <- length(co$A)
  attr(arcs, "capped") <- FALSE          # fixed conjunction: nothing to run out of
  arcs
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
                                       eps = 1e-7, max_cells = 500L,
                                       theta_obs = NULL, lmass = NULL,
                                       mass_tol = 1e-6) {
  labels_at <- function(th) {
    Xt <- X_of_Q(Qtheta(th)); d <- stats::dist(Xt); if (isTRUE(square)) d <- d^2
    as.integer(stats::cutree(stats::hclust(d, method = method), k = K))
  }
  # Point-pair harmonics depend on the geodesic only, so they are built ONCE for
  # the whole traversal rather than rebuilt inside every cell.
  fast  <- method %in% c("average", "centroid")
  pairH <- .hclust_pointpair_harmonics(P, U, W)
  cells_lo <- numeric(0); cells_hi <- numeric(0); cells_pp <- logical(0)
  th <- lo; guard <- 0L; n_floor <- 0L
  use_mass <- !is.null(theta_obs) && !is.null(lmass)
  lN <- -Inf; lD <- -Inf; early <- FALSE; bracket <- NULL
  while (th < hi - 1e-12 && guard < max_cells) {
    guard <- guard + 1L
    probe <- min(th + eps, hi - 1e-12)
    Xp <- X_of_Q(Qtheta(probe))
    # one dist+hclust per cell, shared with .hclust_merge_structure below
    d_cell <- stats::dist(Xp); if (isTRUE(square)) d_cell <- d_cell^2
    hc_cell <- stats::hclust(d_cell, method = method)
    cl_cell <- as.integer(stats::cutree(hc_cell, k = K))
    ms_cell <- .hclust_merge_structure(Xp, K, method, square, hc = hc_cell,
                                       want_comp = !fast)
    nb <- hi
    if (fast) {
      r <- .hclust_traj_next_root(pairH, ms_cell, th, hi, method = method)
      if (is.finite(r)) nb <- min(nb, r)
    } else {
      co <- .hclust_traj_coeffs(P, U, W, ms_cell, method = method, pairH = pairH)
      if (length(co$A)) {
        # |g_j| <= M + |E| with M = sqrt(A^2+B^2) + sqrt(C^2+D^2), so |E| > M means
        # g_j never crosses zero and contributes no breakpoint. Output-preserving:
        # nb is a min over roots and pruned constraints have none.
        Aj <- co$A; Bj <- co$B; Cj <- co$C; Dj <- co$D; Ej <- co$E
        Mj <- sqrt(Aj * Aj + Bj * Bj) + sqrt(Cj * Cj + Dj * Dj)
        cand <- which(abs(Ej) <= Mj + 1e-12)
        if (length(cand)) {
          rb  <- .harmonic_roots_batch(Aj[cand], Bj[cand], Cj[cand], Dj[cand], Ej[cand])
          rts <- rb$t[rb$t > th + 1e-9 & rb$t <= hi + 1e-9]
          if (length(rts)) nb <- min(nb, min(rts))
        }
      }
    }
    if (nb <= th + 1e-12) {                # numeric floor (never stall)
      nb <- min(th + 1e-4, hi); n_floor <- n_floor + 1L
    }
    cells_lo <- c(cells_lo, th); cells_hi <- c(cells_hi, nb)
    pp <- .bc_pair_preserved(cl_obs, cl_cell, k1, k2)
    cells_pp <- c(cells_pp, pp)
    # CERTIFIED EARLY STOP. The p-value is a ratio of measure over the feasible
    # set, mass(S n [0,theta_obs)) / mass(S), and the measure of the UNEXPLORED
    # tail is available in closed form without knowing feasibility there --
    # feasibility can only remove mass. So with N, D accumulated so far,
    # R = mass(nb, hi] and R' = mass(nb, min(theta_obs, hi)],
    #     p >= N / (D + R)          (no future numerator, all future denominator)
    #     p <= (N + R') / (D + R')  (all future numerator feasible; D' >= N' always)
    # Once that bracket is tighter than mass_tol the remaining cells cannot move
    # the answer, so we stop EXACTLY rather than running out of budget. This
    # matters because the measure concentrates near hi while the march starts at
    # lo: without it the march spends its whole budget on negligible mass, caps
    # before reaching theta_obs, and then numerator == denominator returns a FAKE
    # p-value of exactly 1 (measured: 0% of replicates at n=30, 4.5-15.6% at
    # n=50, 100% at n=200).
    if (use_mass && pp) {
      lD <- .lse2(lD, lmass(th, nb))
      if (th < theta_obs) lN <- .lse2(lN, lmass(th, min(nb, theta_obs)))
    }
    th <- nb
    if (use_mass && is.finite(lD)) {
      lR  <- lmass(th, hi)
      lRp <- if (th < theta_obs) lmass(th, min(theta_obs, hi)) else -Inf
      b_lo <- exp(lN - .lse2(lD, lR))
      b_hi <- exp(.lse2(lN, lRp) - .lse2(lD, lRp))
      if (is.finite(b_lo) && is.finite(b_hi) && b_hi - b_lo <= mass_tol) {
        early <- TRUE; bracket <- c(b_lo, b_hi); break
      }
    }
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
  # Hitting max_cells means the march stopped before reaching `hi`, so the arc set
  # is INCOMPLETE and a p-value built from it has the wrong denominator. Warn and
  # record rather than returning a silently truncated answer. n_floor counts the
  # 1e-4 fallback steps, each a non-analytic boundary placement.
  # An EARLY stop is not a truncation: the bracket certifies the remaining mass
  # cannot move the p-value by more than mass_tol. Only running out of budget is.
  capped <- (guard >= max_cells) && (th < hi - 1e-12) && !early
  if (capped)
    warning(sprintf(paste0("hclust arc march hit max_cells=%d at theta=%.6g of %.6g: ",
                           "arc set is TRUNCATED, p-value not exact. Raise max_cells."),
                    max_cells, th, hi), call. = FALSE)
  attr(arcs, "n_cells") <- guard
  attr(arcs, "capped")  <- capped
  attr(arcs, "early")   <- early
  attr(arcs, "bracket") <- bracket
  attr(arcs, "reached") <- th
  attr(arcs, "n_floor") <- n_floor
  arcs
}
