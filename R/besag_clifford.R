# =============================================================================
# besag_clifford.R -- ALGORITHM-AGNOSTIC Besag-Clifford selective Hotelling test.
#
# This is the GENERIC, clustering-agnostic counterpart of the k-means selective
# test. It DROPS every k-means-specific assumption and conditions only on the
# abstract selection event
#
#     E = { Q in St(r-1, p) : cl_fun( X(Q) ) PRESERVES the tested pair (k1,k2) }
#
# for an ARBITRARY clustering routine cl_fun(X) -> integer label vector. The only
# things reused from the Stiefel machinery are:
#   * the polar decomposition X_sub = A Q S^{1/2} (prepare_pair_objects), which
#     under H0 makes Q UNIFORM on St(r-1, p) regardless of Sigma;
#   * the affine reconstruction X(Q), so a draw of Q maps back to a full data
#     matrix on which cl_fun can be re-run;
#   * the contrast statistic u(Q) = ||Q' v||^2, strictly monotone in the pair's
#     two-sample Hotelling T^2.
#
# WHY A METROPOLIS KERNEL.  For a GENERIC cl_fun we have NO closed form for E,
# only the black-box predicate 1{ cl_fun(X(Q)) preserves (k1,k2) }. So we use a
# Metropolis kernel:
#
#     K:  propose Q' = U Q  with U ~ Cayley(random skew, step eps)  (a SYMMETRIC
#         proposal: the law of Q'->Q equals the law of Q->Q', because the skew
#         generator G - G' is sign-symmetric and U, U^{-1} are equally likely);
#         accept Q' iff feasible(Q') := cl_fun(X(Q')) preserves (k1,k2),
#         otherwise stay at Q.
#
# For a SYMMETRIC proposal with a FLAT target restricted to E, the Metropolis
# acceptance ratio is exactly the indicator 1{Q' in E}. Hence K is reversible
# w.r.t. Uniform(E), and K^L (K applied L times) is also reversible w.r.t.
# Uniform(E).
#
# THE BESAG-CLIFFORD GUARANTEE (validity does NOT depend on the kernel).
#   Under H0, Q_obs ~ Uniform(E) = the stationary law of K. Plant Q_obs at a
#   uniform random index J in {0..M}; run K^L forward J..M and backward J..0.
#   Because the seed is stationary and K^L reversible, the rank of u(Q_obs)
#   among the M+1 path marks is EXACTLY Uniform{1..M+1}. So the randomized
#   upper-tail rank p-value is EXACTLY Uniform(0,1) under H0 -- for ANY
#   reversible K, ANY step eps, ANY path length M, ANY clustering cl_fun. Mixing
#   (governed by eps and the acceptance rate) buys POWER only, never validity.
#   (Besag & Clifford 1991.)
# =============================================================================

# -----------------------------------------------------------------------------
# .bc_pair_preserved(cl_obs, cl_new, k1, k2): does cl_new PRESERVE the (k1,k2)
# pair defined on cl_obs? Generic, algorithm-agnostic: the rows that were in
# cluster k1 (resp k2) under cl_obs must ALL share one label under cl_new, all
# rows in k2 must ALL share a DIFFERENT label, and no foreign rows may leak into
# those two labels. Allows label permutations. Returns TRUE/FALSE.
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.bc_pair_preserved <- function(cl_obs, cl_new, k1, k2) {
  i1 <- which(cl_obs == k1); i2 <- which(cl_obs == k2)
  if (!length(i1) || !length(i2)) return(FALSE)
  L1 <- unique(cl_new[i1]); L2 <- unique(cl_new[i2])
  if (length(L1) != 1L || length(L2) != 1L || L1 == L2) return(FALSE)
  if (!all(which(cl_new == L1) %in% i1)) return(FALSE)
  if (!all(which(cl_new == L2) %in% i2)) return(FALSE)
  TRUE
}

# -----------------------------------------------------------------------------
# .bc_setup(X, cl_obs, k1, k2, cl_fun): build the per-dataset polar objects and
# the generic feasibility predicate. cl_obs is the OBSERVED clustering
# (cl_fun(X)). Returns:
#   prep      polar objects (A, S_sqrt, X_fixed, idx1, idx2, m, r, Q_obs)
#   X_of_Q    affine reconstruction Q -> n x p data matrix
#   Q_obs     observed polar factor on St(r-1, p)
#   v         contrast s.t. u(Q) = ||Q' v||^2 is monotone in T^2 of (k1,k2)
#   feasible  predicate Q -> logical: cl_fun(X(Q)) preserves (k1,k2)
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.bc_setup <- function(X, cl_obs, k1, k2, cl_fun) {
  prep  <- prepare_pair_objects(X, cl_obs, k1, k2)
  m <- prep$m; r <- prep$r
  Q_obs <- prep$Q_obs

  # Affine reconstruction X(Q): only the pair rows move; everything else is
  # X_fixed.
  X_of_Q <- function(Q) {
    AQS <- prep$A %*% Q %*% prep$S_sqrt
    XQ <- prep$X_fixed
    XQ[prep$idx1, ] <- XQ[prep$idx1, ] + AQS[1:m, , drop = FALSE]
    XQ[prep$idx2, ] <- XQ[prep$idx2, ] + AQS[(m + 1):r, , drop = FALSE]
    XQ
  }

  # Contrast statistic, monotone in the pair's Hotelling T^2.
  n1 <- m; n2 <- r - m
  contrast <- c(rep(1 / n1, n1), rep(-1 / n2, n2))
  v <- as.numeric(crossprod(prep$A, contrast))

  # Generic feasibility: re-run cl_fun on the reconstructed data and check the
  # tested pair is preserved. cl_fun may return either a bare integer vector or a
  # list carrying $cluster; accept both.
  .as_labels <- function(z) if (is.list(z) && !is.null(z$cluster)) z$cluster else z
  feasible <- function(Q) {
    cl_new <- tryCatch(.as_labels(cl_fun(X_of_Q(Q))), error = function(e) NULL)
    if (is.null(cl_new) || length(cl_new) != length(cl_obs)) return(FALSE)
    .bc_pair_preserved(cl_obs, cl_new, k1, k2)
  }

  list(prep = prep, X_of_Q = X_of_Q, Q_obs = Q_obs, v = v,
       feasible = feasible, m = m, r = r)
}

# -----------------------------------------------------------------------------
# .bc_kernel_step(Q, eps, feasible): ONE symmetric-proposal Metropolis step on
# the Stiefel manifold. Propose Q' = Cayley(eps) Q; accept iff feasible(Q'),
# else stay. Returns list(Q, moved). Reversible w.r.t. Uniform(E).
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.bc_kernel_step <- function(Q, eps, feasible) {
  Qp <- cayley_left_step(Q, eps)          # symmetric proposal on St(r-1, p)
  if (isTRUE(feasible(Qp))) list(Q = Qp, moved = TRUE)
  else                      list(Q = Q,  moved = FALSE)
}

# .bc_kernel_L(Q, L, eps, feasible): K^L. Returns list(Q, n_moves in {0..L}).
#' @keywords internal
#' @noRd
.bc_kernel_L <- function(Q, L, eps, feasible) {
  n_moves <- 0L
  for (s in seq_len(L)) {
    st <- .bc_kernel_step(Q, eps, feasible)
    if (st$moved) n_moves <- n_moves + 1L
    Q <- st$Q
  }
  list(Q = Q, n_moves = n_moves)
}

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
                              min_eps = 1e-4, max_eps = 2.0) {
  eps <- eps0
  for (rd in seq_len(rounds)) {
    acc <- mean(vapply(seq_len(n_probe), function(.)
      isTRUE(feasible(cayley_left_step(Q_obs, eps))), logical(1)))
    if (acc < target_lo)      eps <- max(eps * 0.5, min_eps)
    else if (acc > target_hi) eps <- min(eps * 1.5, max_eps)
    else break
    if (eps <= min_eps || eps >= max_eps) break
  }
  eps
}

#' Algorithm-agnostic Besag-Clifford selective Hotelling test (arbitrary clustering)
#'
#' @description
#' A Monte Carlo selective test of equality of the means of a discovered cluster
#' pair `(k1, k2)`, valid after clustering `X` with an **arbitrary** black-box
#' clustering routine `cl_fun`, treating the covariance \eqn{\Sigma} as unknown.
#' It is the package's entry point for selective inference after *any* clustering
#' algorithm (hierarchical clustering, fixed-init k-means, spectral clustering,
#' ...), in contrast to [selective_hotelling_test()], which is hard-wired to
#' k-means.
#'
#' The construction is exactly that of Besag and Clifford (1991), transported to
#' the Stiefel manifold. Under the null the observed polar factor \eqn{Q_{obs}}
#' is uniform on the feasible region
#' \eqn{E = \{Q \in St(r-1, p) : cl\_fun(X(Q))\ \mathrm{preserves}\ (k_1,k_2)\}}
#' regardless of \eqn{\Sigma}. We build a Metropolis kernel \eqn{K} on \eqn{E}
#' using a **symmetric** Cayley proposal \eqn{Q' = U Q} (with \eqn{U} a random
#' Cayley rotation of step `eps`) that is accepted iff \eqn{Q' \in E}; for a flat
#' target restricted to \eqn{E} this acceptance ratio is the indicator
#' \eqn{1\{Q' \in E\}}, so \eqn{K} (and \eqn{K^L}) is reversible w.r.t.
#' \eqn{\mathrm{Uniform}(E)}. We then plant \eqn{Q_{obs}} at a uniform random
#' index \eqn{J \in \{0,\dots,M\}} and run \eqn{K^L} forward \eqn{J \to M} and
#' backward \eqn{J \to 0}. Because the seed is stationary and \eqn{K^L}
#' reversible, the rank of the contrast statistic \eqn{u(Q) = \|Q'v\|^2} (which
#' is strictly monotone in the pair's two-sample Hotelling \eqn{T^2}) of
#' \eqn{Q_{obs}} among the \eqn{M+1} path marks is exactly
#' \eqn{\mathrm{Uniform}\{1,\dots,M+1\}}. Hence the randomized upper-tail rank
#' p-value is exactly \eqn{\mathrm{Uniform}(0,1)} under \eqn{H_0}.
#'
#' **Validity holds for ANY reversible kernel, ANY step `eps`, ANY path length
#' `M`, and ANY clustering `cl_fun`.** Mixing (governed by `eps`, `L`, and the
#' acceptance rate) buys *power* only, never validity.
#'
#' @param X A numeric \eqn{n \times p} data matrix with observations in rows.
#' @param cl_fun A black-box clustering routine: a function
#'   `function(Xmat) -> integer label vector` (or a list carrying `$cluster`,
#'   e.g. the result of [stats::kmeans()]). It **must be a pure function of its
#'   data argument** (only data-independent randomness is allowed, e.g. a fixed
#'   initialization or a deterministic algorithm such as hierarchical
#'   clustering). It is called once on `X` to obtain the observed clustering and
#'   many times on reconstructed data \eqn{X(Q)} during sampling.
#' @param k1,k2 The two tested cluster labels, in the labelling of `cl_fun(X)`.
#'   Each cluster must have at least three members for a Hotelling \eqn{T^2}.
#' @param M Path length; the path carries `M + 1` marks. Validity holds for any
#'   `M >= 1`; larger `M` reduces Monte Carlo error.
#' @param L Number of kernel steps per path edge (i.e. \eqn{K^L}). Larger `L`
#'   improves spread/mixing and hence power.
#' @param eps Cayley step size. If `NULL` and `calibrate = TRUE` it is
#'   auto-tuned toward a workable acceptance rate; if `NULL` and
#'   `calibrate = FALSE` it defaults to `0.3`.
#' @param randomized If `TRUE` (default) returns the randomized upper-tail rank
#'   p-value, which is exactly \eqn{\mathrm{Uniform}(0,1)} under \eqn{H_0}; if
#'   `FALSE`, ties are broken conservatively.
#' @param calibrate If `TRUE` (default) and `eps` is `NULL`, auto-tune the
#'   Cayley step size for a workable acceptance rate.
#' @param cl_obs Optional precomputed `cl_fun(X)` (saves one call); otherwise it
#'   is computed internally.
#' @param ridge Small non-negative ridge added to the pooled covariance when
#'   computing the reported observed \eqn{T^2} via [hotelling_T2()].
#' @param stiefel_tol Guard tolerance on \eqn{\max|Q'Q - I|} over the path
#'   (GUARD a).
#' @param guard_feas_checks Number of path points to re-verify with the generic
#'   feasibility predicate (GUARD b).
#' @param oracle_Sigma Optional known \eqn{p \times p} covariance matrix. If
#'   supplied, an additional known-\eqn{\Sigma} oracle p-value is computed on the
#'   same path by ranking the Mahalanobis statistic
#'   \eqn{d' \Sigma^{-1} d} of the pair-mean difference; returned as
#'   `p_oracle`. Validity is identical; only the statistic (hence power) differs.
#' @param verbose If `TRUE`, print a one-block diagnostic.
#'
#' @return An S3 object of class `"selectiveHotelling"`: a list whose primary
#'   fields are `pval` (the randomized Besag-Clifford p-value), `se` (its Monte
#'   Carlo standard error), `T_obs` (the observed two-sample Hotelling
#'   \eqn{T^2} of the tested pair), `draws` (the vector of path statistic values
#'   used for ranking), `acc` (kernel move fraction), `n_unique_draws`, `k1`,
#'   `k2`, `K` (number of distinct labels in the observed clustering),
#'   `condition_on` (`"pair-preserved"`), and `sampler` (`"besag-clifford"`).
#'   Additional diagnostic fields are also returned: `p_cons` (conservative
#'   tie-broken p-value), `p_oracle`, `M`, `L`, `eps`, `u_obs`, `path_u`, `J`,
#'   `move_frac`, `stiefel_max`, `feas_ok`, and `n_feasible_path`.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(30 * 2), 30, 2)
#' # A deterministic, data-pure clustering routine (hierarchical clustering).
#' hclust_cl <- function(Xmat) {
#'   hc <- stats::hclust(stats::dist(Xmat), method = "average")
#'   as.integer(stats::cutree(hc, k = 3))
#' }
#' cl <- hclust_cl(X)
#' # Pick a pair with at least 3 members each.
#' sizes <- table(cl)
#' big <- as.integer(names(sizes)[sizes >= 3])
#' res <- selective_hotelling_bc(X, hclust_cl, big[1], big[2], M = 60)
#' res$pval
#' print(res)
#'
#' @seealso [selective_hotelling_test()], [hotelling_T2()]
#' @export
selective_hotelling_bc <- function(X, cl_fun, k1, k2,
                                   M = 300L, L = 2L, eps = NULL,
                                   randomized = TRUE, calibrate = TRUE,
                                   cl_obs = NULL, ridge = 1e-8,
                                   stiefel_tol = 1e-8, guard_feas_checks = 25L,
                                   oracle_Sigma = NULL, verbose = FALSE) {
  stopifnot(is.matrix(X), M >= 1L, L >= 1L)
  M <- as.integer(M); L <- as.integer(L)

  .as_labels <- function(z) if (is.list(z) && !is.null(z$cluster)) z$cluster else z
  if (is.null(cl_obs)) cl_obs <- .as_labels(cl_fun(X))
  cl_obs <- as.integer(cl_obs)
  if (length(cl_obs) != nrow(X)) stop("cl_fun(X) length != nrow(X).")
  if (sum(cl_obs == k1) < 3L || sum(cl_obs == k2) < 3L)
    stop("Tested clusters must each have >= 3 members for a Hotelling T^2.")

  setup    <- .bc_setup(X, cl_obs, k1, k2, cl_fun)
  Q_obs    <- setup$Q_obs
  feasible <- setup$feasible
  v        <- setup$v
  rdim <- nrow(Q_obs); pdim <- ncol(Q_obs); Ip <- diag(pdim)

  # ---- GUARD (c): the observed point must be feasible (it is, by construction:
  #      X(Q_obs) == X and cl_fun(X) == cl_obs preserves its own pair). Assert.
  if (!isTRUE(feasible(Q_obs)))
    stop("GUARD (c) FAILED: Q_obs is not feasible (cl_fun not reproducible on X(Q_obs)?).")

  # ---- step size ----
  if (is.null(eps)) {
    eps <- if (isTRUE(calibrate)) .bc_calibrate_eps(Q_obs, feasible) else 0.3
  }

  # ---- statistic u(Q) = ||Q' v||^2, monotone in the pair's Hotelling T^2 ----
  uval  <- function(Q) sum(crossprod(Q, v)^2)
  u_obs <- uval(Q_obs)

  stiefel_max <- max(abs(crossprod(Q_obs) - Ip))
  path_Q <- vector("list", M + 1L)

  # ---- RANDOM-INDEX PATH: plant Q_obs at J, run K^L forward & backward ----
  J <- sample.int(M + 1L, 1L) - 1L            # J ~ Uniform{0..M}
  path_u <- numeric(M + 1L)
  path_u[J + 1L] <- u_obs
  path_Q[[J + 1L]] <- Q_obs

  n_steps <- 0L; n_moves <- 0L
  record <- function(idx, Q) {
    path_u[idx + 1L] <<- uval(Q)
    path_Q[[idx + 1L]] <<- Q
    e <- max(abs(crossprod(Q) - Ip)); if (e > stiefel_max) stiefel_max <<- e
  }

  Q <- Q_obs
  if (J < M) for (j in (J + 1L):M) {
    kl <- .bc_kernel_L(Q, L, eps, feasible); Q <- kl$Q
    n_steps <- n_steps + L; n_moves <- n_moves + kl$n_moves
    record(j, Q)
  }
  Q <- Q_obs
  if (J > 0L) for (j in (J - 1L):0L) {
    kl <- .bc_kernel_L(Q, L, eps, feasible); Q <- kl$Q
    n_steps <- n_steps + L; n_moves <- n_moves + kl$n_moves
    record(j, Q)
  }
  move_frac <- if (n_steps > 0L) n_moves / n_steps else NA_real_

  # ---- GUARD (a): every path Q is on the Stiefel ----
  if (stiefel_max > stiefel_tol)
    stop(sprintf("GUARD (a) FAILED: path point off Stiefel, max|Q'Q-I| = %.3e > %.1e",
                 stiefel_max, stiefel_tol))

  # ---- GUARD (b): a random subset of path Q's passes the GENERIC feasibility
  #      predicate (cl_fun(X(Q)) preserves the pair). ----
  npath <- M + 1L
  feas_idx <- if (npath <= guard_feas_checks) seq_len(npath)
              else sort(sample.int(npath, guard_feas_checks))
  feas_chk <- vapply(feas_idx, function(ii) isTRUE(feasible(path_Q[[ii]])), logical(1))
  feas_ok  <- all(feas_chk)
  if (!feas_ok) {
    bad <- feas_idx[!feas_chk]
    stop(sprintf("GUARD (b) FAILED: %d/%d sampled path Q's not feasible. first bad path-index=%d",
                 sum(!feas_chk), length(feas_chk), bad[1] - 1L))
  }
  # Every recorded path point is an accepted/stay state in E, so all M+1 are
  # feasible by construction; the random subset above re-confirms it.
  n_feasible_path <- npath

  # ---- randomized upper-tail rank p-value + conservative variant ----
  n_gt <- sum(path_u >  u_obs)
  n_eq <- sum(path_u == u_obs)
  n_ge <- sum(path_u >= u_obs)
  U <- if (isTRUE(randomized)) stats::runif(1L) else 1
  p_bc   <- (n_gt + U * n_eq) / (M + 1L)
  p_cons <- n_ge / (M + 1L)

  # ---- optional KNOWN-Sigma oracle p-value on the SAME path: rank the
  #      Mahalanobis statistic d' Sigma^{-1} d of the pair-mean difference, with
  #      d recomputed from X(Q) at each path point. Same J/plant, independent
  #      tie-break. Validity is identical; only the statistic (hence power)
  #      differs. ----
  p_oracle <- NA_real_
  if (!is.null(oracle_Sigma)) {
    SigInv <- tryCatch(solve(oracle_Sigma), error = function(e) NULL)
    if (!is.null(SigInv)) {
      i1 <- setup$prep$idx1; i2 <- setup$prep$idx2
      uo <- function(Q) { XQ <- setup$X_of_Q(Q)
        d <- colMeans(XQ[i1, , drop = FALSE]) - colMeans(XQ[i2, , drop = FALSE])
        as.numeric(crossprod(d, SigInv %*% d)) }
      path_uo <- vapply(path_Q, uo, numeric(1))
      uo_obs  <- path_uo[J + 1L]
      Uo <- if (isTRUE(randomized)) stats::runif(1L) else 1
      p_oracle <- (sum(path_uo > uo_obs) + Uo * sum(path_uo == uo_obs)) / (M + 1L)
    }
  }

  # ---- observed two-sample Hotelling T^2 (for interpretability) ----
  T_obs <- hotelling_T2(X, cl_obs, k1, k2, ridge)

  if (verbose) {
    cat(sprintf("[bc] n=%d p=%d  Q in St(%d,%d)  M=%d L=%d eps=%.4g\n",
                nrow(X), pdim, rdim, pdim, M, L, eps))
    cat(sprintf("  J=%d u_obs=%.5g  #>=%d #>%d #==%d\n", J, u_obs, n_ge, n_gt, n_eq))
    cat(sprintf("  p_bc=%.4f p_cons=%.4f  acc(move_frac)=%.3f  steps=%d\n",
                p_bc, p_cons, move_frac, n_steps))
    cat(sprintf("  GUARDS: stiefel max|Q'Q-I|=%.2e  feas=%d/%d (%s)\n",
                stiefel_max, sum(feas_chk), length(feas_chk),
                if (feas_ok) "OK" else "FAIL"))
  }

  out <- list(
    pval = p_bc,
    se = sqrt(p_bc * (1 - p_bc) / (M + 1L)),
    T_obs = T_obs,
    draws = path_u,
    acc = move_frac,
    n_unique_draws = length(unique(path_u)),
    k1 = k1,
    k2 = k2,
    K = length(unique(cl_obs)),
    condition_on = "pair-preserved",
    sampler = "besag-clifford",
    # richer diagnostics
    p_cons = p_cons, p_oracle = p_oracle, M = M, L = L, eps = eps,
    u_obs = u_obs, path_u = path_u, J = J, move_frac = move_frac,
    stiefel_max = stiefel_max, feas_ok = feas_ok,
    n_feasible_path = n_feasible_path
  )
  class(out) <- "selectiveHotelling"
  out
}
