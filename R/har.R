# =============================================================================
# har.R -- ANALYTIC-ARC uniform-Q hit-and-run selective Hotelling test (k-means).
#
# This is the CLOSED-FORM, k-means specialization of selective_hotelling_bc: same
# Besag-Clifford random-index validity wrapper and same full Hotelling T^2
# statistic, but the kernel moves are EXACT. Instead of the black-box
# Cayley-Metropolis accept/reject of selective_hotelling_bc (which re-runs the
# clustering at every proposal), each step here draws a two-row Givens plane and
# samples the rotation angle theta UNIFORMLY BY ARC LENGTH inside the CLOSED-FORM
# feasible arcs of the fixed-init Lloyd TRAJECTORY event E_tau -- with NO
# re-clustering in the inner loop. Along a Givens geodesic Q(theta) the
# reconstructed data is a harmonic curve, so each k-means assignment constraint is
# a harmonic A cos2t + B sin2t + C cos t + D sin t + E whose feasible set is a
# union of arcs (harmonic_arc.R). This is the package port of the verified sim
# engine (sim/path_har.R + sim/analytic_arc.R + sim/path_bc.R); it reuses the
# GLOBAL-FREE arc machinery of harmonic_arc.R verbatim.
#
# WHY THIS IS THE ANALYTIC-FIRST k-means DEFAULT.  The move is rejection-free and
# exact (uniform on E_tau restricted to the geodesic), so it mixes far better than
# the demoted black-box MCMC samplers (samplers.R) and needs no step-size tuning.
# The Besag-Clifford wrapper gives finite-sample exactness (Uniform(0,1) under H0)
# regardless of mixing; mixing only buys power. See selective_hotelling_bc for the
# general (any-clustering) black-box counterpart, which is what you fall back to
# when the clustering rule has no closed-form geodesic constraints.
# =============================================================================

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

# -----------------------------------------------------------------------------
# .har_step(Q, ...): ONE analytic-arc hit-and-run step. Pick two random rows
# (i, j), build the closed-form feasible arcs of E_tau on that Givens geodesic,
# and draw theta ~ Uniform-by-arc-length over the union. theta = 0 (the current
# Q) is always feasible, so the union is non-empty; the move stays in E_tau by
# construction. Returns the moved Q. Reversible w.r.t. Uniform(E_tau).
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.har_step <- function(Q, X_of_Q, X_fixed, init_idx, traj, K) {
  rq <- nrow(Q)
  ij <- sample.int(rq, 2L); i <- ij[1]; j <- ij[2]
  puw  <- .har_givens_PUW(X_of_Q, X_fixed, Q, i, j)
  co   <- .kmeans_traj_coeffs(puw$P, puw$U, puw$W, init_idx, traj, K)
  arcs <- .harmonic_feasible_arcs(co, lo = 0, hi = 2 * pi)
  if (!nrow(arcs)) return(Q)                 # (should not happen: t=0 is feasible)
  len <- arcs[, 2] - arcs[, 1]; total <- sum(len)
  if (total <= 0) return(Q)
  u <- stats::runif(1, 0, total); cum <- cumsum(len)
  wa <- which(u <= cum)[1]
  prev <- if (wa == 1L) 0 else cum[wa - 1L]
  theta <- arcs[wa, 1] + (u - prev)
  givens_left(Q, i, j, theta)
}

# .har_kernel_L(Q, L, ...): K^L (L analytic-arc HAR steps).
#' @keywords internal
#' @noRd
.har_kernel_L <- function(Q, L, X_of_Q, X_fixed, init_idx, traj, K) {
  for (s in seq_len(L)) Q <- .har_step(Q, X_of_Q, X_fixed, init_idx, traj, K)
  Q
}

#' Analytic-arc uniform-Q hit-and-run selective Hotelling test (k-means, unknown covariance)
#'
#' @description
#' An **exact** (finite-sample Uniform(0,1) under \eqn{H_0}) selective test of
#' equality of the means of a discovered \eqn{k}-means cluster pair `(k1, k2)`,
#' with the covariance \eqn{\Sigma} treated as **unknown**. It is the closed-form,
#' analytic-arc specialization of [selective_hotelling_bc()]: same Besag--Clifford
#' random-index validity and same full Hotelling \eqn{T^2} statistic, but the
#' Stiefel-manifold moves are **rejection-free and exact**. Conditioning is on the
#' fixed-init Lloyd **trajectory** \eqn{E_\tau}.
#'
#' @details
#' A two-row Givens geodesic \eqn{Q(\theta)} makes the reconstructed data a
#' harmonic curve, so every \eqn{k}-means assignment constraint of \eqn{E_\tau} is
#' a harmonic \eqn{A\cos 2\theta + B\sin 2\theta + C\cos\theta + D\sin\theta + E}
#' whose feasible set is a closed-form union of arcs (see the package engine in
#' `harmonic_arc.R`). Each kernel step draws a random plane and samples \eqn{\theta}
#' uniformly by arc length inside those arcs -- **no re-clustering in the inner
#' loop** -- giving an exact \eqn{\mathrm{Uniform}(E_\tau)}-reversible move. The
#' Besag--Clifford wrapper then makes the rank \eqn{p}-value exactly
#' \eqn{\mathrm{Uniform}(0,1)} under \eqn{H_0} for any path length `M`, kernel
#' repeat `L`, and mixing; mixing only affects power. For an arbitrary
#' (non-\eqn{k}-means) clustering rule with no closed-form geodesic constraints,
#' use the black-box [selective_hotelling_bc()] instead.
#'
#' @param X A numeric \eqn{n \times p} data matrix with observations in rows.
#' @param k1,k2 The two tested cluster labels (in the labelling of the fixed-init
#'   \eqn{k}-means run). Each cluster must have at least three members.
#' @param K Number of clusters.
#' @param init_centers_idx Fixed initial-centre row indices; drawn via
#'   [kmeans_fun()] if `NULL`. Conditioning is on that exact fixed-init run's
#'   Lloyd trajectory.
#' @param M Besag--Clifford path length (the path carries `M + 1` marks). Validity
#'   holds for any `M >= 1`; larger `M` reduces Monte Carlo error.
#' @param L Number of analytic-arc HAR steps per path edge (\eqn{K^L}). Larger `L`
#'   improves mixing and hence power.
#' @param randomized If `TRUE` (default) returns the randomized upper-tail rank
#'   \eqn{p}-value (exactly \eqn{\mathrm{Uniform}(0,1)} under \eqn{H_0}); if
#'   `FALSE`, ties are broken conservatively.
#' @param ridge Small non-negative ridge for the reported observed \eqn{T^2}.
#' @param verbose If `TRUE`, print a one-line diagnostic.
#'
#' @return An S3 object of class `"selectiveHotelling"` with primary fields `pval`,
#'   `se`, `T_obs`, `k1`, `k2`, `K`, `condition_on` (`"kmeans-trajectory"`), and
#'   `sampler` (`"har-analytic"`), plus diagnostics `p_cons`, `M`, `L`, `J`,
#'   `u_obs`, `path_u`, `draws`, `n_unique_draws`, and `stiefel_max`.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40 * 3), 40, 3)
#' km <- kmeans_fun(X, k = 3)
#' sizes <- table(km$cluster)
#' big <- as.integer(names(sizes)[sizes >= 3])
#' res <- selective_hotelling_har(X, big[1], big[2], K = 3,
#'                                init_centers_idx = km$init_centers_idx, M = 100)
#' res$pval
#'
#' @seealso [selective_hotelling_bc()], [selective_hotelling_test()],
#'   [selective_truncated_F()]
#' @export
selective_hotelling_har <- function(X, k1, k2, K = 3L, init_centers_idx = NULL,
                                    M = 300L, L = 2L, randomized = TRUE,
                                    ridge = 1e-8, verbose = FALSE) {
  stopifnot(is.matrix(X), M >= 1L, L >= 1L)
  M <- as.integer(M); L <- as.integer(L)

  if (is.null(init_centers_idx)) init_centers_idx <- kmeans_fun(X, k = K)$init_centers_idx
  init_centers_idx <- as.integer(init_centers_idx)
  ref    <- kmeans_with_fixed_init(X, init_centers_idx)
  cl_obs <- ref$cluster; traj <- ref$trajectory
  if (sum(cl_obs == k1) < 3L || sum(cl_obs == k2) < 3L)
    stop("Tested clusters must each have >= 3 members for a Hotelling T^2.")

  prep <- prepare_pair_objects(X, cl_obs, k1, k2)
  m <- prep$m; r <- prep$r
  X_fixed <- prep$X_fixed
  X_of_Q <- function(Q) {
    AQS <- prep$A %*% Q %*% prep$S_sqrt
    XQ <- prep$X_fixed
    XQ[prep$idx1, ] <- XQ[prep$idx1, ] + AQS[1:m, , drop = FALSE]
    XQ[prep$idx2, ] <- XQ[prep$idx2, ] + AQS[(m + 1):r, , drop = FALSE]
    XQ
  }
  Q_obs <- prep$Q_obs

  # Contrast statistic u(Q) = ||Q' v||^2, strictly monotone in the pair's T^2.
  n1 <- m; n2 <- r - m
  contrast <- c(rep(1 / n1, n1), rep(-1 / n2, n2))
  v <- as.numeric(crossprod(prep$A, contrast))
  uval  <- function(Q) sum(crossprod(Q, v)^2)
  u_obs <- uval(Q_obs)

  Ip <- diag(ncol(Q_obs))
  stiefel_max <- max(abs(crossprod(Q_obs) - Ip))

  # ---- Besag-Clifford random-index path with the analytic-arc HAR kernel ----
  J <- sample.int(M + 1L, 1L) - 1L
  path_u <- numeric(M + 1L); path_u[J + 1L] <- u_obs
  kern <- function(Q) .har_kernel_L(Q, L, X_of_Q, X_fixed, init_centers_idx, traj, K)
  bump <- function(Q) { e <- max(abs(crossprod(Q) - Ip)); if (e > stiefel_max) stiefel_max <<- e }

  Q <- Q_obs
  if (J < M) for (jj in (J + 1L):M) { Q <- kern(Q); path_u[jj + 1L] <- uval(Q); bump(Q) }
  Q <- Q_obs
  if (J > 0L) for (jj in (J - 1L):0L) { Q <- kern(Q); path_u[jj + 1L] <- uval(Q); bump(Q) }

  if (stiefel_max > 1e-6)
    stop(sprintf("GUARD FAILED: path point off Stiefel, max|Q'Q-I| = %.3e", stiefel_max))

  n_gt <- sum(path_u >  u_obs)
  n_eq <- sum(path_u == u_obs)
  n_ge <- sum(path_u >= u_obs)
  U <- if (isTRUE(randomized)) stats::runif(1L) else 1
  p_bc   <- (n_gt + U * n_eq) / (M + 1L)
  p_cons <- n_ge / (M + 1L)

  T_obs <- hotelling_T2(X, cl_obs, k1, k2, ridge)

  if (verbose)
    cat(sprintf("[har-analytic] n=%d p=%d M=%d L=%d  J=%d u_obs=%.5g  p=%.4f\n",
                nrow(X), ncol(X), M, L, J, u_obs, p_bc))

  out <- list(
    pval = p_bc,
    se = sqrt(p_bc * (1 - p_bc) / (M + 1L)),
    T_obs = T_obs,
    draws = path_u,
    acc = NA_real_,
    n_unique_draws = length(unique(path_u)),
    k1 = k1, k2 = k2, K = length(unique(cl_obs)),
    condition_on = "kmeans-trajectory",
    sampler = "har-analytic",
    p_cons = p_cons, M = M, L = L, u_obs = u_obs, path_u = path_u, J = J,
    stiefel_max = stiefel_max
  )
  class(out) <- "selectiveHotelling"
  out
}

# =============================================================================
# HIERARCHICAL analytic-arc HAR (T6): the hclust twin of selective_hotelling_har.
# Same Besag-Clifford wrapper, same full T^2, same two-row Givens geodesic HAR --
# but the feasible arcs come from .hclust_pairpreserved_arcs (the cutree-K
# PAIR-PRESERVED event via event-driven cell traversal), so this conditions on
# the hclust CUT event (what selective_hotelling_bc and Gao et al. use), not a
# trajectory. Closed-form breakpoints are exact for the SINGLE-HARMONIC linkages
# (average, centroid); other linkages have no closed-form arcs -> use
# selective_hotelling_bc (black-box) instead.
# =============================================================================

# ---- one hierarchical analytic-arc HAR step (pair-preserved cut event) --------
#' @keywords internal
#' @noRd
.hclust_har_step <- function(Q, X_of_Q, X_fixed, K, cl_obs, k1, k2,
                             method, square) {
  rq <- nrow(Q)
  ij <- sample.int(rq, 2L); i <- ij[1]; j <- ij[2]
  puw    <- .har_givens_PUW(X_of_Q, X_fixed, Q, i, j)
  Qtheta <- function(th) givens_left(Q, i, j, th)
  arcs <- .hclust_pairpreserved_arcs(puw$P, puw$U, puw$W, X_of_Q, Qtheta, K,
                                     cl_obs, k1, k2, method = method,
                                     square = square, lo = 0, hi = 2 * pi)
  if (!nrow(arcs)) return(Q)
  len <- arcs[, 2] - arcs[, 1]; total <- sum(len)
  if (total <= 0) return(Q)
  u <- stats::runif(1, 0, total); cum <- cumsum(len)
  wa <- which(u <= cum)[1]
  prev <- if (wa == 1L) 0 else cum[wa - 1L]
  theta <- arcs[wa, 1] + (u - prev)
  givens_left(Q, i, j, theta)
}

#' Analytic-arc uniform-Q hit-and-run selective Hotelling test (hierarchical clustering)
#'
#' @description
#' The hierarchical-clustering counterpart of [selective_hotelling_har()]: an exact
#' (finite-sample \eqn{\mathrm{Uniform}(0,1)} under \eqn{H_0}) selective test of a
#' discovered cluster pair `(k1, k2)`, with \eqn{\Sigma} unknown, using the same
#' Besag--Clifford wrapper and full Hotelling \eqn{T^2} but a **closed-form**
#' analytic-arc hit-and-run kernel. Conditioning is on the **cutree-\eqn{K}
#' pair-preserved** event (the same event [selective_hotelling_bc()] uses), whose
#' feasible arcs are computed in closed form per Givens plane via event-driven
#' cell traversal.
#'
#' Closed-form arcs are exact only for the **single-harmonic** linkages
#' (`"average"`, `"centroid"`); for any other linkage there is no closed-form arc,
#' so use the black-box [selective_hotelling_bc()] instead (this function will
#' error via its arc engine).
#'
#' @param X A numeric \eqn{n \times p} data matrix with observations in rows.
#' @param k1,k2 The two tested cluster labels (in the labelling of
#'   `hclust_fun(X, k = K, method = hclust_method, square = hclust_square)`). Each
#'   cluster must have at least three members.
#' @param K Number of clusters to cut the dendrogram into.
#' @param hclust_method Agglomeration method; closed-form arcs require `"average"`
#'   or `"centroid"`.
#' @param hclust_square If `TRUE` (default) the tree is built on squared Euclidean
#'   distances (the `clusterpval` convention; also the convention centroid linkage
#'   requires).
#' @param M,L Besag--Clifford path length and analytic-arc HAR steps per edge.
#' @param randomized If `TRUE` (default) the randomized upper-tail rank p-value.
#' @param ridge Small non-negative ridge for the reported observed \eqn{T^2}.
#' @param cl_obs Optional precomputed observed clustering (saves one call).
#' @param verbose If `TRUE`, print a one-line diagnostic.
#'
#' @return An S3 `"selectiveHotelling"` object with `pval`, `se`, `T_obs`, `k1`,
#'   `k2`, `K`, `condition_on = "pair-preserved"`, `sampler = "har-analytic"`, and
#'   the same diagnostics as [selective_hotelling_har()].
#'
#' @seealso [selective_hotelling_har()], [selective_hotelling_bc()],
#'   [selective_truncated_F()]
#' @export
selective_hclust_har <- function(X, k1, k2, K = 3L,
                                 hclust_method = "average", hclust_square = TRUE,
                                 M = 300L, L = 2L, randomized = TRUE,
                                 ridge = 1e-8, cl_obs = NULL, verbose = FALSE) {
  stopifnot(is.matrix(X), M >= 1L, L >= 1L)
  M <- as.integer(M); L <- as.integer(L)
  if (!(hclust_method %in% c("average", "centroid")))
    stop(sprintf(paste0("closed-form analytic-arc HAR needs a single-harmonic linkage ",
                        "(average, centroid); got '%s'. Use selective_hotelling_bc() ",
                        "for this linkage."), hclust_method))

  if (is.null(cl_obs))
    cl_obs <- hclust_fun(X, k = K, method = hclust_method, square = hclust_square)$cluster
  cl_obs <- as.integer(cl_obs)
  if (sum(cl_obs == k1) < 3L || sum(cl_obs == k2) < 3L)
    stop("Tested clusters must each have >= 3 members for a Hotelling T^2.")

  prep <- prepare_pair_objects(X, cl_obs, k1, k2)
  m <- prep$m; r <- prep$r
  X_fixed <- prep$X_fixed
  X_of_Q <- function(Q) {
    AQS <- prep$A %*% Q %*% prep$S_sqrt
    XQ <- prep$X_fixed
    XQ[prep$idx1, ] <- XQ[prep$idx1, ] + AQS[1:m, , drop = FALSE]
    XQ[prep$idx2, ] <- XQ[prep$idx2, ] + AQS[(m + 1):r, , drop = FALSE]
    XQ
  }
  Q_obs <- prep$Q_obs
  n1 <- m; n2 <- r - m
  contrast <- c(rep(1 / n1, n1), rep(-1 / n2, n2))
  v <- as.numeric(crossprod(prep$A, contrast))
  uval  <- function(Q) sum(crossprod(Q, v)^2)
  u_obs <- uval(Q_obs)

  Ip <- diag(ncol(Q_obs)); stiefel_max <- max(abs(crossprod(Q_obs) - Ip))

  J <- sample.int(M + 1L, 1L) - 1L
  path_u <- numeric(M + 1L); path_u[J + 1L] <- u_obs
  step <- function(Q) .hclust_har_step(Q, X_of_Q, X_fixed, K, cl_obs, k1, k2,
                                       hclust_method, hclust_square)
  kern <- function(Q) { for (s in seq_len(L)) Q <- step(Q); Q }
  bump <- function(Q) { e <- max(abs(crossprod(Q) - Ip)); if (e > stiefel_max) stiefel_max <<- e }

  Q <- Q_obs
  if (J < M) for (jj in (J + 1L):M) { Q <- kern(Q); path_u[jj + 1L] <- uval(Q); bump(Q) }
  Q <- Q_obs
  if (J > 0L) for (jj in (J - 1L):0L) { Q <- kern(Q); path_u[jj + 1L] <- uval(Q); bump(Q) }

  if (stiefel_max > 1e-6)
    stop(sprintf("GUARD FAILED: path point off Stiefel, max|Q'Q-I| = %.3e", stiefel_max))

  n_gt <- sum(path_u >  u_obs); n_eq <- sum(path_u == u_obs); n_ge <- sum(path_u >= u_obs)
  U <- if (isTRUE(randomized)) stats::runif(1L) else 1
  p_bc   <- (n_gt + U * n_eq) / (M + 1L)
  p_cons <- n_ge / (M + 1L)
  T_obs  <- hotelling_T2(X, cl_obs, k1, k2, ridge)

  if (verbose)
    cat(sprintf("[hclust har-analytic] n=%d p=%d M=%d L=%d method=%s  J=%d p=%.4f\n",
                nrow(X), ncol(X), M, L, hclust_method, J, p_bc))

  out <- list(
    pval = p_bc, se = sqrt(p_bc * (1 - p_bc) / (M + 1L)), T_obs = T_obs,
    draws = path_u, acc = NA_real_, n_unique_draws = length(unique(path_u)),
    k1 = k1, k2 = k2, K = length(unique(cl_obs)),
    condition_on = "pair-preserved", sampler = "har-analytic",
    p_cons = p_cons, M = M, L = L, u_obs = u_obs, path_u = path_u, J = J,
    stiefel_max = stiefel_max
  )
  class(out) <- "selectiveHotelling"
  out
}
