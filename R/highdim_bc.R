# =============================================================================
# highdim_bc.R -- selective two-sample mean tests on the Besag-Clifford path,
# with statistics that keep power when p >= |C1| + |C2|.
#
# PARAMETRISATION.  Let A be an r x (r-1) orthonormal basis of 1^perp and
# Z = A' X_sub the (r-1) x p centred coordinates of the tested rows. Under H0'
# (all tested rows share a mean) Z ~ MN(0, I_{r-1}, Sigma), so O Z has the law of
# Z for every O in O(r-1), for every Sigma. The statistic G = Z'Z (the centred
# total scatter, p x p) is invariant, and when rank(Z) = r - 1 (any p >= r - 1)
# or rank(Z) = p (p <= r - 1) its level set is the single orbit O(r-1) Z_obs.
# Haar averaging then gives: given G, the pair mean, the other rows and the
# selection event, Z is uniform on {O Z_obs : O in O(r-1)} restricted to E.
# For p <= r - 2 this is the polar statement Z = Q S^{1/2} with Q uniform on
# St(r-1, p); left Givens moves on Z and on Q coincide. Nothing here needs p < r.
#
# STATISTICS.  d(Z) = Z'v is the difference of the two pair means. G is fixed
# along the path, so any weight matrix built from G is a constant:
#   l2       ||d||^2
#   diag     sum_j d_j^2 / (S0_jj + floor),             S0 = G / (r - 1)
#   max      max_j d_j^2 / (S0_jj + floor)
#   ridge    d'(S0 + lambda I)^{-1} d,                   lambda = tr(S0) / p
#   hotelling d'G^{-1} d  (only for p <= r - 2; for p >= r - 1 it is constant)
#   t2       regular Hotelling T^2, pooled within-pair covariance + ridge 1e-8
#            (ranks like hotelling for p <= r - 2; not constant for p >= r - 1)
#   cq       Chen-Qin ||d||^2 - tr(S1)/n1 - tr(S2)/n2 (S1, S2 vary on the path)
#   diag_w   sum_j d_j^2 / S_w,jj, within-pair pooled variances (vary on the path)
#   proj     random-projection Hotelling, averaged over fixed projections
#   adapt    min-p of the components in `adapt`, calibrated on the same path
# All are multiplied by n_eff = n1 n2 / (n1 + n2) for readability. The ridge and
# Hotelling statistics with the WITHIN-pair pooled covariance are monotone in
# the total-scatter versions (Sherman-Morrison: S_w = (G - n_eff d d')/(r-2)), so
# they are not listed separately.
#
# VALIDITY.  The Besag-Clifford construction makes the index J of the observed
# point uniform given the (unordered) path. Any statistic that depends on a
# point and on the multiset of path points -- in particular each fixed-weight
# statistic above and the min-p combination -- gives an exactly valid rank
# p-value, for any M, L and reversible kernel.
# =============================================================================

#' Centred coordinates of the tested rows and the reconstruction map
#'
#' @keywords internal
#' @noRd
.hd_setup <- function(X, cl_obs, k1, k2, cl_fun) {
  idx1 <- which(cl_obs == k1); idx2 <- which(cl_obs == k2)
  n1 <- length(idx1); n2 <- length(idx2); r <- n1 + n2
  X_sub <- rbind(X[idx1, , drop = FALSE], X[idx2, , drop = FALSE])

  P_perp <- diag(r) - matrix(1 / r, r, r)
  es <- eigen(P_perp, symmetric = TRUE)
  A  <- es$vectors[, 1:(r - 1), drop = FALSE]           # orthonormal basis of 1^perp
  Z_obs <- crossprod(A, X_sub)                           # (r-1) x p

  xbar <- colMeans(X_sub)
  X_fixed <- X
  X_fixed[idx1, ] <- rep(xbar, each = n1)
  X_fixed[idx2, ] <- rep(xbar, each = n2)
  X_of_Z <- function(Z) {
    AZ <- A %*% Z
    XZ <- X_fixed
    XZ[idx1, ] <- XZ[idx1, ] + AZ[1:n1, , drop = FALSE]
    XZ[idx2, ] <- XZ[idx2, ] + AZ[(n1 + 1):r, , drop = FALSE]
    XZ
  }

  contrast <- c(rep(1 / n1, n1), rep(-1 / n2, n2))
  v <- as.numeric(crossprod(A, contrast))                # d(Z) = Z'v

  .as_labels <- function(z) if (is.list(z) && !is.null(z$cluster)) z$cluster else z
  feasible <- function(Z) {
    cl_new <- tryCatch(.as_labels(cl_fun(X_of_Z(Z))), error = function(e) NULL)
    if (is.null(cl_new) || length(cl_new) != length(cl_obs)) return(FALSE)
    .bc_pair_preserved(cl_obs, cl_new, k1, k2)
  }

  list(idx1 = idx1, idx2 = idx2, n1 = n1, n2 = n2, r = r, A = A,
       Z_obs = Z_obs, X_of_Z = X_of_Z, v = v, feasible = feasible)
}

#' Statistic evaluator: Z -> named vector of the requested statistics
#'
#' @keywords internal
#' @noRd
.hd_stat_fun <- function(Z_obs, v, n_eff, stats, A = NULL, n1 = NULL, var_floor = 0,
                         ridge_lambda = NULL, oracle_Sigma = NULL, proj = NULL,
                         t2_ridge = 1e-8) {
  rm1 <- nrow(Z_obs); p <- ncol(Z_obs); r <- rm1 + 1L
  g   <- colSums(Z_obs^2)                                # diag(G), fixed on the path
  s0  <- g / rm1                                         # diag(S0)
  w   <- 1 / (s0 + var_floor)
  sv  <- svd(Z_obs, nu = 0)
  keep <- sv$d > max(sv$d) * 1e-10
  V   <- sv$v[, keep, drop = FALSE]; d2 <- sv$d[keep]^2
  if (is.null(ridge_lambda)) ridge_lambda <- mean(s0)
  w_ridge <- 1 / (d2 / rm1 + ridge_lambda)               # (S0 + lambda I)^{-1} on span(V)
  w_hot   <- 1 / d2                                      # G^{-1} on span(V)
  hot_ok  <- p <= rm1 - 1L
  # Regular T^2 with the pooled within-pair covariance and ridge t2_ridge (as
  # hotelling_T2()): S_w = (G - n_eff d d') / (r - 2), so by Sherman-Morrison
  # T = n_eff q / (1 - c q), q = d'(G/(r-2) + eps I)^{-1} d, c = n_eff / (r - 2).
  # 1 - c q = (1 - U) + sum_k u_k e_k / (1 + e_k), u_k = n_eff a_k^2 / D_k^2,
  # e_k = eps (r - 2) / D_k^2, U = sum u_k; U = 1 exactly when rank(G) = r - 1.
  e_t2 <- t2_ridge * (r - 2) / d2
  SigInv  <- if (!is.null(oracle_Sigma)) solve(oracle_Sigma) else NULL
  # Projected Hotelling: (P'd)'(P'GP)^{-1}(P'd) for each fixed projection P,
  # averaged. Monotone in the within-pair version (Sherman-Morrison) for one P.
  proj_inv <- if (!is.null(proj)) lapply(proj, function(P) {
    ZP <- Z_obs %*% P; solve(crossprod(ZP)) }) else NULL
  n2 <- if (!is.null(n1)) r - n1 else NULL

  function(Z) {
    d  <- as.numeric(crossprod(Z, v))
    out <- numeric(length(stats)); names(out) <- stats
    for (s in stats) out[s] <- switch(s,
      l2        = sum(d^2),
      diag      = sum(w * d^2),
      max       = max(w * d^2),
      ridge     = { a <- as.numeric(crossprod(V, d)); sum(w_ridge * a^2) },
      hotelling = if (hot_ok) { a <- as.numeric(crossprod(V, d)); sum(w_hot * a^2) } else NA_real_,
      t2        = {
        a <- as.numeric(crossprod(V, d)); u <- n_eff * a^2 / d2
        one_m_U <- if (hot_ok) 1 - sum(u) else 0
        q <- sum(a^2 / (d2 / (r - 2) + t2_ridge))
        q / (one_m_U + sum(u * e_t2 / (1 + e_t2))) },
      oracle    = if (!is.null(SigInv)) as.numeric(crossprod(d, SigInv %*% d)) else NA_real_,
      # Chen-Qin: ||d||^2 - tr(S1)/n1 - tr(S2)/n2 (S_k varies along the path).
      cq        = {
        AZ <- A %*% Z
        W1 <- sum(AZ[1:n1, ]^2) - n1 * sum(colMeans(AZ[1:n1, , drop = FALSE])^2)
        W2 <- sum(AZ[(n1 + 1):r, ]^2) - n2 * sum(colMeans(AZ[(n1 + 1):r, , drop = FALSE])^2)
        (sum(d^2) - W1 / ((n1 - 1) * n1) - W2 / ((n2 - 1) * n2)) },
      # Diagonal statistic with the within-pair pooled variances,
      # S_w,jj = (G_jj - n_eff d_j^2) / (r - 2). (The max with these variances
      # ranks exactly like "max": x / (1 - c x) is increasing.)
      diag_w    = sum(d^2 / ((g - n_eff * d^2) / (r - 2) + var_floor)),
      proj      = mean(vapply(seq_along(proj), function(k) {
        a <- as.numeric(crossprod(proj[[k]], d)); sum(a * (proj_inv[[k]] %*% a)) }, numeric(1))),
      stop("unknown statistic: ", s))
    n_eff * out
  }
}

# Data-independent Gaussian projections, drawn under a local seed so the
# global RNG stream (hence the Besag-Clifford path) is unchanged.
.hd_draw_proj <- function(p, k, n_proj, seed) {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    saved <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  lapply(seq_len(n_proj), function(i) qr.Q(qr(matrix(stats::rnorm(p * k), p, k))))
}

# Randomized upper-tail rank p-value of position j0 among the column t.
.hd_rank_p <- function(t, j0, U) {
  if (anyNA(t)) return(NA_real_)
  (sum(t > t[j0]) + U * sum(t == t[j0])) / length(t)
}

#' Selective two-sample mean tests for a cluster pair, valid when p >= n1 + n2
#'
#' @description
#' Exact selective tests of equality of the means of a discovered cluster pair
#' for an arbitrary clustering routine `cl_fun`, with the covariance unknown and
#' **no restriction on the dimension**. The tested rows are reduced to their
#' centred coordinates \eqn{Z = A'X_{sub}}; under the null \eqn{Z} given its
#' total scatter \eqn{G = Z'Z} is uniform on the orbit \eqn{O(r-1) Z_{obs}}, for
#' every \eqn{\Sigma} and every \eqn{p}. The random Givens kernel of
#' [selective_hotelling_bc()] runs on this orbit restricted to the selection
#' event, and several statistics are ranked on the same Besag-Clifford path:
#' squared Euclidean distance (`"l2"`), a diagonal-standardised distance
#' (`"diag"`), the coordinatewise maximum (`"max"`), a ridge Hotelling
#' statistic (`"ridge"`), the Hotelling statistic (`"hotelling"`, only for
#' \eqn{p \le r - 2}) and a min-p combination (`"adapt"`). All weights are
#' built from \eqn{G}, which is constant on the path, so each p-value is exact
#' for any `M`, `L` and reversible kernel.
#'
#' @inheritParams selective_hotelling_bc
#' @param stats Statistics to rank. Any of `"l2"`, `"diag"`, `"max"`,
#'   `"ridge"`, `"hotelling"`, `"t2"` (the regular two-sample \eqn{T^2} with the
#'   pooled covariance plus `t2_ridge`, as [hotelling_T2()]), `"cq"` (Chen-Qin), `"diag_w"` (diagonal with
#'   the within-pair pooled variances), `"proj"` (random-projection Hotelling,
#'   averaged over `proj_n` projections), `"oracle"` (needs `oracle_Sigma`).
#' @param adapt Components of the min-p combination (reported as `"adapt"`), or
#'   a named list of such vectors, one p-value per name; `NULL` skips it. Each
#'   component must be in `stats`.
#' @param var_floor Non-negative floor added to the diagonal variances in
#'   `"diag"` and `"max"`.
#' @param proj_dim,proj_n,proj_seed Dimension, number and seed of the
#'   data-independent Gaussian projections for `"proj"`; `proj_dim` defaults to
#'   \eqn{\lfloor (r-2)/2 \rfloor}. The projections are drawn under a local seed,
#'   so the global random-number stream is unchanged.
#' @param t2_ridge Ridge added to the pooled covariance in `"t2"`.
#' @param ridge_lambda Ridge level for `"ridge"`; default \eqn{tr(S_0)/p}, a
#'   function of \eqn{G} only.
#'
#' @return A list with `pvals` (named p-values), `obs` (observed statistics),
#'   `path` (the \eqn{(M+1)} by statistics matrix), `J`, `move_frac`,
#'   `n_path_calls`, `scatter_max` (largest relative change of the diagonal of
#'   \eqn{G} along the path) and `feas_ok`.
#' @export
selective_meandiff_bc <- function(X, cl_fun, k1, k2,
                                  stats = c("l2", "diag", "max", "ridge"),
                                  adapt = c("diag", "max"),
                                  M = 5000L, L = 2L, eps = "mixture",
                                  eps_range = c(1e-3, 1),
                                  kernel = c("givens", "cayley", "givens_slice"),
                                  slice_max_props = 200L,
                                  cl_obs = NULL, var_floor = 0, ridge_lambda = NULL,
                                  oracle_Sigma = NULL, proj_dim = NULL, proj_n = 10L,
                                  proj_seed = 1L, t2_ridge = 1e-8,
                                  randomized = TRUE,
                                  guard_feas_checks = 25L, scatter_tol = 1e-8,
                                  verbose = FALSE) {
  stopifnot(is.matrix(X), M >= 1L, L >= 1L)
  kernel <- match.arg(kernel)
  M <- as.integer(M); L <- as.integer(L)
  if (!is.null(adapt) && !is.list(adapt)) adapt <- list(adapt = adapt)
  if (!all(unlist(adapt) %in% stats)) stop("`adapt` components must be in `stats`.")
  .as_labels <- function(z) if (is.list(z) && !is.null(z$cluster)) z$cluster else z
  if (is.null(cl_obs)) cl_obs <- .as_labels(cl_fun(X))
  cl_obs <- as.integer(cl_obs)
  if (sum(cl_obs == k1) < 2L || sum(cl_obs == k2) < 2L)
    stop("Tested clusters must each have >= 2 members.")

  st <- .hd_setup(X, cl_obs, k1, k2, cl_fun)
  n_path_calls <- 0L
  feasible <- function(Z) { n_path_calls <<- n_path_calls + 1L; st$feasible(Z) }
  if (!isTRUE(st$feasible(st$Z_obs)))
    stop("Z_obs is not feasible (cl_fun not reproducible on X(Z_obs)?).")
  n_eff <- st$n1 * st$n2 / st$r
  proj <- NULL
  if ("proj" %in% stats) {
    if (is.null(proj_dim)) proj_dim <- max(1L, min(ncol(X), (st$r - 2L) %/% 2L))
    if (proj_dim > st$r - 2L) stop("proj_dim must be at most n1 + n2 - 2.")
    proj <- .hd_draw_proj(ncol(X), proj_dim, proj_n, proj_seed)
  }
  tfun  <- .hd_stat_fun(st$Z_obs, st$v, n_eff, stats, A = st$A, n1 = st$n1,
                        var_floor = var_floor, ridge_lambda = ridge_lambda,
                        oracle_Sigma = oracle_Sigma, proj = proj, t2_ridge = t2_ridge)

  slice <- identical(kernel, "givens_slice")
  if (slice) {
    eps <- NA_real_; draw_eps <- function() NA_real_
  } else if (identical(eps, "mixture")) {
    leps <- log(eps_range)
    draw_eps <- function() exp(stats::runif(1L, leps[1], leps[2]))
  } else {
    if (is.null(eps)) eps <- 0.3
    draw_eps <- function() eps
  }

  # Guard points: evenly spaced indices (no random numbers), their Z kept.
  npath <- M + 1L
  g_idx <- unique(round(seq(0, M, length.out = min(guard_feas_checks, npath))))
  g_Z <- list()
  diagG <- colSums(st$Z_obs^2)
  scatter_max <- 0
  path <- matrix(NA_real_, npath, length(stats), dimnames = list(NULL, stats))

  J <- sample.int(npath, 1L) - 1L
  record <- function(j, Z) {
    path[j + 1L, ] <<- tfun(Z)
    e <- max(abs(colSums(Z^2) - diagG)) / max(diagG)
    if (e > scatter_max) scatter_max <<- e
    if (j %in% g_idx) g_Z[[as.character(j)]] <<- Z
  }
  record(J, st$Z_obs)
  n_moves <- 0L; n_steps <- 0L
  Z <- st$Z_obs
  if (J < M) for (j in (J + 1L):M) {
    kl <- .bc_kernel_L(Z, L, draw_eps(), feasible, kernel, slice_max_props); Z <- kl$Q
    n_moves <- n_moves + kl$n_moves; n_steps <- n_steps + L; record(j, Z)
  }
  Z <- st$Z_obs
  if (J > 0L) for (j in (J - 1L):0L) {
    kl <- .bc_kernel_L(Z, L, draw_eps(), feasible, kernel, slice_max_props); Z <- kl$Q
    n_moves <- n_moves + kl$n_moves; n_steps <- n_steps + L; record(j, Z)
  }
  if (scatter_max > scatter_tol)
    stop(sprintf("GUARD FAILED: total scatter changed along the path (%.2e).", scatter_max))
  feas_ok <- all(vapply(g_Z, function(Zg) isTRUE(st$feasible(Zg)), logical(1)))
  if (!feas_ok) stop("GUARD FAILED: a stored path point is not feasible.")

  j0 <- J + 1L
  pv <- vapply(stats, function(s)
    .hd_rank_p(path[, s], j0, if (randomized) stats::runif(1L) else 1), numeric(1))
  # Min-p combinations: conservative within-path tail probabilities, min over the
  # components; small is extreme. A function of a point and the path multiset,
  # hence exact.
  for (nm in names(adapt)) {
    tail_p <- sapply(adapt[[nm]], function(s) {
      t <- path[, s]; (npath - rank(t, ties.method = "min") + 1) / npath })
    mp <- apply(tail_p, 1, min)
    U  <- if (randomized) stats::runif(1L) else 1
    pv[nm] <- (sum(mp < mp[j0]) + U * sum(mp == mp[j0])) / npath
  }

  out <- list(pvals = pv, obs = path[j0, ], path = path, J = J, M = M, L = L,
              kernel = kernel, n1 = st$n1, n2 = st$n2, p = ncol(X),
              move_frac = if (n_steps) n_moves / n_steps else NA_real_,
              n_path_calls = n_path_calls, scatter_max = scatter_max,
              feas_ok = feas_ok)
  if (verbose) {
    cat(sprintf("[meandiff] n1=%d n2=%d p=%d M=%d L=%d move_frac=%.3f\n",
                st$n1, st$n2, ncol(X), M, L, out$move_frac))
    print(round(pv, 4))
  }
  out
}
