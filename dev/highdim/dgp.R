# Data-generating process for the high-dimensional study. Covariances and the
# "line" / "block" means follow the main grid (sim/common.R, sim/dgp_ext.R in the
# paper repo); "dense" spreads the shift over all p coordinates.
# Truth: K = 3 groups of sizes as equal as possible.
# Signal size s = delta^2 / sqrt(2 tr Sigma^2), so delta is comparable across
# Sigma and p. In "line" and "dense" the centres are -delta u, 0, +delta u
# (adjacent pairs delta apart); in "block" every pair of centres is delta apart.
make_sigma <- function(type, p) switch(type,
  iso    = diag(p),
  ar1    = 0.6^abs(outer(seq_len(p), seq_len(p), "-")),
  cs     = { S <- matrix(0.5, p, p); diag(S) <- 1; S },
  spike5 = diag(p) + 5 * tcrossprod(rep(1 / sqrt(p), p)),
  hetero = diag(exp(seq(-log(10) / 2, log(10) / 2, length.out = p)), p),
  stop("unknown Sigma: ", type))

# tr(Sigma^2) without forming Sigma^2
tr_sigma2 <- function(Sig) sum(Sig^2)

block_centers <- function(p, delta, K = 3, basis_seed = 20260905L) {
  q <- as.integer(min(p, max(K, round(p / 4))))
  had <- exists(".Random.seed", envir = .GlobalEnv)
  if (had) saved <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(basis_seed)
  U <- qr.Q(qr(matrix(stats::rnorm(q * K), q, K)))
  if (had) assign(".Random.seed", saved, envir = .GlobalEnv)
  C <- matrix(0, K, p); C[, seq_len(q)] <- t(U) * (delta / sqrt(2)); C
}

gen_highdim <- function(n, p, sig, signal, s, seed, K = 3) {
  set.seed(seed)
  Sig <- make_sigma(sig, p)
  truth <- sort(rep_len(seq_len(K), n))
  delta <- sqrt(s * sqrt(2 * tr_sigma2(Sig)))
  C <- if (delta == 0) matrix(0, K, p) else if (signal == "block") block_centers(p, delta, K) else {
    u <- switch(signal, line = c(1, rep(0, p - 1)), dense = rep(1 / sqrt(p), p),
                stop("unknown signal: ", signal))
    outer(c(-1, 0, 1) * delta, u) }
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(Sig) + C[truth, , drop = FALSE]
  list(X = X, truth = truth, Sigma = Sig, delta = delta)
}
