# Size/power pilot for selective_meandiff_bc() across p < r and p > r.
# n = 50, p = 100 (pairs of ~33 rows, so p > r), k-means K = 3 from a fixed random init, test estimated labels (1, 2).
# Truth: 3 clusters of 17/17/16 with collinear means -delta u, 0, +delta u (delta = 0: null).
# delta is set from kappa = n_eff delta^2 / sqrt(2 tr Sigma^2), n_eff = n/6.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
out_file <- if (length(args)) args[1] else "dev/highdim/sim_raw.rds"
reps_null <- as.integer(Sys.getenv("REPS_NULL", 1000)); reps_pow <- as.integer(Sys.getenv("REPS_POW", 200))
M <- as.integer(Sys.getenv("M", 1000))

# Same covariances as the main grid (sim/common.R, sim/dgp_ext.R in the paper repo).
make_sigma <- function(type, p) switch(type,
  iso    = diag(p),
  ar1    = 0.6^abs(outer(1:p, 1:p, "-")),
  cs     = { S <- matrix(0.5, p, p); diag(S) <- 1; S },
  spike5 = diag(p) + 5 * tcrossprod(rep(1 / sqrt(p), p)),
  hetero = diag(exp(seq(-log(10) / 2, log(10) / 2, length.out = p))))

# The main grid's quarter-sparse "block" means (sim/dgp_ext.R .blocksparse_centers):
# support q = round(p/4), centres (delta/sqrt 2) u_k with u_1..u_3 orthonormal
# on the block, drawn once from a fixed seed.
block_centers <- function(p, delta, K = 3, basis_seed = 20260905L) {
  q <- as.integer(min(p, max(K, round(p / 4))))
  had <- exists(".Random.seed", envir = .GlobalEnv)
  if (had) saved <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(basis_seed)
  U <- qr.Q(qr(matrix(stats::rnorm(q * K), q, K)))
  if (had) assign(".Random.seed", saved, envir = .GlobalEnv)
  C <- matrix(0, K, p); C[, seq_len(q)] <- t(U) * (delta / sqrt(2)); C
}

one_rep <- function(p, sig, signal, kappa, seed, n = 50, K = 3) {
  set.seed(seed)
  Sig <- make_sigma(sig, p); R <- chol(Sig)
  truth <- sort(rep_len(1:3, n))
  delta <- sqrt(kappa * sqrt(2 * sum(Sig^2)) / (n / 6))
  if (signal == "block") {
    C <- block_centers(p, delta)                    # every pair of centres delta apart
  } else {
    u <- switch(signal, dense = rep(1, p), sparse = c(rep(1, 5), rep(0, p - 5)),
                line = c(1, rep(0, p - 1)))           # line: the main grid's 1-sparse design
    C <- outer(c(-1, 0, 1) * delta, u / sqrt(sum(u^2)))   # adjacent centres delta apart
  }
  X <- matrix(rnorm(n * p), n, p) %*% R + C[truth, , drop = FALSE]
  init <- sample.int(n, K)
  cl_fun <- function(x) kmeans_with_fixed_init(x, init)$cluster
  cl <- cl_fun(X)
  n1 <- sum(cl == 1); n2 <- sum(cl == 2)
  base <- data.frame(p = p, sig = sig, signal = signal, kappa = kappa, seed = seed,
                     n1 = n1, n2 = n2, r = n1 + n2,
                     pure = length(unique(truth[cl == 1])) == 1 && length(unique(truth[cl == 2])) == 1 &&
                            truth[cl == 1][1] != truth[cl == 2][1])
  if (n1 < 2 || n2 < 2) return(cbind(base, ok = FALSE))
  st <- c("l2", "cq", "diag", "diag_w", "max", "ridge", "proj", "oracle", if (p <= n1 + n2 - 3) "hotelling")
  res <- tryCatch(selective_meandiff_bc(X, cl_fun, 1, 2, M = M, stats = st,
                                        adapt = c("diag", "max"), cl_obs = cl, oracle_Sigma = Sig),
                  error = function(e) NULL)
  if (is.null(res)) return(cbind(base, ok = FALSE))
  pv <- as.list(res$pvals); if (is.null(pv$hotelling)) pv$hotelling <- NA_real_
  cbind(base, ok = TRUE, as.data.frame(pv), move = res$move_frac)
}

grid <- rbind(
  expand.grid(p = 100, sig = c("iso", "ar1", "cs", "spike5", "hetero"), signal = "dense",
              kappa = 0, rep = seq_len(reps_null), stringsAsFactors = FALSE),
  expand.grid(p = 100, sig = c("iso", "ar1", "cs", "spike5", "hetero"), signal = c("dense", "sparse", "block"),
              kappa = c(4, 8, 16), rep = seq_len(reps_pow), stringsAsFactors = FALSE))
grid$seed <- 1e6 + seq_len(nrow(grid))
t0 <- Sys.time()
rows <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]; one_rep(g$p, g$sig, g$signal, g$kappa, g$seed) }, mc.cores = 14)
bad <- vapply(rows, inherits, logical(1), "try-error")
cat("errors:", sum(bad), "\n")
res <- do.call(rbind, lapply(rows[!bad], function(z) { for (c in c("l2","cq","diag","diag_w","max","ridge","proj","oracle","adapt","hotelling","move")) if (is.null(z[[c]])) z[[c]] <- NA_real_; z }))
saveRDS(res, out_file)
cat("done in", format(Sys.time() - t0), "\n")
