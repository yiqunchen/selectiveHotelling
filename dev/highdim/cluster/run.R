# =============================================================================
# run.R -- shard runner for the high-dimensional study (dev/highdim/cluster/manifest.R).
# One manifest row = one task = N_shard testable data sets of one cell.
#
# Per data set: gen_highdim() (dev/highdim/dgp.R) -> fast k-means from rows drawn
# from the seed, or fast average linkage (R/fast_cluster.R) -> every pair whose
# clusters both have >= min_size rows -> selective_meandiff_bc() with all
# statistics on one Besag-Clifford path (the package defaults for M, L, kernel, eps).
#
# Output: one row per tested pair with every p-value (columns p_<stat>), the
# observed statistics are not stored. `pair_alt` = both tested clusters pure and
# from different true groups; `recover` = the whole partition equals the truth.
#
# Env: HD_MANIFEST, HD_SHARDS (output dir), PKG_DIR (package checkout),
# SLURM_ARRAY_TASK_ID (+ HD_TASK_OFFSET), SLURM_CPUS_PER_TASK, HD_SMOKE=k (k data sets into
# <shards>/_smoke), HD_OVERDRAW (cap on screened draws = HD_OVERDRAW x N_shard).
# Resume: a task whose shard exists and reads back exits at once.
# =============================================================================
PKG_DIR  <- Sys.getenv("PKG_DIR", "selectiveHotelling")
MANIFEST <- Sys.getenv("HD_MANIFEST", file.path(PKG_DIR, "dev/highdim/cluster/manifest.rds"))
SHARDS   <- Sys.getenv("HD_SHARDS", file.path(PKG_DIR, "dev/highdim/cluster/shards"))
task     <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")) +
            as.integer(Sys.getenv("HD_TASK_OFFSET", "0"))     # chunked arrays
cores    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("CORES", "1")))
SMOKE    <- as.integer(Sys.getenv("HD_SMOKE", "0"))
OVER     <- as.numeric(Sys.getenv("HD_OVERDRAW", "200"))

suppressMessages(pkgload::load_all(PKG_DIR, quiet = TRUE))
source(file.path(PKG_DIR, "dev/highdim/dgp.R"))
man <- readRDS(MANIFEST); row <- man[man$task_id == task, , drop = FALSE]
stopifnot(nrow(row) == 1L)
NS <- if (SMOKE > 0L) min(SMOKE, row$N_shard) else row$N_shard
out_dir  <- if (SMOKE > 0L) file.path(SHARDS, "_smoke") else file.path(SHARDS, row$alg)
out_file <- file.path(out_dir, sprintf("shard_%05d.rds", task))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (SMOKE == 0L && file.exists(out_file)) {
  ok <- tryCatch({ d <- readRDS(out_file); is.data.frame(d) && nrow(d) > 0 &&
    all(d$cell_id == row$cell_id) }, error = function(e) FALSE)
  if (isTRUE(ok)) { cat(sprintf("[hd %05d] SKIP (verified shard)\n", task)); quit(save = "no") }
}
cat(sprintf("[hd %05d] cell=%d alg=%s n=%d p=%d sig=%s signal=%s s=%g N=%d cores=%d BLAS=%s\n",
            task, row$cell_id, row$alg, row$n, row$p, row$sig, row$signal, row$s, NS,
            cores, basename(sessionInfo()$BLAS)))

STATS <- c("l2", "cq", "diag", "diag_w", "max", "ridge", "proj", "t2", "oracle")
ADAPT <- list(adapt_l2_max = c("l2", "max"), adapt_dw_max = c("diag_w", "max"))
K <- as.integer(row$K)

setup_seed <- function(seed) {
  d <- gen_highdim(row$n, row$p, row$sig, if (row$signal == "null") "line" else row$signal,
                   row$s, seed = seed, K = K)
  init <- NULL
  if (row$alg == "kmeans") {
    init <- local({ old <- .Random.seed; on.exit(assign(".Random.seed", old, envir = globalenv()))
                    set.seed(as.integer((seed * 7919) %% 2147483647)); sample.int(row$n, K) })
    cl_fun <- make_fast_kmeans_cl_fun(init)
  } else cl_fun <- make_fast_hclust_cl_fun(K, row$alg)
  cl <- cl_fun(d$X)
  labs <- sort(unique(cl)); sz <- tabulate(cl, max(labs))
  prs <- if (length(labs) >= 2L) utils::combn(labs, 2L) else matrix(integer(0), 2L, 0L)
  prs <- prs[, sz[prs[1, ]] >= row$min_size & sz[prs[2, ]] >= row$min_size, drop = FALSE]
  list(d = d, cl = cl, cl_fun = cl_fun, prs = prs)
}

one <- function(seed) {
  st <- setup_seed(seed); d <- st$d; cl <- st$cl
  if (!ncol(st$prs)) return(NULL)
  rec <- as.integer(all(table(cl, d$truth) %in% c(0, as.numeric(table(d$truth)))) &&
                    length(unique(cl)) == K)
  out <- list()
  for (q in seq_len(ncol(st$prs))) {
    k1 <- st$prs[1, q]; k2 <- st$prs[2, q]
    t1 <- unique(d$truth[cl == k1]); t2 <- unique(d$truth[cl == k2])
    set.seed(as.integer((seed * 4 + q) %% 2147483647))
    t0 <- proc.time()
    res <- selective_meandiff_bc(d$X, st$cl_fun, k1, k2, stats = STATS, adapt = ADAPT,
                                 cl_obs = cl, oracle_Sigma = d$Sigma)
    tt <- proc.time() - t0
    pv <- as.list(res$pvals); names(pv) <- paste0("p_", names(pv))
    out[[q]] <- data.frame(rep = seed, task = task, cell_id = row$cell_id, alg = row$alg,
                           n = row$n, p = row$p, ratio = row$ratio, sig = row$sig,
                           signal = row$signal, s = row$s, delta = d$delta, k1 = k1, k2 = k2,
                           n1 = res$n1, n2 = res$n2,
                           pair_alt = as.integer(length(t1) == 1L && length(t2) == 1L && t1 != t2),
                           recover = rec, pv, move_frac = res$move_frac,
                           n_path_calls = res$n_path_calls,
                           cpu = tt[["user.self"]] + tt[["sys.self"]], secs = tt[["elapsed"]],
                           M = res$M, L = res$L, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# Screen seeds (one clustering each) until NS testable ones, then run the test on
# exactly those. Base block first, then a disjoint per-task extension range.
EXT_BASE <- 5e8 + task * 1e6
t_start <- proc.time()[["elapsed"]]
cap <- ceiling(OVER * NS); n_drawn <- 0L; keep <- numeric(0)
cand <- function(i) if (i <= NS) row$base_seed + i else EXT_BASE + (i - NS)
while (length(keep) < NS && n_drawn < cap) {
  n_drawn <- n_drawn + 1L; s <- cand(n_drawn)
  if (ncol(setup_seed(s)$prs) > 0L) keep <- c(keep, s)
}
cat(sprintf("[hd %05d] screened %d draws: %d testable (yield %.3f)\n",
            task, n_drawn, length(keep), length(keep) / n_drawn))
f <- function(s) tryCatch(one(s), error = function(e) e)
res <- if (cores > 1L) parallel::mclapply(keep, f, mc.cores = cores, mc.preschedule = FALSE) else lapply(keep, f)
bad <- vapply(res, function(z) inherits(z, "error") || inherits(z, "try-error"), logical(1))
if (any(bad)) {
  msg <- unique(vapply(res[bad], function(z) if (inherits(z, "error")) conditionMessage(z) else as.character(z), ""))
  stop(sprintf("task %d: %d/%d data sets errored; no shard written. First: %s", task, sum(bad), length(bad), msg[1]))
}
df <- do.call(rbind, Filter(is.data.frame, res))
if (!is.null(df)) { df$n_drawn <- n_drawn; df$n_testable <- length(keep) }
tmp <- paste0(out_file, ".tmp", Sys.getpid()); saveRDS(df, tmp); file.rename(tmp, out_file)
cat(sprintf("[hd %05d] DONE %d data sets, %d pairs, wall %.0fs, %.1f CPU-s/pair -> %s\n",
            task, length(keep), if (is.null(df)) 0L else nrow(df), proc.time()[["elapsed"]] - t_start,
            if (is.null(df)) NA else mean(df$cpu), out_file))
