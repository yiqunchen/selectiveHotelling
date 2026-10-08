# =============================================================================
# manifest.R -- task list for the high-dimensional study (p > n1 + n2) of
# selective_meandiff_bc(). One row = one array task = N_shard data sets of one
# cell. Run from the package root:  Rscript dev/highdim/cluster/manifest.R
#
# Design (author, 2026-10-07):
#   clustering  k-means (K = 3, initial rows drawn from the seed, then fixed) and
#               average linkage (squared Euclidean, cut at K = 3)
#   n           50, 100, 500;   p = 2n and 5n
#   Sigma       ar1 (0.6), cs (0.5), hetero (diagonal, a decade), spike5
#   signal      line (1-sparse, the main grid's design), block (1/4-sparse),
#               dense; size s = delta^2 / sqrt(2 tr Sigma^2)
#   null        s = 0 (N_NULL data sets per cell)
#   power       s in S_GRID (N_POW data sets per cell)
#   test        every pair whose clusters both have >= MIN_SIZE rows (with
#               K = 3 the pair event fixes the whole partition, so any pair rule
#               has constant weight on it); M = 5000, L = 2, Givens kernel,
#               eps = "mixture"; all statistics ranked on the one path.
#
# Shards: a whole cell per task for n <= 100 (<= ~10 CPU-h), SHARD_500 data sets
# per task for n = 500 (~13 CPU-h); with 4 cores a task takes a few hours.
# Env overrides: N_NULL, N_POW, N_NULL_500, N_POW_500, SHARD_SMALL, SHARD_500, M,
# S_GRID (comma list).
# =============================================================================
env_num <- function(nm, d) as.numeric(Sys.getenv(nm, d))
N_NULL     <- env_num("N_NULL", 500);  N_POW     <- env_num("N_POW", 200)
N_NULL_500 <- env_num("N_NULL_500", 300); N_POW_500 <- env_num("N_POW_500", 100)
SHARD_SMALL <- env_num("SHARD_SMALL", 500); SHARD_500 <- env_num("SHARD_500", 50)
M          <- as.integer(env_num("M", 5000)); L <- 2L; MIN_SIZE <- 3L
S_GRID     <- as.numeric(strsplit(Sys.getenv("S_GRID", "1,2,4,8,16,32"), ",")[[1]])
OUT        <- Sys.getenv("HD_MANIFEST", "dev/highdim/cluster/manifest.rds")

cells <- rbind(
  expand.grid(alg = c("kmeans", "average"), n = c(50, 100, 500), ratio = c(2, 5),
              sig = c("ar1", "cs", "hetero", "spike5"), signal = "null", s = 0,
              stringsAsFactors = FALSE),
  expand.grid(alg = c("kmeans", "average"), n = c(50, 100, 500), ratio = c(2, 5),
              sig = c("ar1", "cs", "hetero", "spike5"), signal = c("line", "block", "dense"),
              s = S_GRID, stringsAsFactors = FALSE))
cells$p <- cells$n * cells$ratio
cells$N <- ifelse(cells$signal == "null", ifelse(cells$n == 500, N_NULL_500, N_NULL),
                                          ifelse(cells$n == 500, N_POW_500, N_POW))
cells$cell_id <- seq_len(nrow(cells))

# Cost model (seconds per tested pair at M = 5000, L = 2, from local timings with
# an optimised BLAS: clustering + reconstruction per kernel step, statistics per
# path point). Two eligible pairs per data set on average. Check it against the
# smoke run: the runner logs the CPU seconds per pair.
t_pair <- function(n, p) ifelse(n <= 50, 25, ifelse(n <= 100, 35, ifelse(p <= 1000, 420, 500)))

man <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  g <- cells[i, ]; SHARD <- if (g$n >= 500) SHARD_500 else SHARD_SMALL; ns <- ceiling(g$N / SHARD)
  data.frame(g[rep(1, ns), c("cell_id", "alg", "n", "p", "ratio", "sig", "signal", "s")],
             shard = seq_len(ns), N_shard = pmin(SHARD, g$N - (seq_len(ns) - 1) * SHARD),
             row.names = NULL, stringsAsFactors = FALSE)
}))
man$task_id   <- seq_len(nrow(man))
# disjoint seed blocks: cell c, shard k -> 1e7 + c * 1e5 + (k - 1) * shard size + 1..
shard_size <- ifelse(man$n >= 500, SHARD_500, SHARD_SMALL)
man$base_seed <- 1e7 + man$cell_id * 1e5 + (man$shard - 1) * shard_size
man$M <- M; man$L <- L; man$min_size <- MIN_SIZE; man$K <- 3L
man$cpu_h <- man$N_shard * 2 * t_pair(man$n, man$p) / 3600
saveRDS(man, OUT)
write.csv(man, sub("\\.rds$", "_tasks.csv", OUT), row.names = FALSE)

cat(sprintf("%d cells, %d tasks -> %s\n", nrow(cells), nrow(man), OUT))
agg <- aggregate(cbind(N_shard, cpu_h) ~ alg + n + p, man, sum)
print(agg[order(agg$alg, agg$n, agg$p), ], row.names = FALSE)
cat(sprintf("TOTAL model ~%.0f CPU-h; largest task ~%.1f CPU-h; array range 1-%d\n",
            sum(man$cpu_h), max(man$cpu_h), nrow(man)))
cat("Arrays above the cluster's MaxArraySize: submit chunks with HD_TASK_OFFSET (see README).\n")
