# Thorough check that the fast clustering routines are the package routines:
# (1) identical labels on many data sets from the cluster design's DGPs;
# (2) identical Besag-Clifford paths and p-values with the same seed.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("dev/highdim/dgp.R")
set.seed(2026); mis <- c(km = 0, hc = 0); N <- 0
cfg <- expand.grid(n = c(50, 100, 500), ratio = c(2, 5), sig = c("ar1", "cs", "hetero", "spike5"),
                   signal = c("null", "line", "block", "dense"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(cfg))) for (rep in 1:(if (cfg$n[i] == 500) 2 else 6)) {
  g <- cfg[i, ]; p <- g$n * g$ratio; s <- if (g$signal == "null") 0 else sample(c(2, 8), 1)
  d <- gen_highdim(g$n, p, g$sig, if (g$signal == "null") "line" else g$signal, s, seed = 1e5 * i + rep)
  init <- sample.int(g$n, 3); N <- N + 1
  if (!identical(make_fast_kmeans_cl_fun(init)(d$X), as.integer(kmeans_with_fixed_init(d$X, init)$cluster))) mis["km"] <- mis["km"] + 1
  if (!identical(make_fast_hclust_cl_fun(3)(d$X), hclust_fun(d$X, 3, square = TRUE)$cluster)) mis["hc"] <- mis["hc"] + 1
}
cat(sprintf("(1) %d data sets: k-means mismatches %d, average-linkage mismatches %d\n", N, mis["km"], mis["hc"]))
st <- c("l2", "cq", "diag", "diag_w", "max", "ridge", "proj", "t2")
for (cf in list(list(50, 250, "km"), list(50, 250, "hc"), list(100, 500, "km"), list(100, 500, "hc"),
                list(500, 1000, "km"), list(500, 1000, "hc"))) {
  n <- cf[[1]]; p <- cf[[2]]; alg <- cf[[3]]
  d <- gen_highdim(n, p, "cs", "line", 8, seed = 7)
  init <- sample.int(n, 3)
  f_new <- if (alg == "km") make_fast_kmeans_cl_fun(init) else make_fast_hclust_cl_fun(3)
  f_old <- if (alg == "km") function(x) kmeans_with_fixed_init(x, init)$cluster else function(x) hclust_fun(x, 3, square = TRUE)$cluster
  cl <- f_new(d$X); tb <- sort(table(cl), decreasing = TRUE); k <- as.integer(names(tb)[1:2])
  M <- if (n == 500) 150 else 1000
  set.seed(3); t1 <- system.time(a <- selective_meandiff_bc(d$X, f_old, k[1], k[2], M = M, stats = st))[["elapsed"]]
  set.seed(3); t2 <- system.time(b <- selective_meandiff_bc(d$X, f_new, k[1], k[2], M = M, stats = st))[["elapsed"]]
  cat(sprintf("(2) n=%d p=%d %s r=%d M=%d: identical path %s, identical p-values %s, move %.2f; %.1fs -> %.1fs\n",
              n, p, alg, b$n1 + b$n2, M, identical(a$path, b$path), identical(a$pvals, b$pvals), b$move_frac, t1, t2))
}
