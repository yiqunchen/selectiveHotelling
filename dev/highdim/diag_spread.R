# Spread of each statistic on the orbit: unrestricted vs selection-restricted,
# and how much a long path moves. Null data, k-means K = 3 with fixed init.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
spread <- function(p, M, sel, seed = 4, n = 90) {
  set.seed(seed); X <- matrix(rnorm(n * p), n, p)
  init <- sample.int(n, 3)
  cl_fun <- if (sel) function(x) kmeans_with_fixed_init(x, init)$cluster else
    { cl0 <- kmeans_with_fixed_init(X, init)$cluster; function(x) cl0 }
  b <- selective_meandiff_bc(X, cl_fun, 1, 2, M = M, stats = c("l2", "diag", "max", "ridge"),
                             cl_obs = kmeans_with_fixed_init(X, init)$cluster)
  q <- apply(b$path, 2, function(z) round(c(cv = sd(z) / mean(z), obs_q = mean(z <= z[b$J + 1])), 3))
  cat(sprintf("p=%d sel=%s M=%d n1=%d n2=%d move=%.2f\n", p, sel, M, b$n1, b$n2, b$move_frac)); print(q)
}
for (p in c(10, 200)) for (sel in c(FALSE, TRUE)) spread(p, 3000, sel)
