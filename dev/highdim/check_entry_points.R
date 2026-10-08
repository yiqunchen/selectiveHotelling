# The package's Monte Carlo entry points give the same p-values with the fast
# clustering routines as with the original ones (same seed).
suppressMessages(pkgload::load_all(".", quiet = TRUE))
old_km <- function(init) function(x) kmeans_with_fixed_init(x, init)$cluster
old_hc <- function(K) function(x) hclust_fun(x, K, square = TRUE)$cluster
# outcome = the p-value, or the error message (both routes must agree on either)
oc <- function(expr) tryCatch(expr, error = function(e) paste("ERR", conditionMessage(e)))
res <- NULL
for (s in 1:12) for (p in c(2, 10, 40)) {
  set.seed(s); X <- matrix(rnorm(90 * p), 90, p) + outer(rep(c(0, 2.5, 5), each = 30), c(1, rep(0, p - 1)))
  # selective_hotelling_test, k-means: replicate its init draw, then the old cl_fun
  set.seed(100 + s); a <- oc(selective_hotelling_test(X, 1, 2, K = 3, M = 200)$pval)
  set.seed(100 + s); km <- kmeans_with_trajectory(X, k = 3)
  b <- oc(selective_hotelling_bc(X, old_km(km$init_centers_idx), 1, 2, M = 200, cl_obs = as.integer(km$cluster))$pval)
  # selective_hotelling_test, hclust
  cl <- old_hc(3)(X); ok <- sum(cl == 1) >= 3 && sum(cl == 2) >= 3
  if (ok) { set.seed(200 + s); c1 <- selective_hotelling_test(X, 1, 2, K = 3, clustering = "hclust", M = 200)$pval
            set.seed(200 + s); c2 <- selective_hotelling_bc(X, old_hc(3), 1, 2, M = 200)$pval } else c1 <- c2 <- NA
  # selective_truncated_F, method = "bc", k-means and hclust
  init <- c(1, 31, 61)
  set.seed(300 + s); d1 <- oc(selective_truncated_F(X, 1, 2, clustering = "kmeans", method = "bc", init_centers_idx = init, M = 200)$pval)
  set.seed(300 + s); d2 <- oc(selective_truncated_F(X, 1, 2, clustering = "custom", method = "bc", cl_fun = old_km(init),
                                                 cl_obs = old_km(init)(X), M = 200)$pval)
  if (ok) { set.seed(400 + s); e1 <- selective_truncated_F(X, 1, 2, clustering = "hclust", method = "bc", M = 200)$pval
            set.seed(400 + s); e2 <- selective_truncated_F(X, 1, 2, clustering = "custom", method = "bc", cl_fun = old_hc(3),
                                                           cl_obs = old_hc(3)(X), M = 200)$pval } else e1 <- e2 <- NA
  res <- rbind(res, data.frame(s, p, test_km = identical(a, b), test_hc = if (ok) identical(c1, c2) else NA,
                               truncF_km = identical(d1, d2), truncF_hc = if (ok) identical(e1, e2) else NA))
}
print(colSums(res[, -(1:2)], na.rm = TRUE)); cat("cases:", nrow(res), " (hclust cases with a testable pair:", sum(!is.na(res$test_hc)), ")\n")
