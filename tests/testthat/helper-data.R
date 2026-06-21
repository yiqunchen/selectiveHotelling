# Shared helpers for the selectiveHotelling test suite.
#
# Several functions exercised here are package-internal (not exported in
# NAMESPACE), so they are reached via the `:::` operator. testthat::test_local()
# / devtools::test() load the package with pkgload, which makes the internal
# namespace available, so `:::` resolves correctly during testing.

# Pull a non-exported object out of the package namespace.
sh_internal <- function(name) {
  get(name, envir = asNamespace("selectiveHotelling"))
}

# A small, reproducible null data set together with a clustering that has a
# usable (>= 2 members each) pair of clusters (1, 2). Returns everything the
# polar / Hotelling tests need.
make_null_clustering <- function(seed = 42, n = 24, p = 2, k = 3,
                                  rho = 0.2) {
  set.seed(seed)
  X <- generate_exchangeable_null_data(n, p, rho = rho)$dat
  km <- kmeans_fun(X, k = k)
  list(X = X, km = km, cl = km$cluster)
}

# Reconstruct the data matrix from a polar factor Q exactly as the samplers do
# (mirrors the X_from_Q closure defined inside samplers.R).
x_from_Q <- function(prep, Q) {
  AQS <- prep$A %*% Q %*% prep$S_sqrt
  XQ <- prep$X_fixed
  m <- prep$m
  r <- prep$r
  XQ[prep$idx1, ] <- XQ[prep$idx1, ] + AQS[1:m, , drop = FALSE]
  XQ[prep$idx2, ] <- XQ[prep$idx2, ] + AQS[(m + 1):r, , drop = FALSE]
  XQ
}
