# =============================================================================
# fast_cluster.R -- fast clustering routines to pass as `cl_fun`.
#
# Both are deterministic functions of the data, so they are valid black boxes
# for the Besag-Clifford tests. They are faster versions of the package
# routines, not new algorithms:
#   * make_fast_kmeans_cl_fun(): Lloyd's algorithm from fixed initial rows, the
#     algorithm of kmeans_with_fixed_init(), with the point-to-centre distances
#     and the centre updates as BLAS products. At n = 500, p = 2500 a call takes
#     about 27 ms, against 140 ms for stats::kmeans(algorithm = "Lloyd") and
#     300 ms for kmeans_with_fixed_init().
#   * make_fast_hclust_cl_fun(): squared Euclidean distances from one BLAS
#     product, ||x_i||^2 + ||x_j||^2 - 2 x_i'x_j, and fastcluster::hclust()
#     (Muellner 2013) when installed; the clustering of hclust_fun(square = TRUE).
# The distance matrix dominates hierarchical clustering when p is large:
# stats::dist() costs about 1.4 s at n = 500, p = 2500 and the BLAS product
# about 0.03 s with an optimised BLAS.
# =============================================================================

#' Fast k-means from fixed initial rows
#'
#' Returns a clustering routine that runs Lloyd's algorithm from the rows
#' `init_centers_idx` of its argument, the algorithm of
#' [kmeans_with_fixed_init()]: assign each row to the nearest centre (ties to the
#' first), stop when the assignment repeats, otherwise move each centre to the
#' mean of its rows (an empty cluster keeps its centre). The distances,
#' \eqn{\|c_j\|^2 - 2 x_i'c_j} up to a term constant in \eqn{j}, and the centre
#' sums are BLAS products, so a call is several times faster than
#' [kmeans_with_fixed_init()] when \eqn{p} is large.
#'
#' @param init_centers_idx Row indices of the initial centres, fixed before the
#'   data are seen.
#' @param max_iter Maximum number of Lloyd iterations.
#' @return A function `function(x) -> integer labels`, for use as `cl_fun`.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 5), 60, 5)
#' cl_fun <- make_fast_kmeans_cl_fun(c(1, 21, 41))
#' table(cl_fun(X))
#' @export
make_fast_kmeans_cl_fun <- function(init_centers_idx, max_iter = 50L) {
  force(init_centers_idx); force(max_iter)
  function(x) {
    x <- as.matrix(x); n <- nrow(x); k <- length(init_centers_idx)
    C <- x[init_centers_idx, , drop = FALSE]; cl_prev <- NULL
    for (t in seq_len(max_iter)) {
      D  <- rep(rowSums(C^2), each = n) - 2 * tcrossprod(x, C)
      cl <- max.col(-D, ties.method = "first")
      if (!is.null(cl_prev) && identical(cl, cl_prev)) break
      cl_prev <- cl
      nk  <- tabulate(cl, k)
      Ind <- matrix(0, n, k); Ind[cbind(seq_len(n), cl)] <- 1
      S   <- crossprod(Ind, x)
      ok  <- nk > 0; C[ok, ] <- S[ok, , drop = FALSE] / nk[ok]
    }
    as.integer(cl)
  }
}

#' Squared Euclidean distances from one matrix product
#'
#' @keywords internal
#' @noRd
.sq_dist_blas <- function(x) {
  s <- rowSums(x^2)
  D <- outer(s, s, "+") - 2 * tcrossprod(x)
  D[D < 0] <- 0
  stats::as.dist(D)
}

#' Fast agglomerative clustering on squared Euclidean distances
#'
#' Returns a clustering routine that cuts an agglomerative clustering into `k`
#' groups, as `hclust_fun(x, k, method, square = TRUE)`, with the squared
#' distances computed from one BLAS product and the tree built by
#' `fastcluster::hclust()` when the fastcluster package is installed
#' ([stats::hclust()] otherwise).
#'
#' @param k Number of clusters.
#' @param method Agglomeration method (e.g. `"average"`).
#' @return A function `function(x) -> integer labels`, for use as `cl_fun`.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 5), 60, 5)
#' cl_fun <- make_fast_hclust_cl_fun(3)
#' table(cl_fun(X))
#' @export
make_fast_hclust_cl_fun <- function(k = 3L, method = "average") {
  force(k); force(method)
  hc <- if (requireNamespace("fastcluster", quietly = TRUE)) fastcluster::hclust else stats::hclust
  function(x) as.integer(stats::cutree(hc(.sq_dist_blas(as.matrix(x)), method = method), k = k))
}
