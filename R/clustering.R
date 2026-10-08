# ---------- Hierarchical clustering helpers (for the Besag-Clifford route) ----
#
# Hierarchical clustering is a DETERMINISTIC, data-pure clustering rule (no Lloyd
# iterations, no RNG), which makes it an ideal black-box `cl_fun` for the
# algorithm-agnostic Besag-Clifford selective test (selective_hotelling_bc()).
#
# Convention note: by default the dendrogram is built on SQUARED Euclidean
# distances (`dist(x)^2`), matching the clusterpval / KmeansInference selective
# conditioning convention. Set `square = FALSE` to use ordinary Euclidean
# distances instead.

#' Hierarchical clustering helper for selective inference
#'
#' @title Hierarchical (agglomerative) clustering of a data matrix
#'
#' @description
#' A thin, deterministic wrapper around [stats::hclust()] + [stats::cutree()]
#' that cuts an agglomerative hierarchical clustering into `k` groups. It is a
#' **pure function of `x`** (no randomness, no global RNG use), which makes it a
#' valid black-box clustering routine to pass to
#' [selective_hotelling_bc()] (the package's preferred, exact selective test for
#' arbitrary clustering algorithms).
#'
#' By default the tree is built on **squared** Euclidean distances
#' (`stats::dist(x)^2`), matching the `clusterpval` / `KmeansInference`
#' selective-conditioning convention; set `square = FALSE` for ordinary
#' Euclidean distances.
#'
#' Unlike [kmeans_fun()], the returned object deliberately carries **no**
#' `init_centers_idx` field: hierarchical clustering has no random initialization
#' to condition on, so it should be treated as a pure black box and never
#' re-seeded.
#'
#' @param x A numeric \eqn{n \times p} data matrix with observations in rows.
#' @param k Number of clusters to cut the dendrogram into.
#' @param method Agglomeration method passed to [stats::hclust()] (e.g.
#'   `"average"`, `"single"`, `"complete"`, `"ward.D2"`).
#' @param square If `TRUE` (default) the dendrogram is built on squared Euclidean
#'   distances (`stats::dist(x)^2`, the `clusterpval` convention); if `FALSE`,
#'   ordinary Euclidean distances (`stats::dist(x)`) are used.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`cluster`}{An integer vector of length `nrow(x)` of cluster labels
#'       in `1:k` (from [stats::cutree()]).}
#'     \item{`trajectory`}{A length-1 list, `list(cluster)`. Hierarchical
#'       clustering has no Lloyd iterations, so the "trajectory" is a single
#'       (final) assignment vector.}
#'   }
#'   The object intentionally does **not** include an `init_centers_idx` field.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(30 * 2), 30, 2)
#' cl <- hclust_fun(X, k = 3)$cluster
#' table(cl)
#'
#' @seealso [selective_hotelling_bc()], [selective_hotelling_test()]
#' @export
hclust_fun <- function(x, k = 3, method = "average", square = TRUE) {
  if (!is.matrix(x)) x <- as.matrix(x)
  d <- stats::dist(x)
  if (isTRUE(square)) d <- d^2
  hc <- stats::hclust(d, method = method)
  cluster <- as.integer(stats::cutree(hc, k = k))
  # NOTE: no `init_centers_idx` -- hierarchical clustering is a pure black box
  # with no random initialization to condition on. The single-element trajectory
  # reflects the absence of Lloyd iterations.
  list(cluster = cluster, trajectory = list(cluster))
}

# -----------------------------------------------------------------------------
# make_hclust_cl_fun(k, method, square): internal factory returning a BARE
# closure `function(Xmat) -> integer labels`, suitable for passing directly to
# selective_hotelling_bc() as its `cl_fun`. Deterministic; captures k/method/
# square by value. Same algorithm as sim/bc_general.R::make_hclust_cl_fun; here
# `square` defaults to TRUE (the clusterpval squared-distance convention used by
# the hierarchical sweep), whereas the sim factory defaults to FALSE.
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
make_hclust_cl_fun <- function(k = 3L, method = "average", square = TRUE) {
  force(k); force(method); force(square)
  # Squared distances (the default): the same clustering from one BLAS product
  # and fastcluster (R/fast_cluster.R), much faster when p is large.
  if (isTRUE(square)) return(make_fast_hclust_cl_fun(k, method))
  function(Xmat) {
    if (!is.matrix(Xmat)) Xmat <- as.matrix(Xmat)
    d <- stats::dist(Xmat)
    if (isTRUE(square)) d <- d^2
    hc <- stats::hclust(d, method = method)
    as.integer(stats::cutree(hc, k = k))
  }
}
