# ---------- User-facing wrapper & S3 methods -------------------------------

#' Selective Hotelling T-squared test after clustering
#'
#' Clusters `X` (k-means or hierarchical clustering, `K` clusters) and tests whether
#' the discovered clusters `k1` and `k2` have the same mean, with Hotelling's
#' \eqn{T^2} and an unknown covariance. The p-value is the Besag--Clifford test
#' [selective_hotelling_bc()] with the random Givens kernel on the event that the
#' clustering returns the tested pair (Algorithm 1 of the paper). It is exactly
#' uniform under the null for any `M` and `L`.
#'
#' For k-means, the initial centres are drawn once (from the current random seed)
#' and then held fixed, so that re-running the clustering on perturbed data is a
#' deterministic function of the data. For the one-dimensional test, use
#' [selective_truncated_F()].
#'
#' @param X A numeric data matrix with observations in rows.
#' @param k1,k2 The two cluster labels (in `1:K`) to compare.
#' @param K Number of clusters.
#' @param clustering `"kmeans"` (default, Lloyd's algorithm from fixed initial
#'   centres) or `"hclust"` (agglomerative clustering cut at `K`).
#' @param hclust_method Agglomeration method for `clustering = "hclust"`, passed
#'   to [hclust_fun()] / [stats::hclust()]. Ignored for `clustering = "kmeans"`.
#' @param hclust_square If `TRUE` (default) the tree is built on squared Euclidean
#'   distances (the `clusterpval` convention). Ignored for `clustering = "kmeans"`.
#' @param M Path length (the path carries `M + 1` points). The p-value cannot be
#'   smaller than \eqn{1/(M+1)}.
#' @param L Number of kernel steps between consecutive path points.
#' @param kernel Kernel passed to [selective_hotelling_bc()]: `"givens"` (default),
#'   `"cayley"`, or `"givens_slice"`.
#' @param ridge Ridge added to the pooled covariance in [hotelling_T2()].
#' @param verbose Whether to print progress messages.
#'
#' @return An S3 object of class `"selectiveHotelling"`; see
#'   [selective_hotelling_bc()] for its elements. It also contains `cluster`, the
#'   cluster labels of the clustering that was tested, so that `k1` and `k2` can be
#'   matched to observations.
#'
#' @examples
#' set.seed(2)
#' Xb <- rbind(
#'   matrix(rnorm(20), 10, 2),
#'   matrix(rnorm(20) + 6, 10, 2),
#'   matrix(rnorm(20) - 6, 10, 2)
#' )
#' res <- selective_hotelling_test(Xb, 1, 2, K = 3, M = 200)
#' res$pval
#' print(res)
#' res_h <- selective_hotelling_test(Xb, 1, 2, K = 3, clustering = "hclust", M = 200)
#' res_h$pval
#'
#' @seealso [selective_hotelling_bc()], [selective_truncated_F()],
#'   [kmeans_fun()], [hclust_fun()]
#' @export
selective_hotelling_test <- function(X, k1, k2, K = 3,
                                     clustering = c("kmeans", "hclust"),
                                     hclust_method = "average",
                                     hclust_square = TRUE,
                                     M = 5000L, L = 2L, kernel = "givens",
                                     ridge = 1e-8, verbose = FALSE) {
  clustering <- match.arg(clustering)
  if (clustering == "hclust") {
    cl_fun <- make_hclust_cl_fun(k = K, method = hclust_method, square = hclust_square)
    cl_obs <- NULL
  } else {
    km <- kmeans_with_trajectory(X, k = K)          # draws the initial centres once
    init <- km$init_centers_idx
    cl_fun <- function(x) kmeans_with_fixed_init(x, init_centers_idx = init)$cluster
    cl_obs <- as.integer(km$cluster)
    if (!identical(as.integer(cl_fun(X)), cl_obs))
      stop("k-means from the fixed initial centres does not reproduce the clustering.")
  }
  res <- selective_hotelling_bc(X, cl_fun, k1, k2, M = M, L = L, eps = "mixture",
                                cl_obs = cl_obs, ridge = ridge, verbose = verbose,
                                kernel = kernel)
  res$cluster <- if (is.null(cl_obs)) as.integer(cl_fun(X)) else cl_obs   # labels of the clustering that was tested
  res
}

#' Print a selectiveHotelling object
#'
#' @param x A `"selectiveHotelling"` object from [selective_hotelling_test()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.selectiveHotelling <- function(x, ...) {
  cat("Selective Hotelling T^2 test (unknown covariance)\n")
  cat(sprintf("  Clusters tested : %s vs %s (of K = %s)\n",
              x$k1, x$k2, x$K))
  cat(sprintf("  Sampler         : %s\n", x$sampler))
  cat(sprintf("  Conditioning    : %s%s\n", x$condition_on,
              if (identical(x$condition_on, "trajectory")) " (EXPERIMENTAL)" else ""))
  cat(sprintf("  Observed T^2    : %s\n",
              if (is.null(x$T_obs) || !is.finite(x$T_obs)) "NA" else formatC(x$T_obs, digits = 4, format = "g")))
  cat(sprintf("  Selective p-val : %s%s\n",
              format(x$pval, digits = 4),
              if (!is.null(x$se) && is.finite(x$se)) sprintf(" (MC s.e. ~ %.3f)", x$se) else ""))
  invisible(x)
}

#' Summarize a selectiveHotelling object
#'
#' @param object A `"selectiveHotelling"` object from
#'   [selective_hotelling_test()].
#' @param ... Ignored.
#'
#' @return `object`, invisibly (after printing a summary).
#'
#' @export
summary.selectiveHotelling <- function(object, ...) {
  print(object)
  n_draws <- sum(is.finite(object$draws))
  cat(sprintf("  Valid draws     : %d\n", n_draws))
  cat(sprintf("  Unique T^2 draws: %s\n",
              if (is.null(object$n_unique_draws) || is.na(object$n_unique_draws)) "NA"
              else as.character(object$n_unique_draws)))
  cat(sprintf("  Acceptance rate : %s\n",
              if (is.null(object$acc) || is.na(object$acc)) "NA (geodesic sampler)"
              else formatC(object$acc, digits = 3, format = "f")))
  invisible(object)
}
