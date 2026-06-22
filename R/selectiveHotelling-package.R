#' selectiveHotelling: Selective Hotelling T-squared Testing After
#' Clustering with Unknown Covariance
#'
#' Implements valid selective inference for a Hotelling \eqn{T^2} test comparing
#' the means of two clusters that were *discovered* by a data-driven clustering
#' algorithm, while treating the covariance matrix \eqn{\Sigma} as unknown. The
#' package computes Monte Carlo p-values that account for the data-driven
#' selection event by sampling the polar factor \eqn{Q} on the Stiefel manifold
#' within the selection region.
#'
#' @section Clustering algorithms:
#' Two entry points cover different clustering routines:
#'
#' * [selective_hotelling_test()] tests a pair of **k-means** clusters via the
#'   Cayley/hit-and-run MCMC samplers (and can also dispatch to the exact path
#'   test below for **hierarchical** clustering, `clustering = "hclust"`).
#' * [selective_hotelling_bc()] is **algorithm-agnostic**: it takes any
#'   black-box clustering function `cl_fun` and returns an *exact* selective
#'   p-value using the Besag & Clifford (1991) serial Monte Carlo (random-index)
#'   construction on a reversible Stiefel kernel. Validity holds for *any*
#'   reversible kernel, *any* path length, and *any* clustering routine; mixing
#'   only buys power. [hclust_fun()] is the bundled hierarchical-clustering
#'   `cl_fun`.
#'
#' @section Conditioning strategies:
#' The user-facing entry point [selective_hotelling_test()] supports two
#' definitions of the selection event via `condition_on`:
#'
#' * `"final"` (default, **validated/supported**): conditions only on the final
#'   cluster assignments of the selected pair. This is the path with empirical
#'   support for valid selective Type I error control.
#' * `"trajectory"` (**EXPERIMENTAL**): conditions on *all* intermediate
#'   k-means assignments. This region is theoretically appealing (see
#'   Chen & Witten, JMLR 2023) but is *much* tighter, and the local Cayley
#'   random walk under-mixes there, so it is not currently recommended.
#'
#' @section Main functions:
#' * [selective_hotelling_test()] - user-facing wrapper returning a
#'   `"selectiveHotelling"` S3 object (k-means MCMC or hierarchical BC path).
#' * [selective_hotelling_bc()] - algorithm-agnostic exact Besag-Clifford
#'   selective test for an arbitrary clustering `cl_fun`.
#' * [selective_hotelling_mcmc_fast()] - fast Cayley random-walk sampler.
#' * [kmeans_fun()] / [kmeans_with_fixed_init()] - k-means with full trajectory.
#' * [hclust_fun()] - hierarchical clustering `cl_fun` (squared-distance default).
#' * [generate_exchangeable_null_data()] - simulate exchangeable null data.
#' * [hotelling_T2()] - classical two-sample Hotelling \eqn{T^2} statistic.
#'
#' @references
#' Chen, Y. T. and Witten, D. M. (2023). Selective inference for k-means
#' clustering. *Journal of Machine Learning Research*.
#'
#' Besag, J. and Clifford, P. (1991). Sequential Monte Carlo p-values.
#' *Biometrika*, 78(2), 301-304.
#'
#' @keywords internal
#' @importFrom stats cov sd rnorm runif
"_PACKAGE"
