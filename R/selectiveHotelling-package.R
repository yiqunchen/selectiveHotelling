#' selectiveHotelling: Selective Hotelling T-squared Testing After k-Means
#' Clustering with Unknown Covariance
#'
#' Implements valid selective inference for a Hotelling \eqn{T^2} test comparing
#' the means of two clusters that were *discovered* by k-means clustering, while
#' treating the covariance matrix \eqn{\Sigma} as unknown. The package computes
#' Monte Carlo p-values that account for the data-driven selection event by
#' sampling uniformly (within the selection region) from the polar factor
#' \eqn{Q} on the Stiefel manifold.
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
#'   `"selectiveHotelling"` S3 object.
#' * [selective_hotelling_mcmc_fast()] - fast Cayley random-walk sampler.
#' * [kmeans_fun()] / [kmeans_with_fixed_init()] - k-means with full trajectory.
#' * [generate_exchangeable_null_data()] - simulate exchangeable null data.
#' * [hotelling_T2()] - classical two-sample Hotelling \eqn{T^2} statistic.
#'
#' @references
#' Chen, Y. T. and Witten, D. M. (2023). Selective inference for k-means
#' clustering. *Journal of Machine Learning Research*.
#'
#' @keywords internal
#' @importFrom stats cov sd rnorm runif
"_PACKAGE"
