#' selectiveHotelling: Selective Hotelling T-squared Testing After
#' Clustering with Unknown Covariance
#'
#' Tests whether two clusters found by a clustering algorithm have the same mean,
#' with Hotelling's \eqn{T^2} and an unknown covariance matrix \eqn{\Sigma}. The
#' p-values condition on the clustering and on the statistics that remove the
#' unknown means and \eqn{\Sigma}, so they control the selective Type I error.
#'
#' @section Main functions:
#' * [selective_hotelling_test()] - clusters the data (k-means or hierarchical
#'   clustering) and runs the main test.
#' * [selective_hotelling_bc()] - the main test for any clustering function
#'   `cl_fun`: the Besag and Clifford (1989) construction with the random Givens
#'   kernel on the event that the clustering returns the tested pair. The p-value
#'   is exactly uniform under the null for any path length and any reversible
#'   kernel; the kernel and the path length affect only power.
#' * [selective_truncated_F()] - the one-dimensional test, which also fixes all of
#'   the polar factor except one angle. It is computed by the same construction on
#'   the angle for any clustering function, and in closed form for k-means when
#'   the conditioning includes the whole Lloyd trajectory.
#' * [kmeans_fun()] / [kmeans_with_fixed_init()] - k-means with its trajectory.
#' * [hclust_fun()] - hierarchical clustering (squared-distance default).
#' * [hotelling_T2()] - the two-sample Hotelling \eqn{T^2} statistic.
#' * [generate_exchangeable_null_data()] - simulate null data.
#'
#' @references
#' Besag, J. and Clifford, P. (1989). Generalized Monte Carlo significance tests.
#' *Biometrika*, 76(4), 633-642.
#'
#' Chen, Y. T. and Witten, D. M. (2023). Selective inference for k-means
#' clustering. *Journal of Machine Learning Research*.
#'
#' Gao, L. L., Bien, J. and Witten, D. (2024). Selective inference for
#' hierarchical clustering. *Journal of the American Statistical Association*.
#'
#' @keywords internal
#' @importFrom stats cov sd rnorm runif
"_PACKAGE"
