# ---------- Utilities ------------------------------------------------------

#' Symmetric positive-semidefinite matrix square root
#'
#' Symmetrizes a matrix, eigendecomposes it, thresholds tiny or negative
#' eigenvalues at `tol`, and returns both the symmetric square root and inverse
#' square root.
#'
#' @param M A (numerically) symmetric matrix.
#' @param tol Lower bound applied to eigenvalues before taking square roots.
#'
#' @return A list with elements `sqrt` (a symmetric square root) and `invsqrt`
#'   (its inverse).
#'
#' @keywords internal
chol_psd_sqrt <- function(M, tol = 1e-10) {
  # Symmetrize, eigendecompose, threshold tiny/neg eigenvalues, return sqrt and invsqrt
  S  <- (M + t(M))/2
  es <- eigen(S, symmetric = TRUE)
  vals <- pmax(es$values, tol)
  U <- es$vectors
  list(
    sqrt    = U %*% (diag(sqrt(vals), length(vals))) %*% t(U),
    invsqrt = U %*%  diag(1/sqrt(vals), length(vals)) %*% t(U)
  )
}

#' Two-sample Hotelling T-squared statistic
#'
#' Computes the classical two-sample Hotelling \eqn{T^2} statistic comparing the
#' means of the rows of `X` assigned to cluster `k1` versus cluster `k2`, using
#' a pooled (ridge-stabilized) covariance estimate.
#'
#' @param X A numeric data matrix with observations in rows.
#' @param cl Integer/character vector of cluster labels, one per row of `X`.
#' @param k1,k2 The two cluster labels to compare.
#' @param ridge Small non-negative ridge added to the diagonal of the pooled
#'   covariance for numerical stability.
#'
#' @return A single numeric \eqn{T^2} value, or `NA_real_` if either cluster has
#'   fewer than two members.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40), 20, 2)
#' cl <- rep(1:2, each = 10)
#' hotelling_T2(X, cl, 1, 2)
#'
#' @export
hotelling_T2 <- function(X, cl, k1, k2, ridge = 1e-8) {
  X1 <- X[cl == k1, , drop = FALSE]
  X2 <- X[cl == k2, , drop = FALSE]
  n1 <- nrow(X1); n2 <- nrow(X2)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  m1 <- colMeans(X1); m2 <- colMeans(X2)
  S1 <- stats::cov(X1); S2 <- stats::cov(X2)
  Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  Sp <- Sp + diag(ridge, ncol(Sp))
  d <- m1 - m2
  as.numeric((n1 * n2)/(n1 + n2) * t(d) %*% solve(Sp, d))
}

#' Generate exchangeable (equicorrelated) Gaussian null data
#'
#' Draws an `n` by `p` matrix of multivariate Gaussian observations with a
#' common mean `mu` and an exchangeable (equicorrelation) covariance: variance
#' `sig^2` on the diagonal and `sig^2 * rho` off-diagonal. Under this model all
#' observations share the same mean, so it serves as a null for cluster-mean
#' comparisons.
#'
#' @param n Number of observations (rows).
#' @param p Number of variables (columns).
#' @param rho Equicorrelation parameter (off-diagonal correlation).
#' @param mu Common mean added to every entry.
#' @param sig Common standard deviation.
#'
#' @return A list with elements `dat` (the `n` by `p` data matrix) and `Sig`
#'   (the `p` by `p` covariance matrix used).
#'
#' @examples
#' set.seed(1)
#' out <- generate_exchangeable_null_data(30, 2, rho = 0.1)
#' dim(out$dat)
#'
#' @export
generate_exchangeable_null_data <- function(n, p, rho = 0, mu = 0, sig = 1) {
  Sig <- matrix(sig^2 * rho, p, p); diag(Sig) <- sig^2
  Z <- matrix(stats::rnorm(n * p), n, p)
  dat <- Z %*% chol(Sig) + mu
  list(dat = dat, Sig = Sig)
}
