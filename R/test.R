# ---------- User-facing wrapper & S3 methods -------------------------------

#' Selective Hotelling T-squared test after k-means clustering
#'
#' Runs the full selective inference pipeline: cluster `X` with k-means
#' (`K` clusters, via [kmeans_fun()]), then test the equality of the means of
#' the discovered clusters `k1` and `k2` with a Hotelling \eqn{T^2} statistic
#' whose null distribution is obtained by MCMC sampling of the polar factor `Q`
#' *within the selection region*, treating the covariance \eqn{\Sigma} as
#' unknown. The returned p-value controls the selective Type I error.
#'
#' @param X A numeric data matrix with observations in rows.
#' @param k1,k2 The two cluster labels (in `1:K`) to compare.
#' @param K Number of clusters for the k-means clustering.
#' @param condition_on Selection-event definition. `"final"` (default,
#'   **validated/supported**) conditions only on the final cluster assignments
#'   of the selected pair. `"trajectory"` (**EXPERIMENTAL**) conditions on all
#'   intermediate k-means assignments; its feasible region is much tighter and
#'   the local Cayley walk under-mixes there, so it is not recommended for
#'   inference.
#' @param sampler MCMC sampler to use: `"fast"` (default, fast Cayley
#'   random walk with exact fixed-init feasibility checks), `"har"` (two-row
#'   geodesic hit-and-run / slice), or `"cayley"` (re-clustering Cayley random
#'   walk).
#' @param burn Burn-in iterations.
#' @param B Number of post-burn-in draws to keep.
#' @param thin Thinning factor.
#' @param eps Initial Cayley step size (used by `"fast"`/`"cayley"`; ignored by
#'   `"har"`, which uses geodesic moves).
#' @param ridge Ridge added to the pooled covariance in [hotelling_T2()].
#' @param show_progress Whether to display an MCMC progress bar (requires the
#'   `progress` package).
#' @param verbose Whether to print sampler progress messages.
#'
#' @return An S3 object of class `"selectiveHotelling"`: a list with elements
#'   `pval` (selective p-value), `se` (Monte Carlo standard error), `T_obs`
#'   (observed \eqn{T^2}), `draws` (sampled \eqn{T^2} values), `acc` (acceptance
#'   rate, `NA` for `"har"`), `n_unique_draws`, `k1`, `k2`, `K`, `condition_on`,
#'   and `sampler`.
#'
#' @examples
#' set.seed(1)
#' X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat
#' res <- selective_hotelling_test(
#'   X, 1, 2, K = 3, sampler = "fast",
#'   burn = 200, B = 400, thin = 2
#' )
#' res$pval
#' print(res)
#'
#' @seealso [selective_hotelling_mcmc_fast()], [kmeans_fun()]
#' @export
selective_hotelling_test <- function(X, k1, k2, K = 3,
                                     condition_on = "final",
                                     sampler = c("fast", "har", "cayley"),
                                     burn = 1000, B = 2000, thin = 5,
                                     eps = 0.05, ridge = 1e-8,
                                     show_progress = FALSE, verbose = FALSE) {
  sampler <- match.arg(sampler)
  condition_on <- match.arg(condition_on, c("final", "trajectory"))

  cl_fun <- function(x) kmeans_fun(x, k = K)

  res <- switch(
    sampler,
    fast = selective_hotelling_mcmc_fast(
      X, cl_fun, k1, k2,
      burn = burn, B = B, thin = thin,
      eps = eps, ridge = ridge,
      verbose = verbose, show_progress = show_progress,
      condition_on = condition_on
    ),
    cayley = selective_hotelling_mcmc(
      X, cl_fun, k1, k2,
      burn = burn, B = B, thin = thin,
      eps = eps, ridge = ridge,
      verbose = verbose, show_progress = show_progress,
      condition_on = condition_on
    ),
    har = selective_hitandrun_mcmc(
      X, cl_fun, k1, k2,
      burn = burn, B = B, thin = thin,
      ridge = ridge,
      verbose = verbose, show_progress = show_progress,
      condition_on = condition_on
    )
  )

  out <- list(
    pval = res$pval,
    se = res$se,
    T_obs = res$T_obs,
    draws = res$draws,
    acc = if (!is.null(res$acc)) res$acc else NA_real_,
    n_unique_draws = if (!is.null(res$n_unique_draws)) res$n_unique_draws else NA_integer_,
    k1 = k1,
    k2 = k2,
    K = K,
    condition_on = condition_on,
    sampler = sampler
  )
  class(out) <- "selectiveHotelling"
  out
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
