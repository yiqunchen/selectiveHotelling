# ---------- User-facing wrapper & S3 methods -------------------------------

#' Selective Hotelling T-squared test after clustering
#'
#' Runs the full selective inference pipeline: cluster `X` (with k-means or
#' hierarchical clustering, `K` clusters), then test the equality of the means
#' of the discovered clusters `k1` and `k2` with a Hotelling \eqn{T^2}
#' statistic, treating the covariance \eqn{\Sigma} as unknown. The returned
#' p-value controls the selective Type I error.
#'
#' Two clustering routes are supported via `clustering`:
#' \itemize{
#'   \item `"kmeans"` (default): clusters with k-means (via [kmeans_fun()]) and
#'     obtains the null distribution of \eqn{T^2} by MCMC sampling of the polar
#'     factor `Q` *within the k-means selection region*. The MCMC controls
#'     (`sampler`, `condition_on`, `burn`, `B`, `thin`, `eps`) apply here.
#'   \item `"hclust"`: clusters with deterministic agglomerative hierarchical
#'     clustering (via [hclust_fun()]) and routes through the preferred,
#'     **exact** algorithm-agnostic Besag-Clifford engine
#'     ([selective_hotelling_bc()]). The Besag-Clifford controls (`M`, `L`)
#'     apply here; the MCMC-only controls (`sampler`, `condition_on`, `burn`,
#'     `B`, `thin`, `eps`) are **ignored**.
#' }
#'
#' @param X A numeric data matrix with observations in rows.
#' @param k1,k2 The two cluster labels (in `1:K`) to compare.
#' @param K Number of clusters.
#' @param clustering Clustering route: `"kmeans"` (default, MCMC sampler) or
#'   `"hclust"` (hierarchical clustering through the exact Besag-Clifford
#'   engine).
#' @param condition_on Selection-event definition for `clustering = "kmeans"`.
#'   `"final"` (default, **validated/supported**) conditions only on the final
#'   cluster assignments of the selected pair. `"trajectory"`
#'   (**EXPERIMENTAL**) conditions on all intermediate k-means assignments; its
#'   feasible region is much tighter and the local Cayley walk under-mixes
#'   there, so it is not recommended for inference. Ignored for
#'   `clustering = "hclust"`.
#' @param sampler MCMC sampler to use when `clustering = "kmeans"`: `"fast"`
#'   (default, fast Cayley random walk with exact fixed-init feasibility
#'   checks), `"har"` (two-row geodesic hit-and-run / slice), or `"cayley"`
#'   (re-clustering Cayley random walk). Ignored for `clustering = "hclust"`.
#' @param burn Burn-in iterations (k-means/MCMC route only; ignored for
#'   `clustering = "hclust"`).
#' @param B Number of post-burn-in draws to keep (k-means/MCMC route only;
#'   ignored for `clustering = "hclust"`).
#' @param thin Thinning factor (k-means/MCMC route only; ignored for
#'   `clustering = "hclust"`).
#' @param eps Initial Cayley step size (k-means/MCMC route only: used by
#'   `"fast"`/`"cayley"`, ignored by `"har"`; ignored entirely for
#'   `clustering = "hclust"`, which auto-calibrates its own step).
#' @param hclust_method Agglomeration method for `clustering = "hclust"`, passed
#'   to [hclust_fun()] / [stats::hclust()] (e.g. `"average"`, `"single"`,
#'   `"complete"`, `"ward.D2"`). Ignored for `clustering = "kmeans"`.
#' @param hclust_square If `TRUE` (default) the hierarchical clustering tree is
#'   built on squared Euclidean distances (the `clusterpval` convention); see
#'   [hclust_fun()]. Ignored for `clustering = "kmeans"`.
#' @param M Besag-Clifford path length for `clustering = "hclust"` (the path
#'   carries `M + 1` marks). Ignored for `clustering = "kmeans"`.
#' @param L Number of Besag-Clifford kernel steps per path edge for
#'   `clustering = "hclust"`. Ignored for `clustering = "kmeans"`.
#' @param ridge Ridge added to the pooled covariance in [hotelling_T2()].
#' @param show_progress Whether to display an MCMC progress bar (requires the
#'   `progress` package; k-means/MCMC route only).
#' @param verbose Whether to print sampler progress messages.
#'
#' @return An S3 object of class `"selectiveHotelling"`: a list with elements
#'   `pval` (selective p-value), `se` (Monte Carlo standard error), `T_obs`
#'   (observed \eqn{T^2}), `draws` (sampled \eqn{T^2} or path statistic values),
#'   `acc` (acceptance / move rate, `NA` for `"har"`), `n_unique_draws`, `k1`,
#'   `k2`, `K`, `condition_on`, and `sampler` (`"fast"`/`"har"`/`"cayley"` for
#'   k-means, or `"besag-clifford"` for `clustering = "hclust"`).
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
#' # Hierarchical clustering routed through the exact Besag-Clifford engine.
#' # Three separated blobs so hclust recovers three sizeable clusters.
#' set.seed(2)
#' Xb <- rbind(
#'   matrix(rnorm(20), 10, 2),
#'   matrix(rnorm(20) + 6, 10, 2),
#'   matrix(rnorm(20) - 6, 10, 2)
#' )
#' res_h <- selective_hotelling_test(
#'   Xb, 1, 2, K = 3, clustering = "hclust", M = 80
#' )
#' res_h$pval
#'
#' @seealso [selective_hotelling_mcmc_fast()], [selective_hotelling_bc()],
#'   [kmeans_fun()], [hclust_fun()]
#' @export
selective_hotelling_test <- function(X, k1, k2, K = 3,
                                     clustering = c("kmeans", "hclust"),
                                     condition_on = "final",
                                     sampler = c("fast", "har", "cayley"),
                                     burn = 1000, B = 2000, thin = 5,
                                     eps = 0.05,
                                     hclust_method = "average",
                                     hclust_square = TRUE,
                                     M = 300L, L = 2L, ridge = 1e-8,
                                     show_progress = FALSE, verbose = FALSE) {
  clustering <- match.arg(clustering)

  # -------------------------------------------------------------------------
  # HIERARCHICAL route: deterministic hclust black box -> EXACT Besag-Clifford
  # engine. The MCMC-only args (sampler, condition_on, burn, B, thin, eps) do
  # NOT apply here and are ignored. Returns the BC "selectiveHotelling" object
  # directly (it already carries the right class and fields for print/summary).
  # -------------------------------------------------------------------------
  if (clustering == "hclust") {
    cl_fun <- make_hclust_cl_fun(k = K, method = hclust_method,
                                 square = hclust_square)
    return(selective_hotelling_bc(
      X, cl_fun, k1, k2,
      M = M, L = L, ridge = ridge, verbose = verbose
    ))
  }

  # -------------------------------------------------------------------------
  # K-MEANS route (DEFAULT): unchanged MCMC pipeline.
  # -------------------------------------------------------------------------
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
