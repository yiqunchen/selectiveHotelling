# ---------- Main MCMC: uniform within selection via rejection --------------
# Supports both FULL TRAJECTORY conditioning and FINAL-ONLY conditioning.
# Set condition_on = "trajectory" (default, theoretically correct per Chen & Witten)
# or condition_on = "final" (faster, but potentially anti-conservative).

#' Cayley random-walk MCMC sampler for the selective Hotelling test
#'
#' Samples the polar factor `Q` uniformly within the selection region using a
#' Cayley random walk with rejection, re-running k-means (from a fixed
#' initialization) at every step to check feasibility. Used internally by
#' [selective_hotelling_test()] when `sampler = "cayley"`.
#'
#' The `condition_on = "trajectory"` mode (default for this low-level function)
#' is **EXPERIMENTAL**: its feasible region is much tighter and the local Cayley
#' walk under-mixes there. The validated path is `condition_on = "final"`.
#'
#' @param X A numeric data matrix with observations in rows.
#' @param cl_fun A clustering function (e.g. [kmeans_fun()]) returning a list
#'   with `cluster`, `trajectory`, and `init_centers_idx`.
#' @param k1,k2 The two cluster labels defining the tested pair.
#' @param burn Burn-in iterations.
#' @param B Number of post-burn-in draws to keep.
#' @param thin Thinning factor.
#' @param eps Initial Cayley step size.
#' @param adapt_burn Whether to calibrate/adapt `eps` during burn-in.
#' @param ridge Ridge for the pooled covariance in [hotelling_T2()].
#' @param verbose Whether to print progress messages.
#' @param show_progress Whether to show a progress bar (needs the `progress`
#'   package).
#' @param condition_on Either `"trajectory"` (experimental) or `"final"`.
#'
#' @return A list with the selective `pval`, `se`, acceptance diagnostics,
#'   `T_obs`, `draws`, and related fields.
#'
#' @keywords internal
selective_hotelling_mcmc <- function(X, cl_fun, k1, k2,
                                     burn = 1500, B = 5000, thin = 5,
                                     eps = 0.08, adapt_burn = TRUE,
                                     ridge = 1e-8, verbose = TRUE,
                                     show_progress = TRUE,
                                     condition_on = c("trajectory", "final")) {
  condition_on <- match.arg(condition_on)

  if (show_progress && !requireNamespace("progress", quietly = TRUE)) {
    show_progress <- FALSE
  }

  # Run clustering and get full trajectory
  cl_result <- cl_fun(X)
  if (is.list(cl_result)) {
    cl0 <- cl_result$cluster
    traj0 <- cl_result$trajectory
  } else {
    cl0 <- cl_result
    traj0 <- list(cl_result)
  }

  if (!(k1 %in% cl0 && k2 %in% cl0)) stop("k1/k2 not present in clustering")
  prep <- prepare_pair_objects(X, cl0, k1, k2)

  # Re-clustering used for EVERY per-step feasibility check. We reuse the fixed
  # initialization recorded from the original clustering and run k-means with
  # kmeans_with_fixed_init (which never touches the global RNG). This keeps the
  # selection event well-defined across steps AND keeps the chain's proposals
  # random (the old code called cl_fun -> kmeans_with_trajectory -> set.seed()
  # on every step, freezing the walk).
  init_centers_idx <- if (is.list(cl_result)) cl_result$init_centers_idx else NULL
  recluster <- if (!is.null(init_centers_idx)) {
    function(Xnew) kmeans_with_fixed_init(Xnew, init_centers_idx)
  } else {
    cl_fun  # fallback: assume cl_fun does not re-seed the global RNG
  }

  # Define feasibility checker based on conditioning mode
  if (condition_on == "trajectory") {
    check_feasible <- function(cl_result_new) {
      if (is.list(cl_result_new)) {
        preserve_trajectory(traj0, cl_result_new$trajectory, k1, k2)
      } else {
        preserve_cl_strict(cl0, cl_result_new, k1, k2)
      }
    }
  } else {
    # Final-only conditioning
    check_feasible <- function(cl_result_new) {
      cl_new <- if (is.list(cl_result_new)) cl_result_new$cluster else cl_result_new
      preserve_cl_strict(cl0, cl_new, k1, k2)
    }
  }

  with(prep, {
    X_from_Q <- function(Q) {
      AQS <- A %*% Q %*% S_sqrt
      XQ <- X_fixed
      XQ[idx1, ] <- XQ[idx1, ] + AQS[1:m, , drop = FALSE]
      XQ[idx2, ] <- XQ[idx2, ] + AQS[(m+1):r, , drop = FALSE]
      XQ
    }

    # Verify observed data reconstruction preserves selection
    X_obs <- X_from_Q(Q_obs)
    cl_obs_result <- recluster(X_obs)

    if (!check_feasible(cl_obs_result)) {
      stop("Observed reconstruction does not preserve selection event.")
    }

    cl_obs <- if (is.list(cl_obs_result)) cl_obs_result$cluster else cl_obs_result
    T_obs <- hotelling_T2_selected(X_obs, cl0, cl_obs, k1, k2, ridge)
    if (!is.finite(T_obs)) stop("Observed T^2 not finite.")
    if (verbose) {
      cat(sprintf("Observed T^2 = %.4f (conditioning: %s, traj iters: %d)\n",
                  T_obs, condition_on, length(traj0)))
    }

    feasible_Q <- function(Q_candidate) {
      X_candidate <- X_from_Q(Q_candidate)
      check_feasible(recluster(X_candidate))
    }
    pilot <- if (adapt_burn) {
      calibrate_cayley_eps(Q_obs, eps, feasible_Q)
    } else {
      list(eps = eps, pilot_acceptance = NA_real_)
    }
    eps <- pilot$eps

    Q <- Q_obs; acc <- 0L; iters <- 0L
    acc_burn <- 0L; prop_burn <- 0L
    acc_post <- 0L; prop_post <- 0L
    acc_window <- 0L; prop_window <- 0L
    draws <- rep(NA_real_, B); kept <- 0L
    target_acc <- 0.4; adapt_every <- 100L
    total_needed <- burn + B * thin

    if (show_progress) {
      pb <- progress::progress_bar$new(
        format = "  MCMC [:bar] :percent eta: :eta",
        total = total_needed, clear = FALSE, width = 60
      )
      pb$tick(0)
    }

    while (iters < total_needed) {
      iters <- iters + 1L
      Q_try <- cayley_left_step(Q, eps)
      X_try <- X_from_Q(Q_try)
      cl_try_result <- recluster(X_try)

      accepted <- FALSE
      if (check_feasible(cl_try_result)) {
        Q <- Q_try
        acc <- acc + 1L
        accepted <- TRUE
      }

      if (iters <= burn) {
        prop_burn <- prop_burn + 1L
        prop_window <- prop_window + 1L
        if (accepted) {
          acc_burn <- acc_burn + 1L
          acc_window <- acc_window + 1L
        }
      } else {
        prop_post <- prop_post + 1L
        if (accepted) acc_post <- acc_post + 1L
      }

      if (adapt_burn && iters %% adapt_every == 0L && iters <= burn) {
        eps <- adapt_cayley_eps(eps, acc_window / prop_window, target_acc)
        acc_window <- 0L
        prop_window <- 0L
      }

      if (iters > burn && ((iters - burn) %% thin == 0L)) {
        kept <- kept + 1L
        X_keep <- X_from_Q(Q)
        cl_keep_result <- recluster(X_keep)
        cl_keep <- if (is.list(cl_keep_result)) cl_keep_result$cluster else cl_keep_result
        draws[kept] <- hotelling_T2_selected(X_keep, cl0, cl_keep, k1, k2, ridge)
      }
      if (show_progress) pb$tick()
    }

    acc_rate <- acc / iters
    acc_burn_rate <- if (prop_burn > 0L) acc_burn / prop_burn else NA_real_
    acc_post_rate <- if (prop_post > 0L) acc_post / prop_post else NA_real_
    ok <- is.finite(draws)
    survives <- sum(ok)
    n_unique_draws <- length(unique(round(draws[ok], 10)))

    if (survives == 0) {
      warning("No valid draws preserved the selection. Returning NA p-value.")
      return(list(pval = NA_real_, se = NA_real_,
                  acc = acc_rate, kept = kept, survives = survives,
                  n_unique_draws = n_unique_draws,
                  T_obs = T_obs, draws = draws, eps = eps,
                  acc_burn = acc_burn_rate, acc_post = acc_post_rate,
                  pilot_acceptance = pilot$pilot_acceptance,
                  n_iter_traj = length(traj0), condition_on = condition_on))
    }
    if (n_unique_draws <= 1L) {
      warning("Saved MCMC draws have one unique T^2 value; sampler did not explore the conditional law.")
    }

    mc <- mc_tail_pvalue(draws, T_obs)
    pval <- mc$pval
    se <- mc$se

    if (verbose) {
      cat(sprintf("MCMC acceptance %.2f (burn %.2f, post %.2f), kept %d draws, valid %d/%d\n",
                  acc_rate, acc_burn_rate, acc_post_rate, kept, survives, kept))
      cat(sprintf("Selective p-value = %s (MC s.e. \u2248 %s)\n",
                  format(pval, digits=4), ifelse(is.na(se),"NA",sprintf("%.3f",se))))
    }
    list(pval = pval, se = se, acc = acc_rate,
         kept = kept, survives = survives, n_unique_draws = n_unique_draws,
         T_obs = T_obs, draws = draws, eps = eps,
         acc_burn = acc_burn_rate, acc_post = acc_post_rate,
         pilot_acceptance = pilot$pilot_acceptance,
         n_iter_traj = length(traj0), condition_on = condition_on)
  })
}

# ---------- Two-row geodesic Hit-and-Run / Slice extensions ----------------
# Updated to support both trajectory and final-only conditioning

#' Givens rotation acting on two rows of Q
#'
#' Applies a Givens rotation of angle `theta` to rows `i` and `j` of `Q` (a left
#' action on those two rows), used by the geodesic hit-and-run / slice samplers.
#'
#' @param Q An `r` by `p` matrix with orthonormal columns.
#' @param i,j Row indices to rotate.
#' @param theta Rotation angle.
#'
#' @return The rotated matrix.
#'
#' @keywords internal
givens_left <- function(Q, i, j, theta) {
  Qi <- Q[i, , drop = FALSE]
  Qj <- Q[j, , drop = FALSE]
  cth <- cos(theta); sth <- sin(theta)
  Q[i, ] <-  cth * Qi + (-sth) * Qj
  Q[j, ] <-  sth * Qi +  cth  * Qj
  Q
}

#' Feasibility of a Givens-rotated point along a geodesic
#'
#' Rotates `Q0` by `theta` in the `(i, j)` plane, reconstructs the data, and
#' checks whether the selection event is preserved.
#'
#' @param theta Rotation angle.
#' @param Q0 Base Stiefel point.
#' @param i,j Row indices being rotated.
#' @param X_from_Q Function mapping `Q` to a reconstructed data matrix.
#' @param cl_fun Re-clustering function.
#' @param check_feasible Feasibility predicate on a clustering result.
#'
#' @return Logical scalar.
#'
#' @keywords internal
feasible_theta_gen <- function(theta, Q0, i, j, X_from_Q, cl_fun, check_feasible) {
  Q <- givens_left(Q0, i, j, theta)
  XQ <- X_from_Q(Q)
  cl_result <- cl_fun(XQ)
  check_feasible(cl_result)
}

#' Find the feasible arc around theta = 0 in the (i, j) plane
#'
#' Expands outward in both angular directions and bisects to locate the boundary
#' of the single feasible arc containing `theta = 0`.
#'
#' @param Q0 Base Stiefel point.
#' @param i,j Row indices being rotated.
#' @param X_from_Q,cl_fun,check_feasible As in [feasible_theta_gen()].
#' @param step Initial angular step for expansion.
#' @param max_expand Maximum number of expansion steps per direction.
#' @param bisect_tol Bisection tolerance for the boundary.
#'
#' @return A length-2 numeric `c(down, up)` giving the arc endpoints.
#'
#' @keywords internal
find_local_arc_gen <- function(Q0, i, j, X_from_Q, cl_fun, check_feasible,
                               step = 0.05, max_expand = 200, bisect_tol = 1e-3) {
  # Expand in + direction
  up <- 0; k <- 0
  repeat {
    k <- k + 1
    cand <- up + step
    if (k > max_expand || !feasible_theta_gen(cand, Q0, i, j, X_from_Q, cl_fun, check_feasible)) break
    up <- cand
  }
  if (k > 1) {
    L <- up; R <- up + step
    while (R - L > bisect_tol) {
      mid <- 0.5 * (L + R)
      if (feasible_theta_gen(mid, Q0, i, j, X_from_Q, cl_fun, check_feasible)) L <- mid else R <- mid
    }
    up <- L
  }
  # Expand in - direction
  down <- 0; k <- 0
  repeat {
    k <- k + 1
    cand <- down - step
    if (k > max_expand || !feasible_theta_gen(cand, Q0, i, j, X_from_Q, cl_fun, check_feasible)) break
    down <- cand
  }
  if (k > 1) {
    R <- down; L <- down - step
    while (R - L > bisect_tol) {
      mid <- 0.5 * (L + R)
      if (feasible_theta_gen(mid, Q0, i, j, X_from_Q, cl_fun, check_feasible)) R <- mid else L <- mid
    }
    down <- R
  }
  c(down, up)  # single arc containing 0
}

#' One geodesic hit-and-run step over two random rows
#'
#' Picks two random rows, finds the feasible arc in their rotation plane, and
#' draws a uniform angle within that arc.
#'
#' @param Q0 Current Stiefel point.
#' @param X_from_Q,cl_fun,check_feasible As in [feasible_theta_gen()].
#' @param step Initial expansion step passed to [find_local_arc_gen()].
#' @param bisect_tol Bisection tolerance.
#' @param ... Ignored (for interface compatibility).
#'
#' @return The new Stiefel point (or `Q0` if the arc is degenerate).
#'
#' @keywords internal
har_two_row_step_gen <- function(Q0, X_from_Q, cl_fun, check_feasible,
                                 step = 0.05, bisect_tol = 1e-3, ...) {
  r <- nrow(Q0)
  idx <- sample(r, 2, replace = FALSE); i <- min(idx); j <- max(idx)
  arc <- find_local_arc_gen(Q0, i, j, X_from_Q, cl_fun, check_feasible,
                            step = step, bisect_tol = bisect_tol)
  if (arc[1] >= arc[2]) return(Q0)  # degenerate; stay
  theta <- stats::runif(1, arc[1], arc[2])
  givens_left(Q0, i, j, theta)
}

#' One geodesic slice-sampling step over two random rows
#'
#' Picks two random rows and performs a stepping-out / shrinkage slice update on
#' the rotation angle within the feasible region.
#'
#' @param Q0 Current Stiefel point.
#' @param X_from_Q,cl_fun,check_feasible As in [feasible_theta_gen()].
#' @param w Initial slice width.
#' @param max_steps Maximum stepping-out steps per side.
#' @param ... Ignored (for interface compatibility).
#'
#' @return The new Stiefel point (or `Q0` if the slice collapses).
#'
#' @keywords internal
slice_two_row_step_gen <- function(Q0, X_from_Q, cl_fun, check_feasible,
                                   w = 0.2, max_steps = 50, ...) {
  r <- nrow(Q0)
  idx <- sample(r, 2, replace = FALSE); i <- min(idx); j <- max(idx)
  L <- -stats::runif(1, 0, w); R <- L + w
  steps_L <- 0; steps_R <- 0
  while (steps_L < max_steps &&
         feasible_theta_gen(L, Q0, i, j, X_from_Q, cl_fun, check_feasible)) {
    L <- L - w; steps_L <- steps_L + 1
  }
  while (steps_R < max_steps &&
         feasible_theta_gen(R, Q0, i, j, X_from_Q, cl_fun, check_feasible)) {
    R <- R + w; steps_R <- steps_R + 1
  }
  repeat {
    theta <- stats::runif(1, L, R)
    if (feasible_theta_gen(theta, Q0, i, j, X_from_Q, cl_fun, check_feasible)) {
      return(givens_left(Q0, i, j, theta))
    } else {
      if (theta < 0) L <- theta else R <- theta
      if (R - L < 1e-6) return(Q0)
    }
  }
}

#' Two-row geodesic hit-and-run / slice MCMC for the selective Hotelling test
#'
#' Samples `Q` within the selection region using two-row geodesic hit-and-run
#' (or slice) moves. Used internally by [selective_hotelling_test()] when
#' `sampler = "har"`.
#'
#' As with the other samplers, `condition_on = "trajectory"` is **EXPERIMENTAL**;
#' the validated path is `condition_on = "final"`.
#'
#' @param X A numeric data matrix with observations in rows.
#' @param cl_fun A clustering function (e.g. [kmeans_fun()]).
#' @param k1,k2 The two cluster labels defining the tested pair.
#' @param burn Burn-in iterations.
#' @param B Number of post-burn-in draws to keep.
#' @param thin Thinning factor.
#' @param use_slice If `TRUE`, use slice steps; otherwise hit-and-run steps.
#' @param step Initial expansion step for hit-and-run.
#' @param w Initial slice width.
#' @param ridge Ridge for the pooled covariance in [hotelling_T2()].
#' @param verbose Whether to print progress messages.
#' @param show_progress Whether to show a progress bar.
#' @param condition_on Either `"trajectory"` (experimental) or `"final"`.
#'
#' @return A list with the selective `pval`, `se`, `T_obs`, `draws`, and move
#'   diagnostics.
#'
#' @keywords internal
selective_hitandrun_mcmc <- function(X, cl_fun, k1, k2,
                                     burn = 1000, B = 5000, thin = 5,
                                     use_slice = TRUE,
                                     step = 0.05, w = 0.2, ridge = 1e-8,
                                     verbose = TRUE, show_progress = TRUE,
                                     condition_on = c("trajectory", "final")) {
  condition_on <- match.arg(condition_on)

  if (show_progress && !requireNamespace("progress", quietly = TRUE)) {
    show_progress <- FALSE
  }

  # Run clustering and get trajectory
  cl_result <- cl_fun(X)
  if (is.list(cl_result)) {
    cl0 <- cl_result$cluster
    traj0 <- cl_result$trajectory
  } else {
    cl0 <- cl_result
    traj0 <- list(cl_result)
  }

  prep <- prepare_pair_objects(X, cl0, k1, k2)

  # Fixed-init re-clustering for per-step feasibility checks (loop-safe, never
  # re-seeds the global RNG). See note in selective_hotelling_mcmc.
  init_centers_idx <- if (is.list(cl_result)) cl_result$init_centers_idx else NULL
  recluster <- if (!is.null(init_centers_idx)) {
    function(Xnew) kmeans_with_fixed_init(Xnew, init_centers_idx)
  } else {
    cl_fun
  }

  # Define feasibility checker
  if (condition_on == "trajectory") {
    check_feasible <- function(cl_result_new) {
      if (is.list(cl_result_new)) {
        preserve_trajectory(traj0, cl_result_new$trajectory, k1, k2)
      } else {
        preserve_cl_strict(cl0, cl_result_new, k1, k2)
      }
    }
  } else {
    check_feasible <- function(cl_result_new) {
      cl_new <- if (is.list(cl_result_new)) cl_result_new$cluster else cl_result_new
      preserve_cl_strict(cl0, cl_new, k1, k2)
    }
  }

  with(prep, {
    X_from_Q <- function(Q) {
      AQS <- A %*% Q %*% S_sqrt
      XQ <- X_fixed
      XQ[idx1, ] <- XQ[idx1, ] + AQS[1:m, , drop = FALSE]
      XQ[idx2, ] <- XQ[idx2, ] + AQS[(m+1):r, , drop = FALSE]
      XQ
    }

    # Verify observed reconstruction
    X_obs <- X_from_Q(Q_obs)
    cl_obs_result <- recluster(X_obs)
    if (!check_feasible(cl_obs_result)) {
      stop("Observed reconstruction does not preserve selection event.")
    }
    cl_obs <- if (is.list(cl_obs_result)) cl_obs_result$cluster else cl_obs_result
    T_obs <- hotelling_T2_selected(X_obs, cl0, cl_obs, k1, k2, ridge)
    if (!is.finite(T_obs)) stop("Observed T^2 not finite.")

    iters <- 0L; kept <- 0L; moves <- 0L
    draws <- rep(NA_real_, B)
    total_needed <- burn + B * thin

    # Use generalized step functions with check_feasible
    step_fun <- if (use_slice) {
      function(Q, ...) slice_two_row_step_gen(Q, X_from_Q, recluster, check_feasible, w = w, ...)
    } else {
      function(Q, ...) har_two_row_step_gen(Q, X_from_Q, recluster, check_feasible, step = step, ...)
    }
    Q <- Q_obs

    # progress bar
    if (show_progress) {
      pb <- progress::progress_bar$new(
        format = "  HAR   [:bar] :percent eta: :eta",
        total = total_needed, clear = FALSE, width = 60
      )
      pb$tick(0)
    }

    # helper to count a "move"
    changed_Q <- function(Q_old, Q_new) {
      max(abs(Q_old - Q_new)) > 1e-12
    }

    if (verbose) {
      cat(sprintf("Starting HAR MCMC (conditioning: %s, traj iters: %d)\n",
                  condition_on, length(traj0)))
    }

    while (iters < total_needed) {
      iters <- iters + 1L
      Q_prev <- Q
      Q <- step_fun(Q)
      if (changed_Q(Q_prev, Q)) moves <- moves + 1L

      if (iters > burn && ((iters - burn) %% thin == 0L)) {
        kept <- kept + 1L
        X_keep <- X_from_Q(Q)
        cl_keep_result <- recluster(X_keep)
        cl_keep <- if (is.list(cl_keep_result)) cl_keep_result$cluster else cl_keep_result
        draws[kept] <- hotelling_T2_selected(X_keep, cl0, cl_keep, k1, k2, ridge)
      }

      if (show_progress) pb$tick()
    }

    valid <- sum(is.finite(draws))
    n_unique_draws <- length(unique(round(draws[is.finite(draws)], 10)))
    if (valid == 0) {
      warning("No valid (finite) T^2 among saved draws. Returning NA p-value.")
      if (verbose) {
        cat(sprintf("HAR summary: moves %d/%d, kept %d, valid %d/%d\n",
                    moves, iters, kept, valid, kept))
      }
      return(list(pval = NA_real_, se = NA_real_,
                  moves = moves, iters = iters,
                  kept = kept, valid = valid, n_unique_draws = n_unique_draws,
                  T_obs = T_obs, draws = draws,
                  n_iter_traj = length(traj0), condition_on = condition_on))
    }
    if (n_unique_draws <= 1L) {
      warning("Saved HAR draws have one unique T^2 value; sampler did not explore the conditional law.")
    }

    mc <- mc_tail_pvalue(draws, T_obs)
    pval <- mc$pval
    se <- mc$se

    if (verbose) {
      cat(sprintf("HAR summary: moves %d/%d, kept %d, valid %d/%d\n",
                  moves, iters, kept, valid, kept))
      cat(sprintf("Selective p-value = %s (MC s.e. \u2248 %s)\n",
                  format(pval, digits=4), ifelse(is.na(se),"NA",sprintf("%.3f",se))))
    }
    list(pval = pval, se = se,
         moves = moves, iters = iters,
         kept = kept, valid = valid, n_unique_draws = n_unique_draws,
         T_obs = T_obs, draws = draws,
         n_iter_traj = length(traj0), condition_on = condition_on)
  })
}

#' Fast Cayley random-walk MCMC for the selective Hotelling test
#'
#' Fast sampler that draws `Q` via a Cayley random walk and checks feasibility
#' with exact fixed-initialization k-means at each step. This is the sampler
#' dispatched by [selective_hotelling_test()] for `sampler = "fast"`.
#'
#' For `condition_on = "final"` (the validated path), feasibility means the
#' selected pair's *final* clusters are preserved, and `T^2` is recomputed from
#' the re-clustered data at each kept draw. For `condition_on = "trajectory"`
#' (**EXPERIMENTAL**), the *entire* pairwise membership trajectory must be
#' preserved; this region is much tighter and the local Cayley walk under-mixes,
#' so it is not recommended for inference.
#'
#' @param X A numeric data matrix with observations in rows.
#' @param cl_fun A clustering function returning `cluster`, `trajectory`, and
#'   `init_centers_idx` (e.g. [kmeans_fun()]). The fast sampler *requires*
#'   `init_centers_idx`.
#' @param k1,k2 The two cluster labels defining the tested pair.
#' @param burn Burn-in iterations.
#' @param B Number of post-burn-in draws to keep.
#' @param thin Thinning factor.
#' @param eps Initial Cayley step size.
#' @param adapt_burn Whether to calibrate/adapt `eps` during burn-in.
#' @param ridge Ridge for the pooled covariance in [hotelling_T2()].
#' @param verbose Whether to print progress messages.
#' @param show_progress Whether to show a progress bar.
#' @param condition_on Either `"final"` (validated) or `"trajectory"`
#'   (experimental).
#'
#' @return A list with the selective `pval`, `se`, acceptance diagnostics,
#'   `T_obs`, `draws`, timing, and related fields.
#'
#' @examples
#' set.seed(1)
#' X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat
#' cl_fun <- function(x) kmeans_fun(x, k = 3)
#' res <- selective_hotelling_mcmc_fast(
#'   X, cl_fun, 1, 2,
#'   burn = 200, B = 400, thin = 2,
#'   condition_on = "final", verbose = FALSE, show_progress = FALSE
#' )
#' res$pval
#'
#' @export
selective_hotelling_mcmc_fast <- function(X, cl_fun, k1, k2,
                                          burn = 1500, B = 5000, thin = 5,
                                          eps = 0.08, adapt_burn = TRUE,
                                          ridge = 1e-8, verbose = TRUE,
                                          show_progress = TRUE,
                                          condition_on = c("trajectory", "final")) {
  condition_on <- match.arg(condition_on)

  # Run clustering to get trajectory
  cl_result <- cl_fun(X)
  if (is.list(cl_result)) {
    cl0 <- cl_result$cluster
    traj0 <- cl_result$trajectory
    init_centers_idx <- cl_result$init_centers_idx
  } else {
    cl0 <- cl_result
    traj0 <- list(cl_result)
    init_centers_idx <- NULL
  }

  K <- length(unique(cl0))
  if (!(k1 %in% cl0 && k2 %in% cl0)) stop("k1/k2 not present in clustering")
  if (is.null(init_centers_idx)) {
    stop("Fast sampler requires cl_fun to return init_centers_idx.")
  }

  prep <- prepare_pair_objects(X, cl0, k1, k2)

  with(prep, {
    X_from_Q <- function(Q) {
      AQS <- A %*% Q %*% S_sqrt
      XQ <- X_fixed
      XQ[idx1, ] <- XQ[idx1, ] + AQS[1:m, , drop = FALSE]
      XQ[idx2, ] <- XQ[idx2, ] + AQS[(m+1):r, , drop = FALSE]
      XQ
    }

    if (verbose) {
      cat("Using exact fixed-init k-means checks for selected-pair conditioning...\n")
    }
    t_precomp <- system.time({
      constraints_fast <- NULL
    })
    if (verbose) {
      cat(sprintf("Precomputed %d constraints in %.2f sec\n",
                  if (is.null(constraints_fast)) 0L else constraints_fast$n_constraints,
                  t_precomp[3]))
    }

    check_feasible_fast <- function(Q) {
      X_candidate <- X_from_Q(Q)
      if (condition_on == "trajectory") {
        cl_candidate <- kmeans_with_fixed_init(X_candidate, init_centers_idx)
        preserve_trajectory(traj0, cl_candidate$trajectory, k1, k2)
      } else {
        cl_candidate <- kmeans_with_fixed_init(X_candidate, init_centers_idx)$cluster
        preserve_cl_strict(cl0, cl_candidate, k1, k2)
      }
    }

    # Verify observed reconstruction preserves selection
    if (!check_feasible_fast(Q_obs)) {
      stop("Observed reconstruction does not preserve selection event.")
    }

    X_obs <- X_from_Q(Q_obs)
    T_obs <- hotelling_T2_selected(X_obs, cl0, cl0, k1, k2, ridge)
    if (!is.finite(T_obs)) stop("Observed T^2 not finite.")
    if (verbose) {
      cat(sprintf("Observed T^2 = %.4f (conditioning: %s, traj iters: %d)\n",
                  T_obs, condition_on, length(traj0)))
    }

    pilot <- if (adapt_burn) {
      calibrate_cayley_eps(Q_obs, eps, check_feasible_fast)
    } else {
      list(eps = eps, pilot_acceptance = NA_real_)
    }
    eps <- pilot$eps

    # MCMC with fast constraint checking
    Q <- Q_obs; acc <- 0L; iters <- 0L
    acc_burn <- 0L; prop_burn <- 0L
    acc_post <- 0L; prop_post <- 0L
    acc_window <- 0L; prop_window <- 0L
    draws <- rep(NA_real_, B); kept <- 0L
    target_acc <- 0.4; adapt_every <- 100L
    total_needed <- burn + B * thin

    if (show_progress && requireNamespace("progress", quietly = TRUE)) {
      pb <- progress::progress_bar$new(
        format = "  MCMC-fast [:bar] :percent eta: :eta",
        total = total_needed, clear = FALSE, width = 60
      )
      pb$tick(0)
    } else {
      pb <- NULL
    }

    t_mcmc <- system.time({
      while (iters < total_needed) {
        iters <- iters + 1L
        Q_try <- cayley_left_step(Q, eps)

        # FAST constraint check (no k-means!)
        accepted <- FALSE
        if (check_feasible_fast(Q_try)) {
          Q <- Q_try
          acc <- acc + 1L
          accepted <- TRUE
        }

        if (iters <= burn) {
          prop_burn <- prop_burn + 1L
          prop_window <- prop_window + 1L
          if (accepted) {
            acc_burn <- acc_burn + 1L
            acc_window <- acc_window + 1L
          }
        } else {
          prop_post <- prop_post + 1L
          if (accepted) acc_post <- acc_post + 1L
        }

        if (adapt_burn && iters %% adapt_every == 0L && iters <= burn) {
          eps <- adapt_cayley_eps(eps, acc_window / prop_window, target_acc)
          acc_window <- 0L
          prop_window <- 0L
        }

        if (iters > burn && ((iters - burn) %% thin == 0L)) {
          kept <- kept + 1L
          X_keep <- X_from_Q(Q)
          cl_keep <- if (condition_on == "trajectory") {
            cl0
          } else {
            kmeans_with_fixed_init(X_keep, init_centers_idx)$cluster
          }
          draws[kept] <- hotelling_T2_selected(X_keep, cl0, cl_keep, k1, k2, ridge)
        }
        if (!is.null(pb)) pb$tick()
      }
    })

    acc_rate <- acc / iters
    acc_burn_rate <- if (prop_burn > 0L) acc_burn / prop_burn else NA_real_
    acc_post_rate <- if (prop_post > 0L) acc_post / prop_post else NA_real_
    ok <- is.finite(draws)
    survives <- sum(ok)
    n_unique_draws <- length(unique(round(draws[ok], 10)))

    if (survives == 0) {
      warning("No valid draws preserved the selection. Returning NA p-value.")
      return(list(pval = NA_real_, se = NA_real_,
                  acc = acc_rate, kept = kept, survives = survives,
                  n_unique_draws = n_unique_draws,
                  T_obs = T_obs, draws = draws, eps = eps,
                  acc_burn = acc_burn_rate, acc_post = acc_post_rate,
                  pilot_acceptance = pilot$pilot_acceptance,
                  n_iter_traj = length(traj0), condition_on = condition_on,
                  time_precomp = t_precomp[3], time_mcmc = t_mcmc[3]))
    }
    if (n_unique_draws <= 1L) {
      warning("Saved fast-MCMC draws have one unique T^2 value; sampler did not explore the conditional law.")
    }

    mc <- mc_tail_pvalue(draws, T_obs)
    pval <- mc$pval
    se <- mc$se

    if (verbose) {
      cat(sprintf("MCMC: %.1f sec (precomp: %.2f sec), acc %.2f (burn %.2f, post %.2f), kept %d, valid %d/%d\n",
                  t_mcmc[3], t_precomp[3], acc_rate, acc_burn_rate,
                  acc_post_rate, kept, survives, kept))
      cat(sprintf("Selective p-value = %s (MC s.e. \u2248 %s)\n",
                  format(pval, digits=4), ifelse(is.na(se),"NA",sprintf("%.3f",se))))
    }

    list(pval = pval, se = se, acc = acc_rate,
         kept = kept, survives = survives, n_unique_draws = n_unique_draws,
         T_obs = T_obs, draws = draws, eps = eps,
         acc_burn = acc_burn_rate, acc_post = acc_post_rate,
         pilot_acceptance = pilot$pilot_acceptance,
         n_iter_traj = length(traj0), condition_on = condition_on,
         time_precomp = t_precomp[3], time_mcmc = t_mcmc[3],
         n_constraints = if (is.null(constraints_fast)) 0L else constraints_fast$n_constraints,
         fast_event = if (condition_on == "trajectory") {
           "exact_fixed_init_pair_membership_trajectory"
         } else {
           "fixed_init_final_selected_clusters"
         })
  })
}
