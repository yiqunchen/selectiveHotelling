# ---------- K-means with full trajectory tracking ---------------------------
# Critical: We must condition on ALL intermediate cluster assignments,
# not just the final clustering. See Chen & Witten (JMLR 2023).

#' K-means recording the full assignment trajectory
#'
#' Runs Lloyd's k-means while recording the cluster assignment vector at *every*
#' iteration (the "trajectory"), which is needed for trajectory-based selective
#' conditioning.
#'
#' This function does **not** silently re-seed the global RNG. Reproducibility is
#' the caller's responsibility: seed once at the call site, or pass an explicit
#' `init_centers_idx`. If `seed` is supplied it is applied *locally* and the
#' global RNG stream is restored on exit, so calling this inside an MCMC loop can
#' never freeze the chain's proposals.
#'
#' @param x A numeric data matrix with observations in rows.
#' @param k Number of clusters.
#' @param max_iter Maximum number of Lloyd iterations.
#' @param init_centers_idx Optional integer vector of length `k` giving row
#'   indices of `x` to use as initial centers. If `NULL`, centers are chosen by
#'   `init_method`.
#' @param seed Optional integer seed applied locally (global RNG restored on
#'   exit).
#' @param init_method Initialization scheme when `init_centers_idx` is `NULL`:
#'   `"random_indices"` (default, data-independent) or `"kmeans++"`.
#'
#' @return A list with elements `cluster` (final assignments), `trajectory` (a
#'   list of per-iteration assignment vectors), `center_history`,
#'   `init_centers_idx`, and `n_iter`.
#'
#' @keywords internal
kmeans_with_trajectory <- function(x, k = 3, max_iter = 50,
                                   init_centers_idx = NULL, seed = NULL,
                                   init_method = c("random_indices", "kmeans++")) {
  # K-means that records cluster assignments at EVERY iteration.
  # Returns list with $cluster (final) and $trajectory (list of assignments).
  #
  # RNG policy (IMPORTANT): this function does NOT silently re-seed the global
  # RNG. Reproducibility is the CALLER's responsibility -- seed once at the call
  # site, or pass an explicit `init_centers_idx`. If `seed` is supplied we seed
  # LOCALLY and RESTORE the global RNG stream on exit, so calling this inside an
  # MCMC loop can never freeze the chain's proposals (the old footgun that broke
  # calibration: a per-step set.seed() pinned every Cayley proposal to one fixed
  # rotation -> deterministic, non-mixing chain).
  init_method <- match.arg(init_method)
  n <- nrow(x); p <- ncol(x)

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      .saved_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .saved_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  # Default initialization is independent of X: condition on these external
  # random indices, not on data-dependent k-means++ probabilities.
  if (is.null(init_centers_idx)) {
    if (init_method == "random_indices") {
      init_centers_idx <- sample.int(n, k, replace = FALSE)
    } else {
      init_centers_idx <- integer(k)
      init_centers_idx[1] <- sample.int(n, 1)
      for (j in 2:k) {
        dists_to_centers <- apply(x[init_centers_idx[1:(j-1)], , drop = FALSE], 1,
                                  function(ctr) rowSums((x - rep(ctr, each = n))^2))
        if (is.matrix(dists_to_centers)) {
          min_dists <- apply(dists_to_centers, 1, min)
        } else {
          min_dists <- dists_to_centers
        }
        probs <- min_dists / sum(min_dists)
        init_centers_idx[j] <- sample.int(n, 1, prob = probs)
      }
    }
  }
  centers_idx <- init_centers_idx
  centers <- x[centers_idx, , drop = FALSE]

  trajectory <- vector("list", max_iter)
  center_history <- vector("list", max_iter)
  cl_prev <- NULL

  for (t in seq_len(max_iter)) {
    center_history[[t]] <- centers

    # Assignment step: assign each point to nearest center
    dists <- matrix(0, n, k)
    for (j in 1:k) {
      dists[, j] <- rowSums((x - rep(centers[j, ], each = n))^2)
    }
    cl <- apply(dists, 1, which.min)
    trajectory[[t]] <- cl

    # Check for convergence (assignments unchanged)
    if (!is.null(cl_prev) && identical(cl, cl_prev)) {
      trajectory <- trajectory[1:t]  # trim unused slots
      center_history <- center_history[1:t]
      break
    }
    cl_prev <- cl

    # Update step: recompute centers
    for (j in 1:k) {
      members <- which(cl == j)
      if (length(members) > 0) {
        centers[j, ] <- colMeans(x[members, , drop = FALSE])
      }
    }
  }

  list(cluster = cl, trajectory = trajectory, center_history = center_history,
       init_centers_idx = centers_idx, n_iter = length(trajectory))
}

#' Trajectory-recording k-means with a simple interface
#'
#' Convenience wrapper around [kmeans_with_trajectory()] that runs k-means with
#' `k` clusters and returns the full trajectory object. This is the clustering
#' function expected by [selective_hotelling_test()] and the MCMC samplers
#' (because it returns `init_centers_idx`, enabling loop-safe fixed-init
#' re-clustering).
#'
#' @param x A numeric data matrix with observations in rows.
#' @param k Number of clusters.
#'
#' @return The list returned by [kmeans_with_trajectory()] (with `cluster`,
#'   `trajectory`, `init_centers_idx`, etc.).
#'
#' @examples
#' set.seed(1)
#' X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat
#' km <- kmeans_fun(X, k = 3)
#' table(km$cluster)
#'
#' @export
kmeans_fun <- function(x, k = 3) {
  kmeans_with_trajectory(x, k = k)
}

#' Trajectory-recording k-means returning only the final labels
#'
#' Backward-compatible wrapper around [kmeans_with_trajectory()] that returns
#' just the final cluster-assignment vector.
#'
#' @param x A numeric data matrix with observations in rows.
#' @param k Number of clusters.
#'
#' @return An integer vector of final cluster labels.
#'
#' @keywords internal
kmeans_fun_simple <- function(x, k = 3) {
  kmeans_with_trajectory(x, k = k)$cluster
}

#' K-means from a fixed initialization
#'
#' Runs Lloyd's k-means from an explicitly supplied set of initial center
#' indices, recording the full assignment trajectory. Because the initialization
#' is fixed and no RNG is used, this is the loop-safe re-clustering routine used
#' for per-step feasibility checks inside the MCMC samplers.
#'
#' @param x A numeric data matrix with observations in rows.
#' @param init_centers_idx Integer vector of row indices of `x` used as initial
#'   centers; its length determines `k`.
#' @param max_iter Maximum number of Lloyd iterations.
#'
#' @return A list with elements `cluster`, `trajectory`, `center_history`,
#'   `init_centers_idx`, and `n_iter`.
#'
#' @examples
#' set.seed(1)
#' X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat
#' km <- kmeans_with_fixed_init(X, init_centers_idx = c(1, 11, 21))
#' km$n_iter
#'
#' @export
kmeans_with_fixed_init <- function(x, init_centers_idx, max_iter = 50) {
  # Pure-R Lloyd's k-means from a fixed initialization with full trajectory
  # tracking, implemented by kmeans_with_fixed_init_R(). Consumes NO randomness,
  # so the loop-safe / RNG-safety contract is preserved (per-step feasibility
  # checks inside the MCMC samplers never touch the global RNG stream).
  #
  # NOTE: an earlier compiled (Rcpp) variant was removed. It was bit-for-bit
  # identical to this R code but gave NO end-to-end speedup for the selective
  # test -- inside the sampler the re-clustering of a tiny perturbation converges
  # in ~1 warm iteration, so k-means is a negligible fraction of a sampler step
  # (measured ~1.00x at p = 2, 10, 50). The package is therefore intentionally
  # pure R and needs no C++ toolchain.
  if (!is.matrix(x)) x <- as.matrix(x)
  kmeans_with_fixed_init_R(x, init_centers_idx, max_iter)
}

#' K-means from a fixed initialization (implementation)
#'
#' Pure-R implementation backing [kmeans_with_fixed_init()]. Runs Lloyd's
#' k-means from a fixed initialization with full trajectory tracking.
#'
#' @param x A numeric data matrix with observations in rows.
#' @param init_centers_idx Integer vector of row indices of `x` used as initial
#'   centers; its length determines `k`.
#' @param max_iter Maximum number of Lloyd iterations.
#'
#' @return A list with elements `cluster`, `trajectory`, `center_history`,
#'   `init_centers_idx`, and `n_iter` (identical in shape to
#'   [kmeans_with_fixed_init()]).
#'
#' @keywords internal
kmeans_with_fixed_init_R <- function(x, init_centers_idx, max_iter = 50) {
  n <- nrow(x)
  k <- length(init_centers_idx)
  centers <- x[init_centers_idx, , drop = FALSE]

  trajectory <- vector("list", max_iter)
  center_history <- vector("list", max_iter)
  cl_prev <- NULL

  for (t in seq_len(max_iter)) {
    center_history[[t]] <- centers

    dists <- matrix(0, n, k)
    for (j in seq_len(k)) {
      dists[, j] <- rowSums((x - rep(centers[j, ], each = n))^2)
    }
    cl <- apply(dists, 1, which.min)
    trajectory[[t]] <- cl

    if (!is.null(cl_prev) && identical(cl, cl_prev)) {
      trajectory <- trajectory[seq_len(t)]
      center_history <- center_history[seq_len(t)]
      break
    }
    cl_prev <- cl

    for (j in seq_len(k)) {
      members <- which(cl == j)
      if (length(members) > 0) {
        centers[j, ] <- colMeans(x[members, , drop = FALSE])
      }
    }
  }

  list(cluster = cl, trajectory = trajectory, center_history = center_history,
       init_centers_idx = init_centers_idx, n_iter = length(trajectory))
}

# ---------- Trajectory-based cluster preservation check ---------------------

#' Check preservation of the full k-means trajectory for a cluster pair
#'
#' Returns `TRUE` only if, at *every* k-means iteration, the observations that
#' were in clusters `k1` and `k2` of the original trajectory `traj0` form exactly
#' two distinct clusters in the candidate trajectory `traj1`. This is the
#' (experimental) trajectory-conditioning feasibility check.
#'
#' @param traj0,traj1 Lists of per-iteration cluster-assignment vectors.
#' @param k1,k2 The two original cluster labels being tested.
#'
#' @return Logical scalar.
#'
#' @keywords internal
preserve_trajectory <- function(traj0, traj1, k1, k2) {
  # Check that ALL intermediate cluster assignments are preserved for clusters k1, k2.
  # traj0, traj1: lists of cluster assignment vectors (one per k-means iteration)
  # Returns TRUE only if every iteration's assignments match for the tested clusters.

  # Must have same number of iterations
  if (length(traj0) != length(traj1)) return(FALSE)

  # Check each iteration
  for (t in seq_along(traj0)) {
    cl0_t <- traj0[[t]]
    cl1_t <- traj1[[t]]

    # Find which observations were in k1, k2 at iteration t in original
    i1 <- which(cl0_t == k1)
    i2 <- which(cl0_t == k2)

    # In the new trajectory, these same observations must form exactly two clusters
    L1 <- unique(cl1_t[i1])
    L2 <- unique(cl1_t[i2])

    # Each set must map to exactly one cluster label
    if (length(L1) != 1 || length(L2) != 1) return(FALSE)

    # The two clusters must be different
    if (L1 == L2) return(FALSE)

    # All obs in new cluster L1 must be exactly i1, and similarly for L2
    if (!all(cl1_t[i1] == L1)) return(FALSE)
    if (!all(cl1_t[i2] == L2)) return(FALSE)
    if (!all(which(cl1_t == L1) %in% i1)) return(FALSE)
    if (!all(which(cl1_t == L2) %in% i2)) return(FALSE)
  }
  TRUE
}

#' Check preservation of the final cluster pair
#'
#' Legacy/final-only feasibility check: requires that the observations in
#' original clusters `k1` and `k2` form exactly those two clusters (up to
#' relabeling) in the candidate final assignment `cl1`.
#'
#' @param cl0,cl1 Original and candidate final cluster-assignment vectors.
#' @param k1,k2 The two original cluster labels being tested.
#'
#' @return Logical scalar.
#'
#' @keywords internal
preserve_cl_strict <- function(cl0, cl1, k1, k2) {
  # Require exact membership preservation for the two tested clusters
  # NOTE: This only checks final clusters - use preserve_trajectory for full conditioning!
  i1 <- which(cl0 == k1); i2 <- which(cl0 == k2)
  L1 <- unique(cl1[i1]); L2 <- unique(cl1[i2])
  if (length(L1) != 1 || length(L2) != 1) return(FALSE)
  all(cl1[i1] == L1) && all(cl1[i2] == L2) &&
    all(which(cl1 == L1) %in% i1) && all(which(cl1 == L2) %in% i2)
}

#' Map originally selected clusters to their current labels
#'
#' Feasibility checks permit label permutations, so downstream \eqn{T^2}
#' calculations must use the *mapped* labels rather than reusing `k1`/`k2`. This
#' returns the current labels of the originally selected clusters, or `NULL` if
#' the pair is not preserved.
#'
#' @param cl0,cl1 Original and candidate cluster-assignment vectors.
#' @param k1,k2 The two original cluster labels being tested.
#'
#' @return A named length-2 vector `c(k1 = L1, k2 = L2)`, or `NULL`.
#'
#' @keywords internal
selected_pair_labels <- function(cl0, cl1, k1, k2) {
  # Return the current labels corresponding to the originally selected clusters.
  # Feasibility checks allow label permutations, so downstream T^2 calculations
  # must use the mapped labels rather than blindly reusing k1/k2.
  i1 <- which(cl0 == k1); i2 <- which(cl0 == k2)
  if (length(i1) == 0 || length(i2) == 0) return(NULL)

  L1 <- unique(cl1[i1]); L2 <- unique(cl1[i2])
  if (length(L1) != 1 || length(L2) != 1 || L1 == L2) return(NULL)
  if (!all(cl1[i1] == L1) || !all(cl1[i2] == L2)) return(NULL)
  if (!all(which(cl1 == L1) %in% i1)) return(NULL)
  if (!all(which(cl1 == L2) %in% i2)) return(NULL)

  c(k1 = L1, k2 = L2)
}
