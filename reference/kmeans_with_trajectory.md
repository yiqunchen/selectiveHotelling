# K-means recording the full assignment trajectory

Runs Lloyd's k-means while recording the cluster assignment vector at
*every* iteration (the "trajectory"), which is needed for
trajectory-based selective conditioning.

## Usage

``` r
kmeans_with_trajectory(
  x,
  k = 3,
  max_iter = 50,
  init_centers_idx = NULL,
  seed = NULL,
  init_method = c("random_indices", "kmeans++")
)
```

## Arguments

- x:

  A numeric data matrix with observations in rows.

- k:

  Number of clusters.

- max_iter:

  Maximum number of Lloyd iterations.

- init_centers_idx:

  Optional integer vector of length `k` giving row indices of `x` to use
  as initial centers. If `NULL`, centers are chosen by `init_method`.

- seed:

  Optional integer seed applied locally (global RNG restored on exit).

- init_method:

  Initialization scheme when `init_centers_idx` is `NULL`:
  `"random_indices"` (default, data-independent) or `"kmeans++"`.

## Value

A list with elements `cluster` (final assignments), `trajectory` (a list
of per-iteration assignment vectors), `center_history`,
`init_centers_idx`, and `n_iter`.

## Details

This function does **not** silently re-seed the global RNG.
Reproducibility is the caller's responsibility: seed once at the call
site, or pass an explicit `init_centers_idx`. If `seed` is supplied it
is applied *locally* and the global RNG stream is restored on exit, so
calling this inside an MCMC loop can never freeze the chain's proposals.
