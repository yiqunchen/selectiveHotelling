# K-means from a fixed initialization

Runs Lloyd's k-means from an explicitly supplied set of initial center
indices, recording the full assignment trajectory. Because the
initialization is fixed and no RNG is used, this is the loop-safe
re-clustering routine used for per-step feasibility checks inside the
MCMC samplers.

## Usage

``` r
kmeans_with_fixed_init(x, init_centers_idx, max_iter = 50)
```

## Arguments

- x:

  A numeric data matrix with observations in rows.

- init_centers_idx:

  Integer vector of row indices of `x` used as initial centers; its
  length determines `k`.

- max_iter:

  Maximum number of Lloyd iterations.

## Value

A list with elements `cluster`, `trajectory`, `center_history`,
`init_centers_idx`, and `n_iter`.

## Examples

``` r
set.seed(1)
X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat
km <- kmeans_with_fixed_init(X, init_centers_idx = c(1, 11, 21))
km$n_iter
#> [1] 7
```
