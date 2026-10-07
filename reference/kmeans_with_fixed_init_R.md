# K-means from a fixed initialization (implementation)

Pure-R implementation backing
[`kmeans_with_fixed_init()`](https://yiqunchen.github.io/selectiveHotelling/reference/kmeans_with_fixed_init.md).
Runs Lloyd's k-means from a fixed initialization with full trajectory
tracking.

## Usage

``` r
kmeans_with_fixed_init_R(x, init_centers_idx, max_iter = 50)
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
`init_centers_idx`, and `n_iter` (identical in shape to
[`kmeans_with_fixed_init()`](https://yiqunchen.github.io/selectiveHotelling/reference/kmeans_with_fixed_init.md)).
