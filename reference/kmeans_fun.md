# Trajectory-recording k-means with a simple interface

Convenience wrapper around
[`kmeans_with_trajectory()`](https://yiqunchen.github.io/selectiveHotelling/reference/kmeans_with_trajectory.md)
that runs k-means with `k` clusters and returns the full trajectory
object. This is the clustering function expected by
[`selective_hotelling_test()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_test.md)
and the MCMC samplers (because it returns `init_centers_idx`, enabling
loop-safe fixed-init re-clustering).

## Usage

``` r
kmeans_fun(x, k = 3)
```

## Arguments

- x:

  A numeric data matrix with observations in rows.

- k:

  Number of clusters.

## Value

The list returned by
[`kmeans_with_trajectory()`](https://yiqunchen.github.io/selectiveHotelling/reference/kmeans_with_trajectory.md)
(with `cluster`, `trajectory`, `init_centers_idx`, etc.).

## Examples

``` r
set.seed(1)
X <- generate_exchangeable_null_data(30, 2, rho = 0.1)$dat
km <- kmeans_fun(X, k = 3)
table(km$cluster)
#> 
#>  1  2  3 
#> 12 14  4 
```
