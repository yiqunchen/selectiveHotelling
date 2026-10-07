# Tutorial: any clustering algorithm

``` r

library(selectiveHotelling)
```

The test needs only one thing from the clustering algorithm: a function
`cl_fun` that takes a data matrix and returns a vector of cluster
labels.
[`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md)
reruns this function on perturbed versions of the data and keeps the
perturbations under which the tested clusters are returned again.

The function must be **deterministic**: run twice on the same data, it
must return the same labels. Fix any random seed or random
initialization inside the function, as below.

## Partitioning around medoids

``` r

library(cluster)
set.seed(123)
X <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100), ncol = 2) + 4,
           matrix(rnorm(100), ncol = 2) - 4)

cl_fun <- function(x) cluster::pam(x, k = 3, cluster.only = TRUE)
cl <- cl_fun(X)
table(cl)
#> cl
#>  1  2  3 
#> 50 50 50
selective_hotelling_bc(X, cl_fun, k1 = 1, k2 = 2, M = 1000)$pval
#> [1] 0.0003431453
```

## A clustering function with a random step

If the algorithm uses random numbers, fix the seed inside the function
and restore the caller’s random number stream afterwards:

``` r

cl_fun_random <- function(x) {
  old <- .Random.seed
  on.exit(assign(".Random.seed", old, envir = globalenv()))
  set.seed(1)
  kmeans(x, centers = 3, nstart = 5)$cluster
}
selective_hotelling_bc(X, cl_fun_random, k1 = 1, k2 = 2, M = 1000)$pval
#> [1] 0.007543708
```

The same approach works for community detection on a nearest-neighbor
graph, as used for single-cell data, or for any other method.

## What `selective_hotelling_bc()` returns

``` r

res <- selective_hotelling_bc(X, cl_fun, k1 = 1, k2 = 3, M = 1000)
c(pval = res$pval, T2 = res$T_obs, move_fraction = res$move_frac)
#>          pval            T2 move_fraction 
#>  4.792621e-04  9.510972e+02  7.510000e-01
```

`move_fraction` is the fraction of steps at which the Monte Carlo path
moved. Values far below about $`0.3`$ suggest that the path barely
explores the set of data sets that reproduce the clusters, which lowers
power but not validity ([the technical
details](https://yiqunchen.github.io/selectiveHotelling/articles/technical_details.md)).
