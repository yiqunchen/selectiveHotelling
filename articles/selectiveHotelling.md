# Getting started with selectiveHotelling

This vignette shows the basic workflow. The tutorial in the package
repository (`tutorial/`) covers each function in more detail.

``` r

library(selectiveHotelling)
```

## Cluster, then test

We simulate $`n = 60`$ observations from three groups and cluster them
with $`k`$-means.

``` r

set.seed(1)
X <- rbind(matrix(rnorm(40), 20, 2),
           matrix(rnorm(40), 20, 2) + 4,
           matrix(rnorm(40), 20, 2) - 4)
res <- selective_hotelling_test(X, k1 = 1, k2 = 2, K = 3, M = 1000)
res
#> Selective Hotelling T^2 test (unknown covariance)
#>   Clusters tested : 1 vs 2 (of K = 3)
#>   Sampler         : besag-clifford
#>   Conditioning    : pair-preserved
#>   Observed T^2    :  1405
#>   Selective p-val : 0.004655 (MC s.e. ~ 0.002)
table(res$cluster)
#> 
#>  1  2  3 
#> 20 20 20
```

`res$pval` is the selective p-value for the null hypothesis that the
observations in clusters 1 and 2 share one mean. The default path length
is `M = 5000`; we use `M = 1000` here so that the vignette builds
quickly. The covariance matrix of the features is never estimated.

## The naive test is not valid

When all observations come from one group, the clusters still look
different, and the classical Hotelling test rejects almost always. The
selective test does not.

``` r

set.seed(2)
X0 <- matrix(rnorm(60 * 2), 60, 2)
res0 <- selective_hotelling_test(X0, k1 = 1, k2 = 2, K = 3, M = 1000)
cl <- res0$cluster
n1 <- sum(cl == 1); n2 <- sum(cl == 2); p <- 2
F_stat <- (n1 + n2 - p - 1) / ((n1 + n2 - 2) * p) * hotelling_T2(X0, cl, 1, 2)
c(naive = pf(F_stat, p, n1 + n2 - p - 1, lower.tail = FALSE), selective = res0$pval)
#>        naive    selective 
#> 6.786752e-13 4.538853e-03
```

## Other clustering algorithms

``` r

selective_hotelling_test(X, k1 = 1, k2 = 2, K = 3, clustering = "hclust", M = 1000)$pval
#> [1] 0.00832417

cl_fun <- function(x) cutree(hclust(dist(x), method = "complete"), k = 3)
selective_hotelling_bc(X, cl_fun, k1 = 1, k2 = 2, M = 1000)$pval
#> [1] 0.0006205319
```

Any function that maps a data matrix to cluster labels can be used as
`cl_fun`, as long as it gives the same labels every time it is run on
the same data.
