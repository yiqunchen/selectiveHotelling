# Tutorial: hierarchical clustering

``` r

library(selectiveHotelling)
```

For hierarchical clustering, set `clustering = "hclust"`. By default the
package uses average linkage on squared Euclidean distances and cuts the
tree into `K` clusters; `hclust_method` chooses another linkage.

``` r

set.seed(123)
Sigma <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
X <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100), ncol = 2) + 4) %*% chol(Sigma)
res <- selective_hotelling_test(X, k1 = 1, k2 = 2, K = 2, clustering = "hclust", M = 1000)
res
#> Selective Hotelling T^2 test (unknown covariance)
#>   Clusters tested : 1 vs 2 (of K = 2)
#>   Sampler         : besag-clifford
#>   Conditioning    : pair-preserved
#>   Observed T^2    : 842.1
#>   Selective p-val : 0.002837 (MC s.e. ~ 0.002)
table(res$cluster)
#> 
#>  1  2 
#> 50 50
```

The two clusters come from two different groups, and the test rejects.
With a single group, it does not:

``` r

set.seed(321)
X0 <- matrix(rnorm(200), ncol = 2) %*% chol(Sigma)
res0 <- selective_hotelling_test(X0, k1 = 1, k2 = 2, K = 2, clustering = "hclust",
                                 M = 1000)
res0$pval
#> [1] 0.4715096
```

With complete linkage:

``` r

selective_hotelling_test(X, k1 = 1, k2 = 2, K = 2, clustering = "hclust",
                         hclust_method = "complete", M = 1000)$pval
#> [1] 0.0002489197
```

Average-linkage hierarchical clustering is sensitive to outliers. The
[technical
details](https://yiqunchen.github.io/selectiveHotelling/articles/technical_details.md)
discuss when the tests rely on Gaussian data.
