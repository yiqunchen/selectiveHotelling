# Selective Hotelling T-squared test after clustering

Clusters `X` (k-means or hierarchical clustering, `K` clusters) and
tests whether the discovered clusters `k1` and `k2` have the same mean,
with Hotelling's \\T^2\\ and an unknown covariance. The p-value is the
Besag–Clifford test
[`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md)
with the random Givens kernel on the event that the clustering returns
the tested pair (Algorithm 1 of the paper). It is exactly uniform under
the null for any `M` and `L`.

## Usage

``` r
selective_hotelling_test(
  X,
  k1,
  k2,
  K = 3,
  clustering = c("kmeans", "hclust"),
  hclust_method = "average",
  hclust_square = TRUE,
  M = 5000L,
  L = 2L,
  kernel = "givens",
  ridge = 1e-08,
  verbose = FALSE
)
```

## Arguments

- X:

  A numeric data matrix with observations in rows.

- k1, k2:

  The two cluster labels (in `1:K`) to compare.

- K:

  Number of clusters.

- clustering:

  `"kmeans"` (default, Lloyd's algorithm from fixed initial centres) or
  `"hclust"` (agglomerative clustering cut at `K`).

- hclust_method:

  Agglomeration method for `clustering = "hclust"`, passed to
  [`hclust_fun()`](https://yiqunchen.github.io/selectiveHotelling/reference/hclust_fun.md)
  / [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html). Ignored
  for `clustering = "kmeans"`.

- hclust_square:

  If `TRUE` (default) the tree is built on squared Euclidean distances
  (the `clusterpval` convention). Ignored for `clustering = "kmeans"`.

- M:

  Path length (the path carries `M + 1` points). The p-value cannot be
  smaller than \\1/(M+1)\\.

- L:

  Number of kernel steps between consecutive path points.

- kernel:

  Kernel passed to
  [`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md):
  `"givens"` (default), `"cayley"`, or `"givens_slice"`.

- ridge:

  Ridge added to the pooled covariance in
  [`hotelling_T2()`](https://yiqunchen.github.io/selectiveHotelling/reference/hotelling_T2.md).

- verbose:

  Whether to print progress messages.

## Value

An S3 object of class `"selectiveHotelling"`; see
[`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md)
for its elements. It also contains `cluster`, the cluster labels of the
clustering that was tested, so that `k1` and `k2` can be matched to
observations.

## Details

For k-means, the initial centres are drawn once (from the current random
seed) and then held fixed, so that re-running the clustering on
perturbed data is a deterministic function of the data. For the
one-dimensional test, use
[`selective_truncated_F()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_truncated_F.md).

## See also

[`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md),
[`selective_truncated_F()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_truncated_F.md),
[`kmeans_fun()`](https://yiqunchen.github.io/selectiveHotelling/reference/kmeans_fun.md),
[`hclust_fun()`](https://yiqunchen.github.io/selectiveHotelling/reference/hclust_fun.md)

## Examples

``` r
set.seed(2)
Xb <- rbind(
  matrix(rnorm(20), 10, 2),
  matrix(rnorm(20) + 6, 10, 2),
  matrix(rnorm(20) - 6, 10, 2)
)
res <- selective_hotelling_test(Xb, 1, 2, K = 3, M = 200)
res$pval
#> [1] 0.02457418
print(res)
#> Selective Hotelling T^2 test (unknown covariance)
#>   Clusters tested : 1 vs 2 (of K = 3)
#>   Sampler         : besag-clifford
#>   Conditioning    : pair-preserved
#>   Observed T^2    : 213.2
#>   Selective p-val : 0.02457 (MC s.e. ~ 0.011)
res_h <- selective_hotelling_test(Xb, 1, 2, K = 3, clustering = "hclust", M = 200)
res_h$pval
#> [1] 0.006160104
```
