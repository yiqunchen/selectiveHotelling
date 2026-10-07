# Hierarchical (agglomerative) clustering of a data matrix

A thin, deterministic wrapper around
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) +
[`stats::cutree()`](https://rdrr.io/r/stats/cutree.html) that cuts an
agglomerative hierarchical clustering into `k` groups. It is a **pure
function of `x`** (no randomness, no global RNG use), which makes it a
valid black-box clustering routine to pass to
[`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md)
(the package's preferred, exact selective test for arbitrary clustering
algorithms).

By default the tree is built on **squared** Euclidean distances
(`stats::dist(x)^2`), matching the `clusterpval` / `KmeansInference`
selective-conditioning convention; set `square = FALSE` for ordinary
Euclidean distances.

Unlike
[`kmeans_fun()`](https://yiqunchen.github.io/selectiveHotelling/reference/kmeans_fun.md),
the returned object deliberately carries **no** `init_centers_idx`
field: hierarchical clustering has no random initialization to condition
on, so it should be treated as a pure black box and never re-seeded.

## Usage

``` r
hclust_fun(x, k = 3, method = "average", square = TRUE)
```

## Arguments

- x:

  A numeric \\n \times p\\ data matrix with observations in rows.

- k:

  Number of clusters to cut the dendrogram into.

- method:

  Agglomeration method passed to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) (e.g.
  `"average"`, `"single"`, `"complete"`, `"ward.D2"`).

- square:

  If `TRUE` (default) the dendrogram is built on squared Euclidean
  distances (`stats::dist(x)^2`, the `clusterpval` convention); if
  `FALSE`, ordinary Euclidean distances (`stats::dist(x)`) are used.

## Value

A list with elements:

- `cluster`:

  An integer vector of length `nrow(x)` of cluster labels in `1:k` (from
  [`stats::cutree()`](https://rdrr.io/r/stats/cutree.html)).

- `trajectory`:

  A length-1 list, `list(cluster)`. Hierarchical clustering has no Lloyd
  iterations, so the "trajectory" is a single (final) assignment vector.

The object intentionally does **not** include an `init_centers_idx`
field.

## Details

Hierarchical clustering helper for selective inference

## See also

[`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md),
[`selective_hotelling_test()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_test.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(30 * 2), 30, 2)
cl <- hclust_fun(X, k = 3)$cluster
table(cl)
#> cl
#>  1  2  3 
#> 27  1  2 
```
