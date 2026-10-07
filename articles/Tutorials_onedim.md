# Tutorial: the one-dimensional test

``` r

library(selectiveHotelling)
```

[`selective_truncated_F()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_truncated_F.md)
computes a second selective $`p`$-value, $`p_{\mathrm{polar,1d}}`$. It
holds fixed everything about the tested observations except one angle,
the one that determines Hotelling’s $`T^2`$. Its reference distribution
is then a truncated $`F`$ distribution in one dimension, which makes it
faster to compute. The price is power: fixing more of the data leaves
less room to detect a difference ([the technical
details](https://yiqunchen.github.io/selectiveHotelling/articles/technical_details.md)).

``` r

set.seed(123)
X <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100), ncol = 2) + 3,
           matrix(rnorm(100), ncol = 2) - 3)
init <- c(1, 51, 101)
cl_fun <- function(x) kmeans_with_fixed_init(x, init_centers_idx = init)$cluster
cl <- cl_fun(X)
```

## Any clustering algorithm

With `clustering = "custom"`, the test takes the same clustering
function as
[`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md)
and computes the $`p`$-value by a Monte Carlo path in one dimension:

``` r

res1d <- selective_truncated_F(X, 1, 2, clustering = "custom", method = "bc",
                               cl_fun = cl_fun, cl_obs = cl, M = 1000)
c(onedim = res1d$pval, classical_F = res1d$p_uncond)
#>      onedim classical_F 
#>  0.02469193  0.00000000
```

`p_uncond` is the classical $`F`$$`p`$-value, which ignores the
clustering.

## A closed form for $`k`$-means

For $`k`$-means, the $`p`$-value can be computed exactly, without Monte
Carlo, if we also condition on the whole sequence of assignments made by
Lloyd’s algorithm:

``` r

res_exact <- selective_truncated_F(X, 1, 2, clustering = "kmeans", method = "exact",
                                   K = 3, init_centers_idx = init)
res_exact$pval
#> [1] 0.04102747
```

This version conditions on the most, and it is the least powerful of the
three tests.

## Which test to use

| Test | Function | Conditions on | Speed | Power |
|----|----|----|----|----|
| $`p_{\mathrm{polar}}`$ | [`selective_hotelling_test()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_test.md), [`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md) | the clustering | slowest | highest |
| $`p_{\mathrm{polar,1d}}`$ | `selective_truncated_F(method = "bc")` | the clustering and all but one angle | faster | lower |
| $`p_{\mathrm{polar,1d}}^{\mathrm{traj}}`$ | `selective_truncated_F(method = "exact")` | also the $`k`$-means trajectory | instant | lowest |

All three control the selective Type I error. We recommend
$`p_{\mathrm{polar}}`$ unless computing time is a constraint.
