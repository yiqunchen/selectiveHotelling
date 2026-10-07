# One-dimensional selective Hotelling test (unknown covariance)

The one-dimensional test \\p\_{\mathrm{polar,1d}}\\ of the paper: it
conditions on every coordinate of the polar factor \\Q\\ except the
Hotelling angle \\\theta\\, so that the selective null law of \\\theta\\
is a truncated Beta (equivalently, a truncated F for the scaled
\\T^2\\).

## Usage

``` r
selective_truncated_F(
  X,
  k1,
  k2,
  clustering = c("kmeans", "hclust", "custom"),
  method = c("bc", "exact", "is"),
  K = 3L,
  init_centers_idx = NULL,
  hclust_method = "average",
  hclust_square = TRUE,
  hclust_event = "trace",
  cl_fun = NULL,
  cl_obs = NULL,
  M = 5000L,
  L = 2L,
  onedim_kernel = c("givens", "slice"),
  n_is = 2000L,
  is_mix = 0.5,
  is_scale = NULL,
  ridge = 1e-08,
  verbose = FALSE
)
```

## Arguments

- X:

  A numeric \\n \times p\\ data matrix with observations in rows.

- k1, k2:

  The two tested cluster labels. Each cluster must have at least three
  members, and \\r - 1 \> p\\.

- clustering:

  `"kmeans"`, `"hclust"`, or `"custom"` (use `cl_fun`).

- method:

  `"bc"` (default), `"exact"` (closed form; `"kmeans"` and `"hclust"`
  only), or `"is"`.

- K:

  Number of clusters (built-in `"kmeans"`/`"hclust"`).

- init_centers_idx:

  Fixed initial-centre row indices for `"kmeans"`; drawn via
  [`kmeans_fun()`](https://yiqunchen.github.io/selectiveHotelling/reference/kmeans_fun.md)
  if `NULL`.

- hclust_method, hclust_square:

  Linkage and squared-distance flag for `"hclust"` (defaults match
  [`hclust_fun()`](https://yiqunchen.github.io/selectiveHotelling/reference/hclust_fun.md)).

- hclust_event:

  Kept for compatibility; the only value is `"trace"`, the merge order
  used by `method = "exact"`.

- cl_fun:

  A clustering routine `function(Xmat) -> labels` (a pure function of
  its data) for `clustering = "custom"`; optional otherwise.

- cl_obs:

  Optional precomputed observed clustering (saves one call).

- M, L:

  Path length and number of kernel steps between path points for
  `method = "bc"`.

- onedim_kernel:

  For `method = "bc"`: `"givens"` (default; random step on the angle,
  step size drawn log-uniformly on `[1e-3, 1]` per block of `L` steps)
  or `"slice"` (slice sampling with shrinkage).

- n_is, is_mix, is_scale:

  Number of draws, weight of the null law in the defensive mixture, and
  the standard deviation of its normal part, for `method = "is"`.

- ridge:

  Small non-negative ridge for the reported observed \\T^2\\.

- verbose:

  If `TRUE`, print a one-line diagnostic.

## Value

An S3 object of class `"selectiveHotelling"` with `pval`, `T_obs`, `k1`,
`k2`, `K`, `condition_on`, `sampler`, `p_uncond` (the classical F
p-value), `theta_obs`, `F_obs`, `df1`, `df2`, and, for
`method = "exact"`, the feasible arcs `intervals`.

## Details

With \\\cos^2\theta = u(Q)\\, \\T^2 = (r-2)\cot^2\theta\\, and under the
null \\\cos^2\theta \sim \mathrm{Beta}(p/2, (r-1-p)/2)\\ independently
of the other coordinates of \\Q\\. Given those coordinates, the
clustering restricts \\\theta\\ to a set \\S\_\theta\\, and the p-value
is the Beta mass of \\S\_\theta \cap (0, \theta\_{obs}\]\\ divided by
the Beta mass of \\S\_\theta\\.

`method = "bc"` (default) computes it for any clustering, on the event
that the clustering returns the tested pair, with the Besag–Clifford
construction on the angle (Algorithm 2 of the paper); it is exact for
any `M`. `method = "exact"` gives a closed form, with no Monte Carlo,
for two stronger events: the whole Lloyd trajectory
(`clustering = "kmeans"`) and the merge order up to the cut
(`clustering = "hclust"`, average or single linkage). These events
condition on more and cost power. `method = "is"` uses importance
sampling (consistent, not exact).

## See also

[`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md),
[`selective_hotelling_test()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_test.md),
[`hotelling_T2()`](https://yiqunchen.github.io/selectiveHotelling/reference/hotelling_T2.md)

## Examples

``` r
set.seed(1)
X <- rbind(matrix(rnorm(40), 20, 2), matrix(rnorm(40) + 5, 20, 2),
           matrix(rnorm(40) - 5, 20, 2))
cl_fun <- function(x) hclust_fun(x, k = 3)$cluster
res <- selective_truncated_F(X, 1, 2, clustering = "custom", cl_fun = cl_fun, M = 200)
res$pval
#> [1] 0.01469148
```
