# selectiveHotelling

Selective Hotelling $T^2$ testing after k-means clustering, with an
**unknown** covariance matrix.

When you cluster data and then test whether two of the discovered clusters
have different means, the clustering has already used the data to maximize
between-cluster separation. A textbook two-sample test that ignores this
selection step is badly anti-conservative. `selectiveHotelling` computes a
**valid selective p-value** for the Hotelling $T^2$ statistic that accounts
for the k-means selection event, while treating the covariance $\Sigma$ as
unknown.

## Installation

The package is not on CRAN. Install it from a local source checkout:

```r
# with remotes
remotes::install_local("path/to/selectiveHotelling")

# or with devtools (also builds the vignette)
devtools::install("path/to/selectiveHotelling", build_vignettes = TRUE)
```

You can also load it for development without installing:

```r
devtools::load_all("path/to/selectiveHotelling")
```

The core test depends only on base R and `stats`. The packages in `Suggests`
(`knitr`, `rmarkdown`, `progress`, `parallel`, `rstiefel`) are optional and
used only for the vignette, a progress bar, and experimental samplers.

## Quickstart

```r
library(selectiveHotelling)

set.seed(1)
X  <- generate_exchangeable_null_data(n = 30, p = 2, rho = 0.2)$dat  # global null
cl <- kmeans_fun(X, k = 3)$cluster                                   # discover clusters
res <- selective_hotelling_test(X, k1 = 1, k2 = 2, K = 3,
                                sampler = "fast",
                                burn = 200, B = 400, thin = 2)        # selective test
res$pval   # valid selective p-value
print(res) # full summary
```

## Why naive Hotelling fails after clustering

K-means assigns observations to clusters so as to *maximize* the separation
between cluster means. If you then feed those clusters into a standard
two-sample Hotelling $T^2$ test and read off the usual $F$ p-value, the null
distribution is wrong: you have already selected the most-separated grouping
of the data, so even under a global null (all observations share one mean)
the test rejects far too often. In simulations on null data, the naive p-value
piles up near zero and rejects at essentially every threshold, whereas a valid
test should produce p-values that are approximately uniform. The vignette
demonstrates this anti-conservatism directly and shows that
`selective_hotelling_test()` restores calibration.

## How it works (one paragraph)

We test equality of the means of two k-means clusters with an unknown
covariance $\Sigma$. Writing the pair's data in terms of a contrast direction,
a **polar decomposition** factors out the scale/covariance information and
leaves a polar factor $Q$ that lives on a Stiefel manifold. The key fact is
that under the null this polar factor $Q$ is **uniform on the Stiefel
manifold regardless of $\Sigma$**, so the nuisance covariance drops out of the
reference distribution. To respect the data-driven selection, we run an MCMC
sampler that draws $Q$ **uniformly over the subset of the Stiefel manifold
that is consistent with the observed selection event** (the cluster pair we
chose to test). Each accepted $Q$ is mapped back to a dataset, its Hotelling
$T^2$ is recomputed, and the Monte Carlo upper-tail probability of the observed
$T^2$ gives a **valid selective p-value**. The supported, validated path
conditions on the **final cluster pair** (`condition_on = "final"`, the
default); conditioning on the full k-means **trajectory**
(`condition_on = "trajectory"`) is **experimental** and not recommended for
inference, because that region is much tighter and the local random walk
under-mixes.

## Beyond k-means: hierarchical and arbitrary clustering

The same Stiefel-manifold argument works for *any* clustering routine, not just
k-means. `selective_hotelling_bc()` implements an **exact** selective test based
on the Besag & Clifford (1991) serial Monte Carlo (random-index) construction:
it runs a reversible Markov chain on the Stiefel manifold restricted to the
selection event, and the rank of the observed statistic among the path is
*exactly* uniform under the null. Validity holds for **any** reversible kernel,
**any** path length `M`, and **any** clustering function — only the *power*
depends on these choices.

```r
# Hierarchical clustering, via the one-line wrapper:
res_h <- selective_hotelling_test(X, k1 = 1, k2 = 2, K = 3,
                                  clustering = "hclust",   # average linkage, squared distances
                                  M = 200)
res_h$pval

# ...or call the algorithm-agnostic engine directly with ANY cl_fun:
cl_fun <- function(x) hclust_fun(x, k = 3, method = "average")$cluster
res_g  <- selective_hotelling_bc(X, cl_fun, k1 = 1, k2 = 2, M = 200)
res_g$pval
```

`hclust_fun()` is the bundled hierarchical-clustering `cl_fun` (using squared
Euclidean distances by default, matching the `clusterpval` convention). Any
deterministic, data-only clustering routine can be passed to
`selective_hotelling_bc()` in the same way.

## Background and references

This package extends selective inference for clustering to the
unknown-covariance setting. For background see:

- Chen, Y. T. and Witten, D. M. (2023). *Selective inference for k-means
  clustering.* Journal of Machine Learning Research.
- Gao, L. L., Bien, J. and Witten, D. (2024). *Selective inference for
  hierarchical clustering.* Journal of the American Statistical Association.
- Besag, J. and Clifford, P. (1991). *Sequential Monte Carlo p-values.*
  Biometrika, 78(2), 301-304.

## License

MIT (see `LICENSE`).
