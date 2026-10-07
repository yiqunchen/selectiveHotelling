# Algorithm-agnostic Besag-Clifford selective Hotelling test (arbitrary clustering)

A Monte Carlo selective test of equality of the means of a discovered
cluster pair `(k1, k2)`, valid after clustering `X` with an
**arbitrary** black-box clustering routine `cl_fun`, treating the
covariance \\\Sigma\\ as unknown. It is the package's entry point for
selective inference after *any* clustering algorithm (hierarchical
clustering, fixed-init k-means, spectral clustering, ...), in contrast
to
[`selective_hotelling_test()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_test.md),
which is hard-wired to k-means.

The construction is exactly that of Besag and Clifford (1991),
transported to the Stiefel manifold. Under the null the observed polar
factor \\Q\_{obs}\\ is uniform on the feasible region \\E = \\Q \in
St(r-1, p) : cl\\fun(X(Q))\\ \mathrm{preserves}\\ (k_1,k_2)\\\\
regardless of \\\Sigma\\. We build a Metropolis kernel \\K\\ on \\E\\
using a **symmetric** Cayley proposal \\Q' = U Q\\ (with \\U\\ a random
Cayley rotation of step `eps`) that is accepted iff \\Q' \in E\\; for a
flat target restricted to \\E\\ this acceptance ratio is the indicator
\\1\\Q' \in E\\\\, so \\K\\ (and \\K^L\\) is reversible w.r.t.
\\\mathrm{Uniform}(E)\\. We then plant \\Q\_{obs}\\ at a uniform random
index \\J \in \\0,\dots,M\\\\ and run \\K^L\\ forward \\J \to M\\ and
backward \\J \to 0\\. Because the seed is stationary and \\K^L\\
reversible, the rank of the contrast statistic \\u(Q) = \\Q'v\\^2\\
(which is strictly monotone in the pair's two-sample Hotelling \\T^2\\)
of \\Q\_{obs}\\ among the \\M+1\\ path marks is exactly
\\\mathrm{Uniform}\\1,\dots,M+1\\\\. Hence the randomized upper-tail
rank p-value is exactly \\\mathrm{Uniform}(0,1)\\ under \\H_0\\.

**Validity holds for ANY reversible kernel, ANY step `eps`, ANY path
length `M`, and ANY clustering `cl_fun`.** Mixing (governed by `eps`,
`L`, and the acceptance rate) buys *power* only, never validity.

## Usage

``` r
selective_hotelling_bc(
  X,
  cl_fun,
  k1,
  k2,
  M = 5000L,
  L = 2L,
  eps = "mixture",
  randomized = TRUE,
  eps_range = c(0.001, 1),
  cl_obs = NULL,
  ridge = 1e-08,
  stiefel_tol = 1e-08,
  guard_feas_checks = 25L,
  oracle_Sigma = NULL,
  verbose = FALSE,
  kernel = c("givens", "cayley", "givens_slice"),
  slice_max_props = 200L
)
```

## Arguments

- X:

  A numeric \\n \times p\\ data matrix with observations in rows.

- cl_fun:

  A black-box clustering routine: a function
  `function(Xmat) -> integer label vector` (or a list carrying
  `$cluster`, e.g. the result of
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html)). It **must
  be a pure function of its data argument** (only data-independent
  randomness is allowed, e.g. a fixed initialization or a deterministic
  algorithm such as hierarchical clustering). It is called once on `X`
  to obtain the observed clustering and many times on reconstructed data
  \\X(Q)\\ during sampling.

- k1, k2:

  The two tested cluster labels, in the labelling of `cl_fun(X)`. Each
  cluster must have at least three members for a Hotelling \\T^2\\.

- M:

  Path length; the path carries `M + 1` marks. Validity holds for any
  `M >= 1`; larger `M` reduces Monte Carlo error.

- L:

  Number of kernel steps per path edge (i.e. \\K^L\\). Larger `L`
  improves spread/mixing and hence power.

- eps:

  `"mixture"` (default) draws a step size independently for each path
  edge from the fixed log-uniform `eps_range`. A positive numeric value
  uses a fixed step size, and `NULL` uses `0.3`. The step size is never
  tuned on the observed point, because that can break exactness.

- randomized:

  If `TRUE` (default) returns the randomized upper-tail rank p-value,
  which is exactly \\\mathrm{Uniform}(0,1)\\ under \\H_0\\; if `FALSE`,
  ties are broken conservatively.

- eps_range:

  Bounds of the fixed log-uniform step-size mixture.

- cl_obs:

  Optional precomputed `cl_fun(X)` (saves one call); otherwise it is
  computed internally.

- ridge:

  Small non-negative ridge added to the pooled covariance when computing
  the reported observed \\T^2\\ via
  [`hotelling_T2()`](https://yiqunchen.github.io/selectiveHotelling/reference/hotelling_T2.md).

- stiefel_tol:

  Guard tolerance on \\\max\|Q'Q - I\|\\ over the path (GUARD a).

- guard_feas_checks:

  Number of path points to re-verify with the generic feasibility
  predicate (GUARD b).

- oracle_Sigma:

  Optional known \\p \times p\\ covariance matrix. If supplied, an
  additional known-\\\Sigma\\ oracle p-value is computed on the same
  path by ranking the Mahalanobis statistic \\d' \Sigma^{-1} d\\ of the
  pair-mean difference; returned as `p_oracle`. Validity is identical;
  only the statistic (hence power) differs.

- verbose:

  If `TRUE`, print a one-block diagnostic.

- kernel:

  Proposal of the Metropolis kernel. `"givens"` (default; the random
  Givens kernel of the paper) picks two rows \\i \< j\\ of \\Q\\
  uniformly at random and proposes \\Q' = G\_{ij}(\omega) Q\\, a
  rotation by angle \\\omega \sim N(0, (\pi\\\mathrm{eps})^2)\\ in that
  coordinate plane (see
  [`givens_left()`](https://yiqunchen.github.io/selectiveHotelling/reference/givens_left.md)).
  Both proposals are symmetric and preserve the uniform law on the
  Stiefel manifold, so both kernels are reversible w.r.t.
  \\\mathrm{Uniform}(E)\\ and the p-value is exact for either; they
  differ only in mixing, hence power. `"cayley"` proposes \\Q' = UQ\\
  with \\U\\ a random Cayley rotation of step `eps`, which moves every
  row of \\Q\\ at once. `eps`, `eps_range` and the `"mixture"` rule
  apply to both. `"givens_slice"` is the shrinking-bracket Givens move
  (the shrinkage step of elliptical slice sampling; Neal 2003; Murray,
  Adams and MacKay 2010): it picks rows \\i \< j\\ uniformly, draws
  \\\omega_0 \sim \mathrm{Uniform}(0, 2\pi)\\, sets the bracket
  \\\[\omega_0 - 2\pi, \omega_0\]\\ (which contains the current point
  \\\omega = 0\\), proposes \\\omega = \omega_0\\ and, while the
  proposal is infeasible, shrinks the bracket toward 0 and redraws
  \\\omega\\ uniformly in it. It is reversible w.r.t.
  \\\mathrm{Uniform}(E)\\, needs no step size (`eps` and `eps_range` are
  ignored) and nearly always moves, at the price of several `cl_fun`
  calls per step (see `n_path_calls`).

- slice_max_props:

  For `kernel = "givens_slice"`: maximum number of proposals (calls of
  `cl_fun`) in one shrinking-bracket step; if all fail the chain stays
  put. The cap is symmetric in the two directions of a move, so the
  kernel stays reversible and the p-value exact.

## Value

An S3 object of class `"selectiveHotelling"`: a list whose primary
fields are `pval` (the randomized Besag-Clifford p-value), `se` (its
Monte Carlo standard error), `T_obs` (the observed two-sample Hotelling
\\T^2\\ of the tested pair), `draws` (the vector of path statistic
values used for ranking), `acc` (kernel move fraction),
`n_unique_draws`, `k1`, `k2`, `K` (number of distinct labels in the
observed clustering), `condition_on` (`"pair-preserved"`), and `sampler`
(`"besag-clifford"`). Additional diagnostic fields are also returned:
`p_cons` (conservative tie-broken p-value), `p_oracle`, `M`, `L`, `eps`,
`u_obs`, `path_u`, `J`, `move_frac`, `stiefel_max`, `feas_ok`,
`n_feasible_path`, `kernel`, `n_path_calls` (calls of `cl_fun` made by
the kernel along the path, the compute measure for comparing kernels)
and `n_cluster_calls` (all calls of `cl_fun` made by this function,
including the observed clustering when `cl_obs` is not supplied,
step-size calibration and the guards).

## See also

[`selective_hotelling_test()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_test.md),
[`hotelling_T2()`](https://yiqunchen.github.io/selectiveHotelling/reference/hotelling_T2.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(30 * 2), 30, 2)
# A deterministic, data-pure clustering routine (hierarchical clustering).
hclust_cl <- function(Xmat) {
  hc <- stats::hclust(stats::dist(Xmat), method = "average")
  as.integer(stats::cutree(hc, k = 3))
}
cl <- hclust_cl(X)
# Pick a pair with at least 3 members each.
sizes <- table(cl)
big <- as.integer(names(sizes)[sizes >= 3])
res <- selective_hotelling_bc(X, hclust_cl, big[1], big[2], M = 60)
res$pval
#> [1] 0.6985158
print(res)
#> Selective Hotelling T^2 test (unknown covariance)
#>   Clusters tested : 2 vs 3 (of K = 3)
#>   Sampler         : besag-clifford
#>   Conditioning    : pair-preserved
#>   Observed T^2    : 38.05
#>   Selective p-val : 0.6985 (MC s.e. ~ 0.059)
```
