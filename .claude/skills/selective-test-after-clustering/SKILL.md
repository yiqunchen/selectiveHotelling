---
name: selective-test-after-clustering
description: Test whether two estimated clusters differ in mean with a valid selective p-value (selectiveHotelling). Use when someone clusters data (k-means, hierarchical, Leiden, or any algorithm) and then wants to test, compare, or report differences between the clusters, or asks why a t-test or Hotelling test after clustering is invalid.
---

# Selective test after clustering

Testing clusters with a textbook test after estimating them from the same data inflates the
Type I error: the naive p-value is tiny even when there are no groups. `selectiveHotelling`
gives p-values that are exact for any clustering algorithm and any unknown covariance matrix,
under Gaussian rows with a common covariance.

## Before running anything

Check these with the user; each one changes whether the answer is valid.

1. **The pair is fixed in advance, or every pair is tested.** Picking the pair that looks most
   different and testing only that one is not covered. If they did, test all pairs and report all.
2. **Enough observations.** The two clusters need `n1 + n2 >= p + 2`. With many features, reduce
   the dimension first, with principal components whose loadings come from data that are not
   tested (held-out rows, as in the paper's data examples). If the loadings come from the
   tested rows themselves, say that the guarantee does not cover this step.
3. **A deterministic clustering.** The clustering is rerun thousands of times; it must return
   the same labels on the same data. For k-means the package fixes the initial centers; for a
   custom `cl_fun`, fix seeds inside it.
4. **`M` and `L` are chosen before looking at results.** Defaults `M = 5000`, `L = 2`. `M = 1000`
   is fine for exploration. Never rerun with a larger `M` because a p-value is near 0.05.
5. **Rows roughly Gaussian.** Heavy tails can inflate the error, especially after average-linkage
   hierarchical clustering. Log-transform skewed features first.

## Running it

If the `selectiveHotelling` MCP server is connected, use its tools:
`cluster_data` (sizes, means, testable pairs), then `selective_test` with the same data,
`clustering`, `K` and `seed`. `onedim_test` gives the one-dimensional tests.

Otherwise, in R:

```r
library(selectiveHotelling)
set.seed(1)                                            # fixes the k-means initial centers
res <- selective_hotelling_test(X, k1 = 1, k2 = 2, K = 3)          # k-means
res <- selective_hotelling_test(X, k1 = 1, k2 = 2, K = 3, clustering = "hclust")
res$pval; res$cluster

# any algorithm: a function from a data matrix to integer labels
cl_fun <- function(x) my_clustering(x)                 # deterministic
selective_hotelling_bc(X, cl_fun, k1 = 1, k2 = 2)
```

For all pairs, loop over `combn(K, 2)` and reset the same seed before each call so the
clustering is identical.

## Which p-value

| p-value | Call | Use when |
|---|---|---|
| `p_polar` | `selective_hotelling_test()` / `selective_hotelling_bc()` | default; most power |
| `p_polar,1d` | `selective_truncated_F(..., clustering = "custom", method = "bc", onedim_kernel = "givens", cl_fun, cl_obs)` | faster; less power |
| `p_polar,1d^traj` | `selective_truncated_F(..., clustering = "kmeans", method = "exact")` | k-means only; instant; least power |

## Reporting

Report, for each tested pair: the cluster sizes, the naive p-value next to the selective one,
`M`, `L`, and the Monte Carlo standard error (`res$se`). Say that the selective p-value accounts
for the clusters having been estimated from the same data and does not require estimating the
covariance matrix. A large selective p-value with a tiny naive one means the apparent
difference can be explained by the clustering itself.

Do not claim the test is valid for non-Gaussian data, for pairs chosen after looking, or when
`M` was tuned on the observed result.
