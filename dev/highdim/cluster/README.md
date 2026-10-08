# High-dimensional study (p > n1 + n2)

This study tests `selective_meandiff_bc()` when the number of features is larger than the
number of tested rows. In that case the Hotelling T² of `selective_hotelling_bc()` is
constant. The test conditions on the centred total scatter of the tested rows, uses the
same random Givens kernel on the pair event, and ranks several statistics on one
Besag-Clifford path. All statistics therefore cost the same clustering calls; each extra
statistic adds only a few milliseconds per path point.

## Design

| Factor | Levels |
|---|---|
| Clustering | k-means (K = 3, initial rows drawn from the seed and then fixed), average linkage on squared distances (cut at K = 3) |
| n | 50, 100, 500 |
| p | 2n and 5n |
| Σ | `ar1` (0.6), `cs` (0.5), `hetero` (diagonal, variances over a factor of 10), `spike5` (I + 5vv') |
| Mean shape | `line` (1-sparse, as in the main grid), `block` (1/4-sparse, as in the main grid), `dense` |
| Signal size | s = δ² / sqrt(2 tr Σ²) in {1, 2, 4, 8, 16, 32}; null cells have s = 0 |
| Data sets per cell | null 500 and power 200 for n ≤ 100; null 300 and power 100 for n = 500 |
| Test | every pair whose clusters both have at least 3 rows; the package defaults (M = 5000, L = 2, Givens kernel, `eps = "mixture"`) |

The statistics are `l2`, `cq`, `diag`, `diag_w`, `max`, `ridge`, `proj`, `t2` (the regular
pooled T² with a 1e-8 ridge), the known-Σ `oracle`, and two min-p combinations,
`adapt_l2_max` and `adapt_dw_max`. Every p-value is exact.

The data-generating process is in `dev/highdim/dgp.R`. The clustering uses the fast routines
`make_fast_kmeans_cl_fun()` and `make_fast_hclust_cl_fun()`, which are on `main`.

## Files

- `manifest.R` writes `manifest.rds`, with one row per array task. It prints a cost model of
  about 10,800 CPU-h. The n = 500 cells account for about 8,600 of them. The sizes can be
  changed with environment variables (see the header of `manifest.R`).
- `run.R` runs one task. It screens seeds until it has `N_shard` testable data sets, then
  writes one row per tested pair with every p-value (`p_<statistic>`), `pair_alt`,
  `recover` and the CPU time. It skips a shard that already exists.
- `highdim.sbatch` is the Slurm wrapper. It uses 4 cores and one BLAS thread per worker.
- `aggregate.R` collects the shards and writes Type I and power tables.

## Running

```bash
cd $REPO_ROOT/selectiveHotelling && git checkout highdim-stats
source $REPO_ROOT/slurm/renv.sh
Rscript -e 'install.packages("fastcluster")'          # once; stats::hclust is the fallback
Rscript dev/highdim/cluster/manifest.R                # 1280 tasks
HD_SMOKE=2 sbatch -p normal --array=1 dev/highdim/cluster/highdim.sbatch
sbatch -p normal --array=1-1000 dev/highdim/cluster/highdim.sbatch
sbatch -p normal --array=1-280 --export=ALL,HD_TASK_OFFSET=1000 dev/highdim/cluster/highdim.sbatch
Rscript dev/highdim/cluster/aggregate.R
```

If the cluster's `MaxArraySize` is below the number of tasks, submit chunks with
`HD_TASK_OFFSET`. The speed of `tcrossprod` depends on the BLAS, so check `sessionInfo()`
on the cluster. The cost model assumes an optimised BLAS (OpenBLAS, MKL or Accelerate).

The runner has not been run end to end yet. Run the smoke task first and compare the CPU
seconds per pair in its log with the cost model.

## Pilot (local, n = 50, p = 100)

`dev/highdim/pilot_typeI_power.png` is drawn by `plot_pilot.R` from `sim*_raw.rds`.

- Type I error was 0.025 to 0.072 over about 195 null reps per Σ (Monte Carlo standard
  error about 0.016).
- When k-means recovered the pair, `l2`, `cq`, `ridge` and `diag_w` matched the oracle.
- `diag` with total-scatter variances lost power for sparse signals.
- `max` lost power for dense signals.
