#!/bin/bash
# High-dimensional study array (p = 2n, 5n): selective_meandiff_bc() after fast
# k-means or fast average linkage. Manifest: dev/highdim/cluster/manifest.R.
# Runner: dev/highdim/cluster/run.R. See dev/highdim/cluster/README.md.
#
#   cd $REPO_ROOT/selectiveHotelling && git checkout highdim-stats
#   source $REPO_ROOT/slurm/renv.sh
#   Rscript dev/highdim/cluster/manifest.R                 # prints the array range
#   HD_SMOKE=2 sbatch -p normal --array=1 dev/highdim/cluster/highdim.sbatch   # smoke
#   sbatch -p normal --array=1-<N>%200 dev/highdim/cluster/highdim.sbatch
#   Rscript dev/highdim/cluster/aggregate.R
#SBATCH --job-name=selhot_hd
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/hd_%A_%a.log
export REPO_ROOT="${REPO_ROOT:-/users/ychen5/selectiveHotelling}"
export PKG_DIR="${PKG_DIR:-$REPO_ROOT/selectiveHotelling}"
source "$REPO_ROOT/slurm/renv.sh"
command -v Rscript >/dev/null 2>&1 || { echo "FATAL: R env not set up"; exit 1; }
export HD_MANIFEST="${HD_MANIFEST:-$PKG_DIR/dev/highdim/cluster/manifest.rds}"
export HD_SHARDS="${HD_SHARDS:-$PKG_DIR/dev/highdim/cluster/shards}"
# one BLAS thread per R worker (the runner forks SLURM_CPUS_PER_TASK workers)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$REPO_ROOT/logs"
cd "$REPO_ROOT"
echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID cpus=$SLURM_CPUS_PER_TASK manifest=$HD_MANIFEST"
Rscript "$PKG_DIR/dev/highdim/cluster/run.R" < /dev/null
