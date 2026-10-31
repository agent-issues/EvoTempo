#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH -t 24:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-cost52/logs/%x-%A_%a.out
# One chain of data-raw/move-cost-52.R (#52): local path-move count vs tree size.
#   Rscript data-raw/move-cost-52.R trees     # once, before the array
#   sbatch --array=1-128 data-raw/move-cost-52-array.sh
#   Rscript data-raw/move-cost-52.R summary
# Build the package with the incremental pruning of #44 (merged, or its branch):
# the per-iteration cost the default must be tuned for.
set +e
set +o pipefail
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=/nobackup/pjjg18/evotempo-cost52
export R_LIBS=$W/lib:/nobackup/pjjg18/evotempo-dev/lib
export EVOTEMPO_SHA="package=$(cat $W/lib/SHA);scripts=$(git -C $W/EvoTempo rev-parse --short HEAD)"
export EVOTEMPO_COST_OUT=${EVOTEMPO_COST_OUT:-$W/out}
cd $W/EvoTempo
Rscript data-raw/move-cost-52.R task $SLURM_ARRAY_TASK_ID
