#!/bin/bash
# Penguin CET pipeline: production chains as a Slurm array on Hamilton.
#
#   sbatch -J evotempo-penguin-cet --array=1-16 run-array.sh
#
# Task t runs cell (t - 1) %% 4 + 1 (cells (0,0), (0,1), (1,0), (1,1)),
# replicate (t - 1) %/% 4 + 1.  Each task runs run-chain.R, which writes a
# checkpoint after every segment and stops cleanly after CET_HOURS (default
# 70); resubmitting the same array continues every unfinished chain from its
# checkpoint (finished chains exit at once).  Settings pass through the
# environment (see run-chain.R), e.g.
#   sbatch -J evotempo-penguin-cet --array=1-16 \
#     --export=ALL,CET_NITER=200000,CET_THIN=100 run-array.sh
# Ladder rungs and power-posterior chains: set CET_CELL (a label such as
# b0.5g1) and/or CET_POWER, and use --array=1-<reps>.
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --time=3-00:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%A_%a.out
set +e; set +o pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
module purge; module load r/4.5.1 gcc/14.2
W=${W:-/nobackup/pjjg18/evotempo-penguin-cet}
SRC=${SRC:-$W/src}
LIB=${LIB:-$W/lib-$(cat $W/out/build-sha.txt)}
export R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib
export CET_OUT=${CET_OUT:-$W/runs}
cd $SRC
t=$SLURM_ARRAY_TASK_ID
if [ -n "$CET_CELL" ]; then
  cell=$CET_CELL; rep=$t
else
  cell=$(( (t - 1) % 4 + 1 )); rep=$(( (t - 1) / 4 + 1 ))
fi
Rscript data-raw/empirical/penguins/cet/run-chain.R $cell $rep \
  >> $W/out/run-$cell-$rep.txt 2>&1
echo "$(date) exit $?" >> $W/out/run-$cell-$rep.txt
