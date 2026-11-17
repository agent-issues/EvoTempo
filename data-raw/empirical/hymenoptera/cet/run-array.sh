#!/bin/bash
# Hymenoptera CET pipeline: production chains as a Slurm array on Hamilton
# (as the penguin pipeline's run-array.sh).
#
#   sbatch -J evotempo-hymcet-prod --array=1-16 run-array.sh
#
# Task t runs cell (t - 1) %% 4 + 1 (cells (0,0), (0,1), (1,0), (1,1)),
# replicate (t - 1) %/% 4 + 1.  Each task runs run-chain.R, which writes a
# checkpoint after every segment and stops cleanly after CET_HOURS (default
# 70); resubmitting the same array continues every unfinished chain from its
# checkpoint (finished chains exit at once).  Settings pass through the
# environment (see run-chain.R).  Ladder rungs and power-posterior chains:
# set CET_CELL (a label such as b0.5g1) and/or CET_POWER, and use
# --array=1-<reps>.
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-hym-cet/logs/%x-%A_%a.out
set +e; set +o pipefail
source /etc/profile > /dev/null 2>&1 || true
module purge || true
module load r/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=${W:-/nobackup/pjjg18/evotempo-hym-cet}
SRC=${SRC:-$W/src}
LIB=${LIB:-$W/lib-$(head -c 7 $W/out/build-sha.txt)}
export R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib
export CET_OUT=${CET_OUT:-$W/runs}
cd $SRC
t=$SLURM_ARRAY_TASK_ID
if [ -n "$CET_CELL" ]; then
  cell=$CET_CELL; rep=$t
else
  cell=$(( (t - 1) % 4 + 1 )); rep=$(( (t - 1) / 4 + 1 ))
fi
mkdir -p $W/out $W/logs
Rscript data-raw/empirical/hymenoptera/cet/run-chain.R $cell $rep \
  >> $W/out/run-$cell-$rep.txt 2>&1
echo "$(date) exit $?" >> $W/out/run-$cell-$rep.txt
