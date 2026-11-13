#!/bin/bash
# Trilobite CET pipeline: production chains as a Slurm array on Hamilton.
#
#   sbatch -J evotempo-trilobite-cet --array=1-16 run-array.sh
#
# Task t runs cell (t - 1) %% 4 + 1 (cells (0,0), (0,1), (1,0), (1,1)),
# replicate (t - 1) %/% 4 + 1, as chain k = t, which starts from draw k of
# prepare.R: one of 16 seeded random draws without replacement from the pooled
# post-burn-in RevBayes posterior (start-tree-<kk>.nwk and its row of
# start-values.csv; README, "Starting trees").  Each task runs run-chain.R, which
# writes a checkpoint after every segment and stops cleanly after CET_HOURS
# (default 70); resubmitting the same array continues every unfinished chain
# from its checkpoint (finished chains exit at once).  Settings pass through
# the environment (see run-chain.R), e.g.
#   sbatch -J evotempo-trilobite-cet --array=1-16 \
#     --export=ALL,CET_NITER=200000,CET_THIN=100,CET_SEG=10000 run-array.sh
# Ladder rungs and power-posterior chains: set CET_CELL (a label such as
# b0.5g1) and/or CET_POWER, and use --array=1-<reps>; add a row for them to
# the launch record $CET_OUT/launch.csv (cell,reps[,power]; README, "Running on
# Hamilton"), or bayes-factors.R combine stops on chains it does not expect
# (#101).  Logs:
# $W/out/run-<label>-<rep>.txt, with the Mkv variant and power in the label.
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --time=3-00:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-trilobite-cet/logs/%x-%A_%a.out
set +e; set +o pipefail
source /etc/profile || true
module purge || true
module load r/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=${W:-/nobackup/pjjg18/evotempo-trilobite-cet}
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
tag=$cell
if [ -n "$CET_MKV" ] && [ "$CET_MKV" != pattern ]; then tag=${tag}_mkv$CET_MKV; fi
if [ -n "$CET_POWER" ] && [ "$CET_POWER" != 1 ]; then tag=${tag}_t$CET_POWER; fi
Rscript data-raw/empirical/trilobites/cet/run-chain.R $cell $rep >> $W/out/run-$tag-$rep.txt 2>&1
echo "$(date) exit $?" >> $W/out/run-$tag-$rep.txt
