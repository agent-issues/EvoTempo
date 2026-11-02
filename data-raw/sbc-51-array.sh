#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH -o /nobackup/pjjg18/evotempo-sbc/logs/%x-%A_%a.out
# Independent calibration of the TempoMoves() schedule (#51): rung-6 SBC
# (data-raw/rung6-sbc-full.R, EVOTEMPO_SBC_SCHEDULE=default) in cell
# beta = 0, gamma = 1 of #21.  Array task i runs replicate g = 4 i - 2.
#   ARM=standard  tasks 1-200: g = 2, 6, ..., 798, i.e. every cell-2 dataset of
#                 #21.  g >= 202 (150 datasets) were never run under the new
#                 schedule; g <= 198 repeat the #19 run on the current code.
#                 Base 20 000 iterations, up to 3 ESS extensions, as #19.
#   ARM=long      tasks 1-50: g = 2, ..., 198, base 200 000 iterations (10x),
#                 no extensions.  Same seeds as the standard arm, so its first
#                 20 000 iterations repeat the standard chain exactly.
# Setup, once: a library of EvoTempo at the commit under test, in $W/lib-<sha>
# with the commit in $W/lib-<sha>/SHA, and this script and rung6-sbc-full.R in
# $W/scripts.  Submit:
#   sbatch --array=1-200 -t 04:00:00 -J sbc51-standard --export=ALL,ARM=standard,SHA=<sha> sbc-51-array.sh
#   sbatch --array=1-50 -t 12:00:00 -J sbc51-long --export=ALL,ARM=long,SHA=<sha> sbc-51-array.sh
# Cost, from the #19 run of g = 2..198 (1.26 core-minutes per 1000 iterations,
# at most 1.64): standard arm about 110 core-hours, longest task 2.2 h with
# three extensions; long arm about 210 core-hours, longest task 5.5 h.
# Summary: Rscript data-raw/sbc-51-summary.R <standard dir> <long dir>
set +e
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=/nobackup/pjjg18/evotempo-sbc
export R_LIBS=$W/lib-$SHA:/nobackup/pjjg18/evotempo-dev/lib
export EVOTEMPO_SHA=$(cat $W/lib-$SHA/SHA)
export EVOTEMPO_SBC_SCHEDULE=default
export EVOTEMPO_SBC_OUT=$W/out-51-$ARM
case "$ARM" in
  standard) ;;
  long)
    export EVOTEMPO_SBC_ITER=10000,200000,10000,200000
    export EVOTEMPO_SBC_MAXEXT=0
    ;;
  *) echo "ARM must be standard or long" >&2; exit 1 ;;
esac
cd $W
Rscript scripts/rung6-sbc-full.R reps $((4 * SLURM_ARRAY_TASK_ID - 2))
