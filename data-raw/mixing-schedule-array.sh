#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 02:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-sbc/logs/%x-%A_%a.out
# Dispersed-start arm of the schedule study (#51; data-raw/mixing-schedule.R):
# 48 tasks = 2 trees x 3 character counts x 2 schedules x 4 replicates, each
# 2000 + 5000 iterations from an overdispersed start.  Needs trees.rds from
# `Rscript dev/mixing/run.R sim` in $W/mixing, and a library as for
# sbc-51-array.sh.  About 10 core-hours (README section 3 timings); the
# longest task, large tree at 300 characters under TempoMoves(), about 35 min.
#   sbatch --array=1-48 -J mixing51 --export=ALL,SHA=<sha> mixing-schedule-array.sh
# Then: EVOTEMPO_MIXING_OUT=$W/mixing Rscript scripts/mixing-schedule.R summary
set +e
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=/nobackup/pjjg18/evotempo-sbc
export R_LIBS=$W/lib-$SHA:/nobackup/pjjg18/evotempo-dev/lib
export EVOTEMPO_SHA=$(cat $W/lib-$SHA/SHA)
export EVOTEMPO_MIXING_OUT=$W/mixing EVOTEMPO_MIXING_TREES=$W/mixing/trees.rds
i=$((SLURM_ARRAY_TASK_ID - 1))
tree=$((i / 24 + 1))
chars=(30 100 300)
nChar=${chars[$((i / 8 % 3))]}
configs=(new old)
config=${configs[$((i / 4 % 2))]}
rep=$((i % 4 + 1))
cd $W
Rscript scripts/mixing-schedule.R chain $tree $nChar $config $rep 5000 2000 dispersed
