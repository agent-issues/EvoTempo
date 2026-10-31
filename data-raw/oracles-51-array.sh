#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 01:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-sbc/logs/%x-%A_%a.out
# Long runs of the `s` oracle (#51; data-raw/oracles-51.R) with and without
# `shift`, to settle the lean of its upper quartiles seen locally (README
# section 4).  40 tasks: seeds 301-320 of `s` (odd tasks) and of `sNoShift`
# (even tasks), each 2 chains x 50 000 iterations after 2000 burn-in, about
# 6 core-minutes each, 4 core-hours in all.  Run from a checkout of the branch
# in $W/src, with the library of sbc-51-array.sh:
#   sbatch --array=1-40 -J oracles51 --export=ALL,SHA=<sha> data-raw/oracles-51-array.sh
# Rows are appended to data-raw/mixing-results/oracles-51-{s,sNoShift}.csv.
set +e
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=/nobackup/pjjg18/evotempo-sbc
export R_LIBS=$W/lib-$SHA:/nobackup/pjjg18/evotempo-dev/lib
oracles=(sNoShift s)
oracle=${oracles[$((SLURM_ARRAY_TASK_ID % 2))]}
seed=$((300 + (SLURM_ARRAY_TASK_ID + 1) / 2))
cd $W/src
Rscript data-raw/oracles-51.R $oracle none $seed 50000 2 2000
