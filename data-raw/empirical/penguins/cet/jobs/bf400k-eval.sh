#!/bin/bash
# eval with CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20 on the 40 chains of the 400k
# results, via the link root bf400k-root (Slurm 19472716; code: a scratch
# worktree at 363af64 with the bayes-factors.R and common.R committed to #105
# as ee4b364).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH --time=0-02:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%A_%a.out
source /etc/profile || true; module purge || true; module load r/4.5.1 gcc/14.2
export R_LIBS=/nobackup/pjjg18/evotempo-penguin-cet/lib-48e6bf7:/nobackup/pjjg18/evotempo-dev/lib OMP_NUM_THREADS=1 CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20
cd /nobackup/pjjg18/evotempo-penguin-cet/tmp-cap
labels=(b0g0_run1 b0g1_run1 b1g0_run1 b1g1_run1 b0g0_run2 b0g1_run2 b1g0_run2 b1g1_run2 b0g0_run3 b0g1_run3 b1g0_run3 b1g1_run3 b0g0_run4 b0g1_run4 b1g0_run4 b1g1_run4 b0g0.25_run1 b0g0.25_run2 b0g0.25_run3 b0g0.25_run4 b0g0.5_run1 b0g0.5_run2 b0g0.5_run3 b0g0.5_run4 b0g0.75_run1 b0g0.75_run2 b0g0.75_run3 b0g0.75_run4 b1g0.25_run1 b1g0.25_run2 b1g0.25_run3 b1g0.25_run4 b1g0.5_run1 b1g0.5_run2 b1g0.5_run3 b1g0.5_run4 b1g0.75_run1 b1g0.75_run2 b1g0.75_run3 b1g0.75_run4)
lab=${labels[$((SLURM_ARRAY_TASK_ID-1))]}
Rscript data-raw/empirical/penguins/cet/bayes-factors.R eval /nobackup/pjjg18/evotempo-penguin-cet/bf400k-root/$lab; echo "exit $?"
