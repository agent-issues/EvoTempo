#!/bin/bash
# Production chains (16 tasks): Hamilton batch header, then run-array.sh
# (used for Slurm 19431480, 19448311 and 19466805).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --time=3-00:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%A_%a.out
# Wrapper: batch header required on Hamilton, then the canonical run-array.sh
set +e; set +o pipefail; source /etc/profile || true; module purge || true
module load r/4.5.1
echo "$(date) src $(git -C /nobackup/pjjg18/evotempo-penguin-cet/src rev-parse HEAD) task $SLURM_ARRAY_TASK_ID host $(hostname)"
source /nobackup/pjjg18/evotempo-penguin-cet/src/data-raw/empirical/penguins/cet/run-array.sh
