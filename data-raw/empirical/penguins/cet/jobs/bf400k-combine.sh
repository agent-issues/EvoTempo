#!/bin/bash
# combine on bf400k-root with the same caps (Slurm 19472717).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=0-02:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%j.out
source /etc/profile || true; module purge || true; module load r/4.5.1 gcc/14.2
export R_LIBS=/nobackup/pjjg18/evotempo-penguin-cet/lib-48e6bf7:/nobackup/pjjg18/evotempo-dev/lib OMP_NUM_THREADS=1 CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20
cd /nobackup/pjjg18/evotempo-penguin-cet/tmp-cap
Rscript data-raw/empirical/penguins/cet/check-segcap.R > /nobackup/pjjg18/evotempo-penguin-cet/out/bf400k-check-segcap.txt 2>&1; echo "segcap $?"
Rscript data-raw/empirical/penguins/cet/check-labels.R > /nobackup/pjjg18/evotempo-penguin-cet/out/bf400k-check-labels.txt 2>&1; echo "check-labels $?"; tail -1 /nobackup/pjjg18/evotempo-penguin-cet/out/bf400k-check-labels.txt
Rscript data-raw/empirical/penguins/cet/bayes-factors.R combine /nobackup/pjjg18/evotempo-penguin-cet/bf400k-root > /nobackup/pjjg18/evotempo-penguin-cet/out/bf400k-combine.txt 2>&1; echo "combine $?"
