#!/bin/bash
# stage1-results.R (this commit) at CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20 on
# bf400k-root (Slurm 19472732).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=0-04:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%j.out
source /etc/profile || true; module purge || true; module load r/4.5.1 gcc/14.2
export R_LIBS=/nobackup/pjjg18/evotempo-penguin-cet/lib-48e6bf7:/nobackup/pjjg18/evotempo-penguin-cet/lib-tools:/nobackup/pjjg18/evotempo-dev/lib W=/nobackup/pjjg18/evotempo-penguin-cet CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20
cd /nobackup/pjjg18/evotempo-penguin-cet/src
O=/nobackup/pjjg18/evotempo-penguin-cet/stage1-results-400k-v3; rm -rf $O
/usr/bin/time -f '%e s %M KB' Rscript /nobackup/pjjg18/evotempo-penguin-cet/jobs/stage1-results-v3.R /nobackup/pjjg18/evotempo-penguin-cet/bf400k-root /nobackup/pjjg18/evotempo-penguin-cet/rbmorph-fd807d10 $O 0.25 > $O.log 2>&1; echo "results exit $?" >> $O.log
Rscript /nobackup/pjjg18/evotempo-penguin-cet/jobs/stage1-figures-v3.R $O $O/figures >> $O.log 2>&1; echo "figures exit $?" >> $O.log
