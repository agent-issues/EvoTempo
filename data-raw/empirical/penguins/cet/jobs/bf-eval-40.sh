#!/bin/bash
# eval on the 16 production + 24 rung chains (Slurm 19449553).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH --time=0-04:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%A_%a.out
set +e; set +o pipefail; source /etc/profile || true; module purge || true
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=/nobackup/pjjg18/evotempo-penguin-cet
export R_LIBS=$W/lib-48e6bf7:$W/lib-tools:/nobackup/pjjg18/evotempo-dev/lib
cd $W/src-105
labels=(b0g0_run1 b0g1_run1 b1g0_run1 b1g1_run1 b0g0_run2 b0g1_run2 b1g0_run2 b1g1_run2 b0g0_run3 b0g1_run3 b1g0_run3 b1g1_run3 b0g0_run4 b0g1_run4 b1g0_run4 b1g1_run4 b0g0.25_run1 b0g0.25_run2 b0g0.25_run3 b0g0.25_run4 b0g0.5_run1 b0g0.5_run2 b0g0.5_run3 b0g0.5_run4 b0g0.75_run1 b0g0.75_run2 b0g0.75_run3 b0g0.75_run4 b1g0.25_run1 b1g0.25_run2 b1g0.25_run3 b1g0.25_run4 b1g0.5_run1 b1g0.5_run2 b1g0.5_run3 b1g0.5_run4 b1g0.75_run1 b1g0.75_run2 b1g0.75_run3 b1g0.75_run4)
lab=${labels[$((SLURM_ARRAY_TASK_ID-1))]}
test -f $W/runs/$lab/seg0040_samples.csv -o -f $W/runs/$lab/seg0020_samples.csv || echo "WARNING $lab looks incomplete"
echo "$(date) src $(git rev-parse HEAD) eval $lab"
/usr/bin/time -v Rscript data-raw/empirical/penguins/cet/bayes-factors.R eval $W/runs/$lab
echo "$(date) exit $?"
