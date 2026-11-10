#!/bin/bash
#SBATCH -p shared
#SBATCH -c 16
#SBATCH --mem=12G
#SBATCH -t 04:00:00
#SBATCH -J evotempo-bias72-chains
#SBATCH -o /nobackup/pjjg18/evotempo-bias72/logs/slurm-%A_%a.out
set +e
source /etc/profile; module purge; module load r/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS=/nobackup/pjjg18/evotempo-bias72/lib:/nobackup/pjjg18/evotempo-dev/lib
export EVOTEMPO_SHA=dc47aedffd36ee4e7dbba64da6fcfff487db31ee
read arm case offset <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" /nobackup/pjjg18/evotempo-bias72/tasks.txt)"
export EVOTEMPO_RUNG5_SEED_OFFSET=$offset
export EVOTEMPO_RUNG5_OUT=/nobackup/pjjg18/evotempo-bias72/out/$arm
if [ "$arm" = B ]; then export EVOTEMPO_VERIFY_PRIOR=gamma; script=oracle-gamma.R
else unset EVOTEMPO_VERIFY_PRIOR; script=rung5-age-oracle.R; fi
cd /nobackup/pjjg18/evotempo-bias72
echo "task $SLURM_ARRAY_TASK_ID arm $arm case $case offset $offset script $script host $(hostname) start $(date -u +%FT%TZ)"
Rscript /nobackup/pjjg18/evotempo-bias72/scripts/$script $case 60000 32 extant power > /nobackup/pjjg18/evotempo-bias72/logs/run-$arm-$offset.txt 2>&1
echo "exit $? end $(date -u +%FT%TZ)"
