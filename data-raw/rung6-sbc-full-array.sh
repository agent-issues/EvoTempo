#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH -t 04:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-sbc/logs/%x-%A_%a.out
# One Slurm array task = NPER consecutive replicates of data-raw/rung6-sbc-full.R
set +e
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=/nobackup/pjjg18/evotempo-sbc
export R_LIBS=$W/lib-a2c3aa1:/nobackup/pjjg18/evotempo-dev/lib
export EVOTEMPO_SHA=$(cat $W/lib-a2c3aa1/SHA)
cd $W
Rscript scripts/rung6-sbc-full.R task $SLURM_ARRAY_TASK_ID ${NPER:-1}
