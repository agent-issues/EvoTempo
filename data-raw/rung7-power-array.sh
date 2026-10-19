#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH -t 12:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-power/logs/%x-%A_%a.out
# One Slurm array task of data-raw/rung7-power.R.
#   MODE=pilot sbatch --array=1-7 data-raw/rung7-power-array.sh
#   MODE=task NPER=1 sbatch --array=1-N data-raw/rung7-power-array.sh
set +e
set +o pipefail
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=/nobackup/pjjg18/evotempo-power
L=/nobackup/pjjg18/evotempo-sbc/lib-a2c3aa1
export R_LIBS=$L:/nobackup/pjjg18/evotempo-dev/lib
export EVOTEMPO_SHA=$(cat $L/SHA)
export EVOTEMPO_POWER_OUT=${EVOTEMPO_POWER_OUT:-$W/out}
cd $W/EvoTempo
if [ "${MODE:-task}" = "pilot" ]; then
  Rscript data-raw/rung7-power.R pilot $SLURM_ARRAY_TASK_ID
else
  Rscript data-raw/rung7-power.R task $SLURM_ARRAY_TASK_ID ${NPER:-1}
fi
