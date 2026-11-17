#!/bin/bash
# Trilobite CET pipeline: free-parameter SBC-lite (check-sbc.R) on Hamilton.
#
#   sbatch -J evotempo-trilobite-cet-sbc data-raw/empirical/trilobites/cet/check-sbc.sh
#
# SBC_ARGS passes [nRep] [nIter] (default 200 19800).
#SBATCH -p shared
#SBATCH -c 16
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-trilobite-cet/logs/%x-%j.out
set +e; set +o pipefail
source /etc/profile || true
module purge || true
module load r/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=${W:-/nobackup/pjjg18/evotempo-trilobite-cet}
SRC=${SRC:-$W/src}
LIB=${LIB:-$W/lib-$(cat $W/out/build-sha.txt)}
export R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib
export EVOTEMPO_CORES=${SLURM_CPUS_PER_TASK:-16}
export CET_SBC_OUT=${CET_SBC_OUT:-$SRC/data-raw/empirical/trilobites/cet/sbc-check}
cd $SRC
/usr/bin/time -v Rscript data-raw/empirical/trilobites/cet/check-sbc.R $SBC_ARGS \
  > $W/out/sbc-${SLURM_JOB_ID}.txt 2>&1
echo "exit $?" >> $W/out/sbc-${SLURM_JOB_ID}.txt
