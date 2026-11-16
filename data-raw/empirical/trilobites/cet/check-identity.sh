#!/bin/bash
# Trilobite CET pipeline: check-identity.R on Hamilton: the character-free
# path of TempoMCMC() (check-window.R, check-sbc.R) under the build that ran
# PR #91's checks ($LIB_OLD, 48e6bf7) and the build the chains now use ($LIB).
#
#   sbatch -J evotempo-trilobite-cet-identity data-raw/empirical/trilobites/cet/check-identity.sh
#
# check-window.R runs at 300 iterations (after its 1000 burn-in), 2 chains and
# 32 simulated trees per conditioning (CET_WINDOW_SEED=7, 2 cores) under each.
# Writes $OUT/identity.txt (default cet/smoke-results).
#SBATCH -p shared
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
set +e; set +o pipefail
source /etc/profile || true
module purge || true
module load r/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=${W:-/nobackup/pjjg18/evotempo-trilobite-gamma}
SRC=${SRC:-$W/src}
LIB=${LIB:-$W/lib-$(cat $W/out/build-sha.txt)}
LIB_OLD=${LIB_OLD:-/nobackup/pjjg18/evotempo-trilobite-cet/lib-48e6bf7}
CET=data-raw/empirical/trilobites/cet
OUT=${OUT:-$SRC/$CET/smoke-results}
I=$W/identity
rm -rf $I; mkdir -p $I/old/window $I/new/window $OUT
cd $SRC
for b in old new; do
  if [ $b = old ]; then L=$LIB_OLD; else L=$LIB; fi
  export R_LIBS=$L:/nobackup/pjjg18/evotempo-dev/lib
  Rscript $CET/check-identity.R chain $I/$b/chain.rds > $I/$b/chain.txt 2>&1
  EVOTEMPO_CORES=2 CET_WINDOW_SEED=7 CET_WINDOW_OUT=$I/$b/window \
    Rscript $CET/check-window.R 300 2 32 > $I/$b/window.txt 2>&1
done
Rscript $CET/check-identity.R compare $I/old $I/new > $OUT/identity.txt 2>&1
echo "exit $?" >> $OUT/identity.txt
