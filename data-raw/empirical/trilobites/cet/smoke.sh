#!/bin/bash
# Trilobite CET pipeline: smoke run on Hamilton's `test` partition.
#
#   sbatch -J evotempo-trilobite-cet-smoke --array=1-5 smoke.sh
#
# From a checkout at $SRC with EvoTempo installed in $LIB and cet/data/ written
# by prepare.R.  Tasks 1-4: cells (0,0), (0,1), (1,0), (1,1), replicate 1,
# CET_NITER=200 after 40 burn-in iterations, in two segments of 100 (so the
# second resumes from the checkpoint).  Task 5: cell (1,1) again with one
# segment of 200; its last 10 samples must equal task 4's segment 2 (exact
# continuation; check-smoke.R).  Logs: $W/out/<smoke dir>-<label>-<task>.txt.
#SBATCH -p test
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --time=00:15:00
#SBATCH -o /nobackup/pjjg18/evotempo-trilobite-cet/logs/%x-%A_%a.out
set +e; set +o pipefail
source /etc/profile || true
module purge || true
module load r/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=${W:-/nobackup/pjjg18/evotempo-trilobite-cet}
SRC=${SRC:-$W/src}
LIB=${LIB:-$W/lib-$(cat $W/out/build-sha.txt)}
export R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib
cd $SRC
# SMOKE_* override the defaults, e.g. a longer pilot (tasks 1-4, `shared`):
#   -p shared --time=02:00:00 --export=ALL,SMOKE_NITER=3000,SMOKE_BURNIN=500,SMOKE_SEG=1000,SMOKE_OUT=$W/pilot
export CET_NITER=${SMOKE_NITER:-200} CET_BURNIN=${SMOKE_BURNIN:-40} CET_THIN=10
export CET_SEG=${SMOKE_SEG:-100} CET_OUT=${SMOKE_OUT:-$W/smoke}
cell=$SLURM_ARRAY_TASK_ID
if [ "$cell" = 5 ]; then cell=4; export CET_SEG=$CET_NITER CET_OUT=$CET_OUT-identity; fi
labels=(x b0g0 b0g1 b1g0 b1g1)
sfx=""; if [ -n "$CET_MKV" ] && [ "$CET_MKV" != pattern ]; then sfx=_mkv$CET_MKV; fi
rm -rf $CET_OUT/${labels[$cell]}${sfx}_run1
log=$W/out/$(basename $CET_OUT)-${labels[$cell]}${sfx}-$SLURM_ARRAY_TASK_ID
/usr/bin/time -v Rscript data-raw/empirical/trilobites/cet/run-chain.R $cell 1 > $log.txt 2> $log.err
echo "exit $?" >> $log.txt
