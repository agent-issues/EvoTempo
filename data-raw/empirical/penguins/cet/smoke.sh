#!/bin/bash
# Penguin CET pipeline: smoke run on Hamilton's `test` partition.
#
#   sbatch -J evotempo-penguin-cet-smoke --array=1-5 smoke.sh
#
# From a checkout of the branch at $SRC with EvoTempo installed in $LIB.
# Tasks 1-4: cells (0,0), (0,1), (1,0), (1,1), CET_NITER=200 after 40 burn-in
# iterations, in two segments of 100 (so the second resumes from the
# checkpoint).  Task 5: cell (1,1) again with one segment of 200; its last 10
# samples must equal task 4's segment 2 (exact continuation; check-smoke.R).
#SBATCH -p test
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --time=00:15:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%A_%a.out
set +e; set +o pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
module purge; module load r/4.5.1 gcc/14.2
W=${W:-/nobackup/pjjg18/evotempo-penguin-cet}
SRC=${SRC:-$W/src}
LIB=${LIB:-$W/lib-$(cat $W/out/build-sha.txt)}
export R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib
cd $SRC
# SMOKE_* override the defaults, e.g. the 12-minute mixing pilot (tasks 1-4):
#   --export=ALL,SMOKE_NITER=3000,SMOKE_BURNIN=500,SMOKE_SEG=1000,SMOKE_OUT=$W/smoke-long
export CET_NITER=${SMOKE_NITER:-200} CET_BURNIN=${SMOKE_BURNIN:-40} CET_THIN=10
export CET_SEG=${SMOKE_SEG:-100} CET_OUT=${SMOKE_OUT:-$W/smoke}
cell=$SLURM_ARRAY_TASK_ID
if [ "$cell" = 5 ]; then cell=4; export CET_SEG=$CET_NITER CET_OUT=$CET_OUT-identity; fi
labels=(x b0g0 b0g1 b1g0 b1g1)
rm -rf $CET_OUT/${labels[$cell]}_run1
log=$W/out/$(basename $CET_OUT)-$SLURM_ARRAY_TASK_ID
/usr/bin/time -v Rscript data-raw/empirical/penguins/cet/run-chain.R $cell 1 > $log.txt 2> $log.err
echo "exit $?" >> $log.txt
