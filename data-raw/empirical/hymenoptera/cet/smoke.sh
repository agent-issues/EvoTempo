#!/bin/bash
# Hymenoptera CET pipeline: smoke run on Hamilton (as the penguin pipeline's
# smoke.sh).
#
#   sbatch -J evotempo-hymcet-smoke --array=1-5 smoke.sh
#
# From a checkout at $SRC with EvoTempo installed in $LIB.  Tasks 1-4: cells
# (0,0), (0,1), (1,0), (1,1), CET_NITER (default 40) sampling iterations after
# CET_BURNIN (default 10), thin 2, in two segments (so the second resumes
# from the checkpoint).  Task 5: cell (1,1) again in one segment; its last
# samples must equal task 4's segment 2 (exact continuation; check-smoke.R).
# The defaults are smaller than the penguin smoke (200 + 40) because an
# iteration here costs several seconds; SMOKE_* override them, e.g. a pilot:
#   --export=ALL,SMOKE_NITER=400,SMOKE_BURNIN=100,SMOKE_SEG=200,SMOKE_THIN=4,SMOKE_OUT=$W/pilot
# Partition and time: `test` (15 min) fits the default smoke; use
# `-p shared -t <hours>` for a pilot.
#SBATCH -p test
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH -o /nobackup/pjjg18/evotempo-hym-cet/logs/%x-%A_%a.out
set +e; set +o pipefail
source /etc/profile > /dev/null 2>&1 || true
module purge || true
module load r/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=${W:-/nobackup/pjjg18/evotempo-hym-cet}
SRC=${SRC:-$W/src}
LIB=${LIB:-$W/lib-$(head -c 7 $W/out/build-sha.txt)}
export R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib
cd $SRC
export CET_NITER=${SMOKE_NITER:-40} CET_BURNIN=${SMOKE_BURNIN:-10} CET_THIN=${SMOKE_THIN:-2}
export CET_SEG=${SMOKE_SEG:-20} CET_OUT=${SMOKE_OUT:-$W/smoke}
cell=$SLURM_ARRAY_TASK_ID
if [ "$cell" = 5 ]; then cell=4; export CET_SEG=$CET_NITER CET_OUT=$CET_OUT-identity; fi
labels=(x b0g0 b0g1 b1g0 b1g1)
sfx=""; if [ -n "$CET_MKV" ] && [ "$CET_MKV" != pattern ]; then sfx=_mkv$CET_MKV; fi
rm -rf $CET_OUT/${labels[$cell]}${sfx}_run1
mkdir -p $W/out $W/logs
log=$W/out/$(basename $CET_OUT)-$SLURM_ARRAY_TASK_ID
/usr/bin/time -v Rscript data-raw/empirical/hymenoptera/cet/run-chain.R $cell 1 > $log.txt 2> $log.err
echo "exit $?" >> $log.txt
