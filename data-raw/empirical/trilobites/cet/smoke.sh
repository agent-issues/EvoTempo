#!/bin/bash
# Trilobite CET pipeline: smoke run on Hamilton's `test` partition.
#
#   sbatch -J evotempo-trilobite-cet-smoke --array=1-10 smoke.sh
#
# From a checkout at $SRC with EvoTempo installed in $LIB and cet/data/ written
# by prepare.R.  Tasks 1-4: cells (0,0), (0,1), (1,0), (1,1), replicate 1,
# CET_NITER=200 after 40 burn-in iterations, in two segments of 100 (so the
# second resumes from the checkpoint).  Task 5: cell (1,1) again with one
# segment of 200; its last 10 samples must equal task 4's segment 2 (exact
# continuation; check-smoke.R).  Tasks 6 and 7: cells (1,1) and (0,0),
# replicate 2, as tasks 1-4: their start trees (draws 8 and 5) have other tip
# orders than replicate 1's, so check-smoke.R's eval of them covers the
# start-tree path of bayes-factors.R (review of #91, A3).  Tasks 8-10: resume
# in a second process, the predictive stop, and a kill before a checkpoint
# (below; review of #114).
# Logs: $W/out/<smoke dir>-<label>_run<rep>-<task>.txt.
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
task=$SLURM_ARRAY_TASK_ID; cell=$task; rep=1; runs=1
if [ "$task" = 5 ]; then cell=4; export CET_SEG=$CET_NITER CET_OUT=$CET_OUT-identity; fi
if [ "$task" = 6 ]; then cell=4; rep=2; fi
if [ "$task" = 7 ]; then cell=1; rep=2; fi
# Tasks 8-10 (review of #114, R114-1, R114-2): each runs run-chain.R in two or
# three Rscript processes, so each later one resumes from the checkpoint on disk.
# 8: cell (1, 1) as task 4 with CET_HOURS=0: the first process stops after
#    segment 1, the second resumes; check-smoke.R compares it with task 4 (one
#    process) and task 5 (one segment).
# 9: cell (0, 0), segments of 50, in three processes.  The first stops after
#    segment 1 (CET_HOURS=0) and its checkpoint records that segment's time,
#    segSecs.  The second resumes with a fake wall limit from that measured
#    time, CET_WALL_HOURS = (1.25 x segSecs + 5 s) / 3600 (margin 0): at its
#    start the time elapsed is a few seconds, so elapsed + 1.5 x segSecs passes
#    the limit and elapsed + 1.0 x segSecs does not (while elapsed < 0.25 x
#    segSecs + 5 s), whatever the node's speed.  It must therefore stop before
#    its first segment, by the predictive term with the segSecs carried over
#    in the checkpoint (a resubmitted job's first segment).  The third resumes
#    with the default limit and completes.  (#114 follow-up check: a fixed fake
#    limit of 0.043 h falsely failed unmodified code on slower nodes.)
# 10: as 9, but the first process is killed after writing segment 2's
#    acceptance.csv rows and before its checkpoint (CET_TEST_KILL_SEGMENT=2):
#    the second redoes segment 2; acceptance.csv must hold each segment once,
#    and the samples must equal task 9's.
case $task in
  8) cell=4; runs=2; export CET_OUT=$CET_OUT-twoproc ;;
  9) cell=1; runs=3; export CET_SEG=50 CET_OUT=$CET_OUT-predict ;;
  10) cell=1; runs=2; export CET_SEG=50 CET_OUT=$CET_OUT-kill ;;
esac
labels=(x b0g0 b0g1 b1g0 b1g1)
sfx=""; if [ -n "$CET_MKV" ] && [ "$CET_MKV" != pattern ]; then sfx=_mkv$CET_MKV; fi
rm -rf $CET_OUT/${labels[$cell]}${sfx}_run$rep
log=$W/out/$(basename $CET_OUT)-${labels[$cell]}${sfx}_run$rep-$task
rm -f $log.txt $log.err
for run in $(seq 1 $runs); do
  first=""
  if [ "$run" = 1 ]; then
    case $task in
      8|9) first="CET_HOURS=0" ;;
      10) first="CET_TEST_KILL_SEGMENT=2" ;;
    esac
  fi
  if [ "$task" = 9 ] && [ "$run" = 2 ]; then
    secs=$(Rscript -e "cat(readRDS('$CET_OUT/${labels[$cell]}${sfx}_run$rep/checkpoint.rds')\$segSecs)")
    wall=$(awk -v s="$secs" 'BEGIN { printf "%.8f", (1.25 * s + 5) / 3600 }')
    echo "segSecs from the checkpoint: $secs; CET_WALL_HOURS = (1.25 x $secs + 5) / 3600 = $wall" >> $log.txt
    first="CET_WALL_HOURS=$wall CET_WALL_MARGIN=0"
  fi
  echo "process $run of $runs: $first" >> $log.txt
  env $first /usr/bin/time -v Rscript data-raw/empirical/trilobites/cet/run-chain.R $cell $rep \
    >> $log.txt 2>> $log.err
  echo "exit $?" >> $log.txt
done
