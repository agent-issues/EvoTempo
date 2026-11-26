#!/bin/bash
# Trilobite CET pipeline: the checkpoint's build and data guards (review of
# #114, R114-3, R114-6), on Hamilton's `test` partition.
#
#   sbatch -J evotempo-trilobite-cet-guards data-raw/empirical/trilobites/cet/check-guards.sh
#
# Cell (0, 0), three segments of 20 iterations, runs one segment
# (CET_HOURS=0 stops it there), then resumes:
#   A. under a copy of the install with another BUILD_SHA: must stop (evotempoSha);
#   B. from a copy of the trilobite data whose characters.csv has a changed
#      byte (a trailing blank line, which read.csv skips): must stop
#      (charactersMd5);
#   C. with the alpha prior's rate changed in common.R (copy): must stop
#      (alphaPrior; the code md5 alphaPriorCode catches edits to the
#      functions themselves);
#   F. under a reinstall of the same commit (a copy of the install whose
#      R/EvoTempo.rdb has one byte appended past the last object, so its md5
#      changes and nothing it reads does): must stop (evotempoMd5) without
#      CET_ALLOW_REINSTALL, and with CET_ALLOW_REINSTALL=1 resume, say so in
#      progress.txt, run segment 2 and stop (CET_HOURS=0) (#114 delta review,
#      N2);
#   D. under an identical copy of the install that wrote the last checkpoint
#      (F's): must resume and complete;
#   E. a chain written before Gamma rates (PR #91's smoke chain b0g0_run1 at
#      $OLD_CHAIN, copied): must stop, naming nRateCat and the new fields as
#      "(none)".
# Writes $OUT/guards.txt (default cet/smoke-results).
#SBATCH -p test
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --time=00:15:00
set +e; set +o pipefail
source /etc/profile || true
module purge || true
module load r/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=${W:-/nobackup/pjjg18/evotempo-trilobite-cet}
SRC=${SRC:-$W/src}
LIB=${LIB:-$W/lib-$(cat $W/out/build-sha.txt)}
DEV=/nobackup/pjjg18/evotempo-dev/lib
CET=data-raw/empirical/trilobites/cet
OUT=${OUT:-$SRC/$CET/smoke-results}
G=$W/guards
rm -rf $G; mkdir -p $G/runs
export CET_OUT=$G/runs CET_NITER=60 CET_BURNIN=10 CET_THIN=10 CET_SEG=20
F=$OUT/guards.txt
Run() {  # label, R_LIBS, source dir, [env]
  echo "== $1" >> $F
  (cd $3 && env R_LIBS=$2:$DEV $4 Rscript $CET/run-chain.R 1 1 > $G/run.txt 2>&1; echo "exit $?" >> $G/run.txt)
  grep -E "^Error|^exit" $G/run.txt | cut -c1-2000 >> $F
  tail -1 $G/runs/b0g0_run1/progress.txt >> $F
}
: > $F
Run "first process: one segment, CET_HOURS=0" $LIB $SRC CET_HOURS=0
cp -r $LIB $G/libA; echo 0000000 > $G/libA/EvoTempo/BUILD_SHA
Run "A: install with another BUILD_SHA (must stop: evotempoSha)" $G/libA $SRC
mkdir -p $G/srcB/data-raw/empirical $G/srcC/data-raw/empirical
cp -r $SRC/data-raw/empirical/trilobites $G/srcB/data-raw/empirical/
cp -r $SRC/data-raw/empirical/trilobites $G/srcC/data-raw/empirical/
echo "" >> $G/srcB/$CET/data/characters.csv
Run "B: characters.csv changed (must stop: charactersMd5)" $LIB $G/srcB
sed -i 's/^cetAlphaRate <- 2$/cetAlphaRate <- 2.5/' $G/srcC/$CET/common.R
grep -n "^cetAlphaRate" $G/srcC/$CET/common.R >> $F
Run "C: alpha prior rate 2.5 (must stop: alphaPrior)" $LIB $G/srcC
cp -r $LIB $G/libF; printf x >> $G/libF/EvoTempo/R/EvoTempo.rdb
Run "F1: reinstall of the same commit (must stop: evotempoMd5)" $G/libF $SRC
Run "F2: the same with CET_ALLOW_REINSTALL=1 (must resume)" $G/libF $SRC "CET_ALLOW_REINSTALL=1 CET_HOURS=0"
grep "CET_ALLOW_REINSTALL" $G/runs/b0g0_run1/progress.txt | cut -c1-400 >> $F
cp -r $G/libF $G/libD
Run "D: identical copy of the install that wrote the last checkpoint (must resume and complete)" $G/libD $SRC
OLD_CHAIN=${OLD_CHAIN:-/nobackup/pjjg18/evotempo-trilobite-cet/smoke/b0g0_run1}
mkdir -p $G/old; cp -r $OLD_CHAIN $G/old/
echo "== E: checkpoint from before Gamma rates (must stop: nRateCat (none))" >> $F
(cd $SRC && env R_LIBS=$LIB:$DEV CET_OUT=$G/old CET_NITER=300 CET_BURNIN=40 CET_THIN=10 CET_SEG=100 \
   Rscript $CET/run-chain.R 1 1 > $G/run.txt 2>&1; echo "exit $?" >> $G/run.txt)
grep -E "^Error|^exit" $G/run.txt | cut -c1-2000 >> $F
pass=0
# Each stop must name exactly the field(s) changed
Line() { grep -A1 "^== $1" $F | grep "^Error"; }
Line A | grep -q "differs in: evotempoSha (checkpoint [0-9a-f]*; this run 0000000)$" && pass=$((pass+1))
Line B | grep -q "differs in: charactersMd5 (checkpoint md5:[0-9a-f]*; this run md5:[0-9a-f]*)$" && pass=$((pass+1))
Line C | grep -q "differs in: alphaPrior (checkpoint alpha ~ Exponential(rate 2); 4 Gamma categories; this run alpha ~ Exponential(rate 2.5); 4 Gamma categories)$" && pass=$((pass+1))
grep -A3 "^== D" $F | grep -q "complete" && pass=$((pass+1))
Line F1 | grep -q "differs in: evotempoMd5 (checkpoint md5:[0-9a-f]*; this run md5:[0-9a-f]*)$" && pass=$((pass+1))
grep -A3 "^== F2" $F | grep -q "^exit 0" && grep -q "CET_ALLOW_REINSTALL: same version" $G/runs/b0g0_run1/progress.txt &&
  grep -q "segment 2/3 done" $G/runs/b0g0_run1/progress.txt && pass=$((pass+1))
Line E | grep -q "differs in: nRateCat (checkpoint (none); this run 4)" && pass=$((pass+1))
echo "check-guards.sh: $pass of 7 guards behave as required" >> $F
[ $pass = 7 ] && echo "PASS" >> $F || echo "FAIL" >> $F
rm -rf $G/libA $G/libD $G/libF $G/srcB $G/srcC $G/old
