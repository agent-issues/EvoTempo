#!/bin/bash
# Validation of #97 on Hamilton (Slurm), from the package root, with EvoTempo
# installed in $LIB and its dependencies in $DEV:
#
#   bash data-raw/gamma-acrv-hamilton.sh rb <scratch>          # validation 1
#   bash data-raw/gamma-acrv-hamilton.sh hymenoptera <scratch> # validation 1, real matrix
#   bash data-raw/gamma-acrv-hamilton.sh oracle <scratch> <gamma>  # validation 2 (-c 16)
#   bash data-raw/gamma-acrv-hamilton.sh sbc <scratch> <task> <nTask>  # validation 3
#   bash data-raw/gamma-acrv-hamilton.sh sbc-summary <scratch>
#
# Results go to data-raw/gamma-acrv-results/.  RevBayes: the development
# build fd807d10 (gcc/14.2 boost/1.88.0); RB=<path to rb> overrides it.
set +e; set +o pipefail
source /etc/profile > /dev/null 2>&1 || true
module purge || true
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
DEV=${DEV:-/nobackup/pjjg18/evotempo-dev/lib}
export R_LIBS=$LIB:$DEV
RB=${RB:-/nobackup/pjjg18/revbayes-dev-fd807d10/projects/cmake/build/rb}
RES=data-raw/gamma-acrv-results
mkdir -p $RES
what=$1
SCR=$2
mkdir -p $SCR
case $what in
  rb)
    Rscript data-raw/gamma-acrv-revbayes.R prepare $SCR/data > $RES/revbayes-prepare.txt 2>&1
    cp data-raw/gamma-acrv-revbayes.Rev $SCR/data/
    (module purge || true; module load gcc/14.2 boost/1.88.0
     cd $SCR/data && $RB --version > $SCR/rb-version.txt 2>&1
     $RB gamma-acrv-revbayes.Rev > $SCR/revbayes-output.txt 2>&1)
    cp $SCR/revbayes-output.txt $RES/revbayes-output.txt
    head -3 $SCR/rb-version.txt > $RES/revbayes-version.txt
    Rscript data-raw/gamma-acrv-revbayes.R compare $SCR/data $RES/revbayes-output.txt \
      $RES/revbayes.csv > $RES/revbayes-compare.txt 2>&1
    ;;
  hymenoptera)
    # Matrix, start tree and RevBayes oracle of the Hymenoptera CET pipeline
    rm -rf $SCR/hym
    git clone -q https://github.com/agent-issues/EvoTempo $SCR/hym
    (cd $SCR/hym && git checkout -q ${HYM_SHA:-e0c1384} && git rev-parse HEAD) > $RES/hymenoptera-sha.txt
    CET=$SCR/hym/data-raw/empirical/hymenoptera/cet
    Rscript data-raw/gamma-acrv-revbayes.R hymenoptera $CET \
      $CET/smoke-results/revbayes-check-mkv-fd807d10.txt $RES/hymenoptera.csv \
      > $RES/hymenoptera-compare.txt 2>&1
    rm -rf $SCR/hym
    ;;
  oracle)
    gamma=$3
    [ -d $SCR/data ] || Rscript data-raw/gamma-acrv-revbayes.R prepare $SCR/data > /dev/null 2>&1
    Rscript data-raw/gamma-acrv-alpha-oracle.R $SCR/data $gamma $RES/alpha-oracle-gamma$gamma \
      > $RES/alpha-oracle-gamma$gamma.txt 2>&1
    ;;
  sbc)
    Rscript data-raw/gamma-acrv-sbc.R run $3 $4 $SCR ${SBC_REPS:-100} ${SBC_ITER:-3000} \
      > $SCR/sbc-task$3.txt 2>&1
    ;;
  sbc-summary)
    Rscript data-raw/gamma-acrv-sbc.R summary $SCR > $RES/sbc-summary.txt 2>&1
    cp $SCR/ranks.csv $SCR/sbc-tests.csv $RES/
    ;;
esac
