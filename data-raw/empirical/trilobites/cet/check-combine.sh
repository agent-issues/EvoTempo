#!/bin/bash
# Trilobite CET pipeline: bayes-factors.R eval and combine on tiny real chains
# with a ladder rung present (#100), combine's expected-chain-set check
# (#101), and the start-tree resume guard (#102), on Hamilton's `test`
# partition.  Three stages, each a Slurm job:
#
#   export CHECK_ROOT=<scratch dir>          # default $W/check-combine
#   a=$(sbatch --parsable -J evotempo-check-combine-chains --array=1-10 check-combine.sh chain)
#   e=$(sbatch --parsable -J evotempo-check-combine-eval --array=1-10 -d afterany:$a check-combine.sh eval)
#   sbatch -J evotempo-check-combine-check -d afterany:$e check-combine.sh check
#
# Tasks 1-10: b0g0, b0g1, b1g0, b1g1 and the rung b0.5g0, replicates 1 and 2
# (task t: label (t - 1) %% 5 + 1, replicate (t - 1) %/% 5 + 1), each 40
# iterations after 10 burn-in, thinned to 4 states, in $CHECK_ROOT/runs.
# `eval` runs bayes-factors.R eval on its chain (betas and gammas found under
# the root, as in production).  If OLD_BF names an earlier bayes-factors.R
# (e.g. main before #100: git show <sha>:<path> > $CHECK_ROOT/bayes-factors-old.R),
# it also copies the chain to $CHECK_ROOT/old-runs and runs that script's eval
# there (CET_BF_BETAS=0,0.5,1 so that it need not wait for the other copies),
# so that check-combine.R can show the old failure.  `check` runs
# check-combine.R, which writes $CHECK_ROOT/check-combine.txt.
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
export R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib
ROOT=${CHECK_ROOT:-$W/check-combine}
OLD_BF=${OLD_BF:-$ROOT/bayes-factors-old.R}
cd $SRC
export CET_NITER=40 CET_BURNIN=10 CET_THIN=10 CET_SEG=40 CET_OUT=$ROOT/runs
labels=(x b0g0 b0g1 b1g0 b1g1 b0.5g0)
t=${SLURM_ARRAY_TASK_ID:-1}
i=$(( (t - 1) % 5 + 1 )); rep=$(( (t - 1) / 5 + 1 ))
chain=${labels[$i]}_run$rep
mkdir -p $ROOT/logs $CET_OUT
case "$1" in
  chain)
    if [ "$i" = 5 ]; then cell=b0.5g0; else cell=$i; fi
    rm -rf $CET_OUT/$chain
    Rscript data-raw/empirical/trilobites/cet/run-chain.R $cell $rep > $ROOT/logs/chain-$chain.txt 2>&1
    echo "exit $?" >> $ROOT/logs/chain-$chain.txt
    ;;
  eval)
    Rscript data-raw/empirical/trilobites/cet/bayes-factors.R eval $CET_OUT/$chain > $ROOT/logs/eval-$chain.txt 2>&1
    echo "exit $?" >> $ROOT/logs/eval-$chain.txt
    if [ -f "$OLD_BF" ]; then
      mkdir -p $ROOT/old-runs
      rm -rf $ROOT/old-runs/$chain
      cp -r $CET_OUT/$chain $ROOT/old-runs/$chain
      rm -f $ROOT/old-runs/$chain/cross.csv
      CET_BF_BETAS=0,0.5,1 CET_BF_GAMMAS=0,1 Rscript $OLD_BF eval $ROOT/old-runs/$chain \
        > $ROOT/logs/old-eval-$chain.txt 2>&1
      echo "exit $?" >> $ROOT/logs/old-eval-$chain.txt
    fi
    ;;
  check)
    Rscript data-raw/empirical/trilobites/cet/check-combine.R $ROOT $OLD_BF > $ROOT/check-combine.log 2>&1
    echo "exit $?" >> $ROOT/check-combine.log
    ;;
  *)
    echo "usage: check-combine.sh chain|eval|check" >&2
    exit 2
    ;;
esac
