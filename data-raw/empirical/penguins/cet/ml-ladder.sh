#!/bin/bash
# Penguin CET pipeline: character-power ladder for the common-footing marginal
# likelihood (README, "Common-footing marginal likelihoods").
#
#   sbatch -J evotempo-penguin-cet-ml-b0g0 --array=1-128 \
#     --export=ALL,CET_CELL=b0g0,ML_SEED=700000 ml-ladder.sh
#
# Task i runs power index j = (i - 1) % ML_NPOW, replicate (i - 1) / ML_NPOW + 1,
# at character power t_j = (j / ML_NPOW)^5, j = 0 .. ML_NPOW - 1 (default
# ML_NPOW = 64).  These are exactly the powers of the RevBayes comparator's
# 65-stone schedule (cats = 65, alpha = 0.2: ((65 - k) / 64)^5), so the two
# programs can be compared on the same ladder (and on its every second power).  t = 1 is the cell's
# production chains, which `bayes-factors.R combine` finds under the same run
# root (symbolic links; ml-root.sh).
# Start tree: the production start tree of the same cell and replicate
# (start-tree-<k>, k = 4 (rep - 1) + cell index); CET_SEED = ML_SEED + 100 j + rep.
# Settings as the production chains except CET_NITER (default 100000).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --time=3-00:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%A_%a.out
set +e; set +o pipefail; source /etc/profile || true; module purge || true
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=${W:-/nobackup/pjjg18/evotempo-penguin-cet}
SRC=${SRC:-$W/src}
LIB=${LIB:-$W/lib-$(cat $W/out/build-sha.txt)}
export R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib
export CET_OUT=${CET_OUT:-$W/ml-runs}
export CET_NITER=${CET_NITER:-100000} CET_THIN=${CET_THIN:-100} CET_SEG=${CET_SEG:-10000}
NPOW=${ML_NPOW:-64}
case $CET_CELL in
  b0g0) k0=1;; b0g1) k0=2;; b1g0) k0=3;; b1g1) k0=4;;
  *) echo "CET_CELL must be one of the production cells, not '$CET_CELL'"; exit 2;;
esac
i=$SLURM_ARRAY_TASK_ID
j=$(( (i - 1) % NPOW )); rep=$(( (i - 1) / NPOW + 1 ))
export CET_POWER=$(Rscript -e "cat(format(($j / $NPOW)^5, digits = 17, scientific = FALSE))")
export CET_SEED=$(( ML_SEED + 100 * j + rep ))
export CET_START_TREE=data-raw/empirical/penguins/cet/data/start-trees/start-tree-$(printf %02d $(( 4 * (rep - 1) + k0 ))).nwk
cd $SRC
echo "$(date) src $(git rev-parse HEAD) cell $CET_CELL rep $rep power $CET_POWER seed $CET_SEED start $CET_START_TREE host $(hostname)"
Rscript data-raw/empirical/penguins/cet/run-chain.R $CET_CELL $rep \
  >> $W/out/ml-$CET_CELL-t$j-$rep.txt 2>&1
echo "$(date) exit $?"
