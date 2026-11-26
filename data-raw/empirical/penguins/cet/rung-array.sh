#!/bin/bash
# Model-index ladder rungs (CET_CELL = b<beta>g<gamma>), replicate = array task.
#
#   sbatch -J evotempo-penguin-cet-rung-b0g0.125 --array=1-4 \
#     --export=ALL,CET_CELL=b0g0.125,RUNG_SEED=540000,RUNG_START=dispersed,CET_NITER=600000 \
#     rung-array.sh
#
# RUNG_START = production (the rungs launched 2026-11-13): the start tree of the
#   gamma = 1 production chain with the same beta and replicate
#   (data/start-trees/start-tree-<k>.nwk, k = 4 (rep - 1) + 2 for beta = 0,
#   + 4 for beta = 1), and CetStart(rep), i.e. the start parameters of
#   production chain k = rep (#107 review, R107-15);
# RUNG_START = dispersed: the rung's own draw from the round-4 posterior,
#   data/start-trees/rungs/start-tree-<k>.nwk with k from rungs/draws.csv
#   (cell, rep), and CetStart(k) (k = 17-32, not used by any other chain).
# CET_SEED = RUNG_SEED + rep.  Other settings pass through to run-chain.R
# (CET_NITER, CET_THIN, CET_SEG).  Resubmitting continues from the checkpoint.
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
export CET_OUT=${CET_OUT:-$W/runs}
export CET_THIN=${CET_THIN:-100} CET_SEG=${CET_SEG:-10000}
cd $SRC
rep=$SLURM_ARRAY_TASK_ID
beta=$(echo $CET_CELL | sed 's/^b\([0-9.]*\)g.*$/\1/')
starts=data-raw/empirical/penguins/cet/data/start-trees
case ${RUNG_START:-production} in
  production)
    cell=$([ "$beta" = 0 ] && echo 2 || echo 4)
    k=$(( 4 * (rep - 1) + cell ))
    export CET_START_TREE=$starts/start-tree-$(printf %02d $k).nwk
    unset CET_START_INDEX;;
  dispersed)
    k=$(awk -F, -v c="\"$CET_CELL\"" -v r=$rep 'NR > 1 && ($2 == c || $2 == substr(c, 2, length(c) - 2)) && $3 == r {print $1}' $starts/rungs/draws.csv)
    [ -n "$k" ] || { echo "no draw for $CET_CELL rep $rep in $starts/rungs/draws.csv"; exit 2; }
    export CET_START_TREE=$starts/rungs/start-tree-$(printf %02d $k).nwk
    export CET_START_INDEX=$k;;
  *) echo "RUNG_START must be production or dispersed"; exit 2;;
esac
export CET_SEED=$(( RUNG_SEED + rep ))
echo "$(date) src $(git rev-parse HEAD) cell $CET_CELL rep $rep start $CET_START_TREE index ${CET_START_INDEX:-rep} seed $CET_SEED host $(hostname)"
Rscript data-raw/empirical/penguins/cet/run-chain.R $CET_CELL $rep >> $W/out/run-$CET_CELL-$rep.txt 2>&1
echo "$(date) exit $?"
