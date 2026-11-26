#!/bin/bash
# Superseded by ../rung-array.sh (RUNG_START=production does the same);
# kept as the record of the 2026-11-13 rung launch (Slurm 19448246-19448251).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --time=3-00:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%A_%a.out
# Model-index ladder rung (CET_CELL = b<beta>g<gamma>), replicate = array task.
# Start tree: the dispersed start of the gamma = 1 production chain with the
# same beta and replicate; seed distinct per rung and replicate (#68, #93).
set +e; set +o pipefail; source /etc/profile || true; module purge || true
module load r/4.5.1
rep=$SLURM_ARRAY_TASK_ID
beta=$(echo $CET_CELL | sed 's/^b\([0-9.]*\)g.*$/\1/')
cell=$([ "$beta" = 0 ] && echo 2 || echo 4)
k=$(( 4 * (rep - 1) + cell ))
export CET_START_TREE=data-raw/empirical/penguins/cet/data/start-trees/start-tree-$(printf %02d $k).nwk
export CET_SEED=$(( RUNG_SEED + rep ))
echo "$(date) src $(git -C /nobackup/pjjg18/evotempo-penguin-cet/src rev-parse HEAD) cell $CET_CELL rep $rep start $CET_START_TREE seed $CET_SEED host $(hostname)"
source /nobackup/pjjg18/evotempo-penguin-cet/src/data-raw/empirical/penguins/cet/run-array.sh
