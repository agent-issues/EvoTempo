#!/bin/bash
# #110 pilot of Metropolis-coupled chains: penguin cell (0,0), 4 replicates
# per arm as array tasks 1-4, at equal numbers of sampling level-iterations
# (iterations x chains: 40 000).  The arms are not at equal cost: burn-in runs
# at every level, so n1 does 42 000 level-iterations in all and the coupled
# arms 48 000 (+14%), plus two likelihood evaluations per accepted swap.
# Compare the arms by the per-hour columns of mc3-pilot-summary.R.
#
#   sbatch -J mc3-pilot-<arm> --array=1-4 --export=ALL,ARM=<arm> mc3-pilot.sh
#
#   ARM=n1     one chain,    40 000 iterations, thin 40
#   ARM=n4     four chains,  10 000 iterations, thin 10, deltaHeat 0.1
#   ARM=n4d3   four chains,  10 000 iterations, thin 10, deltaHeat 0.3
# Each arm has 2000 burn-in iterations and stores 1000 samples.  Replicate r
# starts from the dispersed start tree of production chain 4 (r - 1) + 1.
# W, SRC and LIB name the work directory, the checkout and the library.
# Summary: mc3-pilot-summary.R.
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH --time=1-00:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-mc3/logs/%x-%A_%a.out
set +e; set +o pipefail; source /etc/profile || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
module purge; module load r/4.5.1 gcc/14.2
W=${W:-/nobackup/pjjg18/evotempo-mc3}
SRC=${SRC:-$W/src}
export R_LIBS=${LIB:-$W/lib}:/nobackup/pjjg18/evotempo-dev/lib
rep=$SLURM_ARRAY_TASK_ID
case "$ARM" in
  n1)   export CET_NCHAIN=1 CET_NITER=40000 CET_THIN=40 CET_SEG=4000; arm=1 ;;
  n4)   export CET_NCHAIN=4 CET_DELTA_HEAT=0.1 CET_NITER=10000 CET_THIN=10 CET_SEG=1000; arm=2 ;;
  n4d3) export CET_NCHAIN=4 CET_DELTA_HEAT=0.3 CET_NITER=10000 CET_THIN=10 CET_SEG=1000; arm=3 ;;
  *) echo "unknown ARM '$ARM'"; exit 1 ;;
esac
export CET_BURNIN=2000 CET_HOURS=22 CET_STATES=false CET_OUT=$W/pilot/$ARM
export CET_SEED=$(( 120000 + 1000 * arm + rep ))
mkdir -p $CET_OUT $W/out
cd $SRC
echo "$(date) src $(git rev-parse HEAD) arm $ARM rep $rep seed $CET_SEED host $(hostname)"
Rscript data-raw/empirical/penguins/cet/run-chain.R 1 $rep >> $W/out/pilot-$ARM-$rep.txt 2>&1
status=$?
echo "$(date) exit $status"
exit $status
