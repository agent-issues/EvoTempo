#!/bin/bash
# stage1-results.R and stage1-figures.R on the 400k chains (Slurm 19449555).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH --time=0-04:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%j.out
set +e; set +o pipefail; source /etc/profile || true; module purge || true
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=/nobackup/pjjg18/evotempo-penguin-cet
export R_LIBS=$W/lib-48e6bf7:$W/lib-tools:/nobackup/pjjg18/evotempo-dev/lib
cd $W/src
O=$W/stage1-results-400k
for d in $W/runs/b[01]g[01]_run[1-4]; do test -f $d/seg0040_samples.csv || { echo "$d not at 400k"; exit 2; }; done
/usr/bin/time -v Rscript $W/jobs/stage1-results-v2.R $W/runs $W/rbmorph-fd807d10 $O 0.25 > $O.log 2>&1
echo "results exit $?" >> $O.log
Rscript $W/jobs/stage1-figures.R $O $O/figures >> $O.log 2>&1
echo "figures exit $?" >> $O.log
