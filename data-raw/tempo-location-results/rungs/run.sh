#!/bin/bash
#SBATCH -p shared
#SBATCH -o /nobackup/pjjg18/evotempo-meanrate/logs/%x-%A_%a.out
# Runs $CMD in a fresh copy of src ($W/run-$RUN; RUNDIR=keep reuses it), with
# the lib in $W/libsha.  CLEAN: cached result files to delete first.
set +e
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W=/nobackup/pjjg18/evotempo-meanrate
LIB=$W/lib-$(cat $W/libsha)
export R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib
export EVOTEMPO_SHA=$(cat $LIB/SHA)
D=$W/run-$RUN
if [ "$RUNDIR" != keep ] || [ ! -d $D ]; then rm -rf $D; cp -r $W/src $D; fi
cd $D
[ -n "$CLEAN" ] && rm -f $CLEAN
echo "RUN=$RUN SHA(lib)=$EVOTEMPO_SHA src=$(cat SHA) task=$SLURM_ARRAY_TASK_ID start $(date)"
eval "$CMD"
echo "exit $? $(date)"
