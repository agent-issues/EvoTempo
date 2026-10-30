#!/bin/bash
# One MCMC chain of the hymenoptera RevBayes baseline on Hamilton (Durham HPC).
#
#   sbatch -J evotempo-hyme-rb141-<config>-<rep> --export=ALL,CONFIG=<config>,REP=<rep> run-hamilton.sh
#
# CONFIG = te-shared | te-separate | morph; REP = replicate (seeds the chain).  Optional:
#   N_GEN, SAMPLE_EVERY, BURNIN_GEN  chain length, thinning, tuning iterations
#   START_FROM   output/ directory of an earlier chain: start from its last
#                tree and scalar parameter values
#   RESUME=true  continue this chain from output/hymenoptera_<CONFIG>_run<REP>.ckp,
#                appending to its output files
#   CKP_EVERY    iterations between checkpoints (default 500)
# Run from a copy of this directory under /nobackup.  RB is RevBayes v1.4.1
# built from the release tag (gcc/14.2 + boost/1.88.0); the cluster's
# revbayes/1.2.5 module is not used.
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

set +e; set +o pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
module purge
module load gcc/14.2 boost/1.88.0
RB=${RB:-/nobackup/pjjg18/revbayes-1.4.1/projects/cmake/build/rb}

N_GEN=${N_GEN:-100000}
SAMPLE_EVERY=${SAMPLE_EVERY:-50}
BURNIN_GEN=${BURNIN_GEN:-2000}
RESUME=${RESUME:-false}
CKP_EVERY=${CKP_EVERY:-500}

mkdir -p output out
extra=""
if [ -n "$START_FROM" ] && [ "$RESUME" != true ]; then
  prev_trees=$(ls "$START_FROM"/hymenoptera_${CONFIG}_run*.trees | grep -v _rates | head -1)
  prev_log=${prev_trees%.trees}.log
  last_tree=$(grep ';' "$prev_trees" | tail -n 1)
  start_iter=$(printf '%s\n' "$last_tree" | cut -f1)
  printf '%s\n' "$last_tree" | awk -F'\t' '{print $NF}' | sed 's/\[[^]]*\]//g' > out/start_tree_${CONFIG}_${REP}.nwk
  awk -F'\t' -v it="$start_iter" 'NR == 1 { for (i = 1; i <= NF; i++) name[i] = $i; next }
    $1 == it { for (i = 1; i <= NF; i++)
                 if (name[i] ~ /^(diversification|turnover|sampling_prop|clock_rate|clock_morph|clock_dna|ucld_sd|ucld_sd_morph|ucld_sd_dna|alpha_morph|alpha_dna\[[0-9]+\])$/)
                   printf "%s.setValue(%s)\n", name[i], $i }' "$prev_log" > out/start_values_${CONFIG}_${REP}.Rev
  echo "starting from $prev_trees at iteration $start_iter"
  extra="start_tree_file = \"out/start_tree_${CONFIG}_${REP}.nwk\"
start_values_file = \"out/start_values_${CONFIG}_${REP}.Rev\""
fi

driver=out/driver_${CONFIG}_${REP}.Rev
cat > "$driver" <<EOF
config = "${CONFIG}"
rep = ${REP}
n_gen = ${N_GEN}
sample_every = ${SAMPLE_EVERY}
burnin_gen = ${BURNIN_GEN}
resume = ${RESUME}
checkpoint_every = ${CKP_EVERY}
${extra}
source("hymenoptera_fbd.Rev")
q()
EOF
suffix=""; [ "$RESUME" = true ] && suffix="_resume$(date +%s)"
/usr/bin/time -v "$RB" "$driver" > out/screen_${CONFIG}_${REP}${suffix}.txt 2> out/time_${CONFIG}_${REP}${suffix}.txt
echo "rb exit status $?"
