#!/bin/bash
# One RevBayes (fd807d10) chain of the penguin FBD model for the PR #99
# young-crown test.  Env: CHAIN_DIR CLOCK REP N_GEN SAMPLE_EVERY BURNIN_GEN and
# either START_STEM (young-crown start: <stem>.nwk + <stem>.Rev) or START_FROM
# (output/ dir of a round-4 chain: last tree + scalars, as run-hamilton.sh).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00
set +e; set +o pipefail; source /etc/profile || true; module purge || true
module load gcc/14.2 boost/1.88.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
RB=/nobackup/pjjg18/revbayes-dev-fd807d10/projects/cmake/build/rb
cd "$CHAIN_DIR" || exit 1
mkdir -p output out
if [ -n "$START_STEM" ]; then
  cp "$START_STEM.nwk" out/start_tree.nwk; cp "$START_STEM.Rev" out/start_values.Rev
else
  prev_trees=$(ls "$START_FROM"/penguins_${CLOCK}_run*.trees | grep -v _rates | head -1)
  prev_log=${prev_trees%.trees}.log
  last_tree=$(grep ';' "$prev_trees" | tail -n 1)
  start_iter=$(printf '%s\n' "$last_tree" | cut -f1)
  printf '%s\n' "$last_tree" | awk -F'\t' '{print $NF}' | sed 's/\[[^]]*\]//g' > out/start_tree.nwk
  awk -F'\t' -v it="$start_iter" 'NR == 1 { for (i = 1; i <= NF; i++) name[i] = $i; next }
    $1 == it { for (i = 1; i <= NF; i++)
                 if (name[i] ~ /^(diversification|turnover|sampling_prop|origin_time|clock_morph|clock_dna|ucld_sd|ucld_sd_morph|ucld_sd_dna|alpha_morph|alpha_dna\[[0-9]+\]|r_(ac|at|cg|ct|gt)\[[0-9]+\])$/)
                   printf "%s.setValue(%s)\n", name[i], $i }' "$prev_log" > out/start_values.Rev
  echo "control start: $prev_trees iteration $start_iter"
fi
cat > out/driver.Rev <<EOF
clock_model = "${CLOCK}"
rep = ${REP}
n_gen = ${N_GEN}
sample_every = ${SAMPLE_EVERY}
burnin_gen = ${BURNIN_GEN}
ucln_param = "noncentred"
resume = false
checkpoint_every = 1000
use_up_down = false
under_prior = false
start_tree_file = "out/start_tree.nwk"
start_values_file = "out/start_values.Rev"
source("penguins_fbd.Rev")
q()
EOF
/usr/bin/time -v "$RB" out/driver.Rev > out/screen.txt 2> out/time.txt
echo "rb exit status $?"
