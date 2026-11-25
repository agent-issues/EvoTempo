#!/bin/bash
# One chain of the morphology-only RevBayes comparator (penguins_morph_fbd.Rev).
#
#   sbatch -J evotempo-penguin-rbmorph-ucln-1 --export=ALL,CLOCK=ucln,REP=1 run-hamilton.sh
#
# CLOCK = ucln | strict; REP seeds the chain.  Optional: N_GEN, SAMPLE_EVERY,
# BURNIN_GEN, CKP_EVERY, RESUME=true (continue from the checkpoint),
# ML=true (stepping-stone instead of MCMC; CATS, ML_GEN), START_TREE,
# START_ORIGIN.  One stone per Slurm task (README, "Common-footing marginal
# likelihoods"):
#   sbatch -J evotempo-penguin-rbml-strict-1 --array=1-65 \
#     --export=ALL,CLOCK=strict,REP=1,ML=true,ML_MODE=stone,CATS=65 run-hamilton.sh
# runs stone $SLURM_ARRAY_TASK_ID (ML_STONE overrides it); ML_PREBURN
# (pre-burn-in generations, default ML_GEN), ML_ALPHA (default 0.2),
# ML_MONITOR (true: parameter and tree samples of the stone, every
# ML_SAMPLE_EVERY generations; default true for stones 1 and CATS, the powers
# 1 and 0); then ML_MODE=summarize merges the stones and writes
# output/penguins_morph_<CLOCK>_run<REP>_ml.txt.  Run from a copy of this directory (with data/) under
# /nobackup.  RB defaults to the project's pinned RevBayes development build
# a2a8e2db (gcc/14.2 + boost/1.88.0; `rb --version` prints 1.4.2-preview), not
# v1.4.1: in v1.4.1 mvCollapseExpandFossilBranch is biased when fossil ages
# are uncertain (fixed upstream in c389174b), and every penguin fossil has an
# age range (#66).  Never the cluster's revbayes/1.2.5 module.  The binary
# path and its `--version` output go to out/rb_version_<CLOCK>_<REP>*.txt.
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=3-00:00:00
set +e; set +o pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
module purge
module load gcc/14.2 boost/1.88.0
RB=${RB:-/nobackup/pjjg18/revbayes-dev-a2a8e2db/projects/cmake/build/rb}
N_GEN=${N_GEN:-200000}
SAMPLE_EVERY=${SAMPLE_EVERY:-100}
BURNIN_GEN=${BURNIN_GEN:-2000}
CKP_EVERY=${CKP_EVERY:-5000}
RESUME=${RESUME:-false}
ML=${ML:-false}
CATS=${CATS:-32}
ML_GEN=${ML_GEN:-10000}
ML_MODE=${ML_MODE:-all}
ML_STONE=${ML_STONE:-${SLURM_ARRAY_TASK_ID:-0}}
ML_PREBURN=${ML_PREBURN:-$ML_GEN}
ML_ALPHA=${ML_ALPHA:-0.2}
ML_SAMPLE_EVERY=${ML_SAMPLE_EVERY:-100}
if [ -z "$ML_MONITOR" ]; then
  ML_MONITOR=false
  [ "$ML_STONE" = 1 -o "$ML_STONE" = "$CATS" ] && ML_MONITOR=true
fi
tag=${CLOCK}_${REP}
[ "$ML_MODE" = stone ] && tag=${tag}_stone${ML_STONE}
[ "$ML_MODE" = summarize ] && tag=${tag}_summarize
START_TREE=${START_TREE:-data/start_tree.nwk}
START_ORIGIN=${START_ORIGIN:-$(cat data/start_origin.txt)}
mkdir -p output out
driver=out/driver_${tag}.Rev
cat > "$driver" <<EOF
clock_model = "${CLOCK}"
rep = ${REP}
n_gen = ${N_GEN}
sample_every = ${SAMPLE_EVERY}
burnin_gen = ${BURNIN_GEN}
checkpoint_every = ${CKP_EVERY}
resume = ${RESUME}
ml = ${ML}
cats = ${CATS}
ml_gen = ${ML_GEN}
ml_mode = "${ML_MODE}"
ml_stone = ${ML_STONE}
ml_preburn = ${ML_PREBURN}
ml_alpha = ${ML_ALPHA}
ml_monitor = ${ML_MONITOR}
ml_sample_every = ${ML_SAMPLE_EVERY}
start_tree_file = "${START_TREE}"
start_origin = ${START_ORIGIN}
source("penguins_morph_fbd.Rev")
q()
EOF
suffix=""; [ "$RESUME" = true ] && suffix="_resume$(date +%s)"
{ echo "RB=$RB"; "$RB" --version 2>&1; } > out/rb_version_${tag}${suffix}.txt
/usr/bin/time -v "$RB" "$driver" > out/screen_${tag}${suffix}.txt 2> out/time_${tag}${suffix}.txt
echo "rb exit status $?"
