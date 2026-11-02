#!/bin/bash
# One chain of the morphology-only RevBayes comparator (penguins_morph_fbd.Rev).
#
#   sbatch -J evotempo-penguin-rbmorph-ucln-1 --export=ALL,CLOCK=ucln,REP=1 run-hamilton.sh
#
# CLOCK = ucln | strict; REP seeds the chain.  Optional: N_GEN, SAMPLE_EVERY,
# BURNIN_GEN, CKP_EVERY, RESUME=true (continue from the checkpoint),
# ML=true (stepping-stone instead of MCMC; CATS, ML_GEN), START_TREE,
# START_ORIGIN.  Run from a copy of this directory (with data/) under
# /nobackup; RB is our RevBayes v1.4.1 build (never the revbayes/1.2.5 module).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=3-00:00:00
set +e; set +o pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
module purge
module load gcc/14.2 boost/1.88.0
RB=${RB:-/nobackup/pjjg18/revbayes-1.4.1/projects/cmake/build/rb}
N_GEN=${N_GEN:-200000}
SAMPLE_EVERY=${SAMPLE_EVERY:-100}
BURNIN_GEN=${BURNIN_GEN:-2000}
CKP_EVERY=${CKP_EVERY:-5000}
RESUME=${RESUME:-false}
ML=${ML:-false}
CATS=${CATS:-32}
ML_GEN=${ML_GEN:-10000}
START_TREE=${START_TREE:-data/start_tree.nwk}
START_ORIGIN=${START_ORIGIN:-$(cat data/start_origin.txt)}
mkdir -p output out
driver=out/driver_${CLOCK}_${REP}.Rev
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
start_tree_file = "${START_TREE}"
start_origin = ${START_ORIGIN}
source("penguins_morph_fbd.Rev")
q()
EOF
suffix=""; [ "$RESUME" = true ] && suffix="_resume$(date +%s)"
/usr/bin/time -v "$RB" "$driver" > out/screen_${CLOCK}_${REP}${suffix}.txt 2> out/time_${CLOCK}_${REP}${suffix}.txt
echo "rb exit status $?"
