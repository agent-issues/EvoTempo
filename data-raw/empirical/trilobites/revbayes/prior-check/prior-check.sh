#!/bin/bash
# Build check for the trilobite baseline, as it was run on Hamilton on 2 November 2026 (job
# 19329128; paths as on Hamilton).  Before this script ran, prior-check-tri/{v141,dev}/ each
# held copies of trilobites_fbd.Rev (identical to trilobites_fbd.Rev at commit a3ad184, md5
# e3ba8720ba13eb5da630a665c35cc638), mcc.Rev and data/.  The drivers set under_prior = true
# (mymodel.ignoreAllData()).  The same runs can now be made with run-hamilton.sh:
#   sbatch --export=ALL,CONFIG=morph,REP=<1|2>,UNDER_PRIOR=true,N_GEN=200000,SAMPLE_EVERY=100,BURNIN_GEN=1000,CKP_EVERY=10000,RB=<rb> run-hamilton.sh
# Summary: Rscript prior-summary.R prior-check-tri trilobites num_samp_anc,sampling_prop,root_age,ingroup_age
#SBATCH -J evotempo-tril-prior-check
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-baselines-rb141/prior-check-tri/slurm-%j.txt
set +e; set +o pipefail
export OMP_NUM_THREADS=1
P=/nobackup/pjjg18/evotempo-baselines-rb141/prior-check-tri
module purge; module load gcc/14.2 boost/1.88.0
declare -A RBS=( [v141]=/nobackup/pjjg18/revbayes-1.4.1/projects/cmake/build/rb [dev]=/nobackup/pjjg18/revbayes-dev-a2a8e2db/projects/cmake/build/rb )
for b in v141 dev; do
  ${RBS[$b]} --version > $P/version_$b.txt 2>&1
  for s in 1 2; do
    d=$P/$b; mkdir -p $d/out $d/output
    printf 'config = "morph"\nrep = %d\nn_gen = %d\nsample_every = 100\nburnin_gen = 1000\nunder_prior = true\ncheckpoint_every = 10000\nsource("trilobites_fbd.Rev")\nq()\n' $s ${N_GEN:-200000} > $d/out/driver_$s.Rev
    (cd $d && ${RBS[$b]} out/driver_$s.Rev > out/screen_$s.txt 2>&1) &
  done
done
wait

