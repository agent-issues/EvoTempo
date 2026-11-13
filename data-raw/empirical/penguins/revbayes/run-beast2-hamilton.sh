#!/bin/bash
# BEAST2 reference run of the source analysis (Gavryushkina et al. 2017,
# xmls/8+dna_R.xml in ../raw/xmls.zip) on Hamilton.  Run from this directory:
#
#   sbatch -J evotempo-penguin-beast2-post-1 --export=ALL,SEED=201 run-beast2-hamilton.sh <run dir>
#   PRIOR=true  samples from the prior (fossil ages and tip counts only)
#
# BEAST 2.4.8 with its bundled JRE, SA 1.1.7 and MM 1.0.5, installed with
#   beast/bin/addonmanager -add SA; beast/bin/addonmanager -add MM
# One chain of 1e8 states took 22-28 h on one core.
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
set +e; set +o pipefail
BEAST=${BEAST:-/nobackup/pjjg18/beast-2.4.8/beast}
export JAVA_HOME=$(ls -d "$BEAST"/jre*)
XMLZIP=${XMLZIP:-$(pwd)/../raw/xmls.zip}
RUN=$1; mkdir -p "$RUN"; cd "$RUN"
unzip -o -j -q "$XMLZIP" 'xmls/8+dna_R.xml'
xml=8+dna_R.xml
if [ "${PRIOR:-false}" = true ]; then
  sed -e 's/<run chainLength="100000000" id="mcmc" spec="MCMC">/<run chainLength="50000000" id="mcmc" spec="MCMC" sampleFromPrior="true">/' \
      -e 's/logEvery="10000"/logEvery="5000"/g' 8+dna_R.xml > penguins_prior.xml
  xml=penguins_prior.xml
fi
"$BEAST"/bin/beast -java -threads 1 -seed ${SEED:-201} -overwrite $xml > screen.txt 2>&1
