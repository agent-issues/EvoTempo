#!/bin/bash
# Summaries of the penguin FBD chains on Hamilton.
#
#   bash summarise-hamilton.sh <run root> [burn-in fraction]
#
# <run root> holds one directory per chain, <clock>_run<rep>/output/.  Writes
# ./summary/: the CSVs of summarise.R, and the rank-normalised R-hat / ESS
# table and quartet consensus tree of rhat-ess-quartet.R.  No MCC tree is
# made.  The R library must provide ape, TreeTools, ConsTree and MkPrime
# (R_LIBS below).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH --time=0-08:00:00

set +e; set +o pipefail
# batch jobs need the module system initialised before module purge/load
[ -f /etc/profile ] && source /etc/profile
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ROOT=$1
BURNIN=${2:-0.25}
export R_LIBS=${R_LIBS:-/nobackup/pjjg18/evotempo-penguin-cet/lib-tools:/nobackup/pjjg18/evotempo-dev/lib}
module purge
module load r/4.5.1
mkdir -p summary flat
for clock in shared separate; do
  # flat copies of every chain's log and tree file, without the last line,
  # which may be incomplete if the chain is still running or was stopped
  for d in "$ROOT"/${clock}_run*/output; do
    for f in "$d"/penguins_${clock}_run*.log "$d"/penguins_${clock}_run*.trees; do
      case "$f" in *_rates.trees) continue ;; esac
      [ -f "$f" ] && head -n -1 "$f" > flat/$(basename "$f")
    done
  done
  Rscript summarise.R flat ${clock} ${BURNIN}
  Rscript rhat-ess-quartet.R "$ROOT" ${clock} ${BURNIN} summary 4000 penguins_ > summary/rhat_ess_${clock}.txt 2>&1
done
rm -rf flat
ls -lh summary
