#!/bin/bash
# Launch the RevBayes baseline chains on Hamilton: one directory per chain
# under <run root>, each a copy of this directory (scripts + data/).
#
#   bash launch-hamilton.sh <run root> <config> <first rep> <last rep> [partition] [n resumes]
#
# On partition "shared" (3-day limit) each chain is followed by [n resumes]
# jobs that continue it from its checkpoint (--dependency=afterany).  Job
# names: evotempo-<dataset>-rb141-<config>-<rep>.  Extra sbatch --export
# variables (N_GEN, SAMPLE_EVERY, BURNIN_GEN, CKP_EVERY) come from the caller's
# environment.
set +e; set +o pipefail
ROOT=$1; CONFIG=$2; FIRST=$3; LAST=$4; PART=${5:-shared}; NRES=${6:-0}
HERE=$(cd "$(dirname "$0")" && pwd)
DATASET=$(basename "$(dirname "$HERE")")
SHORT=${DATASET:0:4}
mkdir -p "$ROOT"
for rep in $(seq "$FIRST" "$LAST"); do
  d="$ROOT/${CONFIG}_run${rep}"
  if [ -e "$d" ]; then echo "exists: $d"; continue; fi
  mkdir -p "$d/out" "$d/output"
  cp -r "$HERE"/*.Rev "$HERE"/*.sh "$HERE"/data "$d"/
  cd "$d"
  name="evotempo-${SHORT}-rb141-${CONFIG}-${rep}"
  time=3-00:00:00; [ "$PART" = long ] && time=7-00:00:00
  jid=$(sbatch --parsable -p "$PART" --time="$time" -J "$name" -o out/slurm-%j.txt \
        --export=ALL,CONFIG="$CONFIG",REP="$rep" run-hamilton.sh)
  echo "$name $jid"
  prev=$jid
  for k in $(seq 1 "$NRES"); do
    jid=$(sbatch --parsable -p "$PART" --time="$time" -J "${name}-resume${k}" -o out/slurm-%j.txt \
          --dependency=afterany:"$prev" --export=ALL,CONFIG="$CONFIG",REP="$rep",RESUME=true run-hamilton.sh)
    echo "${name}-resume${k} $jid"
    prev=$jid
  done
  cd - >/dev/null
done
