#!/bin/bash
# Cut the round-4 separate-clock chains at a common iteration after their
# checkpoint extension overshot (RevBayes appends N_GEN iterations on resume;
# see run-hamilton.sh).  This is the command that was used, on Hamilton, after
# the extension jobs had been cancelled:
#
#   bash truncate-chains.sh <run root> 220000
#
# Keeps the header line and every sample with iteration <= the limit in each
# chain's .log, .trees and _rates.trees, and prints the last kept iteration.
ROOT=$1
LIMIT=${2:-220000}
for d in "$ROOT"/separate_run*/output; do
  k=$(basename "$(dirname "$d")" | sed 's/separate_run//')
  for f in "$d"/penguins_separate_run$k.log "$d"/penguins_separate_run$k.trees "$d"/penguins_separate_run${k}_rates.trees; do
    awk -F'\t' -v lim="$LIMIT" '$1 !~ /^[0-9]+$/ || $1 + 0 <= lim' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
    echo "$f $(tail -n 1 "$f" | cut -f1)"
  done
done
