#!/bin/bash
# Link root for reproducing the 400 000-iteration results of #107 from runs/
# (#107 delta review, N1): the 16 production chains and the 24 gamma rungs
# launched on 2026-11-13 (not the gamma = 0.125 / 0.375 rungs), with the
# 400 000-iteration launch record.  Use with CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20.
#
#   bash data-raw/empirical/penguins/cet/jobs/make-bf400k-root.sh <runs> <new root>
#
# Run from the package root (it copies stage1-results/launch.csv).  Note that
# eval writes cross.csv into the linked chain directories, replacing any
# cross.csv evaluated under another cap.
set -eu
runs=$1; root=$2
C=data-raw/empirical/penguins/cet
[ -e "$root" ] && { echo "$root exists; remove it first"; exit 2; }
mkdir -p "$root"
cp $C/stage1-results/launch.csv "$root/launch.csv"
n=0
for cell in b0g0 b0g1 b1g0 b1g1 b0g0.25 b0g0.5 b0g0.75 b1g0.25 b1g0.5 b1g0.75; do
  case $cell in b[01]g[01]) m=40;; *) m=20;; esac
  for rep in 1 2 3 4; do
    d="$runs/${cell}_run$rep"
    for k in $(seq 1 $m); do
      [ -f "$d/$(printf 'seg%04d_samples.csv' $k)" ] || { echo "$d lacks segment $k of $m"; exit 2; }
    done
    ln -s "$(cd "$d" && pwd)" "$root/${cell}_run$rep"
    n=$((n + 1))
  done
done
echo "$n chains linked in $root, launch record $(($(wc -l < "$root/launch.csv") - 1)) groups"
