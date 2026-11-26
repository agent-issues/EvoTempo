#!/bin/bash
# Draws 17-32 of draw-start-trees.R for the new rungs (Slurm 19472362).
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%j.out
set +e; source /etc/profile || true; module purge || true; module load r/4.5.1 gcc/14.2
export R_LIBS=/nobackup/pjjg18/evotempo-penguin-cet/lib-48e6bf7:/nobackup/pjjg18/evotempo-penguin-cet/lib-tools:/nobackup/pjjg18/evotempo-dev/lib
T=/nobackup/pjjg18/evotempo-penguin-cet/tmp-draw32; rm -rf $T; mkdir -p $T/out
cd /nobackup/pjjg18/evotempo-penguin-cet/src
Rscript data-raw/empirical/penguins/cet/draw-start-trees.R /nobackup/pjjg18/evotempo-penguins-rbdev-r4 $T/out 32 20261111 0.25 shared
Rscript -e 'a <- read.csv("data-raw/empirical/penguins/cet/data/start-trees/draws.csv"); b <- read.csv("'$T'/out/draws.csv"); cat("first 16 identical:", identical(a[, c("chain","iteration")], b[1:16, c("chain","iteration")]), "; 32 distinct:", !anyDuplicated(paste(b$chain, b$iteration)), "\n")'
for k in $(seq 17 32); do kk=$(printf %02d $k)
  cp -r data-raw/empirical/penguins/cet $T/cetcopy-$kk 2>/dev/null || true
done
rm -rf $T/cetcopy-*
for k in $(seq 17 32); do kk=$(printf %02d $k)
  Rscript data-raw/empirical/penguins/cet/prepare.R $T/out/start-tree-$kk.nex $(cat $T/out/start-tree-$kk-origin.txt) $T/out/start-tree-$kk.nwk > $T/out/prepare-$kk.txt 2>&1; echo "k $k prepare exit $?"
done
git status --short | head
ls $T/out | wc -l
