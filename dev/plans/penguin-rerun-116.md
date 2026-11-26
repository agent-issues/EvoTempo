# Penguin CET full rerun under the mean-rate tempo location (#116)

**Not submitted.** This is the launch plan for after the tempo-location PR merges. The maintainer decided on #116 to
switch penguins and rerun everything: production, γ-rungs and Bayes factors.

## Preconditions

1. **What must be merged.** This PR, and the penguin pipeline PRs whose scripts the rerun uses:
   - #103, the Bayes-factor safeguards;
   - #105, eval gated on its checks;
   - #107, `rung-array.sh`, `stage1-results.R`, `data/start-trees/rungs/` and the launch record.

   If #110's new topology moves are about to land, wait for them as well. Penguins need a rerun once they exist
   anyway, and one rerun costs about 3,300 core-hours (below).
2. **Use a fresh tree, `W2`.** Never resume in `/nobackup/pjjg18/evotempo-penguin-cet/runs`. The old checkpoints hold
   no `location`: `run-chain.R` now refuses to resume them, and `bayes-factors.R eval` refuses to score them, because
   they were run under `"bm"` (review of #123, M3). Nothing checks the library itself.

## Commands

```bash
W2=/nobackup/pjjg18/evotempo-penguin-cet-v2
SHA=<short sha of main after the merges>
mkdir -p $W2/logs $W2/out $W2/runs
git clone https://github.com/agent-issues/EvoTempo.git $W2/src   # or git archive + scp, as for evotempo-meanrate
cd $W2/src && git checkout $SHA

# Build (Slurm test partition, about 10 min)
module load r/4.5.1 gcc/14.2
mkdir -p $W2/lib-$SHA
R_LIBS=$W2/lib-$SHA:/nobackup/pjjg18/evotempo-dev/lib R CMD INSTALL --library=$W2/lib-$SHA $W2/src
echo $SHA > $W2/out/build-sha.txt
cp $W2/out/build-sha.txt $W2/lib-$SHA/SHA

export W=$W2 SRC=$W2/src LIB=$W2/lib-$SHA

# Pre-flight checks
R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib Rscript -e \
  'source("data-raw/empirical/penguins/cet/common.R"); m <- CetModel(1, 1); stopifnot(m$location == "mean", m$scheme == "sqra")'
R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib Rscript data-raw/empirical/penguins/cet/check-mkv.R
sbatch -J evotempo-penguin-v2-smoke -o $W2/logs/%x-%A_%a.out --export=ALL,W,SRC,LIB \
  --array=1-5 data-raw/empirical/penguins/cet/smoke.sh
R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib Rscript data-raw/empirical/penguins/cet/check-smoke.R $W2

# Production: 4 cells x 4 replicates, dispersed starts (data/start-trees, draws 1-16), 600 000 iterations
sbatch -J evotempo-penguin-v2 -o $W2/logs/%x-%A_%a.out --array=1-16 \
  --export=ALL,W,SRC,LIB,CET_NITER=600000,CET_THIN=100,CET_SEG=10000 \
  data-raw/empirical/penguins/cet/run-array.sh

# gamma rungs, 4 replicates each, 600 000 iterations; the same seeds and starts as stage 1
for spec in b0g0.25:510000:production b0g0.5:520000:production b0g0.75:530000:production \
            b1g0.25:610000:production b1g0.5:620000:production b1g0.75:630000:production \
            b0g0.125:540000:dispersed b0g0.375:550000:dispersed \
            b1g0.125:640000:dispersed b1g0.375:650000:dispersed; do
  IFS=: read cell seed start <<< "$spec"
  sbatch -J evotempo-penguin-v2-rung-$cell -o $W2/logs/%x-%A_%a.out --array=1-4 \
    --export=ALL,W,SRC,LIB,CET_CELL=$cell,RUNG_SEED=$seed,RUNG_START=$start,CET_NITER=600000 \
    data-raw/empirical/penguins/cet/rung-array.sh
done
```

**Wall time.** `CET_HOURS` (default 70) stops a chain cleanly. A chain of 600 000 iterations at about 0.33 s per
iteration needs about 55 h, so it fits in one job. If a chain stops early, resubmit the same command: it continues from
its checkpoint.

**Launch record.** Write `$W2/runs/launch.csv` with the 14 groups (production plus 10 rungs), in the format of
stage 1's `runs/launch.csv`. Record the array IDs in the penguin README.

## After the chains finish

```bash
for d in $W2/runs/*_run*; do
  sbatch -J evotempo-penguin-v2-eval -o $W2/logs/%x-%j.out -p shared -c 1 --mem=4G -t 06:00:00 \
    --wrap "module load r/4.5.1 gcc/14.2; cd $SRC; R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib Rscript data-raw/empirical/penguins/cet/bayes-factors.R eval $d"
done
# then, once every eval has finished:
R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib Rscript data-raw/empirical/penguins/cet/bayes-factors.R combine $W2/runs
R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib Rscript data-raw/empirical/penguins/cet/stage1-results.R $W2/runs ...
```

Apply the #107 acceptance gate before reporting a γ Bayes factor:
- every rung has R-hat ≤ 1.01 and bulk ESS ≥ 400, meanX included;
- every bridge step has Kish ESS ≥ 100 in both directions;
- the chain-half SE is < 0.5;
- the squares close within 3 SE.

## Cost (from stage 1's speeds)

| Item | Core-hours |
|---|---|
| Production: 16 × 600 000 × 0.33 s | ≈ 880 |
| Rungs: 40 × 600 000 × 0.35 s | ≈ 2,330 |
| Eval, combine, summaries | < 20 |
| **Total** | **≈ 3,250** |

## What changes for interpretation

Rate parameters are now lineage-mean rates. Under the stage-1 model they were rates at log-tempo 0. In the γ = 1 cells,
stage-1 c and θ were higher by about e^{s²/2} ≈ 1.3–1.4 at s ≈ 0.75–0.83. The priors (comparator-matched) are unchanged
in value and now describe the same quantity as the comparators' (`ucld_mean`, FBD rates). Stage-1 results are not
comparable parameter by parameter.
