# How the reported chains were run on Hamilton

This file records the commands that produced the chains reported in `README.md`.
`launch-hamilton.sh` and `run-hamilton.sh` cover these steps, but the reported chains also used
an exported `RB`, a per-chain `START_FROM` and manual cancellation. An exact rerun is not
possible:
- the `START_FROM` directories (the v1.4.1 chains) have been deleted;
- each resume reseeds the random number generator with the same `seed(rep * 7919 + 211)`,
  because RevBayes checkpoints do not store its state;
- the segment boundaries were set by Slurm time limits.

`W=/nobackup/pjjg18/evotempo-baselines-rb141`; `H=$W/src-hym/data-raw/empirical/hymenoptera/revbayes`
(a clone of this branch, with `data/` written by `Rscript prepare-data.R; Rscript start-tree.R`).

## 1. Superseded v1.4.1 chains (31 October – 3 November 2026)

At that time `run-hamilton.sh` defaulted to `RB=/nobackup/pjjg18/revbayes-1.4.1/projects/cmake/build/rb`
and `launch-hamilton.sh` named jobs `evotempo-hyme-rb141-<config>-<rep>`.

```bash
export N_GEN=2000000 BURNIN_GEN=2000 CKP_EVERY=500
SAMPLE_EVERY=50  bash $H/launch-hamilton.sh $W/runs/hymenoptera te-shared   1 6 shared 2
SAMPLE_EVERY=50  bash $H/launch-hamilton.sh $W/runs/hymenoptera te-separate 1 6 shared 2
SAMPLE_EVERY=100 bash $H/launch-hamilton.sh $W/runs/hymenoptera morph       1 6 shared 2
# 31 Oct ~03:00: two more chains per te-* configuration
export N_GEN=2000000 BURNIN_GEN=2000 CKP_EVERY=500 SAMPLE_EVERY=50
bash $H/launch-hamilton.sh $W/runs/hymenoptera te-shared   7 8 shared 2
bash $H/launch-hamilton.sh $W/runs/hymenoptera te-separate 7 8 shared 2
```

`N_GEN = 2e6` is a ceiling, never reached; chains ran until stopped. The 22 first segments ended
at their 3-day limit (3 November, 00:40 and 02:54). The `morph` resume jobs had been cancelled
on 2 November, and the running `te-*` resume jobs were cancelled on 3 November at 08:15, once the
development-build chains below were running. These chains are not reported, and their output
(`$W/runs`) has been deleted.

## 2. Reported chains, development build `a2a8e2db` (2–10 November 2026)

The build was selected by exporting `RB`. The run-directory copies of `run-hamilton.sh` are the
older ones, which still default to v1.4.1; the committed version now defaults to `a2a8e2db`. All
54 chain segments print the development build's start-up warning `Unknown user setting with key
'moduledir'`, which v1.4.1 does not print.

```bash
export RB=/nobackup/pjjg18/revbayes-dev-a2a8e2db/projects/cmake/build/rb BUILD_TAG=dev \
       N_GEN=2000000 BURNIN_GEN=2000 CKP_EVERY=500
for cfg in te-shared te-separate; do for k in 1 2 3 4 5 6 7 8; do
  START_FROM=$W/runs/hymenoptera/${cfg}_run$k/output SAMPLE_EVERY=50 \
    bash $H/launch-hamilton.sh $W/runs-dev/hymenoptera $cfg $k $k shared 2
done; done
for k in 1 2 3 4 5 6; do
  START_FROM=$W/runs/hymenoptera/morph_run$k/output SAMPLE_EVERY=100 \
    bash $H/launch-hamilton.sh $W/runs-dev/hymenoptera morph $k $k shared 1
done
```

Each chain *k* started from the last tree logged by v1.4.1 chain *k* at launch (2 November
18:33). It kept that topology, node ages and fossil ages, and the scalar parameters below
(`run-hamilton.sh`, `START_FROM`):
- `te-shared`: 12 values — `diversification`, `turnover`, `sampling_prop`, `clock_rate`,
  `ucld_sd`, `alpha_morph`, `alpha_dna[1..6]`.
- `te-separate`: 14 values — the same, with `clock_morph`, `clock_dna`, `ucld_sd_morph` and
  `ucld_sd_dna` in place of `clock_rate` and `ucld_sd`.
- `morph`: 6 values.

All other variables restarted:
- the UCLN branch multipliers `z` were set to 0;
- the DNA exchangeabilities and frequencies (`er`, `pi`) and the partition-rate proportions
  (`part_props`) took their initial random draws.

| Chain | v1.4.1 iteration started from | Slurm jobs (first, resumes) |
|---|---|---|
| `te-shared_run1` | 86,850 | 19329130, 19329131, 19329189 |
| `te-shared_run2` | 93,400 | 19329132, 19329133, 19329188 |
| `te-shared_run3` | 93,150 | 19329134, 19329135, 19329187 |
| `te-shared_run4` | 83,550 | 19329136, 19329137, 19329186 |
| `te-shared_run5` | 83,250 | 19329138, 19329139, 19329185 |
| `te-shared_run6` | 82,900 | 19329140, 19329141, 19329184 |
| `te-shared_run7` | 83,300 | 19329142, 19329143, 19329183 |
| `te-shared_run8` | 92,300 | 19329144, 19329145, 19329182 |
| `te-separate_run1` | 99,000 | 19329146, 19329147, 19329181 |
| `te-separate_run2` | 99,750 | 19329148, 19329149, 19329180 |
| `te-separate_run3` | 92,400 | 19329150, 19329151, 19329179 |
| `te-separate_run4` | 91,550 | 19329152, 19329153, 19329178 |
| `te-separate_run5` | 92,500 | 19329154, 19329155, 19329177 |
| `te-separate_run6` | 101,700 | 19329156, 19329157, 19329176 |
| `te-separate_run7` | 93,900 | 19329158, 19329159, 19329175 |
| `te-separate_run8` | 106,000 | 19329160, 19329161, 19329174 |
| `morph_run1` | 249,200 | 19329162 |
| `morph_run2` | 248,600 | 19329164 |
| `morph_run3` | 242,200 | 19329166 |
| `morph_run4` | 270,500 | 19329168 |
| `morph_run5` | 275,000 | 19329170 |
| `morph_run6` | 275,400 | 19329172 |

Segments and stopping:
- `te-*`: the first segment and `resume1` each ended at the 3-day limit (5 and 8 November,
  18:33). `resume2` was cancelled by hand on 10 November at 16:33, after a summary showed pooled
  ESS ≥ 200, PSRF < 1.05 and ASDSF < 0.05 for both configurations.
- `morph`: the `resume1` jobs were cancelled before they started, because the chains already met
  those criteria. The first segments ended at the 3-day limit on 5 November, 18:33.

## 3. Summaries

```bash
# copy of this directory with data/; one job per configuration (MCC needs memory)
cd $W/summary-v2-hymenoptera
export RB=/nobackup/pjjg18/revbayes-dev-a2a8e2db/projects/cmake/build/rb
CONFIGS=te-shared   bash summarise-hamilton.sh $W/runs-dev/hymenoptera 0.25   # sbatch --mem=96G
CONFIGS=te-separate bash summarise-hamilton.sh $W/runs-dev/hymenoptera 0.25   # sbatch --mem=96G
CONFIGS=morph       bash summarise-hamilton.sh $W/runs-dev/hymenoptera 0.25   # sbatch --mem=32G
# fossil ages in every logged tree (8 cores)
ls $W/runs-dev/hymenoptera/*_run*/output/hymenoptera_*_run?.trees | xargs -P 8 -I{} \
  sh -c 'Rscript check-fossil-ages.R data/taxa.tsv fossilcheck/$(basename {} .trees).csv {}'
```
