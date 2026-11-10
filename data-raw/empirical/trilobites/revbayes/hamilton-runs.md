# How the reported chains were run on Hamilton

This file records the commands that produced the six chains reported in `README.md`.
`launch-hamilton.sh` and `run-hamilton.sh` cover these steps, but the reported chains also used
an exported `RB` and a per-chain `START_FROM`. An exact rerun is not possible: the
`START_FROM` directories (the v1.4.1 chains) have been deleted. The starting states those
directories gave are committed in `start-states/`.

`W=/nobackup/pjjg18/evotempo-baselines-rb141`; `T=$W/src-tri/data-raw/empirical/trilobites/revbayes`
(a clone of this branch, at `88bf486` for the v1.4.1 launch and `aad8f23` for the reported chains (the clone's reflog shows HEAD at aad8f23 from 2 Nov 18:32 to 3 Nov 08:17; all six chains started 2 Nov 18:44–19:20; the model file is identical, md5 e3ba8720, to the one at a3ad184), with `data/`
written by `Rscript prepare-data.R; Rscript start-tree.R`).

## 1. Superseded v1.4.1 chains (31 October – 3 November 2026)

At that time `run-hamilton.sh` defaulted to `RB=/nobackup/pjjg18/revbayes-1.4.1/projects/cmake/build/rb`
(the v1.4.1 release) and job names were `evotempo-tril-rb141-morph-<rep>`.

```bash
export N_GEN=2000000 BURNIN_GEN=2000 CKP_EVERY=500
SAMPLE_EVERY=200 bash $T/launch-hamilton.sh $W/runs/trilobites morph 1 6 long 0
```

The chains started from the BEAST2 consensus tree (`start-tree.R`). They ran from 31 October
00:46 until they were cancelled on 3 November at 08:15, after the build check. They are not
reported, and their output (`$W/runs`) has been deleted.

## 2. Reported chains, development build `a2a8e2db` (2–9 November 2026)

```bash
export RB=/nobackup/pjjg18/revbayes-dev-a2a8e2db/projects/cmake/build/rb BUILD_TAG=dev \
       N_GEN=3000000 BURNIN_GEN=2000 CKP_EVERY=500 SAMPLE_EVERY=200
for k in 1 2 3 4 5 6; do
  START_FROM=$W/runs/trilobites/morph_run$k/output \
    bash $T/launch-hamilton.sh $W/runs-dev/trilobites morph $k $k long 0
done
```

The build was selected by exporting `RB`. The run-directory copies of `run-hamilton.sh` are the
older ones, which default to v1.4.1, and their `/usr/bin/time` files are empty because the jobs
ended at the time limit. The build is shown by the start-up warnings in every chain's screen
output: `Unknown user setting with key 'moduledir'` is printed only by the development build. The
committed `run-hamilton.sh` now defaults to `a2a8e2db` and writes `rb --version` to
`out/rb_version_*.txt`.

Each chain *k* started from the last tree that v1.4.1 chain *k*, still running, had logged when
the job started (table below). The chain kept that topology and node ages, and took seven scalar values from the
same iteration: `diversification`, `turnover`, `sampling_prop`, `origin_time`, `clock_morph`,
`ucld_cv` and `alpha_morph`. The UCLN branch multipliers took their initial draws.
`start-states/` holds, per chain, the start tree (`start_tree_morph_<k>.nwk`), the values
(`start_values_morph_<k>.Rev`) and the driver (`driver_morph_<k>.Rev`) that `run-hamilton.sh`
wrote.

| Chain | v1.4.1 iteration started from | Slurm job | Ran |
|---|---|---|---|
| `morph_run1` | 217,600 | 19329207 | 2 Nov 18:45 – 9 Nov 18:45 |
| `morph_run2` | 241,800 | 19329208 | 2 Nov 18:45 – 9 Nov 18:45 |
| `morph_run3` | 243,200 | 19329209 | 2 Nov 19:06 – 9 Nov 19:06 |
| `morph_run4` | 240,200 | 19329210 | 2 Nov 19:20 – 9 Nov 19:20 |
| `morph_run5` | 203,000 | 19329211 | 2 Nov 19:20 – 9 Nov 19:20 |
| `morph_run6` | 217,600 | 19329212 | 2 Nov 19:20 – 9 Nov 19:20 |

`N_GEN = 3e6` is a ceiling, never reached. Each chain ran in one 7-day job on `long` and ended
at the time limit. There were no resumes: see **Checkpoints** in `README.md`.

## 3. Summary

```bash
cd $W/summary-dev-trilobites      # a copy of this directory with data/
CONFIGS=morph bash summarise-hamilton.sh $W/runs-dev/trilobites 0.25
```
