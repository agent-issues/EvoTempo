# Cambrian trilobites: conventional tip-dating baseline in RevBayes

UCLN + FBD comparator for the CET analyses of the Paterson, Edgecombe & Lee 2019 trilobite data
(inputs in `..`). The data are 107 fossil species with 106 analysed discrete characters, so there
is a single configuration (`morph`), run as six independent chains. The reported runs use the RevBayes development
build at commit `a2a8e2db` (`rb --version`: `1.4.2-preview`, `/nobackup/pjjg18/revbayes-dev-a2a8e2db` on
Hamilton), not the v1.4.1 release and not the cluster's `revbayes/1.2.5` module. A data-free check
(see **Build check**) gave very different sampled-ancestor priors for v1.4.1 and the development
build, so the first set of 6 v1.4.1 chains was superseded. The development-build chains started
from their last states.

## Files

| File | Content |
|---|---|
| `prepare-data.R` | Writes `data/`: one NEXUS matrix per Mkv partition (unordered K = 2, 3, 4, 5, 6, 7, 9, 16; ordered K = 4, 5); `taxa.tsv` (ages in Myr before 485.4 Ma); `partitions.Rev` (partition names and sizes, the ingroup) |
| `start-tree.R` | Starting tree: the source's BEAST2 FBD median consensus tree, see below (`data/start_tree.nwk`) |
| `trilobites_fbd.Rev` | The model. A driver sets `config = "morph"`, `rep`, `n_gen`, `sample_every`, `burnin_gen` |
| `run-hamilton.sh`, `launch-hamilton.sh` | One chain per Slurm job; chains run on `long` (7 days), see **Checkpoints** |
| `summarise-hamilton.sh`, `summarise.R`, `mcc.Rev` | Pooled post-burn-in trees, MCC tree and the CSVs in `summary/` |
| `check-diagnostics.R` | Checks that the coda-free ESS, PSRF and HPD in `summarise.R` reproduce coda |
| `published.csv` | Root and ingroup ages from the source's consensus trees (`../raw/Dataset_4_BEAST_and_BEAST2.zip`) |
| `summary/` | Results (below). Pooled tree files and traces are not committed |

## Model and priors

| Component | This analysis (RevBayes) | Paterson et al. 2019 |
|---|---|---|
| Tree prior | `dnFBDP(originAge, lambda, mu, psi, rho = 0)` with sampled ancestors; time 0 = 485.4 Ma; `condition = "sampling"` | BEAST2: `SABirthDeathModel` (FBD with sampled ancestors, strict clock, discrete data only). BEAST 1.8 main analyses: birth-death serial sampling without sampled ancestors |
| Net diversification, turnover, sampling proportion | LogNormal(−3.5, 1.5); U(0, 1); U(0, 1) | BEAST2: the same |
| Origin | U(0, 1000) Myr above the youngest tip (486.0 Ma) | BEAST2: origin ~ U(0, 1000) above the youngest tip |
| Tip ages | fixed point ages (`../trilobites_tips.csv`) | fixed |
| Topology | ingroup (all but the three fallotaspidoids) monophyletic | the same constraint |
| Clock | UCLN: mean rate ~ U(0, 100) (real-space mean); coefficient of variation ~ Exponential(mean 1/3), log-SD = √ln(1 + CV²); non-centred multipliers | BEAST 1.8 UCLN analysis: the same priors. BEAST2: strict clock |
| Morphology | Mkv (`coding = "variable"`) + Γ4, one shape ~ Exponential(mean 0.5); partitions by number of states K; characters 45 (K = 5) and 47 (K = 4) ordered; character 17 excluded | BEAST 1.8: Mkv partitioned by number of states, partition rates scaled, Γ shape ~ Exponential(mean 0.5), the same two characters ordered. BEAST2: 106 characters in unordered and ordered partitions with relative rates |
| Other traits | not used (outside Mk) | BEAST 1.8 main analyses: 2 meristic and 6 continuous traits with separate clocks |

### FBD conditioning for an all-fossil tree

All 107 tips are Cambrian (486-519 Ma) and none is extant. The process runs from the origin to
the end of the sampling window, the Cambrian-Ordovician boundary at 485.4 Ma (GTS 2012). That
boundary is time 0, with no sampling there (ρ = 0). Lineages alive at 485.4 Ma are simply
unobserved, which is the right treatment of a sampling window: trilobites persisted to the end of
the Permian, but the data sample only the Cambrian.

The tree is conditioned on at least one sampled lineage (`condition = "sampling"`). Conditioning
on survival to time 0 (`"survival"`) is not possible: with ρ = 0 the probability of a sample at
time 0 is zero. The BEAST 1.8 analyses stop the process at the youngest tip (486.0 Ma) with no
extant sample. The BEAST2 analyses also stop there, and tie ρ to the sampling proportion so that
the youngest tip can be a ρ-sample. Here the youngest tip is 0.6 Myr before the end of the window
and is an ordinary ψ-sample.

### Differences from the source analysis

1. **Clock and tree prior combine two source analyses.** The tree prior is the BEAST2 FBD analysis
   (the only source analysis with sampled ancestors). The clock and Γ priors are the BEAST 1.8 UCLN
   analysis (the only source analysis with a relaxed clock).
2. **Discrete characters only, 106 of 107.** Character 17 shows one state in all taxa; it carries no
   information under Mkv and cannot be scored with `coding = "variable"`. The BEAST2 files also
   have 106 characters, which is consistent with this being the character that analysis dropped
   (the curated README leaves this open). That has not been checked column by column: the BEAST2
   matrices recode some characters. The 2 meristic and 6 continuous traits of the BEAST 1.8 analyses are not used.
3. **One relative rate for all Mkv partitions.** The partitions share one rate; the source scales
   partition rates. K is the highest observed state + 1, so character 76 (states 0-1 and 3-8) has
   K = 9, as in the BEAST 1.8 9-state partition.
4. **No tree-scaling moves.** In both v1.4.1 and the development build, `mvUpDownScale`,
   `mvSubtreeScale` and `mvLayeredScaleProposal` give biased results on FBD trees with sampled
   ancestors; see `../../hymenoptera/revbayes/README.md` (**RevBayes up-down move**). Node ages move by
   `mvNodeTimeSlideUniform` and `mvRootTimeSlideUniform`; topology by `mvFNPR`, `mvNarrow` and
   `mvCollapseExpandFossilBranch`. There is no Metropolis coupling.

## Starting tree

The BEAST2 FBD median consensus tree (`CambrianTrilos_106_DisGamStr_SABD_time_Median_con.tre`).
Heights are converted to ages before 485.4 Ma, and tips are set to their fixed ages. Zero-length
sampled-ancestor branches become 0.01-Myr side branches, and every node is made at least 0.01 Myr
older than its oldest child. The start has root 541.4 Ma and ingroup 539.9 Ma. `readTrees()` puts
the youngest tip of a clock tree at age 0, so the model shifts the tree by the youngest tip's age
(`start_tree.offset(youngest_tip_age)`).

## Checkpoints

`initializeFromCheckpoint()` cannot restore an all-fossil tree with point ages: the tree is re-read
with its youngest tip at age 0, which lies outside that tip's fixed age ("The age of tip ... in the
initial tree is outside of specified bounds"). The chains therefore run on Hamilton's `long`
partition (7 days) in one job each, rather than as 3-day jobs resumed from checkpoints.

## Build check

This is a data-free check: the same model with `mymodel.ignoreAllData()`, 2 chains per build,
200 000 iterations sampled every 100, and 10 % burn-in. Cells give mean ± MCSE, with the two chain
means in brackets. Logs are on Hamilton under
`/nobackup/pjjg18/evotempo-baselines-rb141/prior-check-tri/`.

| Quantity | v1.4.1 | development `a2a8e2db` |
|---|---|---|
| sampled ancestors | 100.6 ± 0.3 [99.6, 101.7] | 6.49 ± 0.13 [6.36, 6.62] |
| sampling proportion s | 0.953 ± 0.002 | 0.0111 ± 0.0004 |
| root age (Ma) | 519.40 ± 0.01 | 519.62 ± 0.01 |
| ingroup age (Ma) | 518.48 ± 0.01 | 518.87 ± 0.01 |

Under v1.4.1, about 101 of the 107 tips become sampled ancestors under the prior, even though all
ages are points. That rules out the interval-age `mvCollapseExpandFossilBranch` bug that affects
the penguin and Hymenoptera data. The development build also contains sampled-ancestor fixes,
including `mvNarrow` reversibility with sampled-ancestor tips. The v1.4.1 chains were therefore
superseded. Their posterior had about 3 sampled ancestors, so the data dominate the prior there,
but the v1.4.1 results are not reported.

## Results

Six chains, run on Hamilton (`long` partition, one 7-day job each, development build `a2a8e2db`).
They started from the last states of the superseded v1.4.1 chains and ran 778 400–849 400
iterations, sampled every 200. `summarise.R` discards the first 25 % of each chain; the shortest
chain keeps 2919 samples. All files are in
`summary/`: `headline_morph.txt`, `convergence_morph.csv`, `parameters_morph.csv`,
`published_morph.csv`, `fossils_morph.csv` (per-fossil sampled-ancestor probability),
`clades_morph.csv` (clade frequencies per chain) and `mcc_morph.nex`.

### Convergence

The criteria were met: ASDSF < 0.05, PSRF < 1.05 and pooled ESS ≥ 200 for every logged parameter.

| Quantity | min ESS per chain | pooled ESS | PSRF |
|---|---|---|---|
| Posterior | 1419 | 9384 | 1.000 |
| root age | 1373 | 9280 | 1.000 |
| ingroup age | 1365 | 8729 | 1.000 |
| sampled ancestors | 2115 | 14347 | 1.000 |
| sampling proportion | 1396 | 11887 | 1.000 |
| clock rate | 1465 | 9524 | 1.000 |
| UCLN CV | 2141 | 14433 | 1.000 |
| λ, μ (lowest ESS) | 452 | 6542 | 1.010 |

Topology: ASDSF 0.0079 over 186 clades with frequency ≥ 0.1 in any chain (max SDSF 0.021).

### Estimates

Median (95 % HPD); ages in Ma.

| Parameter | This analysis |
|---|---|
| root (first split) | 528.53 (524.53–533.57) |
| ingroup crown | 527.52 (523.78–531.94) |
| origin | 531.05 (525.28–539.62) |
| sampled ancestors | 3 (0–7); mean 2.92 |
| sampling proportion s | 0.0026 (< 0.0001–0.0132) |
| net diversification λ − μ (per Myr) | 0.031 (0.003–0.068) |
| turnover μ/λ | 0.992 (0.972–1.000) |
| clock rate (substitutions per character per Myr) | 0.0236 (0.0165–0.0321) |
| UCLN coefficient of variation | 0.79 (0.49–1.13) |
| Γ shape | 0.52 (0.37–0.70) |

No fossil is a sampled ancestor with probability above 0.3. The highest are
*Kuetsingocephalus qiannanensis* (0.28), *Myopsolenites altus* (0.23) and *Xingrenaspis
xingrenensis* (0.23).

### Comparison with the published analyses

Published values are the median heights in the source's consensus trees plus 486.0 Ma (the
youngest tip), with 95 % HPD (`published.csv`).

| Analysis (Paterson et al. 2019) | Root | Ingroup |
|---|---|---|
| **This analysis** (FBD + SA, UCLN, discrete only) | **528.53 (524.53–533.57)** | **527.52 (523.78–531.94)** |
| BEAST 1.8 UCLN, discrete + meristic + continuous | 529.50 (525.59–534.30) | 528.73 (524.99–533.13) |
| BEAST 1.8 epoch clock, all data (their Fig. 1) | 535.23 (530.02–541.30) | 534.40 (529.40–539.95) |
| BEAST 1.8 strict clock, all data | 537.58 (532.26–543.78) | 536.69 (531.88–542.82) |
| BEAST2 FBD + SA, strict clock, discrete only | 541.42 (534.26–548.33) | 539.92 (533.81–546.74) |
| BEAST2 BDSS, strict clock, discrete only | 540.87 (534.53–547.89) | 539.26 (533.73–545.88) |

The ages agree with the source's one relaxed-clock analysis (BEAST 1.8 UCLN): the medians differ by
about 1 Myr and the HPDs overlap almost entirely. They are 7–13 Myr younger than the strict-clock
and epoch-clock analyses. The clock model, not the tree prior, accounts for the difference.
Our tree prior is that of the BEAST2 FBD + SA analysis, yet our ages match the BEAST 1.8 UCLN
analysis, which has no sampled ancestors. Among the strict-clock analyses, swapping FBD + SA for
BDSS moves the root by only 0.5 Myr.

`mcc_morph.nex` is the maximum-clade-credibility tree with RevBayes' node annotations. The tree
files store branch lengths only, so its node ages are in Myr before 486.0 Ma (the youngest tip).
Add 486.0 for ages in Ma.

## Compute

Core-hours on Hamilton (Slurm `CPUTimeRAW`, one core per chain job):

| Jobs | Core-hours |
|---|---|
| Reported chains: 6 × 7 days, development build (`evotempo-tril-dev-morph-*`) | 1008.0 |
| Superseded v1.4.1 chains (`evotempo-tril-rb141-morph-*`) | 476.9 |
| Build check, both builds (`evotempo-tril-prior-check`) | 30.2 |
| Smoke test | 0.1 |
| **Total** | **1515.2** |

The summary jobs ran on 4 cores, about 2 core-hours in all, shared with the Hymenoptera baseline.
