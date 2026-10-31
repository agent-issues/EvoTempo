# Cambrian trilobites: conventional tip-dating baseline in RevBayes 1.4.1

UCLN + FBD comparator for the CET analyses of the Paterson, Edgecombe & Lee 2019 trilobite data
(inputs in `..`). The data are 107 fossil species with 106 analysed discrete characters, so there
is a single configuration (`morph`), run as six independent chains. The run uses RevBayes v1.4.1
(tag `v1.4.1`, commit `3a504eb`), built on Hamilton (`/nobackup/pjjg18/revbayes-1.4.1`), not the
cluster's `revbayes/1.2.5` module.

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
4. **No tree-scaling moves.** RevBayes 1.4.1's `mvUpDownScale` with a tree has a wrong Hastings
   ratio; see `../../hymenoptera/revbayes/README.md` (**RevBayes up-down move**). Node ages move by
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

## Results

{{RESULTS}}

## Compute

{{COMPUTE}}
