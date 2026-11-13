# Penguins: conventional total-evidence tip dating in RevBayes

This is the comparator for the CET analyses of the penguin pilot dataset, and the source of their starting trees. It is a fossilized birth–death (FBD) analysis with sampled ancestors and an uncorrelated lognormal (UCLN) relaxed clock. It follows Gavryushkina et al. (2017, *Syst. Biol.* 66:57–73, doi:10.1093/sysbio/syw060), analysis "8+dna R" (`../raw/xmls.zip`, `xmls/8+dna_R.xml`), as closely as RevBayes allows.

Two clock models were run:

- **shared** (the EvoTempo comparator): one set of UCLN branch-rate multipliers scales both morphology and DNA, each with its own mean rate. This is the conventional analogue of CET, in which one tempo per lineage drives every character.
- **separate** (the source analysis): independent UCLN clocks for morphology and DNA, as in the BEAST2 XML.

**Software.** The headline results (round 4) come from RevBayes built from the `development` branch at commit `a2a8e2db8a1b57eee123c6b8e644335ccc268692` (29 Oct 2026, reports itself as 1.4.2-preview), not from the v1.4.1 release. v1.4.1 has two moves that are biased on trees with sampled ancestors, and rounds 1–3 used them; see [RevBayes faults](#revbayes-faults). v1.4.2 (tagged 6 Nov 2026) is a2a8e2db plus 101 commits, none of which touch the FBD distribution or the tree moves, so these results stand for 1.4.2 too. Later development builds (fd807d10, used for new runs from 11 Nov) fix `mvTreeScale` and `mvSubtreeScale` with sampled ancestors, neither of which is used here. As an independent reference, the source XML was also rerun in BEAST2 (see [BEAST2 reference](#beast2-reference)).

## Headline results (round 4)

Convergence: 4 chains per clock model, the first 25 % of each discarded. R̂ is the rank-normalised split R̂ (the larger of the bulk and folded-tail values) and ESS the rank-normalised bulk and tail ESS of Vehtari et al. (2021, *Bayesian Anal.* 16:667–718), as implemented in MkPrime (`MkPrime:::.Rhat`, `.EssBulk`, `.EssTail`; `rhat-ess-quartet.R`). They are over all logged continuous parameters except Posterior, Likelihood and Prior. ASDSF is over clades with frequency ≥ 0.1 in any chain. `convergence_<clock>.csv` also gives the classical per-chain ESS and Gelman–Rubin PSRF, as a secondary check.

| Clock | Chains | Iterations per chain | Samples per chain after burn-in | Max R̂ | Min bulk ESS | Min tail ESS | Min per-chain bulk ESS | ASDSF (max SDSF) | Run time per chain (h) |
|---|---|---|---|---|---|---|---|---|---|
| shared | 4 | 150,000 | 1,876 | 1.020 (`node_11`) | 404 (`node_11`) | 540 (`node_11`) | 37 (`node_4`) | 0.0098 (110 clades; max 0.034) | 55–67 |
| separate | 4 | 220,000 | 2,751 | 1.031 (`node_4`) | 243 (`node_4`) | 566 (`node_3`) | 58 (`node_4`) | 0.0162 (114 clades; max 0.143) | 62–72 + 22–29 (extension) |

Every R̂ is below 1.05, and every bulk and tail ESS is above 240. Neither model meets the stricter R̂ < 1.01 of Vehtari et al. for every parameter. In the shared model, 7 of 86 parameters exceed 1.01, all of them extant node ages (nodes 4–7, 11, 13, 14), and every bulk ESS is above 400. In the separate model, 12 of 87 exceed 1.01: the crown age, nine extant node ages (nodes 3–5, 7–9, 11–13), the molecular clock mean and its log-SD. Nine of these have a bulk ESS between 243 and 400. The separate-clock chains were therefore extended from 150 000 to 220 000 iterations from their checkpoints, and truncated at 220 000 (`summary/trunc_log.txt`). The slowest feature is the basal split of the crown. In the separate model, the clade of all extant penguins except *Eudyptula* and *Spheniscus*, with their fossil relatives, has frequency 0.34–0.68 among chains (max SDSF 0.143). The node ages at the base of the crown (nodes 3–5) are also the slowest-mixing ages. The per-chain bulk ESS is a split-chain value, so it also penalises drift between a chain's two halves. Its minimum in the shared model (37, node 4 in one chain) reflects such slow drift. Node 4's pooled bulk ESS is 584 and its R̂ is 1.010.

### Comparison with the published analysis and with the BEAST2 rerun

The published values are from the arXiv version of the paper and its SM (arXiv:1506.04797; SM Table 4, Fig. 5 node numbers, main text). "BEAST2 rerun" is the source XML run again here (5 chains; see [BEAST2 reference](#beast2-reference)). The table shows mean and 95 % HPD interval. Ages are in Ma and are common-ancestor ages, as in the paper: the age of the MRCA of the listed extant taxa in every tree, whether or not the group is monophyletic in that tree. `crown_age` and `node_*` in `penguins_fbd.Rev` are the same quantity.

| Node | Published | BEAST2 rerun | RevBayes shared | RevBayes separate |
|---|---|---|---|---|
| Crown age (node 1) | 12.7 [9.91, 15.8] | 12.5 [10.0, 15.8] | 12.1 [10.4, 14.0] | 12.5 [10.3, 15.7] |
| Node 2: Megadyptes + Eudyptes + Eudyptula + Spheniscus | 11.2 [7.92, 13.0] | 11.2 [8.34, 13.5] | 12.1 [10.4, 14.1] | 11.5 [9.84, 13.5] |
| Node 3: Eudyptula + Spheniscus | 10.3 [7.14, 12.0] | 10.4 [7.32, 12.9] | 11.9 [9.89, 14.1] | 10.8 [8.69, 13.2] |
| Node 4: Aptenodytes + Pygoscelis | 9.83 [5.84, 14.0] | 9.59 [5.91, 14.0] | 6.62 [4.02, 9.16] | 9.27 [5.57, 14.2] |
| Node 5: Pygoscelis | 6.12 [3.30, 9.04] | 5.96 [3.22, 8.94] | 3.74 [1.81, 5.70] | 5.73 [3.11, 9.17] |
| Node 6: Megadyptes + Eudyptes | 5.08 [3.27, 6.86] | 5.06 [3.21, 6.89] | 4.29 [1.77, 7.06] | 5.05 [3.12, 7.03] |
| Node 7: P. antarctica + P. papua | 3.47 [1.68, 5.27] | 3.38 [1.71, 5.24] | 1.96 [0.753, 3.48] | 3.25 [1.52, 5.25] |
| Node 8: Eudyptes | 2.42 [1.49, 3.33] | 2.39 [1.49, 3.35] | 1.96 [0.999, 3.04] | 2.34 [1.36, 3.29] |
| Node 9: E. chrysocome, filholi, moseleyi, chrysolophus, schlegeli | 2.20 [1.11, 2.85] | 2.17 [1.19, 3.22] | 1.54 [0.744, 2.54] | 2.08 [1.09, 3.18] |
| Node 10: E. pachyrhynchus, robustus, sclateri | 2.10 [1.05, 2.60] | 2.08 [1.06, 3.05] | 1.85 [0.785, 3.01] | 2.08 [1.09, 3.09] |
| Node 11: Spheniscus | 1.65 [0.990, 2.32] | 1.62 [0.977, 2.32] | 1.40 [0.517, 2.42] | 1.61 [0.929, 2.34] |
| Node 12: Aptenodytes | 1.63 [0.740, 2.49] | 1.60 [0.746, 2.50] | 1.24 [0.263, 2.47] | 1.57 [0.661, 2.52] |
| Node 13: S. humboldti + S. mendiculus | 1.06 [0.630, 1.88] | 1.23 [0.605, 1.88] | 0.908 [0.237, 1.74] | 1.22 [0.564, 1.88] |
| Node 14: S. demersus + S. magellanicus | 1.06 [0.540, 1.60] | 1.04 [0.508, 1.58] | 0.808 [0.227, 1.56] | 1.02 [0.466, 1.59] |
| Node 15: E. chrysocome, filholi, moseleyi | 0.770 [0.350, 1.23] | 0.758 [0.345, 1.22] | 0.579 [0.184, 1.05] | 0.737 [0.305, 1.21] |
| Node 16: E. pachyrhynchus + E. robustus | 0.760 [0.380, 1.16] | 0.754 [0.368, 1.16] | 0.554 [0.149, 1.05] | 0.731 [0.322, 1.14] |
| Node 17: E. chrysolophus + E. schlegeli | 0.440 [0.190, 0.710] | 0.437 [0.189, 0.713] | 0.335 [0.077, 0.715] | 0.404 [0.157, 0.689] |
| Node 18: E. chrysocome + E. filholi | 0.320 [0.090, 0.580] | 0.316 [0.095, 0.582] | 0.233 [0.036, 0.488] | 0.307 [0.078, 0.573] |

| Parameter | Published | BEAST2 rerun | RevBayes shared | RevBayes separate |
|---|---|---|---|---|
| Net diversification d | 0.039 [0.002, 0.089] | 0.038 [0.0016, 0.090] | 0.037 [0.00089, 0.089] | 0.037 [0.0018, 0.085] |
| Turnover ν | 0.880 [0.720, 1.00] | 0.885 [0.720, 0.997] | 0.904 [0.761, 0.997] | 0.892 [0.740, 0.997] |
| Sampling proportion s | 0.230 [0.060, 0.430] | 0.231 [0.061, 0.428] | 0.192 [0.043, 0.371] | 0.207 [0.051, 0.398] |
| λ − μ − ψ | 0.019 [−0.061, 0.099]¹ | -0.050 [-0.119, 0.020] | -0.044 [-0.108, 0.027] | -0.042 [-0.108, 0.024] |
| Sampled ancestors (count) | not reported (≈ 17–18 drawn in Fig. 3) | 16.7 [9.00, 23.0] | 15.4 [9.00, 22.0] | 15.5 [8.00, 22.0] |
| *S. muizoni* ancestor of the four extant *Spheniscus* | 0.61 | 0.616 | 0.589 | 0.641 |
| Origin (Ma) | — | 64.7 [60.6, 71.4] | 65.0 [60.6, 72.0] | 65.3 [60.6, 72.7] |

¹ Quoted in the paper's discussion of the birth–death analysis without sampled ancestors, so it is probably not the FBD value; the BEAST2 rerun of the FBD analysis gives the value in the next column.

**Crown placement** (posterior probability that the fossil lies inside the crown; Gavryushkina et al. Table 3):

| Fossil | Published | BEAST2 rerun | RevBayes shared | RevBayes separate |
|---|---|---|---|---|
| *Spheniscus megaramphus* | 0.9992 | 0.9991 | 0.9836 | 0.9934 |
| *Spheniscus urbinai* | 0.9991 | 0.9991 | 0.9833 | 0.9933 |
| *Pygoscelis grandis* | 0.9928 | 0.9945 | 0.9987 | 0.9924 |
| *Spheniscus muizoni* | 0.9201 | 0.9207 | 0.9892 | 0.9777 |
| *Madrynornis mirandus* | 0.9007 | 0.8946 | 0.9941 | 0.9820 |
| *Marplesornis novaezealandiae* | 0.1652 | 0.1717 | 0.0300 | 0.0965 |

**Fossils with positive evidence** in Gavryushkina et al. (Fig. 4; log Bayes factor > 1 or < −1). The columns give the posterior probability of being a sampled ancestor:

| Fossil | Published evidence | BEAST2 rerun | RevBayes shared | RevBayes separate |
|---|---|---|---|---|
| *Palaeospheniscus patagonicus* | sampled ancestor | 0.707 | 0.668 | 0.652 |
| *Icadyptes salasi* | sampled ancestor | 0.741 | 0.748 | 0.746 |
| *Eretiscus tonnii* | terminal | 0.169 | 0.133 | 0.132 |
| *Marplesornis novaezealandiae* | terminal | 0.206 | 0.227 | 0.195 |
| *Paraptenodytes antarcticus* | terminal | 0.023 | 0.019 | 0.020 |
| *Pygoscelis grandis* | terminal | 0.296 | 0.443 | 0.321 |

What the comparison shows:

1. **The separate-clock run reproduces the source analysis.** Its crown age, extant node ages, diversification parameters, sampled-ancestor count and per-fossil probabilities agree with both the published values and the BEAST2 rerun to within Monte Carlo error. One difference stands out: the near-crown fossils *Madrynornis mirandus* and *S. muizoni* are placed inside the crown more often (0.98–0.99 against 0.90–0.92), and *Marplesornis novaezealandiae* less often (0.03–0.10 against 0.17). Candidate causes are among the model differences listed below (character 81, Gamma discretisation, continuous UCLN).
2. **Linking the clocks (shared) changes the extant topology at the base of the crown.** The published basal split is *Aptenodytes* + *Pygoscelis* against the rest. With linked clocks, the common ancestor of node 2 (*Megadyptes*, *Eudyptes*, *Eudyptula*, *Spheniscus*) is the crown root in 0.985 of samples, against 0.563 (separate) and 0.404 (BEAST2 rerun). Nodes 4 and 5 (*Aptenodytes* + *Pygoscelis*, *Pygoscelis*) are 2.7 and 2.0 Myr younger and the molecular clock mean is 70 % faster. The crown age itself is unchanged (12.1 against 12.5 Ma).
3. **The sampling proportion.** After the fix, s (0.19–0.21) agrees with the published value (0.23). The definitions agree (below); the earlier gap (s ≈ 0.50 in rounds 1–3) came from the v1.4.1 bias.

**Summary trees.** `summary/quartet_consensus_<clock>.nwk` is the quartet consensus (`ConsTree::Quartet()`, default settings) of 4000 pooled post-burn-in samples, thinned evenly. It resolves 42 (shared) and 41 (separate) of the 52 possible splits. Its node labels are the splits' frequencies in all pooled samples, and `split_support_<clock>.csv` gives them per chain. The consensus is unrooted and carries no node ages. Node ages are the common-ancestor ages in the tables above (`clade_ages_<clock>.csv`). No MCC tree is used.

### Why the sampling-proportion definitions match

Both programs use the same reparameterisation and the same tree density:

- **Reparameterisation.** BEAST2 SA (`SABirthDeathModel.transformParameters`, SA 1.1.7) sets λ = d/(1 − ν), μ = νλ, ψ = sμ/(1 − s), so s = ψ/(μ + ψ). `penguins_fbd.Rev` defines `lambda`, `mu` and `psi` in the same way from `diversification`, `turnover` and `sampling_prop`.
- **Constant rates; ρ and r fixed.** Rates are constant through time in both. ρ = 1 and r = 0 are fixed in the XML: `rho` and `removalProbability` are not in the state and have no operator. `dnFBDP` has no removal parameter (r = 0).
- **Conditioning.** With an origin, RevBayes `condition = "survival"` divides by the probability that the origin lineage leaves at least one sample *at the present* (`BirthDeathSamplingTreatmentProcess::pSurvival`, computed without fossil sampling). This is BEAST2's `conditionOnRhoSampling = true`. RevBayes's `condition = "sampling"` would be BEAST2's `conditionOnSampling`.
- **Tree shape.** Both apply a factor of 2 per bifurcation that is not a sampled ancestor (BEAST2: `(leaves − SA − 1) log 2`; RevBayes `lnProbTreeShape`: `(tips − SA − 1) ln 2 − ln(tips!)`, which differs only by a constant).
- **Fossil ages.** Both treat the fossil ages as points inside uniform intervals.

So s needs no conversion. The paper's own definition (its eq. 2) is s = ψ/(μ + ψ). Within every posterior, s rises steeply with the number of sampled ancestors (`summary/sa_vs_s_*.csv`). This is why the v1.4.1 excess of sampled ancestors carried s up to 0.50.

### Clock rates

Units are changes per character per Myr (morphology; normalised Mk matrices, Γ4 with mean 1) and substitutions per site per Myr (DNA), in both programs. The tree is in Myr.

- **Locus rates.** In both programs the five locus rate multipliers sum to 1 (BEAST2: unweighted `DeltaExchangeOperator` on `mutationRate`s that start at 0.2; RevBayes: `locus_rates ~ dnDirichlet(1, …, 1)`). A locus's rate is therefore the clock mean × its multiplier (about 0.2 on average), and the molecular "clock mean" is a rate summed over loci, not a per-site average.
- **Branch multipliers.** The UCLN multipliers have mean 1 in real space in both programs, so the clock mean is the mean branch rate. For the shared model, `clock_morph` and `clock_dna` are the two mean rates on the one set of multipliers, so they compare directly with the separate model's means.
- **What the paper reports.** The paper gives no clock-rate estimates. Its simulation study (SM) used values "close to the parameters estimated in the total-evidence analysis": μ_morph = 0.0138, σ_morph = 0.728, μ_mol = 0.0436, σ_mol = 0.147 (²). The main text reports a coefficient of variation of the morphological rates of 1.15.
- **Coefficient of variation (³).** BEAST2's `RateStatistic` computes it as the SD over the mean of the branch rates, unweighted, over all branches including the zero-length sampled-ancestor branches. `derived-quantities.py` computes the same statistic from the `_rates.trees` files (`cv_branch_all_*`). Its "mean" is length-weighted (`mean_rate_lenweighted_*`).

| Quantity | Published | BEAST2 rerun | RevBayes shared | RevBayes separate |
|---|---|---|---|---|
| Morphological clock mean (changes / character / Myr) | ≈ 0.0138² | 0.033 [0.018, 0.051] | 0.023 [0.013, 0.033] | 0.027 [0.013, 0.043] |
| Molecular clock mean (substitutions / site / Myr, summed over the five loci) | ≈ 0.0436² | 0.035 [0.023, 0.048] | 0.062 [0.042, 0.087] | 0.036 [0.023, 0.051] |
| Branch-length-weighted mean morphological rate | — | 0.034 [0.020, 0.051] | 0.025 [0.015, 0.036] | 0.031 [0.017, 0.048] |
| Branch-length-weighted mean molecular rate | — | 0.034 [0.023, 0.045] | 0.067 [0.045, 0.093] | 0.035 [0.023, 0.047] |
| UCLN log-SD, morphology | ≈ 0.728² | 0.916 [0.582, 1.27] | 0.715 [0.466, 0.971] (shared) | 1.09 [0.717, 1.48] |
| UCLN log-SD, DNA | bimodal, modes ≈ 0.4 and 0 | 0.278 [3.3e-09, 0.521] | (shared) | 0.312 [2.7e-07, 0.552] |
| Coefficient of variation of morphological branch rates³ | 1.15 | 1.13 [0.625, 1.71] | 0.875 [0.488, 1.35] | 1.53 [0.796, 2.38] |

² Simulation settings, not posterior estimates.
³ Unweighted SD/mean over branches, as BEAST2's `RateStatistic`.

In the separate model, morphological rates vary more among branches than in the BEAST2 rerun: log-SD 1.09 against 0.92, CV 1.53 against 1.13. A candidate cause is the continuous UCLN: BEAST2's 110 discretised rate categories cap the most extreme branch rates. Linking the clocks (shared) pulls the morphological rate variation down to that of the combined data (log-SD 0.72).

## Files

| File | Content |
|---|---|
| `prepare-data.R` | Writes `data/` from `../penguins_morph.nex`, `../penguins_dna.nex` and `../penguins_tips.csv` |
| `start-tree.R` | Writes `data/start_tree.nwk` from the BEAST2 XML's starting tree (needs ape) |
| `penguins_fbd.Rev` | The model and MCMC; sourced by a driver that sets `clock_model`, `rep`, `n_gen`, `sample_every`, `burnin_gen` and options |
| `run-hamilton.sh` | Slurm script for one chain: writes the driver and runs `rb` |
| `progress.sh` | Prints the last logged state of every chain under a run root |
| `summarise-hamilton.sh` | Runs `summarise.R` and `rhat-ess-quartet.R` for both clock models |
| `summarise.R` | Burn-in, posterior summaries (mean, median, 95 % HPD), clade frequencies and ASDSF, per-fossil probabilities; classical ESS and PSRF as a secondary check |
| `rhat-ess-quartet.R` | Rank-normalised R̂ and bulk/tail ESS (MkPrime), and the quartet consensus tree (`ConsTree::Quartet()`) of the pooled samples with split frequencies |
| `check-diagnostics.R` | Checks that `summarise.R`'s coda-free ESS, PSRF and HPD equal coda's |
| `derived-quantities.py` | ψ, λ − μ − ψ, branch-rate statistics, s by sampled-ancestor count, basal-split frequency |
| `check-fossil-ages.py` | Checks every sampled tree for fossil ages outside their ranges, and counts fossil sampled ancestors at the root |
| `check-sa-count.py` | Checks that the logged `num_samp_anc` equals the number of zero-length fossil branches in each sampled tree |
| `run-beast2-hamilton.sh` | BEAST2 reference run of the source XML (posterior or `PRIOR=true`) |
| `beast-posterior-summary.py`, `beast-sa-summary.py` | Summaries of the BEAST2 posterior and prior runs |
| `published-gavryushkina2017.csv` | Published estimates compared against |
| `summary/` | Round-4 results: `rhat_ess_<clock>.csv` (posterior summaries with R̂ and ESS), `quartet_consensus_<clock>.nwk` (quartet consensus, unrooted, node labels = split frequencies), `split_support_<clock>.csv` (its splits per chain), `parameters_`, `clade_ages_`, `clades_`, `convergence_`, `fossils_`, `derived_`, `sa_vs_s_<clock>.csv`, `headline_<clock>.txt`, `fossil_age_check.txt`, `trunc_log.txt` |
| `summary/beast2/` | BEAST2 rerun: posterior summaries (5 chains, and per chain) and the data-free run (`beast2_prior_summary.txt`) |
| `summary/data-free/` | Data-free runs of RevBayes v1.4.1 and of a2a8e2db |
| `summary/penguins_r4_<clock>_md5.txt` | md5 of every file in the round-4 chain archive (below) |
| `summary/round3-v1.4.1/` | Round 3 (v1.4.1, biased), kept as evidence of the bias |

`data/`, the trace files and the tree files are not committed (25–40 MB per tree file). `prepare-data.R` and `start-tree.R` regenerate `data/`.

**Round-4 chain archive.** The round-4 parameter logs, tree files and rate-tree files of all eight chains (shared 150 000 iterations; separate truncated at 220 000) are archived as two xz tarballs with md5 manifests: `penguins_r4_shared_chains.tar.xz` and `penguins_r4_separate_chains.tar.xz` (28.3 MB and 41.2 MB; SHA-256 `3b775337…484062d` and `39d7632f…c12125b`; md5 of every member file in `summary/penguins_r4_<clock>_md5.txt`). The canonical copy is the EvoTempo project's artifact store (Claude Science project `proj_07b3619548a5`, artifacts `aadac10f-152b-43de-a186-b0d15d7a3b9b` (shared) and `55110c6f-814d-4c45-88a0-63bdf1bf80f1` (separate)); a working copy stays on Hamilton in `/nobackup/pjjg18/evotempo-penguins-summaries/round4-chains-archive/`, which is not backed up. The CET starting trees (`../cet/data/start-tree-selection.txt`) were drawn from the shared-clock chains in this archive. The BEAST2 traces stay on Hamilton (`/nobackup/pjjg18/evotempo-penguins-beast2/`) until this PR is merged.

## Reproduce

```sh
Rscript prepare-data.R          # writes data/
Rscript start-tree.R            # writes data/start_tree.nwk
# RevBayes development a2a8e2db, built on Hamilton with gcc/14.2 + boost/1.88.0 + cmake/3.31:
#   git clone https://github.com/revbayes/revbayes && git checkout a2a8e2db8a1b57eee123c6b8e644335ccc268692
#   cd projects/cmake && ./regenerate.sh && ./build.sh -boost_root $BOOST_ROOT -j 8
# One copy of this directory per chain (<clock>_run<rep>/), one job each:
sbatch -J evotempo-penguin-r4-shared-13 \
  --export=ALL,CLOCK=shared,REP=13,N_GEN=150000,SAMPLE_EVERY=60,BURNIN_GEN=2000 run-hamilton.sh
#   CLOCK = shared | separate.  Round 4 also set START_FROM=<round-3 chain>/output
#   (start from that chain's last tree and parameters); without it a chain starts from
#   data/start_tree.nwk.  A chain stopped by the wall clock, or extended, continues with
#   RESUME=true; N_GEN is then the number of ADDITIONAL iterations (RevBayes appends
#   N_GEN iterations to the checkpointed run).  The separate-clock extension was
#   launched with N_GEN=220000, ran past 220 000 in total, was stopped and truncated
#   at iteration 220 000 (summary/trunc_log.txt).
sbatch summarise-hamilton.sh <run root holding <clock>_run<rep>/> 0.25   # needs ape, TreeTools, ConsTree, MkPrime
python3 derived-quantities.py <run root> <clock> 0.25 5
python3 check-fossil-ages.py data/taxa.tsv <run root>/*_run*/output/penguins_*_run*[0-9].trees
```

## Model and priors

| Component | Gavryushkina et al. 2017 (BEAST2 2.4, SA package) | Here (RevBayes) |
|---|---|---|
| Tree prior | `SABirthDeathModel` with an origin, ρ = 1, r = 0, `conditionOnRhoSampling` | `dnFBDP(originAge = origin_time, …, rho = 1, condition = "survival")`, r = 0 (see above) |
| Parameterisation | d = λ − μ ~ LogNormal(−3.5, 1.5); ν = μ/λ ~ U(0, 1); s = ψ/(μ + ψ) ~ U(0, 1) | same priors; λ, μ, ψ deterministic |
| Origin | U(0, 160) Ma | same |
| Fossil ages | uniform on the interval (`SamplingDate`) | taxon age ranges in `readTaxonData()` (the density is zero outside each interval); `mvFossilTipTimeUniform` per fossil |
| Morphology | Lewis Mkv, six partitions by number of states (2–7), Γ4 shared, shape ~ U(0, 10), unordered | same (`fnJC(K)`, `coding = "variable"`, `fnDiscretizeGamma(α, α, 4)`) |
| Morphological clock | UCLN: mean ~ LogNormal(−5.5, 2); log-SD ~ Gamma(shape 0.5396, scale 0.3819) | same priors; non-centred multipliers, z ~ N(0, 1), m = exp(σz − σ²/2) |
| DNA | GTR + Γ4 per locus (5 loci, unlinked); r_AG = 1, other exchangeabilities U(0, 100); flat frequencies; shape ~ U(0, 10) | same |
| Locus rates | `mutationRate`s summing to 1, flat on the simplex | `locus_rates ~ dnDirichlet(1, 1, 1, 1, 1)` (same prior) |
| Molecular clock | UCLN: mean ~ LogNormal(−3.5, 1); log-SD ~ Gamma(0.5396, 0.3819) | same priors |
| MCMC | one chain, 10⁸ states | 4 chains per clock model; 309 (shared) or 422 (separate) proposals per iteration; sampled every 60; 25 % burn-in |

### Differences from the source analysis

1. **Clock linkage (shared runs only).** One set of branch-rate multipliers with one log-SD (prior Gamma(0.5396, scale 0.3819)) and two mean rates. The separate runs follow the source.
2. **Continuous UCLN.** BEAST2's `UCRelaxedClockModel` uses 110 discretised rate categories; RevBayes uses a continuous multiplier per branch, in non-centred form (same prior). Centred multipliers (round 1) mixed slowly in the log-SDs.
3. **One morphological character fewer (201 of 202).** RevBayes's `coding = "variable"` drops character 81 of `../penguins_morph.nex` (2-state partition: state 1 in every taxon except two `{01}` polymorphisms). BEAST2 kept it and divided its likelihood by P(variable).
4. **Gamma discretisation.** RevBayes uses the category means (`fnDiscretizeGamma` default); BEAST2's `SiteModel` uses category medians rescaled to mean 1.
5. **Operators.** The moves differ from BEAST2's operators (see `penguins_fbd.Rev`). There is no clock-up/tree-down move (BEAST2's `UpDownOperator`), because RevBayes's `mvUpDownScale` is biased with sampled ancestors (below).
6. **Starting trees.** Chains descend from the XML's starting tree (`start-tree.R`), not from random FBD trees. Round 4 started from the last state of the round-3 chains.
7. **Software.** RevBayes development a2a8e2db (above).

## RevBayes faults

Three faults affect FBD analyses with sampled ancestors and uncertain fossil ages. Rounds 1–3 were affected; round 4 avoids or tests each one.

1. **`mvUpDownScale` with a tree** (v1.4.1 and later). The Hastings ratio counts `getNumberOfInteriorNodes()`. `TreeUtilities::rescaleSubtree()`, however, also scales the root, and it skips every sampled-ancestor parent together with the subtree below it. The Jacobian is wrong whenever the tree has sampled ancestors (prior-only checks on the Hymenoptera baseline showed large distortions; `mvTreeScale` has the same mismatch). Rounds 1–2 used this move. It is off by default from round 3 (`use_up_down`), and is not fixed upstream.
2. **`mvCollapseExpandFossilBranch` irreversibility** (v1.4.1). The expand move drew the new parent age on (fossil's *maximum* range age, grandparent), but collapse acted from any parent age above the fossil's sampled age. Collapses from parent ages between the fossil's age and the top of its range could never be undone, so chains drifted towards sampled ancestors whenever fossil ages are uncertain (all 36 here). Fixed upstream after v1.4.1 (c389174b; also 94a3a7ed, and mvNarrow fixes 98424aa8 and 443ecdbc, all in a2a8e2db). Rounds 1–3 were affected.
3. **Fossil sampled ancestor at a stochastic root age** (found by the Hymenoptera baseline). With `dnFBDP(rootAge = <stochastic>)`, the collapse move can make a fossil the sampled ancestor at the root, and root-age moves can then drag it outside its age range. **The penguin analysis is not exposed.** `penguins_fbd.Rev` (the only penguin script with a tree prior) uses origin mode, `dnFBDP(originAge = origin_time, …)`. Root sampled ancestors do occur: a fossil is a zero-length child of the root in 41 % of round-4 trees (mostly *Waimanu*). But in a2a8e2db `mvRootTimeSlideUniform` refuses to move a root that is a sampled-ancestor parent, and `check-fossil-ages.py` found no fossil outside its age range in any of the 24 672 round-4 trees (`summary/fossil_age_check.txt`), nor in the 20 008 round-3 trees (`summary/fossil_age_check_round3.txt`). In origin mode the logged `num_samp_anc` includes these root sampled ancestors: it equals the number of zero-length fossil branches in every sampled tree checked (`check-sa-count.py`).

### Data-free validation of the fix

The decisive test of fault 2 samples the tree prior given only the tip counts and fossil age ranges (no characters). Gavryushkina et al. used the same quantity as the "prior" for their Bayes factors, and give it as 0.42 for both *Paraptenodytes antarcticus* and *Palaeospheniscus patagonicus*.

| Run | Sampled ancestors | s | P(SA) *P. antarcticus* | P(SA) *P. patagonicus* | Mean / max abs. difference in per-fossil P(SA) from BEAST2 |
|---|---|---|---|---|---|
| BEAST 2.4.8 + SA 1.1.7 (`sampleFromPrior`, 5×10⁷ states) | 16.6 [3.00, 29.0] | 0.254 [0.0025, 0.632] | 0.40 | 0.40 | — |
| RevBayes v1.4.1 (`ignoreAllData`, 2 × 200 000 iterations) | 32.0 [26.0, 36.0] | 0.735 [0.458, 0.982] | 0.81 | 0.86 | 0.428 / 0.472 |
| RevBayes development a2a8e2db (same) | 17.2 [4.00, 29.0] | 0.266 [0.0072, 0.629] | 0.41 | 0.41 | 0.017 / 0.039 |

v1.4.1 nearly doubles the number of sampled ancestors before any data enter. a2a8e2db matches BEAST2 within Monte Carlo error (per-fossil P(SA) has a Monte Carlo SE of about 0.02–0.03 at these ESSs).

### Effect on the posterior: round 3 (v1.4.1) against round 4 (a2a8e2db)

Both rounds use the same model and moves (no up-down move) and the same chain lengths. Round 4 started from round 3's final states.

| Clock | Quantity | Round 3 (v1.4.1) | Round 4 (a2a8e2db) |
|---|---|---|---|
| shared | crown age | 12.0 [10.4, 14.0] | 12.1 [10.4, 14.0] |
| shared | node 4 (*Aptenodytes* + *Pygoscelis*) | 6.83 [4.28, 9.40] | 6.62 [4.02, 9.16] |
| shared | sampled ancestors | 27.4 [23.0, 31.0] | 15.4 [9.00, 22.0] |
| shared | sampling proportion s | 0.504 [0.266, 0.739] | 0.192 [0.043, 0.371] |
| separate | crown age | 12.4 [10.2, 15.5] | 12.5 [10.3, 15.7] |
| separate | node 4 (*Aptenodytes* + *Pygoscelis*) | 9.46 [5.84, 13.9] | 9.27 [5.57, 14.2] |
| separate | sampled ancestors | 27.2 [22.0, 31.0] | 15.5 [8.00, 22.0] |
| separate | sampling proportion s | 0.507 [0.267, 0.760] | 0.207 [0.051, 0.398] |

Node ages hardly move, but the sampled-ancestor count and s are strongly biased in v1.4.1. Per-fossil sampled-ancestor and placement probabilities (`summary/round3-v1.4.1/fossils_*.csv`) are affected accordingly.

## BEAST2 reference

The source XML was run unchanged in BEAST 2.4.8 (bundled JRE) with SA 1.1.7 and MM 1.0.5 (`run-beast2-hamilton.sh`; md5 of `xmls.zip` as in `../raw/MANIFEST.csv`). There were 6 chains of 10⁸ states (22–28 h each on one core), with 10 % burn-in. One chain (seed 202) settled in the young-crown region described below (crown 7.2 Ma). It is excluded from the pooled reference values; the per-chain values are in `summary/beast2/perchain.txt`. The other five reproduce the paper:

- crown age 12.54 [10.02, 15.80] against 12.66 [9.91, 15.75];
- s 0.231 against 0.23;
- *S. muizoni* ancestry 0.616 against 0.61;
- Table 3 crown-placement probabilities to within 0.006.

The posterior number of sampled ancestors, which the paper does not report, is 16.7 [9, 23]. The molecular clock mixes slowly in BEAST2 too (`ucldMean.dna`: summed ESS 390, worst chain 44).

## MCMC history

| Round | RevBayes | Multipliers | Up-down move | Chains × iterations | Start | Status |
|---|---|---|---|---|---|---|
| — | 1.2.5 (Hamilton module) | centred | yes | 16 × ≈ 25 000 (cancelled) | random / XML tree | superseded (old, buggy build); summaries kept on Hamilton |
| 1 | v1.4.1 | centred | yes | 8 × 120 000 | XML tree | biased (faults 1, 2); slow log-SD mixing |
| 2 | v1.4.1 | non-centred | yes | 8 × 120 000 | round-1 final states | biased (faults 1, 2) |
| 3 | v1.4.1 | non-centred | no | 8 × 150 000 | round-2 final states | biased (fault 2); `summary/round3-v1.4.1/` |
| 4 | development a2a8e2db | non-centred | no | 4 × 150 000 (shared), 4 × 220 000 (separate) | round-3 final states | headline |

**Random starting trees get trapped.** Chains started from random FBD trees settled with a crown age of 3–7 Ma and the crown fossils (*Spheniscus urbinai*, *S. megaramphus*, *S. muizoni*, *Pygoscelis grandis*, *Madrynornis mirandus*) outside the crown. A fossil cannot enter a crown younger than itself, and the crown cannot grow older past the fossils on its stem. One of six BEAST2 chains started from the XML's tree fell into the same region, so this is a property of the posterior surface, not of RevBayes. Every RevBayes chain here therefore descends from the XML's starting tree. CET runs on this dataset should also start from a sensible tree and use several chains.

**Round 1 output defect.** Two round-1 output files ended early (a tree file at iteration 15 000, a rate-tree file at 15 300) for no identified I/O reason, while the chains ran on. `run-hamilton.sh`'s `START_FROM` therefore takes the last complete tree and the parameter values logged at the same iteration.
