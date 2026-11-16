# Penguins: conventional total-evidence tip dating in RevBayes

This is the comparator for the CET analyses of the penguin pilot dataset, and the source of their starting trees. It is a fossilized birth–death (FBD) analysis with sampled ancestors and an uncorrelated lognormal (UCLN) relaxed clock. It follows Gavryushkina et al. (2017, *Syst. Biol.* 66:57–73, doi:10.1093/sysbio/syw060), analysis "8+dna R" (`../raw/xmls.zip`, `xmls/8+dna_R.xml`), as closely as RevBayes allows.

Two clock models were run:

- **shared** (the EvoTempo comparator): one set of UCLN branch-rate multipliers scales both morphology and DNA, each with its own mean rate. This is the conventional analogue of CET, in which one tempo per lineage drives every character.
- **separate** (the source analysis): independent UCLN clocks for morphology and DNA, as in the BEAST2 XML.

**Software.** The headline results (round 4) come from RevBayes built from the `development` branch at commit `a2a8e2db8a1b57eee123c6b8e644335ccc268692` (29 Oct 2026, reports itself as 1.4.2-preview), not from the v1.4.1 release. v1.4.1 has two moves that are biased on trees with sampled ancestors, and rounds 1–3 used them; see [RevBayes faults](#revbayes-faults). v1.4.2 (tagged 6 Nov 2026) is a2a8e2db plus 101 commits. None of them touches the FBD distribution or the tree moves, but they do change some code this model uses: the CTMC likelihood cache (RevBayes PR #1136, `65f103b4`, `AbstractPhyloCTMCSiteHomogeneous.h`; described as removing over-invalidation, i.e. a speed fix), the argument cast in `Move_Slide.cpp` (used here for every `z_*` branch multiplier), and the default `alpha` types of `Move_BetaSimplex` and `Move_DirichletSimplex`. None of these is expected to change the posterior, but the analysis has not been rerun on 1.4.2. Later development builds (fd807d10, used for new runs from 11 Nov) fix `mvTreeScale` and `mvSubtreeScale` with sampled ancestors, neither of which is used here. As an independent reference, the source XML was also rerun in BEAST2 (see [BEAST2 reference](#beast2-reference)).

## Headline results (round 4)

Convergence: 4 chains per clock model, the first 25 % of each discarded. R̂ is the rank-normalised split R̂ (the larger of the bulk and folded-tail values) and ESS the rank-normalised bulk and tail ESS of Vehtari et al. (2021, *Bayesian Anal.* 16:667–718), as implemented in MkPrime (`MkPrime:::.Rhat`, `.EssBulk`, `.EssTail`; `rhat-ess-quartet.R`). The counts and minima below are over all logged continuous parameters except Posterior, Likelihood and Prior. `rhat_ess_<clock>.csv` lists these three as well; the separate model's Posterior has R̂ 1.010. ASDSF is over clades with frequency ≥ 0.1 in any chain. `convergence_<clock>.csv` also gives the classical per-chain ESS and Gelman–Rubin PSRF, as a secondary check. Versions: MkPrime 0.0.0.9000 as installed in the project library on Hamilton (built 16 Oct 2026; its source commit was not recorded, so `summary/mkprime-diagnostics-functions.R` records the deparsed definitions of `.Rhat`, `.EssBulk`, `.EssTail` and the 9 helpers they call); ConsTree 1.0.1 (ms609/ConsTree `6124244`, as in `../cet/README.md`); TreeTools 2.4.1.9000. Sample windows: `rhat-ess-quartet.R` and `perchain-summary.py` use every sample after burn-in (1 876 per shared chain, 2 751 per separate chain). `summarise.R` reads copies of the files without their last line, which can be incomplete in a running chain, so its tables use 1 875 and 2 750. The difference is one sample per chain.

| Clock | Chains | Iterations per chain | Samples per chain after burn-in | Max R̂ | Min bulk ESS | Min tail ESS | Min per-chain bulk ESS | ASDSF (max SDSF) | Run time per chain (h) |
|---|---|---|---|---|---|---|---|---|---|
| shared | 4 | 150,000 | 1,876 | 1.020 (`node_11`) | 404 (`node_11`) | 540 (`node_11`) | 37 (`node_4`) | 0.0098 (110 clades; max 0.034) | 55–67 |
| separate | 4 | 220,000 | 2,751 | 1.031 (`node_4`) | 243 (`node_4`) | 566 (`node_3`) | 58 (`node_4`) | 0.0162 (114 clades; max 0.143) | 62–72 + 22.4 / 28.8 (extension; Slurm elapsed) |

Every R̂ is below 1.05, and every bulk and tail ESS is above 240. Neither model meets the stricter R̂ < 1.01 of Vehtari et al. for every parameter. In the shared model, 7 of 86 parameters exceed 1.01, all of them extant node ages (nodes 4–7, 11, 13, 14), and every bulk ESS is above 400. In the separate model, 12 of 87 exceed 1.01: the crown age, nine extant node ages (nodes 3–5, 7–9, 11–13), the molecular clock mean and its log-SD. Nine of these have a bulk ESS between 243 and 400. The separate-clock chains were therefore extended from 150 000 to 220 000 iterations from their checkpoints. The extension jobs (Slurm 19425609, 19425610, 19425612, 19425613) were cancelled after 22 h 24 min (runs 13–15) and 28 h 48 min (run 16), so their `time_*_resume*.txt` files are empty and are not committed. All three output files of each chain were then cut at iteration 220 000 with `truncate-chains.sh` (end iterations in `summary/trunc_log.txt`; the last kept sample is 219 960). The slowest feature is the basal split of the crown. In the separate model, the clade of all extant penguins except *Eudyptula* and *Spheniscus*, with their fossil relatives, has frequency 0.34–0.68 among chains (max SDSF 0.143). The node ages at the base of the crown (nodes 3–5) are also the slowest-mixing ages. The per-chain bulk ESS is a split-chain value, so it also penalises drift between a chain's two halves. Its minimum in the shared model (37, node 4 in one chain) reflects such slow drift. Node 4's pooled bulk ESS is 584 and its R̂ is 1.010.

### Comparison with the published analysis and with the BEAST2 rerun

The published values are from the arXiv version of the paper and its SM (arXiv:1506.04797; SM Table 4, Fig. 5 node numbers, main text). "BEAST2 rerun" is the source XML run again here (5 chains; see [BEAST2 reference](#beast2-reference)). The table shows mean and 95 % HPD interval. Ages are in Ma and are common-ancestor ages, as in the paper: the age of the MRCA of the listed extant taxa in every tree, whether or not the group is monophyletic in that tree. `crown_age` and `node_*` in `penguins_fbd.Rev` are the same quantity.

| Node | Published | BEAST2 rerun | RevBayes shared | RevBayes separate |
|---|---|---|---|---|
| Crown age (node 1) | 12.7 [9.91, 15.8] | 12.5 [10.0, 15.8] | 12.1 [10.4, 14.0] | 12.5 [10.3, 15.7] † |
| Node 2: Megadyptes + Eudyptes + Eudyptula + Spheniscus | 11.2 [7.92, 13.0] | 11.2 [8.34, 13.5] | 12.1 [10.4, 14.1] | 11.5 [9.84, 13.5] |
| Node 3: Eudyptula + Spheniscus | 10.3 [7.14, 12.0] | 10.4 [7.32, 12.9] | 11.9 [9.89, 14.1] | 10.8 [8.69, 13.2] † |
| Node 4: Aptenodytes + Pygoscelis | 9.83 [5.84, 14.0] | 9.59 [5.91, 14.0] | 6.62 [4.02, 9.16] | 9.27 [5.57, 14.2] † |
| Node 5: Pygoscelis | 6.12 [3.30, 9.04] | 5.96 [3.22, 8.94] | 3.74 [1.81, 5.70] | 5.73 [3.11, 9.17] † |
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

1. **Node ages and diversification parameters agree; the sampled-ancestor count, s and per-fossil P(SA) do not.** In the separate model, the crown age, the extant node ages, d, ν and the origin agree with both the published values and the BEAST2 rerun within their intervals. The number of sampled ancestors and s are systematically lower than in BEAST2, beyond Monte Carlo error: the per-chain means do not overlap (SA count 15.35–15.68 in the four RevBayes separate chains against 16.57–16.71 in the five BEAST2 chains; s 0.203–0.209 against 0.229–0.232; `summary/perchain_round4.csv`, `summary/beast2/perchain/`). Per fossil, P(SA) differs by more than 0.05 for 18 of 36 fossils, and the RevBayes and BEAST2 per-chain ranges do not overlap for 30 of 36. 31 of 36 are lower in RevBayes (by up to 0.11, *Perudyptes devriesi*), as a lower s would make them; two are much higher: *Pachydyptes ponderosus* (0.51 against 0.22) and *Madrynornis mirandus* (0.68 against 0.45). The crown placements of *Madrynornis* and *S. muizoni* are also higher (0.98–0.99 against 0.89–0.92) and that of *Marplesornis novaezealandiae* lower (0.03–0.10 against 0.17). The data-free runs agree (below), so the tree prior, conditioning, fossil-age treatment and FBD moves are not the cause; the difference arises only with data. Candidate causes are the documented model differences (continuous UCLN against BEAST2's 108 rate categories, which lets morphological rates vary more; character 81; Gamma category means against medians) or an unidentified RevBayes behaviour with data. The cause is open: issue #109.
2. **Linking the clocks (shared) changes the extant topology at the base of the crown.** The published basal split is *Aptenodytes* + *Pygoscelis* against the rest. With linked clocks, the common ancestor of node 2 (*Megadyptes*, *Eudyptes*, *Eudyptula*, *Spheniscus*) is the crown root in 0.985 of samples (per chain 0.978–0.988). This contrast is robust. The separate model's pooled 0.563 is not distinguishable from BEAST2's 0.404: its per-chain values are 0.362, 0.561, 0.700 and 0.628, and BEAST2's are 0.355–0.475 over the five retained chains (`summary/beast2/beast2_basal_split.txt`, `summary/beast2/perchain/`). With linked clocks, nodes 4 and 5 (*Aptenodytes* + *Pygoscelis*, *Pygoscelis*) are 2.7 and 2.0 Myr younger than in the separate model, and the molecular clock mean is 70 % faster. The crown age is unchanged (12.1 against 12.5 Ma).
3. **The sampling proportion.** The earlier gap (s ≈ 0.50 in rounds 1–3) came from the v1.4.1 bias. After the fix, s is 0.19 (shared) and 0.21 (separate), against 0.23 published and in the BEAST2 rerun; the definitions agree (below), and the remaining difference is the one in item 1.
4. **Both posteriors are conditional on the old-crown mode (separate model).** See [Two crown-age modes](#two-crown-age-modes). Under separate clocks a young-crown mode (crown ≈ 5–6 Ma, four of the five crown fossils on the stem) has the same posterior density as the old-crown mode, in BEAST2 and in RevBayes. Every separate-clock value above, the BEAST2 rerun and the published analysis describe the old-crown mode only. Under the shared clock the young region is transient and the shared results are not materially affected.

**Provisional values (†).** In the separate model, the basal split of the crown, the ages of nodes 3–5, the crown age and the molecular clock have R̂ 1.017–1.031 (crown_age 1.017, node_3 1.020, node_4 1.031, node_5 1.029, clock_dna 1.022, ucld_sd_dna 1.026; bulk ESS 243–430). They are marked † in the tables above and below. ucld_sd_dna also differs between chains (bimodal posterior, modes near 0.4 and 0).

**Summary trees.** `summary/quartet_consensus_<clock>.nwk` is the quartet consensus (`ConsTree::Quartet()`, default settings) of 4000 pooled post-burn-in samples, thinned evenly. It resolves 42 (shared) and 41 (separate) of the 52 possible splits. Its node labels are the splits' frequencies in all pooled samples, and `split_support_<clock>.csv` gives them per chain. The consensus is unrooted and carries no node ages; its top level is an unrooted trifurcation, and its least-supported retained splits have frequencies 0.287 and 0.271 (separate), so it is a summary, not a resolved estimate of the basal relationships. Node ages are the common-ancestor ages in the tables above (`clade_ages_<clock>.csv`). No MCC tree is used.

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
| Molecular clock mean (substitutions / site / Myr, summed over the five loci) | ≈ 0.0436² | 0.035 [0.023, 0.048] | 0.062 [0.042, 0.087] | 0.036 [0.023, 0.051] † |
| Branch-length-weighted mean morphological rate | — | 0.034 [0.020, 0.051] | 0.025 [0.015, 0.036] | 0.031 [0.017, 0.048] |
| Branch-length-weighted mean molecular rate | — | 0.034 [0.023, 0.045] | 0.067 [0.045, 0.093] | 0.035 [0.023, 0.047] † |
| UCLN log-SD, morphology | ≈ 0.728² | 0.916 [0.582, 1.27] | 0.715 [0.466, 0.971] (shared) | 1.09 [0.717, 1.48] |
| UCLN log-SD, DNA | bimodal, modes ≈ 0.4 and 0 | 0.278 [3.3e-09, 0.521] | (shared) | 0.312 [2.7e-07, 0.552] † |
| Coefficient of variation of morphological branch rates³ | 1.15 | 1.13 [0.625, 1.71] | 0.875 [0.488, 1.35] | 1.53 [0.796, 2.38] |

² Simulation settings, not posterior estimates.
³ Unweighted SD/mean over branches, as BEAST2's `RateStatistic`.

In the separate model, morphological rates vary more among branches than in the BEAST2 rerun: log-SD 1.09 against 0.92, CV 1.53 against 1.13. A candidate cause is the continuous UCLN: BEAST2's 108 discretised rate categories (the XML declares 110; the state file holds 108, `rateCategories[108 1]`) cap the most extreme branch rates. Linking the clocks (shared) pulls the morphological rate variation down to that of the combined data (log-SD 0.72).

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
| `perchain-summary.py` | Per-chain sampled-ancestor count, s, basal-split frequency and per-fossil P(SA) from the archived round-4 chains (`summary/perchain_round4.csv`) |
| `truncate-chains.sh` | The command that cut the extended separate-clock chains at iteration 220 000 |
| `check-sa-count.py` | Checks that the logged `num_samp_anc` equals the number of zero-length fossil branches in each sampled tree |
| `run-beast2-hamilton.sh` | BEAST2 reference run of the source XML (posterior or `PRIOR=true`) |
| `beast-posterior-summary.py`, `beast-sa-summary.py` | Summaries of the BEAST2 posterior and prior runs |
| `published-gavryushkina2017.csv` | Published estimates compared against |
| `summary/` | Round-4 results: `rhat_ess_<clock>.csv` (posterior summaries with R̂ and ESS), `quartet_consensus_<clock>.nwk` (quartet consensus, unrooted, node labels = split frequencies), `split_support_<clock>.csv` (its splits per chain), `parameters_`, `clade_ages_`, `clades_`, `convergence_`, `fossils_`, `derived_`, `sa_vs_s_<clock>.csv`, `headline_<clock>.txt`, `fossil_age_check.txt`, `trunc_log.txt` |
| `summary/beast2/` | BEAST2 rerun: posterior summaries (5 chains), per-chain summaries for post1–6 (`perchain/`), the basal-split frequencies (`beast2_basal_split.txt`), the post2 regime table (`beast2_crown_mode_regimes.csv`), the data-free run (`beast2_prior_summary.txt`) and the md5 manifest of the BEAST2 archive |
| `summary/young-crown/` | RevBayes young-crown test: start states with md5s, chain and summary scripts, and summaries (`agg_post25.txt`, `blocks5000.txt`) |
| `summary/perchain_round4.csv`, `summary/mkprime-diagnostics-functions.R` | Per-chain round-4 summaries; the deparsed MkPrime diagnostic functions used |
| `summary/data-free/` | Data-free runs of RevBayes v1.4.1 and of a2a8e2db |
| `summary/penguins_r4_<clock>_md5.txt` | md5 of every file in the round-4 chain archive (below) |
| `summary/round3-v1.4.1/` | Round 3 (v1.4.1, biased), kept as evidence of the bias |

`data/`, the trace files and the tree files are not committed (25–40 MB per tree file). `prepare-data.R` and `start-tree.R` regenerate `data/`.

**Round-4 chain archive.** The round-4 parameter logs, tree files and rate-tree files of all eight chains (shared 150 000 iterations; separate truncated at 220 000) are archived as two xz tarballs with md5 manifests: `penguins_r4_shared_chains.tar.xz` and `penguins_r4_separate_chains.tar.xz` (28.3 MB and 41.2 MB; SHA-256 `3b775337…484062d` and `39d7632f…c12125b`; md5 of every member file in `summary/penguins_r4_<clock>_md5.txt`). The canonical copy is the EvoTempo project's artifact store (Claude Science project `proj_07b3619548a5`, artifacts `aadac10f-152b-43de-a186-b0d15d7a3b9b` (shared) and `55110c6f-814d-4c45-88a0-63bdf1bf80f1` (separate)); a working copy stays on Hamilton in `/nobackup/pjjg18/evotempo-penguins-summaries/round4-chains-archive/`, which is not backed up. The CET starting trees (`../cet/data/start-tree-selection.txt`) were drawn from the shared-clock chains in this archive. **BEAST2 archive.** The parameter logs and tree files of all six BEAST2 posterior chains (post1–6, including the young-crown chain post2) and of the data-free chain are archived as `penguins_beast2_chains.tar.xz` (90.2 MB; SHA-256 `73f953b8af079863654e6a69cc06076f906b006ddb4d18c06ddd738dfa18ac37`; md5 of all 49 member files in `summary/beast2/penguins_beast2_chains_md5.txt`). The canonical copy is the project's artifact store (artifact `e5e18b72-77ea-42c9-ac35-6cc79f55f7c3`); a working copy stays on Hamilton in `/nobackup/pjjg18/evotempo-penguins-beast2/`. Per-chain BEAST2 summaries are committed in `summary/beast2/perchain/`.

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
2. **Continuous UCLN.** BEAST2's `UCRelaxedClockModel` uses 108 discretised rate categories (the XML declares 110; its state file holds `rateCategories[108 1]`); RevBayes uses a continuous multiplier per branch, in non-centred form (same prior). Centred multipliers (round 1) mixed slowly in the log-SDs.
3. **One morphological character fewer (201 of 202).** RevBayes's `coding = "variable"` drops character 81 of `../penguins_morph.nex` (2-state partition: "1" in 20 taxa, `{01}` in 2, "?" in 31 and "-" in 2, so every scored cell is compatible with state 1 and the character counts as constant). BEAST2 kept it and divided its likelihood by P(variable).
4. **Gamma discretisation.** RevBayes uses the category means (`fnDiscretizeGamma` default); BEAST2's `SiteModel` uses category medians rescaled to mean 1.
5. **Operators.** The moves differ from BEAST2's operators (see `penguins_fbd.Rev`). There is no clock-up/tree-down move (BEAST2's `UpDownOperator`), because RevBayes's `mvUpDownScale` is biased with sampled ancestors (below).
6. **Starting trees.** Chains descend from the XML's starting tree (`start-tree.R`), not from random FBD trees. Round 4 started from the last state of the round-3 chains. `START_FROM` restores the tree, the FBD parameters, the clock means and log-SDs, the substitution-model parameters and the Gamma shapes, but not `locus_rates`, the base frequencies `pi` or the non-centred branch multipliers `z_*`, which restart from their initial values. The chains therefore re-equilibrate at the start of each round: in round 4 the sampled-ancestor count fell 24 → 16 → 12 → 9 over the first 2 400 iterations, well inside the 25 % burn-in.
7. **Software.** RevBayes development a2a8e2db (above).

## RevBayes faults

Three faults affect FBD analyses with sampled ancestors and uncertain fossil ages. Rounds 1–3 were affected; round 4 avoids or tests each one.

1. **`mvUpDownScale` with a tree** (v1.4.1 and later). The Hastings ratio counts `getNumberOfInteriorNodes()`. `TreeUtilities::rescaleSubtree()`, however, also scales the root, and it skips every sampled-ancestor parent together with the subtree below it. The Jacobian is wrong whenever the tree has sampled ancestors (prior-only checks on the Hymenoptera baseline showed large distortions; `mvTreeScale` has the same mismatch). Rounds 1–2 used this move. It is off by default from round 3 (`use_up_down`), and is not fixed upstream.
2. **`mvCollapseExpandFossilBranch` irreversibility** (v1.4.1). The expand move drew the new parent age on (fossil's *maximum* range age, grandparent), but collapse acted from any parent age above the fossil's sampled age. Collapses from parent ages between the fossil's age and the top of its range could never be undone, so chains drifted towards sampled ancestors whenever fossil ages are uncertain (all 36 here). Fixed upstream after v1.4.1 (c389174b; also 94a3a7ed, and mvNarrow fixes 98424aa8 and 443ecdbc, all in a2a8e2db). Rounds 1–3 were affected.
3. **Fossil sampled ancestor at a stochastic root age** (found by the Hymenoptera baseline). With `dnFBDP(rootAge = <stochastic>)`, the collapse move can make a fossil the sampled ancestor at the root, and root-age moves can then drag it outside its age range. **The penguin analysis is not exposed.** `penguins_fbd.Rev` (the only penguin script with a tree prior) uses origin mode, `dnFBDP(originAge = origin_time, …)`. Root sampled ancestors do occur: a fossil is a zero-length child of the root in 41 % of round-4 trees (10 118 of 24 672, burn-in included; mostly *Waimanu*). But in a2a8e2db `mvRootTimeSlideUniform` refuses to move a root that is a sampled-ancestor parent, and `check-fossil-ages.py` found no fossil outside its age range in any of the 24 672 round-4 trees (`summary/fossil_age_check.txt`), nor in the 20 008 round-3 trees (`summary/fossil_age_check_round3.txt`). In origin mode the logged `num_samp_anc` includes these root sampled ancestors: it equals the number of zero-length fossil branches in every sampled tree checked (`check-sa-count.py`) but one: in shared_run15 at iteration 132 480 the tree file shows 12 zero-length fossil tips against 11 logged sampled ancestors, because the tree files round branch lengths to 6 decimals and one very short branch prints as zero.

### Data-free validation of the fix

The decisive test of fault 2 samples the tree prior given only the tip counts and fossil age ranges (no characters). Gavryushkina et al. used the same quantity as the "prior" for their Bayes factors, and give it as 0.42 for both *Paraptenodytes antarcticus* and *Palaeospheniscus patagonicus*.

| Run | Sampled ancestors | s | P(SA) *P. antarcticus* | P(SA) *P. patagonicus* | Mean / max abs. difference in per-fossil P(SA) from BEAST2 |
|---|---|---|---|---|---|
| BEAST 2.4.8 + SA 1.1.7 (`sampleFromPrior`, one chain of 5×10⁷ states) | 16.6 [3.00, 29.0] (batch SE 0.43) | 0.254 [0.0025, 0.632] (batch SE 0.012) | 0.40 | 0.40 | — |
| RevBayes v1.4.1 (`ignoreAllData`, 2 × 200 000 iterations) | 32.0 [26.0, 36.0] | 0.735 [0.458, 0.982] | 0.81 | 0.86 | 0.428 / 0.472 |
| RevBayes development a2a8e2db (same) | 17.2 [4.00, 29.0] | 0.266 [0.0072, 0.629] | 0.41 | 0.41 | 0.017 / 0.039 |

v1.4.1 nearly doubles the number of sampled ancestors before any data enter. The BEAST2 reference is a single chain, so its own Monte Carlo error (batch-means SE 0.43 in the sampled-ancestor count and 0.012 in s) limits this comparison. a2a8e2db matches BEAST2 within Monte Carlo error (per-fossil P(SA) has a Monte Carlo SE of about 0.02–0.03 at these ESSs).

### Effect on the posterior: round 3 (v1.4.1) against round 4 (a2a8e2db)

Both rounds use the same model and moves (no up-down move). Round 3 ran 150 000 iterations for both clock models; round 4 ran 150 000 (shared) and 220 000 (separate). Round 4 started from round 3's final states.

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

The source XML was run unchanged in BEAST 2.4.8 (bundled JRE) with SA 1.1.7 and MM 1.0.5 (`run-beast2-hamilton.sh`; md5 of `xmls.zip` as in `../raw/MANIFEST.csv`). There were 6 chains of 10⁸ states (22–28 h each on one core), with 10 % burn-in. One chain (post2, seed 202) moved into the young-crown mode described in [Two crown-age modes](#two-crown-age-modes) at about 2.6×10⁷ states and stayed there (its whole-run crown age averages 7.2 Ma). It is excluded from the pooled reference values, which therefore describe the old-crown mode only; the per-chain values are in `summary/beast2/perchain.txt` and `summary/beast2/perchain/`. The other five reproduce the paper:

- crown age 12.54 [10.02, 15.80] against 12.66 [9.91, 15.75];
- s 0.231 against 0.23;
- *S. muizoni* ancestry 0.616 against 0.61;
- Table 3 crown-placement probabilities to within 0.0065 (largest difference, *Marplesornis*).

The posterior number of sampled ancestors, which the paper does not report, is 16.7 [9, 23]. The molecular clock mixes slowly in BEAST2 too (`ucldMean.dna`: summed ESS 390, worst chain 44).

## Two crown-age modes

Random FBD starting trees often settled in a region with a crown age of 3–7 Ma and four of the five crown fossils (*Spheniscus urbinai*, *S. megaramphus*, *S. muizoni*, *Pygoscelis grandis*, *Madrynornis mirandus*) on the stem. An earlier version of this README called this a trap. It is better described as a second posterior mode, at least under separate clocks (#99 review, R99-2; issue #108).

**BEAST2.** Chain post2 (seed 202) sat in the old-crown mode until about 2.6×10⁷ states, then moved to the young-crown mode and stayed there for the remaining 7.3×10⁷ states (`summary/beast2/beast2_crown_mode_regimes.csv`; per-sample trace regenerated from the archived chain):

| Regime | Samples | Mean log posterior | Max log posterior | Log likelihood (morphology) | Tree prior (birth–death) | Crown age (Ma) | Crown fossils in crown | Sampled ancestors | s |
|---|---|---|---|---|---|---|---|---|---|
| post1, post3–6 (states ≥ 10⁷) | 9001 each | −19556.4 to −19559.5 | −19497.6 to −19509.2 | −1699.5 to −1700.5 | −191.1 to −193.5 | 12.3–12.7 | 4.76–4.87 | 16.6–16.7 | 0.229–0.232 |
| post2, old (1–2.6×10⁷) | 1600 | −19556.7 | −19512.1 | −1699.4 | −191.1 | 12.2 | 4.79 | 16.8 | 0.232 |
| post2, young (≥ 3×10⁷) | 7001 | −19555.4 | −19503.9 | −1717.5 | −170.4 | 5.96 | 0.96 | 15.5 | 0.183 |

The young mode fits the morphology about 18 log units worse and the tree prior about 21 better, so its posterior density is the same as the old mode's. Posterior values are comparable across these chains (same XML, same data, identical constant terms). With one switch in about 5.3×10⁸ old-mode states and none back in 7.3×10⁷ young-mode states, a flux balance puts at least 0.2–4 % of the posterior mass in the young mode, with no upper bound.

**RevBayes test** (fd807d10; 50 000 iterations per chain, 25 % burn-in; Slurm 19449375–82; start states from post2's young regime in `summary/young-crown/starts/`, with md5s; summary `summary/young-crown/agg_post25.txt`):

| Clock | Start | Chains | Mean log posterior | Mean log likelihood | Crown age (Ma) | Crown fossils in crown | Fraction of samples young |
|---|---|---|---|---|---|---|---|
| separate | young crown | 2 | −20011.9, −20010.3 | −19245.5, −19245.3 | 4.89, 4.60 | 0.81, 0.85 | 1.000, 1.000 |
| separate | old crown (control) | 1 | −20013.8 | −19222.6 | 12.25 | 4.97 | 0.002 |
| shared | young crown | 4 | −19877.5 to −19879.2 | −19246.0 to −19247.0 | 11.98–12.34 | 4.94–4.98 | 0.001–0.005 |
| shared | old crown (control) | 1 | −19878.3 | −19247.0 | 12.00 | 4.91 | 0.009 |

- **Separate clocks.** Chains started young stay young for all 50 000 iterations, at the same posterior density as the old-crown control. RevBayes therefore reproduces the BEAST2 mode independently. Every separate-clock value in this README, the BEAST2 reference values and the published analysis describe the old-crown mode only. The mode masses are not yet estimated (#108). The young mode has 13.1–13.6 sampled ancestors, so including it would move RevBayes further from BEAST2 in the sampled-ancestor count (#109).
- **Shared clock (the CET comparator).** All four young starts returned to the old crown within 10 000 iterations, and old-crown chains visit the young region only briefly (fraction of samples young ≤ 0.009). The shared-clock results, and the CET starting trees drawn from them, are not materially affected.

Every RevBayes chain in rounds 1–4 descends from the XML's starting tree. CET runs on this dataset should start from trees of the shared-clock posterior and use several chains; whether CET itself has a young-crown mode is checked separately (#108).

## MCMC history

| Round | RevBayes | Multipliers | Up-down move | Chains × iterations | Start | Status |
|---|---|---|---|---|---|---|
| — | 1.2.5 (Hamilton module) | centred | yes | 16 × ≈ 25 000 (cancelled) | random / XML tree | superseded (old, buggy build); summaries kept on Hamilton |
| 1 | v1.4.1 | centred | yes | 8 × 120 000 | XML tree | biased (faults 1, 2); slow log-SD mixing |
| 2 | v1.4.1 | non-centred | yes | 8 × 120 000 | round-1 final states | biased (faults 1, 2) |
| 3 | v1.4.1 | non-centred | no | 8 × 150 000 | round-2 final states | biased (fault 2); `summary/round3-v1.4.1/` |
| 4 | development a2a8e2db | non-centred | no | 4 × 150 000 (shared), 4 × 220 000 (separate) | round-3 final states | headline |

**Round 1 output defect.** Two round-1 output files ended early (a tree file at iteration 15 000, a rate-tree file at 15 300) for no identified I/O reason, while the chains ran on. `run-hamilton.sh`'s `START_FROM` therefore takes the last complete tree and the parameter values logged at the same iteration.
