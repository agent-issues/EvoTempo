# Gamma rates among characters and per-mask Mkv (#97): validation

`MkCharacters()` (R/Characters.R, src/characters.cpp) adds to the CET character model:

- a number of states per character;
- ordered characters (RevBayes `fnOrderedRateMatrix(maxState = K)`);
- ambiguity sets;
- the Mkv correction per missing-data mask (RevBayes `dnPhyloCTMC(coding = "variable")`; #47, #65);
- discrete Gamma(alpha, alpha) rates among characters (`GammaRates()`: the mean rates of K equiprobable categories, computed exactly with `qgamma()` and `pgamma()`; RevBayes's `fnDiscretizeGamma(alpha, alpha, K)` approximates the same discretisation).

In `TempoMCMC()` the shape `alpha` is `start$alpha`. It is free when `charPriors$alpha` gives its prior, and then has its own move (`.AlphaMove()`).

**Mixture and correction.** With rate categories c (weights w_c = 1/K), a character's likelihood is sum_c w_c L_c. The correction for mask m is -n_m log(1 - sum_c w_c P_c,m), where P_c,m is the probability that every scored tip of m is in one state. This is RevBayes's `PhyloCTMCSiteHomogeneousConditional::sumRootLikelihood()`: the correction is computed within each category and combined over categories before the log is taken. As in RevBayes, a missing cell and a set of all K states are gaps of the mask. A character whose scored cells all hold the same state or set is an error (RevBayes excludes it with a warning).

Scripts (from the package root, with EvoTempo installed; Hamilton driver `data-raw/gamma-acrv-hamilton.sh`):

| script | check |
|---|---|
| `data-raw/gamma-acrv-identity.R` | Gamma off: seeded chains on a character matrix, main vs this branch |
| `data-raw/gamma-acrv-revbayes.R`, `.Rev` | 1. fixed-tree lnL vs RevBayes, 3 rates x 3 shapes (+ no Gamma), ordered characters, ambiguity, both engines; and on the Hymenoptera matrix |
| `data-raw/gamma-acrv-alpha-oracle.R` | 2. alpha move vs the exact 1-D posterior on a fixed tree with fixed paths |
| `data-raw/gamma-acrv-sbc.R` | 3. SBC-lite, cells (beta, gamma) = (1, 0) and (1, 1) |
| `data-raw/gamma-acrv-pt-check.R` | 4. transition probabilities at the empirical settings (trilobites, Hymenoptera) |
| `tests/testthat/test-MkCharacters.R` | unit tests (below) |

## Results

### Gamma off: chains unchanged

Characters given as a matrix take the code path of main. With rate categories, `TempoMCMC()` passes `alpha` to `MkLikelihood()`, so pipelines that rebind `MkLikelihood()` with the old signature are unaffected. Without them, no random numbers are drawn and no columns are added.

The check is `gamma-acrv-identity.R`. Four seeded chains (C++ and R engine x gamma 0 and 1; 20 burn-in + 30 sampled iterations with every move: `tree`, `node`, `tip`, `branch`, `joint`, `shift`, `nodeAge`, `fossilAge`, `slide`, `saJump`) on the test fixture (sampled ancestor, fossil age ranges). They ran under main fb24ed6 and under this branch (Hamilton, Slurm 19431188). `samples`, `occupancy`, `acceptance`, `burninAcceptance` and `moves` are `identical()` in all four (`identity.txt`). The check was repeated after the review fixes, at dd15973 against the same main chains (Slurm 19467760): all four identical (`identity-dd15973.txt`).

### Unit tests (`test-MkCharacters.R`)

- `GammaRates()` against Yang (1994) Table 1 (alpha 0.5) to 1e-7, with mean 1; errors for an invalid `nCat`.
- Parsing of states, multi-digit states, sets (`01`, `0/12`) and unknowns. A full set equals a gap. Input errors. Tips are matched to the tree by row name.
- Transition probabilities. Ordered characters use uniformization (P(0) = I exactly). Unordered characters use the closed form of the JC model, as in MkPrime. The tests compare both with a direct Poisson sum for K = 3 to 10 and t = 1e-300 to 100, and check that a sampled ancestor (zero-length branch) gives the same value as its parent node.
- `MkCharacters()` with one state count, no correction and no Gamma equals `MkLikelihood()` on the matrix to 1e-12, in both engines.
- An independent reference likelihood (`RefLogLik()`: its own pruning, mixture and correction) agrees with both engines on a 30-tip random tree for K = 3, 5, 7, 10, alpha 1e-3, 1 and 1e3, and rates 1e-6 to 10. The C++ engine is compared only where every site lnL is above -690, because it does not rescale (Limitations).
- An ordered K = 5 character with a sampled ancestor, at rates 1e-4 to 1.
- C++ and R engines agree to 1e-10 (relative) for coding all and variable, 1 or 4 categories, and alpha 0.2, 1.1 and 8.
- The Mkv correction: over every variable pattern of a mask (binary with all or 4 tips scored, ordered K = 3, unordered K = 3; with Gamma), the conditional probabilities sum to 1 (to 1e-12). The value equals RevBayes's mixture formula log mean_c L_c - log(1 - mean_c P_c), built from Gamma-free pieces, to 1e-12. A data set whose constant patterns have probability 1 gives -Inf, not NaN.
- `TempoMCMC()` with alpha free: identical seeded chains in both engines; the `alpha` column and acceptance; `logPost` equals the recomputed augmented density + character lnL + priors + alpha prior; fixed alpha; `updateParameters = FALSE`; an error when the starting character lnL is not finite; errors for an invalid `charPriors$alpha`.

At a277d3a (Hamilton, r/4.5.1, Slurm 19467785; fresh clone): `roxygenise()` leaves no uncommitted change. The C++ engine gives 115 test blocks and 5193 expectations, 0 failed and 3 skipped (`skip_on_cran()` oracles); `test-MkCharacters.R` has 13 blocks and 571 expectations. The R engine gives the same counts, and its only failures are the 6 in the 2 stiff-solver tests of `test-Paths.R`. These 6 also fail on main fb24ed6 (Slurm 19431685); #45 addresses them. `R CMD check --no-manual --ignore-vignettes`: `Status: OK`, with no NOTE or WARNING. (Slurm marks the job FAILED. The cause is the last line of its script, a `grep` for WARNING or NOTE, which exits with status 1 when it finds none.)

### 1. Fixed-tree log-likelihood vs RevBayes

`gamma-acrv-revbayes.R prepare` simulates the data. The tree is a coalescent time tree of 24 extant tips with root age 1; 4 tips are cut back to fossils, and one sampled ancestor sits on an internal branch. There are 74 characters: unordered K = 2, 3, 4 (30, 12, 8) and ordered K = 3, 4, 5 (10, 8, 6). Of the cells, 247 are missing (`?` or `-`) and 129 are ambiguity sets, including full binary `{01}`. RevBayes is the development build fd807d10 (`1.5.0-preview`) running `gamma-acrv-revbayes.Rev`: `dnPhyloCTMC(Q = fnJC(K) | fnOrderedRateMatrix(maxState = K), branchRates = rate, siteRates = fnDiscretizeGamma(a, a, 4), type = "Standard", coding = "variable")`, plus a run without `siteRates`. CET is `MkLikelihood(MkCharacters(..., coding = "variable", nRateCat = 4), lengths = rate x duration, alpha = a)` in both engines. The table gives totals over the 6 partitions, and maxima of |difference| per partition (`revbayes.csv`, `revbayes-compare.txt`, `revbayes-output.txt`). The values were computed on Hamilton (Slurm 19448619, branch e743555). They were recomputed after the review fixes (Slurm 19467761, dd15973), with the same values and with RevBayes's rates printed.

| rate | alpha | RevBayes total lnL | CET (C++) - RevBayes, total | max per partition, C++ | max per partition, R | max C++ - R |
|---|---|---|---|---|---|---|
| 0.1 | 0.3 | -569.0939762 | 1.15e-09 | 5.02e-10 | 5.00e-10 | 1.36e-12 |
| 0.1 | 1.1 | -604.0945475 | -6.86e-07 | 2.93e-07 | 2.93e-07 | 1.08e-12 |
| 0.1 | 5 | -633.5497585 | -1.55e-09 | 7.21e-10 | 7.22e-10 | 1.36e-12 |
| 0.1 | none | -649.7540704 | 3.41e-13 | 2.34e-13 | 1.42e-12 | 1.42e-12 |
| 0.5 | 0.3 | -538.0519635 | 1.00e-09 | 2.90e-10 | 2.90e-10 | 7.11e-13 |
| 0.5 | 1.1 | -551.9600503 | -1.05e-07 | 5.08e-08 | 5.08e-08 | 4.55e-13 |
| 0.5 | 5 | -566.6901362 | -1.73e-09 | 3.95e-10 | 3.95e-10 | 3.98e-13 |
| 0.5 | none | -574.7532632 | 1.14e-13 | 4.26e-14 | 3.98e-13 | 4.26e-13 |
| 2 | 0.3 | -564.3943075 | 1.15e-10 | 9.58e-10 | 9.58e-10 | 5.97e-13 |
| 2 | 1.1 | -580.9501423 | -1.54e-07 | 4.60e-08 | 4.60e-08 | 4.26e-13 |
| 2 | 5 | -617.8112790 | -9.10e-09 | 2.19e-09 | 2.19e-09 | 3.41e-13 |
| 2 | none | -647.8422985 | 0.00e+00 | 5.68e-14 | 3.13e-13 | 3.41e-13 |

The largest relative difference is 1.5e-9, against a tolerance of 1e-6. The differences come from the discretised rates. `GammaRates()` computes the category means exactly with `qgamma()` and `pgamma()`. RevBayes's `fnDiscretizeGamma()` takes its quantiles from the AS 91 chi-square approximation, whose tolerance is 5e-7. The printed RevBayes rates (14 significant digits) differ from `GammaRates()` by up to 2.2e-8 at alpha 1.1 (for example 1.01477710109970 against 1.014777123099664), by at most 1.2e-10 at alpha 0.3 and 5. This is the source of the 2.9e-7 per-partition differences at alpha 1.1. When RevBayes's printed rates are given to CET in place of `GammaRates()`, every total agrees with RevBayes to 4.55e-13, in both engines. The two engines agree to 8.0e-16 (relative).

**Hymenoptera matrix** (supplementary). Data are from the CET pipeline at e0c1384 (`hymenoptera-sha.txt`): 113 tips, 343 characters in 9 (ordering, K) partitions, ordered characters and ambiguity, on its start tree. The oracle is its RevBayes `check-mkv.Rev` output (fd807d10), at clock rates 0.0007, 0.002 and 0.006 with `coding = "variable"` (with and without Gamma(4) at alpha 1.10) and `coding = "all"`. The largest difference per partition is 1.1e-10 without Gamma and 8.5e-7 with alpha 1.10 (relative 1.4e-9). With RevBayes's rates at alpha 1.10 the largest difference is 2.17e-11. C++ and R agree to 1.5e-15 (relative) (`hymenoptera.csv`, `hymenoptera-compare.txt`; Slurm 19467761).

### 2. Exact one-dimensional oracle for alpha

`gamma-acrv-alpha-oracle.R` runs on the validation-1 tree and characters (`coding = "variable"`, 4 categories) with every parameter, the tree and the paths fixed, so only alpha moves, under the prior log alpha ~ N(0, 1). It uses 16 chains of 20 000 iterations after 1000 of burn-in, with gamma = 0 and gamma = 1. With gamma = 1 each chain has its own path draw, and so its own target. For each chain the exact posterior CDF of log alpha is computed on a 2001-point grid over +-8 posterior sd, and the probability integral transforms u = F(alpha_i) are pooled; batch-means z-scores use 25 batches per chain (Hamilton, Slurm 19448620, 19448621; repeated after the review fixes at dd15973, Slurm 19467762, 19467763, with the same results).

- gamma = 0: mean(u): z = 0.91, P(u < 0.05): z = -0.74, P(u < 0.25): z = -0.55, P(u < 0.50): z = -0.66, P(u < 0.75): z = -0.21, P(u < 0.95): z = -0.71, mean(log alpha - exact mean): z = 0.76. KS test on every 50th u (n = 6400): KS D 0.0091 p 0.6598; max|z| 0.912.
- gamma = 1: mean(u): z = 1.39, P(u < 0.05): z = 0.53, P(u < 0.25): z = -1.20, P(u < 0.50): z = -1.52, P(u < 0.75): z = 0.52, P(u < 0.95): z = -0.96, mean(log alpha - exact mean): z = 1.21. KS: KS D 0.0126 p 0.2655; max|z| 1.523.

Alpha acceptance is 0.22-0.40 per chain (`alpha-oracle-gamma*-chains.csv`).

### 3. SBC-lite

`gamma-acrv-sbc.R` is set up as rung 6 (`rung6-sbc.R`): prior draws, exact dCET simulation on a K = 16 grid, a fixed true tree with 15-40 tips, and a start from an independent prior draw. On top of that:
- alpha ~ LogNormal(0, 0.5);
- 45 variable characters (30 binary, 15 three-state), each with a Gamma category drawn uniformly, redrawn until variable, scored with `MkCharacters(coding = "variable", nRateCat = 4)`;
- every parameter and alpha free;
- moves `tree = 1, node = 3, tip = 3, branch = 3, shift = 1`, with 600 burn-in;
- ranks on 50 equally spaced draws, and a chi-square test with 10 bins, Holm-adjusted over the quantities of a cell.

The ranks take 51 values (0 to 50), so the 10 bins are not equally likely: the chi-square test uses each bin's exact probability (`chisq.test(p = )`). The p-values below are from that test (Slurm 19467765). The first version of this README used equal bin probabilities; this changes p-values but no conclusion. Per-arm tables are in `sbc-<arm>-summary.txt` and `sbc-<arm>-tests.csv`. Per-replicate ranks remain on Hamilton (`val-fbbf149/scratch/`, archived in `val-fbbf149/diag/sbcres.tgz`).

**3000 iterations, both cells** (100 replicates per cell; Slurm 19448196; `sbc-3k-gamma-*`):
- **Cell (1, 0):** every quantity passes. The smallest raw p is 0.046 (psi; Holm 0.79). Theta p is 0.39, alpha p 0.79.
- **Cell (1, 1):** theta fails (chi^2 30.5, p 0.00035, Holm 0.0064, mean rank 0.537). The next smallest raw p is for s (0.050, Holm 0.81); alpha p is 0.73.
- **Mixing:** median theta ESS over the 50 draws is 13 in cell (1, 1) and 33 in cell (1, 0). This run predates whole-chain ESS, so its summary's "whole chain" label is wrong for it.

**Is the theta failure caused by the Gamma code?** I ran matched arms in cell (1, 1). They share seeds, so they have the same trees and parameter draws and differ only in the characters. The plain arm uses main's character path: 45 binary characters as a matrix, with no `MkCharacters()` and no Mkv.

| arm (cell 1, 1) | iterations | replicates | theta chi^2 p | Holm p (min over quantities) | theta mean rank (z) | median theta ESS |
|---|---|---|---|---|---|---|
| Gamma, Mkv | 3000 | 100 | 0.00035 | 0.0064 (theta) | 0.537 (1.26) | 13 (of 50 draws) |
| no Gamma, Mkv | 3000 | 100 | 0.00085 | 0.0068 (theta) | 0.577 (2.62) | 15 |
| main's path, no Mkv | 20 000 | 200 | 0.034 | 0.27 (theta) | 0.550 (2.38) | 54 |
| no Gamma, Mkv | 20 000 | 200 | 0.34 | 1.00 | 0.547 (2.27) | 51 |
| Gamma, Mkv | 20 000 | 200 | 0.097 | 0.35 (mu) | 0.548 (2.30) | 58 |

Slurm 19448622-4 and 19449104-6 ran the arms; 19449529 and 19467765 summarised them. ESS is Geyer's estimate over the whole chain unless marked. z is the mean rank's distance from 1/2 in units of 0.294 / sqrt(n), from the unrounded mean ranks in `sbc-<arm>-tests.csv`.

- **3000 iterations is too short for theta at gamma = 1.** Rung 6 full (#21) needed a median of 40 000 iterations and ESS >= 50 to pass this cell.
- **At 20 000 iterations no arm fails after Holm adjustment.** Theta's mean rank is about 0.55 in all three arms, including main's path. Split by theta ESS, the shift is at least as large among replicates with ESS >= 50 (`sbc-theta-by-ess.txt`):
  - main's path: mean rank 0.557, KS p 0.017, n = 108;
  - Mkv without Gamma: 0.565, KS p 0.080, n = 101;
  - Mkv with Gamma: 0.568, KS p 0.029, n = 121;
  - ESS < 50: 0.517-0.540.
- **Conclusion.** The theta shift is not introduced by this PR: it appears equally on main's character path, without the Gamma categories or the Mkv path. It was reported against main as #106, which the lead closed as chance. The arms share their truths, so the three agreeing means are one experiment (z about 2.3), not three.
- **The new parameter alpha passes.** Ranks are centred (mean rank 0.497, z -0.15; 3000-iteration arm p 0.73), and it mixes well (median ESS 3322 per 20 000 iterations). With 200 replicates its histogram is uneven (chi^2 15.6, p 0.075, bins 13-35) and passes after Holm adjustment (p 0.60). The exact oracle (validation 2) checks its move directly. mu (p 0.039, Holm 0.35) is the only raw p below 0.05 in the 20 000-iteration arms apart from theta on main's path.

**Not covered by SBC.** The prior LogNormal(0, 0.5) on alpha rarely gives small alpha (alpha < 0.3 has prior probability 0.8%), so strongly heterogeneous rates were not calibrated. The SBC characters are unordered, with no missing cells, no ambiguity and one mask. Ordered characters, sets and several masks are checked only by validations 1 and 4 and by the unit tests. Cells other than (1, 0) and (1, 1) were not run.

### 4. Transition probabilities at the empirical settings

The review (R98-13) asked whether P(t) is accurate at the rates that real analyses visit. `gamma-acrv-pt-check.R` evaluates single-character likelihoods in both engines and compares them with an independent reference. The reference has its own pruning with rescaling, and computes P(t) by a direct Poisson sum without halving or squaring (scaling and squaring of the Taylor series once mu t > 200). The settings:
- trilobites: the 16 start trees and the matrix on main (`data-raw/empirical/trilobites/cet`), K up to 16, two ordered characters;
- Hymenoptera: the ordered partitions (o3, o4, o5) on the start tree of the CET pipeline (e0c1384);
- clock rates 10^-9, 10^-8, ..., 10 and 30 per Myr, each also multiplied by e^-5 and e^5 (gamma = 1 branch rates at z = +-5 with s = 1);
- no Gamma, and Gamma with 4 categories at alpha 1e-3, 1.1 and 1e3.

In total there were 23 472 evaluations over 17 trees (Slurm 19467764, dd15973; `pt-check.txt`, `pt-check-by-rate.txt`). Unordered characters use the closed form in both engines; they were compared between engines and checked for finiteness, not against the reference.

- **R engine:** finite everywhere. It agrees with the reference to 5.8e-10 (relative) or better at every effective rate (rate x scale) from 6.7e-12 to 4.45e3.
- **Hymenoptera ordered (432 evaluations):** both engines finite at every rate. C++ and R agree to 1.4e-13 at effective rates from 1e-8 to 148 per Myr. They agree to 8.7e-12 at 1.5e3 to 4.5e3 per Myr, and to 7.3e-11 at 6.7e-12 per Myr. Both engines agree with the reference to 5.8e-10.
- **Trilobites (4608 ordered, 18 432 unordered):** at effective rates of 1e-7 per Myr and above, C++ is finite everywhere and agrees with R to 1.6e-12. At effective rates of 1e-8 per Myr and below, 309 C++ values are -Inf (281 ordered, 28 unordered). Below 1e-7, C++ also differs from R by up to 8.3e-5. In these cases a single character's lnL is about -740 to -1130, which is below the double-precision range of the unscaled C++ Pruner (Limitations). Such rates are far below any trilobite clock rate a chain would visit. R remains exact there.

## Status

Validations 1, 2 and 4 pass. Validation 3 passes for every quantity in cell (1, 0). In cell (1, 1) every quantity passes after Holm adjustment at 20 000 iterations. Theta's upward rank shift in cell (1, 1) also appears on main's character path (#106, closed as chance).

The lead ran the first validation jobs on Hamilton at fbbf149. After the review, the package changed (P(t), guards, input checks), so validations 1, 2 and 4 and the Gamma-off identity were rerun at dd15973 (Slurm 19467760-19467764). The tests and R CMD check were rerun at a277d3a (Slurm 19467785). Commits after dd15973 change only tests, documentation and these results. The SBC (validation 3) ran at fbbf149. Between fbbf149 and dd15973 the C++ engine changed P(t) only for ordered characters. The SBC characters are all unordered, and their closed-form JC model is unchanged. Only the SBC summaries were recomputed (Slurm 19467765).

To rerun, from a checkout of this branch with `LIB` set to its installed library:

```bash
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh rb <scratch>
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh hymenoptera <scratch>
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh oracle <scratch> 0   # -c 16
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh oracle <scratch> 1   # -c 16
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh ptcheck <scratch>
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh sbc <scratch> <task> <nTask>   # array; SBC_CELLS, SBC_NCAT, SBC_ITER, SBC_REPS
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh sbc-summary <scratch>
```

`pt-check.csv` (3.2 MB, one row per evaluation) is not committed; it is on Hamilton in `val-fbbf149/src/data-raw/gamma-acrv-results/`.

## Incremental pruning (#44)

#44 is a draft and is not on main, so this branch does not touch it. Combining them needs:

1. `CharacterCache` holds one `Pruner` per `MkCharacters` partition, including its constant columns.
2. `set_edge_length()` / `set_tree()` goes to each Pruner.
3. `set_rates()` is called when alpha changes; it marks every node dirty, so an alpha move costs one full pruning.
4. `.AlphaMove()` reports its decision through `.ResolveCharacters()`.
5. The correction is re-summed from the cached site likelihoods.

The bitwise identity of #44 (cache vs full pruning) holds as long as the per-partition pruning order is fixed.

## Limitations

- The C++ Pruner does not rescale partial likelihoods. When a single character's likelihood is below about 1e-308 (lnL below about -709), C++ returns -Inf, or loses accuracy just above that bound, while the R reference stays finite. The limit is per character, not per data set. It is reached with many tips, or at very low rates on long trees: in the trilobite check (validation 4), only at effective rates of 1e-8 per Myr and below.
- As in RevBayes, conditioning on variability is exact only for characters without partial sets: a character such as `0, 0, {01}` is kept although some of the histories it sums over are constant.
- Rate categories are equiprobable, with mean rates (RevBayes default `median = FALSE`).
