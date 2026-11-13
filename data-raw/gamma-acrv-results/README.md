# Gamma rates among characters and per-mask Mkv (#97): validation

`MkCharacters()` (R/Characters.R, src/characters.cpp) adds to the CET character model:

- a number of states per character;
- ordered characters (RevBayes `fnOrderedRateMatrix(maxState = K)`);
- ambiguity sets;
- the Mkv correction per missing-data mask (RevBayes `dnPhyloCTMC(coding = "variable")`; #47, #65);
- discrete Gamma(alpha, alpha) rates among characters (`GammaRates()`, as `fnDiscretizeGamma(alpha, alpha, K)`).

In `TempoMCMC()` the shape `alpha` is `start$alpha`. It is free when `charPriors$alpha` gives its prior, and then has its own move (`.AlphaMove()`).

**Mixture and correction.** With rate categories c (weights w_c = 1/K), a character's likelihood is sum_c w_c L_c. The correction for mask m is -n_m log(1 - sum_c w_c P_c,m), where P_c,m is the probability that every scored tip of m is in one state. This is RevBayes's `PhyloCTMCSiteHomogeneousConditional::sumRootLikelihood()`: the correction is computed within each category and combined over categories before the log is taken. As in RevBayes, a missing cell and a set of all K states are gaps of the mask. A character whose scored cells all hold the same state or set is an error (RevBayes excludes it with a warning).

Scripts (from the package root, with EvoTempo installed; Hamilton driver `data-raw/gamma-acrv-hamilton.sh`):

| script | check |
|---|---|
| `data-raw/gamma-acrv-identity.R` | Gamma off: seeded chains on a character matrix, main vs this branch |
| `data-raw/gamma-acrv-revbayes.R`, `.Rev` | 1. fixed-tree lnL vs RevBayes, 3 rates x 3 shapes (+ no Gamma), ordered characters, ambiguity, both engines; and on the Hymenoptera matrix |
| `data-raw/gamma-acrv-alpha-oracle.R` | 2. alpha move vs the exact 1-D posterior on a fixed tree with fixed paths |
| `data-raw/gamma-acrv-sbc.R` | 3. SBC-lite, cells (beta, gamma) = (1, 0) and (1, 1) |
| `tests/testthat/test-MkCharacters.R` | unit tests (below) |

## Results

### Gamma off: chains unchanged

Characters given as a matrix take the code path of main. With rate categories, `TempoMCMC()` passes `alpha` to `MkLikelihood()`, so pipelines that rebind `MkLikelihood()` with the old signature are unaffected. Without them, no random numbers are drawn and no columns are added.

The check is `gamma-acrv-identity.R`. Four seeded chains (C++ and R engine x gamma 0 and 1; 20 burn-in + 30 sampled iterations with every move: `tree`, `node`, `tip`, `branch`, `joint`, `shift`, `nodeAge`, `fossilAge`, `slide`, `saJump`) on the test fixture (sampled ancestor, fossil age ranges). They ran under main fb24ed6 and under this branch (Hamilton, Slurm 19431188). `samples`, `occupancy`, `acceptance`, `burninAcceptance` and `moves` are `identical()` in all four (`identity.txt`).

### Unit tests (`test-MkCharacters.R`)

- `GammaRates()` against Yang (1994) Table 1 (alpha 0.5), with mean 1.
- Parsing of states, sets and unknowns. A full set equals a gap. Input errors.
- `MkCharacters()` with one state count, no correction and no Gamma equals `MkLikelihood()` on the matrix to 1e-12, in both engines.
- Ordered characters and sets on a two-tip tree equal `expm(Q t)` (scaling and squaring), in both engines.
- C++ and R engines agree to 1e-10 (relative) for coding all and variable, 1 or 4 categories, and alpha 0.2, 1.1 and 8.
- The Mkv correction: over every variable pattern of a mask (binary with all or 4 tips scored, ordered K = 3, unordered K = 3; with Gamma), the conditional probabilities sum to 1 (to 1e-12). The value equals RevBayes's mixture formula log mean_c L_c - log(1 - mean_c P_c), built from Gamma-free pieces, to 1e-12.
- `TempoMCMC()` with alpha free: identical seeded chains in both engines; the `alpha` column and acceptance; `logPost` equals the recomputed augmented density + character lnL + priors + alpha prior; fixed alpha; `updateParameters = FALSE`; error paths.

On Hamilton (Slurm 19431558, 19431559, build of 7e7d0f5): the C++ engine gives 110 test blocks and 4729 expectations, 0 failed and 3 skipped (`skip_on_cran()` oracles). The R engine has 6 failures in the 2 stiff-solver tests of `test-Paths.R`. The same 6 failures occur on main fb24ed6 (Slurm 19431685); #45 addresses them.

`R CMD check --no-manual --ignore-vignettes` of 3e6a4da (Hamilton, r/4.5.1, Slurm 19432046): `Status: OK`, with no NOTE or WARNING. The new tests at 3e6a4da (Slurm 19432045) give 0 failures in both engines.

### 1. Fixed-tree log-likelihood vs RevBayes

`gamma-acrv-revbayes.R prepare` simulates the data. The tree is a coalescent time tree of 24 extant tips with root age 1; 4 tips are cut back to fossils, and one sampled ancestor sits on an internal branch. There are 74 characters: unordered K = 2, 3, 4 (30, 12, 8) and ordered K = 3, 4, 5 (10, 8, 6). Of the cells, 247 are missing (`?` or `-`) and 129 are ambiguity sets, including full binary `{01}`. RevBayes is the development build fd807d10 (`1.5.0-preview`) running `gamma-acrv-revbayes.Rev`: `dnPhyloCTMC(Q = fnJC(K) | fnOrderedRateMatrix(maxState = K), branchRates = rate, siteRates = fnDiscretizeGamma(a, a, 4), type = "Standard", coding = "variable")`, plus a run without `siteRates`. CET is `MkLikelihood(MkCharacters(..., coding = "variable", nRateCat = 4), lengths = rate x duration, alpha = a)` in both engines. The table gives totals over the 6 partitions, and maxima of |difference| per partition (`revbayes.csv`, `revbayes-compare.txt`; Hamilton, Slurm 19448619, branch e743555).

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

The largest relative difference is 1.5e-9, against a tolerance of 1e-6. At alpha 0.3 and 5 the differences (<= 2.2e-9) are at the level of RevBayes's printed precision. At alpha 1.1 they are up to 2.9e-7 per partition, in every partition and at every rate, which points to a difference in the discretised rates at that shape between R's `qgamma()` and RevBayes's chi-square quantile. That is not verified (RevBayes's rates were not printed).

**Hymenoptera matrix** (supplementary). Data are from the CET pipeline at e0c1384 (`hymenoptera-sha.txt`): 113 tips, 343 characters in 9 (ordering, K) partitions, ordered characters and ambiguity, on its start tree. The oracle is its RevBayes `check-mkv.Rev` output (fd807d10), at clock rates 0.0007, 0.002 and 0.006 with `coding = "variable"` (with and without Gamma(4) at alpha 1.10) and `coding = "all"`. The largest difference per partition is 1.1e-10 without Gamma and 8.5e-7 with alpha 1.10 (relative 1.4e-9). C++ and R agree to 6.9e-13 (`hymenoptera.csv`, `hymenoptera-compare.txt`).

### 2. Exact one-dimensional oracle for alpha

`gamma-acrv-alpha-oracle.R` runs on the validation-1 tree and characters (`coding = "variable"`, 4 categories) with every parameter, the tree and the paths fixed, so only alpha moves, under the prior log alpha ~ N(0, 1). It uses 16 chains of 20 000 iterations after 1000 of burn-in, with gamma = 0 and gamma = 1. With gamma = 1 each chain has its own path draw, and so its own target. For each chain the exact posterior CDF of log alpha is computed on a 2001-point grid over +-8 posterior sd, and the probability integral transforms u = F(alpha_i) are pooled; batch-means z-scores use 25 batches per chain (Hamilton, Slurm 19448620, 19448621).

- gamma = 0: mean(u): z = 0.91, P(u < 0.05): z = -0.74, P(u < 0.25): z = -0.55, P(u < 0.50): z = -0.66, P(u < 0.75): z = -0.21, P(u < 0.95): z = -0.71, mean(log alpha - exact mean): z = 0.76. KS test on every 50th u (n = 6400): KS D 0.0091 p 0.6598; max|z| 0.912.
- gamma = 1: mean(u): z = 1.39, P(u < 0.05): z = 0.53, P(u < 0.25): z = -1.20, P(u < 0.50): z = -1.52, P(u < 0.75): z = 0.52, P(u < 0.95): z = -0.96, mean(log alpha - exact mean): z = 1.21. KS: KS D 0.0126 p 0.2655; max|z| 1.523.

Alpha acceptance is 0.22-0.40 per chain (`alpha-oracle-gamma*-chains.csv`).

### 3. SBC-lite

Not run yet (see Status).

## Status

The Hamilton job script that runs validations 1-3 was held by a safety-approval check on 2026-11-11, and has not been run. To run them on Hamilton, from a checkout of this branch with `LIB` set to its installed library:

```bash
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh rb <scratch>
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh hymenoptera <scratch>
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh oracle <scratch> 0   # sbatch -c 16
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh oracle <scratch> 1   # sbatch -c 16
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh sbc <scratch> <task> 25   # array 1-25, -c 8
LIB=<lib> bash data-raw/gamma-acrv-hamilton.sh sbc-summary <scratch>
```

## Incremental pruning (#44)

#44 is a draft and is not on main, so this branch does not touch it. Combining them needs:

1. `CharacterCache` holds one `Pruner` per `MkCharacters` partition, including its constant columns.
2. `set_edge_length()` / `set_tree()` goes to each Pruner.
3. `set_rates()` is called when alpha changes; it marks every node dirty, so an alpha move costs one full pruning.
4. `.AlphaMove()` reports its decision through `.ResolveCharacters()`.
5. The correction is re-summed from the cached site likelihoods.

The bitwise identity of #44 (cache vs full pruning) holds as long as the per-partition pruning order is fixed.

## Limitations

- The C++ Pruner does not rescale partial likelihoods. A single character likelihood below about 1e-308 (thousands of tips) gives -Inf in C++, but a finite value in the R reference.
- Rate categories are equiprobable, with mean rates (RevBayes default `median = FALSE`).
