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

Not run yet (see Status).

### 2. Exact one-dimensional oracle for alpha

Not run yet (see Status).

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
