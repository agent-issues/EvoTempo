# Metropolis-coupled chains in `TempoMCMC()` (#110)

## Target

Level r = 1..N of the ladder has heat h_r = 1 / (1 + deltaHeat (r - 1)); level 1 is cold.
Its target is

    p(parameters) x augmented density(tree, paths | parameters) x L_char^(h_r).

Only the character likelihood is heated.
MrBayes and RevBayes `mcmcmc` heat the whole posterior; here that would also heat the
augmented (tree and path) density, which every path and tree proposal is built to match
(FFBS and frozen bridges), so heated levels would lose their near-exact proposals.
The topology modes of #110 come from the characters, and h -> 0 gives the FBD-CET prior,
which the sampler mixes well on (rung 5).
It is also the power that the stepping-stone ladders use, so the two compose (below).

## Implementation

- `CharLL()` returns h x log L_char, with h the heat of the level being updated.
  Every move already uses `CharLL()` differences, so no move changes.
- State is held per *slot* (`mc3$chain[[c]]`: tree, parameters, paths, frozen set,
  messages, and later a persistent `MkPrime::Pruner`); tuning is held per *level*
  (`mc3$tune[[r]]`: step sizes, adaptation counters, move counts after the burn-in drop,
  acceptance counters). `mc3$slot[r]` names the slot at level r.
- A swap exchanges `slot[a]` and `slot[a + 1]`: states never move, and step sizes stay
  with the temperature they were tuned at. The stored `charLL` of the two states is
  re-evaluated at its new heat.
- Swap of adjacent levels a, a + 1 holding states A, B:
  log ratio = (h_a - h_{a+1}) (log L_B - log L_A). The pair is uniform among the N - 1
  adjacent pairs, independently of the state, so the proposal is symmetric.
- Schedule: every `swapInterval` iterations (default 1), one attempt, during burn-in too.
- Conventions: ladder as MrBayes and RevBayes; `deltaHeat` default 0.1 (MrBayes `Temp`);
  one attempt per interval (both); adjacent pair (RevBayes `swapMethod = "neighbor"`;
  MrBayes picks any pair by default); argument names from RevBayes. An iteration here is
  a full sweep of moves, so `swapInterval = 1` follows MrBayes `Swapfreq = 1`, not
  RevBayes's 10.

## Random numbers

One interleaved stream (R's). Start-up: the paths of level 1, then of levels 2..N.
Each iteration: levels 1..N in order, each consuming what a single chain does; then, if a
swap is due, one uniform for the pair and one for acceptance (both always drawn).
N = 1 draws nothing extra, so it reproduces earlier output exactly; both engines consume
the same numbers.

## Output

Samples, summaries, `acceptance`, `moves` and `state` are those of the cold level.
With N > 1: `heats`, and `swaps` (per adjacent pair: attempts, acceptances, rate, over the
sampling phase).

## Power posteriors

A stepping-stone rung powers the character likelihood by t (pipelines: rebound
`MkLikelihood()`). The heat multiplies it: level r targets L^(t h_r), the cold level is
the rung. Nothing forbids the combination.

## Parallelism

The loop is in R, so levels run in series: N levels cost N times one chain.
Threads would need per-level random streams (breaking the single-stream scheme) and a
loop in C++; follow-up.

## Tests

- N = 1: identical to `main` for the same seed, both engines (checked against saved output).
- N = 3: R and C++ engines agree.
- Exact oracle (`helper-oracles.R`): steep ladder, cold-level posterior of `rate` against
  the grid posterior; a sign-flipped swap ratio must fail it.
