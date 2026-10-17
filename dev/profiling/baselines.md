# Profiling baselines

Current timings for `/profile regress` (flag > 10 % slowdowns). Measure A/B
interleaved, at least three repeats; the Windows machine's load varies by
2x over a day, so compare only within one session.

Workload: `drivers/tempo-mcmc.R` on `drivers/fixture-rung6.R` (25 tips, 5
fossils, 3 sampled ancestors, 30 binary characters; K = 16, z in [-3, 3],
nSlice = 16, beta = gamma = 1, six free parameters, `thin = 1`).

| Date | Platform | Build | `all` ms/iter | `rung6` ms/iter | Top hotspots |
|------|----------|-------|---------------|-----------------|--------------|
| 2026-10-16 | Windows, MinGW (Rtools45), R-devel | before round 1 | 129 | 41.6 | `.SampleTreeCpp` 27 %, `.EGridCpp` 22 %, `.MessagesCpp` 7 %, `.LogWeightCpp` 5 %, `eigen` 4 % (Rprof, `all`) |
| 2026-10-16 | Windows, MinGW (Rtools45), R-devel | after round 1 | 89 | 29 | see log.md round 1 |
| 2026-10-16 | Linux, gcc 14.2, R 4.5.1 (Hamilton `test`) | before round 1 | 83.8 | 29.0 | `.EGridCpp` ~30 %, `.SampleTreeCpp` ~18 %, `.AddSlices` ~12 % (component timings) |
| 2026-10-16 | Linux, gcc 14.2, R 4.5.1 (Hamilton `test`) | after round 1 | 54.3 | 18.8 | `joint` ~44 %, tree moves ~17 %; tree draw 1.9 ms per call |

`all` = moves `tree = 1, node = 3, tip = 3, branch = 3, joint = 1, nodeAge = 1,
fossilAge = 1, slide = 1, saJump = 1` (the full-size SBC design of #21);
`rung6` = `tree = 1, node = 3, tip = 3, branch = 3` (the current
`data-raw/rung6-sbc.R`).

## Components (`drivers/tempo-components.R`, microseconds per call)

Windows "before" was measured under heavier machine load than the A/B runs;
compare Linux columns for ratios.

| Component | Win before | Linux before | Linux after |
|-----------|-----------:|-------------:|------------:|
| `.Precompute` | 4541 | 2563 | 1234 |
| `.EGrid` | 3627 | 2125 | 773 |
| `.AddSlices` | 2767 | 1469 | 477 |
| `.Messages` | 1760 | 828 | 727 |
| `.SampleTree` | 7117 | 2188 | 1906 |
| `.LogWeight` | 551 | 262 | 262 |
| `.LogAugmented` | 191 | 95 | 98 |
| `MkLikelihood` | 86 | 44 | 44 |
| `tree` move | 9495 | 2844 | 2563 |
| `node` move | 1768 | 680 | 656 |
| `nodeAge` move | 5611 | 2688 | 2563 |
| `slide` move | 8042 | 4000 | 3813 |
| `saJump` move | 7145 | 2688 | 2594 |
