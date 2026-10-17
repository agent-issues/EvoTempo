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
| 2026-10-16 | Windows, MinGW (Rtools45), R-devel | after round 2 | 111 | 39 | loaded machine: round-1 build 135 / 44 in the same session |
| 2026-10-16 | Linux, gcc 14.2, R 4.5.1 (Hamilton `test`) | after round 2 | 50.1 | 18.1 | tree draw 1.1 ms per call; see log.md round 2 |
| 2026-10-16 | Windows, MinGW (Rtools45), R-devel | after round 3 (vs round 1, without round 2) | 134 -> 116 (same session) | -- | tree moves now ~5 % |
| 2026-10-16 | Linux, gcc 14.2, R 4.5.1 (Hamilton `shared`) | after round 3 (vs round 1, without round 2) | 58.0 -> 51.7 (same job) | -- | tree moves 17 % -> 8 % |

`all` = moves `tree = 1, node = 3, tip = 3, branch = 3, joint = 1, nodeAge = 1,
fossilAge = 1, slide = 1, saJump = 1` (the full-size SBC design of #21);
`rung6` = `tree = 1, node = 3, tip = 3, branch = 3` (the current
`data-raw/rung6-sbc.R`).

## Components (`drivers/tempo-components.R`, microseconds per call)

Windows "before" was measured under heavier machine load than the A/B runs;
compare Linux columns for ratios.  Columns come from different Hamilton
sessions (node and load differ by a few per cent); same-session A/B ratios are
in log.md.  "r2" is round 2 (bridge jump times); "r3" is round 3 (tree moves),
measured against the round-1 build without round 2, in a `shared`-partition job
where the round-1 build timed 3-10 % above the "r1" column.  Against the
round-1 build in the same session, round 2 left `.Messages` (unchanged code in
the same translation unit as `Bridge()`) 828 -> 859-875 us on Linux, likely a
code-layout effect.

| Component | Win before | Linux before | Linux after r1 | Linux after r2 | Linux after r3 |
|-----------|-----------:|-------------:|---------------:|---------------:|---------------:|
| `.Precompute` | 4541 | 2563 | 1234 | 1219 | 1281 |
| `.EGrid` | 3627 | 2125 | 773 | 766 | 813 |
| `.AddSlices` | 2767 | 1469 | 477 | 469 | 500 |
| `.Messages` | 1760 | 828 | 727 | 875 | 797 |
| `.SampleTree` | 7117 | 2188 | 1906 | 1094 | 2000 |
| `.LogWeight` | 551 | 262 | 262 | 262 | 281 |
| `.LogAugmented` | 191 | 95 | 98 | 96 | 102 |
| `MkLikelihood` | 86 | 44 | 44 | 44 | 45 |
| `tree` move | 9495 | 2844 | 2563 | 1813 | 2688 |
| `node` move | 1768 | 680 | 656 | 594 | 703 |
| `nodeAge` move | 5611 | 2688 | 2563 | 2656 | 1297 |
| `slide` move | 8042 | 4000 | 3813 | 3688 | 1656 |
| `saJump` move | 7145 | 2688 | 2594 | 2531 | 1094 |
