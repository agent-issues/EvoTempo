# Hymenoptera: Ronquist et al. 2012, with interval ages from Zhang et al. 2016

Main total-evidence-dating test set for CET.

- **Study:** Ronquist F., Klopfstein S., Vilhelmsen L., Schulmeister S., Murray D.L., Rasnitsyn A.P. 2012. A total-evidence approach to dating with fossils, applied to the early radiation of the Hymenoptera. *Systematic Biology* 61:973–999. doi:10.1093/sysbio/sys058
- **Reanalysis (age intervals, FBD settings):** Zhang C., Stadler T., Klopfstein S., Heath T.A., Ronquist F. 2016. Total-evidence dating under the fossilized birth–death process. *Systematic Biology* 65:228–249. doi:10.1093/sysbio/syv080
- **Data:** Dryad doi:10.5061/dryad.j2r64 (Ronquist et al.; version 1). Zhang et al. supplementary tables: Dryad doi:10.5061/dryad.26820. Cross-check: `examples/hymfossil.nex` from MrBayes (NBISweden/MrBayes, commit `77da419ecb2a335fb185a79c0d7d42d413562f19`).
- **Licence:** CC0-1.0 for both Dryad packages; the processed files are released under the same terms. `raw/hymfossil.nex` is distributed with MrBayes under GPL-3.

## Files

| File | Content |
|---|---|
| `hymenoptera_morph.nex` | 113 taxa × 353 morphological characters as deposited. `EXSET * published` lists the 10 characters the analyses exclude (343 analysed); `TYPESET * published` lists the 44 characters treated as ordered. `CHARSET`s give the Vilhelmsen (1–236) and Schulmeister (237–353) blocks |
| `hymenoptera_characters.csv` | Per character: block, ordered, excluded, states observed |
| `hymenoptera_dna.nex` | 113 taxa × 5096 sites (7 markers; fossils all-missing), columns 354–5449 of the source renumbered from 1, with `CHARSET`s per marker and codon position |
| `hymenoptera_tips.csv` | Taxon, extant/fossil, Hymenoptera/outgroup, `min_age`/`max_age` (Ma; Zhang et al. 2016 stratum intervals), the point ages used by Ronquist et al. 2012 and by `hymfossil.nex`, stratum, notes |
| `hymenoptera_constraints.csv` | The 10 MrBayes constraints, resolved to taxon names; `used_in_ted` marks the two used in total-evidence dating |
| `curated/fossil-ages.csv` | Per-fossil stratum and species (Ronquist et al. 2012 Table 2, transcribed from the PMC full text) and stratum interval (Zhang et al. 2016 Table 3) |
| `raw/` | Byte-exact downloads; checksums in `raw/MANIFEST.csv` |

## Raw downloads

| File | Bytes | MD5 | SHA-256 | Committed | Source |
|---|---|---|---|---|---|
| `Appendix-1_MorphoCharacterDescriptions.pdf` | 701554 | `4185c2c23080f72d70948c3e800f9fc0` | `71616541c7563655138ef9dd3ab18d84b316ca1aa21aa22483fa22e016c33fac` | yes | doi:10.5061/dryad.j2r64 |
| `MrBayes_Commands.nex` | 659786 | `f4fb160d2cc16d07a5173dff77199577` | `f289112f8ae099a2d5b50dc672a2c4011955486ddabbcdeefb4341e5b7dce464` | yes | doi:10.5061/dryad.j2r64 |
| `R-Scripts_PriorFitting-RelaxedClockModels.r` | 16979 | `b941a11f4d068de130eb44b6678c25f6` | `9b7c480ecc392defa095d72889ae08bb6f20cc11a8d198391a35c920b298846e` | yes | doi:10.5061/dryad.j2r64 |
| `Supplementary-Figure-1_Morphological-versus-molecular-BranchLengths.eps` | 2451664 | `6e1ec597667f032e184a0faa1b3de6fe` | `7cac74c27e864cd681d4be027492b690b0dea3a13dc552ae4125a08550f95c94` | no | doi:10.5061/dryad.j2r64 |
| `Supplementary-Figure-2_DatedTree-ErrorBars.eps` | 1378818 | `ca0ded4f24d0360f50fe7b4ecc24e566` | `3ee175ca3a95720358581db5577220b33dbd7ecc804312fb0b9359ebc05bb4d8` | no | doi:10.5061/dryad.j2r64 |
| `FBD_Suppl.pdf` | 836057 | `0a2bc505285ac77dcdd21d4ef96a85d6` | `1abbc996c5453f5cb6c5883dadb10715d333c1fe737297ec9cdf98d14aaee603` | yes | doi:10.5061/dryad.26820 |
| `hymfossil.nex` | 665462 | `340f0980426d00a535631cb53225efbc` | `817e252b0f04b9040a7818908818747bdc3399c0141c823daac7b4db89e1995c` | yes | NBISweden/MrBayes@77da419ecb2a335fb185a79c0d7d42d413562f19 |

The two EPS figure files are not committed; they carry no data.

## Transformations (all in `../build.R`)

1. Parse the interleaved `MATRIX` of `MrBayes_Commands.nex` (113 taxa × 5449 = 353 standard + 5096 DNA). Split at column 353.
2. Polymorphisms are written compactly: `(0,1)` → `(01)`, `{0,1}` → `{01}`. The source mixes both bracket types; both mean ambiguity to MrBayes. No other recoding.
3. `morph_ordered` and `morph_excluded` are copied from the MrBayes block into `TYPESET`/`EXSET`. Character 218 is in both sets.
4. DNA `CHARSET`s are shifted by −353 so that they index `hymenoptera_dna.nex`.
5. Tip ages. Point ages are parsed from the commented-out `calibrate ... Fixed()` statements of the total-evidence block. Intervals: each fossil's stratum (Ronquist et al. 2012 Table 2) is matched to the stratum interval of Zhang et al. 2016 Table 3. The per-class counts reproduce Table 3 exactly (2, 1, 4, 4, 3, 17, 1, 1, 6, 2, 2, 1, 1 fossils). Extant taxa are 0.
6. Constraints. The numeric taxon ranges in the MrBayes block are resolved against the 114-taxon order of MrBayes's `hymfossil.nex`, then *Lepid_Micropterigidae* (absent from this file) is dropped. `literal_113_matches` records whether a literal reading against this 113-taxon file gives the same clade. See Oddities.

## Source analysis settings

**Ronquist et al. 2012, total-evidence dating** (MrBayes 3.2; block in `MrBayes_Commands.nex`):

| Component | Setting |
|---|---|
| Tree prior | Uniform clock-tree prior (`brlenspr = clock:uniform`); fossils are tips only |
| Root | `calibrate root = offsetexp(315, 0.01234568)`: minimum 315 Ma, mean 396 Ma |
| Other calibration | `holometabola_withFossils = offsetexp(302, 0.0106383)`: minimum 302 Ma, mean 396 Ma |
| Topology | `constraints(root, holometabola_withFossils)`; outgroup Orthoptera |
| Fossil ages | Fixed point ages (`point_age_ronquist2012`) |
| Clock | IGR relaxed clock, `igrvarpr = exp(37.12)`; `clockratepr = lognorm(-7.08069, 2.458582)` |
| Morphology | Mkv (`coding = variable`) + Γ; 343 characters (10 excluded), 44 ordered (one of them excluded) |
| DNA | 7 partitions in all: morphology; 12S+16S; 18S; 28S (SYM+I+Γ, equal frequencies); CO1 positions 1–2; EF1α positions 1–2; EF1α position 3; each GTR+I+Γ except 28S. CO1 position 3 excluded. Parameters unlinked; partition rates variable |
| MCMC | 4 runs × 4 chains, temp 0.1, 2 × 10⁷ generations, sample every 1000 |

The node-dating, strict-clock, TK02 and CPP alternatives are commented out in the same block.

**Zhang et al. 2016, FBD reanalysis** (MrBayes 3.2.6; paper Appendix, Tables 3–5): fossil ages uniform on the intervals in `hymenoptera_tips.csv`; `brlenspr = clock:fossilization`; `treeagepr = offsetexp(315, 396)`; `clockvarpr = igr`; `clockratepr = lognorm(-7.1, 0.5)`; turnover and fossil-sampling proportion ~ Beta(1, 1); extant sampling proportion ρ = 0.0005 with diversified (`samplestrat = diversity`) or random sampling; optional piecewise-constant shifts at 252 Ma (and 66 Ma). The Appendix prints `speciationpr = exp(100)` for net diversification; the paper reports results under several d priors (Table 4) and uses an "intermediate" prior for the full analysis (Table 5). Same partitions as Ronquist et al. 2012 (CO1 position 3 defined but not used), but GTR+Γ / SYM+Γ without +I.

## Oddities

- **Constraint numbering is shifted by one.** The numeric ranges (e.g. `hymenoptera = 10-69`, `Xyelidae = 10-11`) fit the 114-taxon order of `hymfossil.nex`, which has *Lepid_Micropterigidae* at position 7. In this 113-taxon file they are off by one: read literally, "Hymenoptera" would exclude *Xyela* and include the fossil *Triassoxyela*. The two constraints used for total-evidence dating are unaffected (`root = 1-.`; `holometabola_withFossils = 3-114`, which overruns the 113 taxa and so takes everything from *Chrysopidae* on). The eight node-dating constraints are affected; `hymenoptera_constraints.csv` gives the intended clades.
- **Point ages differ from Table 2.** `MrBayes_Commands.nex` fixes *Eoxyela* and *Liadoxyela* at 180 Ma; Ronquist et al. 2012 Table 2 gives 176 Ma (Ichetuy Formation). This does not affect the interval ages.
- **Revised ages.** Zhang et al. moved several fossils well away from the 2012 point ages. *Prosyntexis* went from 140 Ma to 80–86 Ma (Ola Formation, Late Cretaceous), and the Turga fossils from 130 Ma to 133–135 Ma.
- **Provisional assignment.** Zhang et al. Table 3 has one fossil at 138.5 (125–152) Ma and six at 135 (125–145) Ma. Five Baissa fossils certainly take 125–145 Ma. Whether *Onochoius* (Onokhoy) or *Palaeathalia* (Laiyang) takes the wider 125–152 Ma interval is not stated in the paper, its Appendix (which elides most calibrations) or its Dryad supplement. *Palaeathalia* is given 125–152 Ma here because the age of the Laiyang Formation is the more uncertain. Both rows are flagged `assignment = inferred`. `hymfossil.nex` gives both a point age of 135 Ma, which cannot be squared with Table 3's counts.
- **Comment ranges.** The bracketed comments above each DNA block in the matrix (e.g. `[12S 354-559, 203 chars]`) are offset by +3 from the `charset` statements. The `charset`s are used: they tile 354–5449 exactly.
- **Structure.** 60 extant Hymenoptera, 8 holometabolan/paraneopteran outgroups, 45 fossils. All fossils are stem or crown Hymenoptera. They have no DNA and are 88% missing in the morphological matrix (4–20% scored); extant taxa are 30% missing.
