# Cambrian trilobites: Paterson, Edgecombe & Lee 2019

Tip-dated test of the Cambrian case. All 107 terminals are fossils.

- **Study:** Paterson J.R., Edgecombe G.D., Lee M.S.Y. 2019. Trilobite evolutionary rates constrain the duration of the Cambrian explosion. *PNAS* 116:4394–4399. doi:10.1073/pnas.1819366116
- **Data:** Dryad doi:10.5061/dryad.v7q827k (version 1).
- **Licence:** CC0-1.0 (Dryad). The processed files are released under the same terms.

## Files

| File | Content |
|---|---|
| `trilobites_morph.nex` | 107 species × 107 discrete characters as deposited (unordered, `DEFTYPE=UNORD`), with character labels; taxon names use underscores |
| `trilobites_characters.csv` | Per character: label, states observed, whether the BEAST 1.8 analyses treat it as ordered |
| `trilobites_quantitative.csv` | The 2 meristic and 6 continuous traits (characters 108–115) and the raw measurements behind them, as in `Dataset_3`; not used by Mk models |
| `trilobites_tips.csv` | Taxon, status, `min_age` = `max_age` = point age (Ma), order, classification, GTS 2012 stage, biozone, formation |
| `trilobites_constraints.csv` | Ingroup monophyly constraint and the fallotaspidoid outgroup |
| `raw/` | Byte-exact Dryad downloads; checksums in `raw/MANIFEST.csv` |

## Raw downloads

| File | Bytes | MD5 | SHA-256 | Committed | Source |
|---|---|---|---|---|---|
| `Dataset_1_Taxa,Dates,CodingRefs.xlsx` | 53164 | `b998c63b45be99a2bb3ee00c9cc7808f` | `ec9b5d8f197600da6a77d75019972901ff66d7a42a1b936935de69eb73ea4cef` | yes | doi:10.5061/dryad.v7q827k |
| `Dataset_2_Discrete_characters.nex` | 53090 | `ede9a9e048a229204b6bf759bc83bebd` | `9aa2a0a7b1eb9896db7729d731fee1e786973201606e00f052323c1fcef28c29` | yes | doi:10.5061/dryad.v7q827k |
| `Dataset_3_Meristic_Continuous_characters.xlsx` | 29236 | `2f0937ff3af54797dd2ebde87b53f809` | `f8454d65a1b7dd25a541579b8fb0b3331989a727de1c56e93df5c488ee4d14e4` | yes | doi:10.5061/dryad.v7q827k |
| `Dataset_4_BEAST_and_BEAST2.zip` | 723251 | `ee18ce8194f89cd6e4fe80617e00fc90` | `de8fb0122238c5b60274b54a6052f6eec1f1139286b3e573bb16b2d1d69185e5` | yes | doi:10.5061/dryad.v7q827k |
| `Dataset_5_TNT.zip` | 86469 | `499830a858b29d2f49f5b6f382e68d4b` | `b85839b8564fa89bd85ef5792c1c8464f2cf1a3cfc46245981ede19c3f3db204` | yes | doi:10.5061/dryad.v7q827k |

## Transformations (all in `../build.R`)

1. Parse `Dataset_2_Discrete_characters.nex` (written by NDE 0.5.0); replace spaces in taxon names by underscores. Polymorphisms are kept (`{01}`, `(23)`, ...). No recoding. The character labels are taken from `CHARLABELS`; the state labels that follow them in the source are not carried over (see the raw file).
2. Tip ages: `Age (Ma)` column of `Dataset_1`, taken as point ages (`min_age` = `max_age`). The script checks that they equal the `date` values in the BEAST 1.8 XMLs, which use them as fixed tip dates.
3. Constraints: the `Ingroup_node` (104 taxa) and `Fallotaspidoids` (3 taxa) taxon sets of the BEAST 1.8 XMLs.
4. `Dataset_3` is flattened to one header row (two header rows joined with ": "); group rows (EODISCIDA etc.) are dropped. Values are not transformed. The BEAST 1.8 XMLs carry transformed values of the six continuous traits as taxon attributes; that transformation was not reconstructed here.
5. Ordered characters: the BEAST 1.8 XMLs split the discrete data into partitions by number of states (2, 3, 4, 5, 6, 7, 9 and 16 states unordered; one 4-state and one 5-state partition ordered). Their columns were matched to the NEXUS by comparing state patterns. The ordered 4-state character is 47 (exact match) and the ordered 5-state character is 45 (consistent on all 86 jointly scored taxa). This matching was done once during curation (not in `build.R`); `ordered_in_beast` records the result.

## Source analysis settings

**BEAST 1.8 (main analyses, `Dataset_4/BEAST_1.8_executables_UPDATED`)**: discrete + meristic + continuous data, 107 taxa.

| Component | Setting |
|---|---|
| Tree prior | Birth–death serial sampling (`birthDeathSerialSampling`, BEAST 1.8; no extant sample): birth rate ~ U(0, 1000), relative death rate ~ U(0, 1), ψ ~ U(0, 1000), origin ~ U(0, 10⁸) |
| Root | Unconstrained: root height ~ U(0, 10⁶). Alternative `Epoch_..._Root522.xml`: root height ~ Normal(36, 0.001) above the youngest tip (486 Ma), i.e. root fixed at 522 Ma |
| Tip ages | Fixed point ages (tip-date sampling is commented out) |
| Topology | `Ingroup_node` (all but the three fallotaspidoids) monophyletic |
| Discrete model | Mkv (ascertained patterns) partitioned by number of states, partition rates scaled, Γ shape ~ Exp(mean 0.5); two characters ordered |
| Clocks (separate for discrete, meristic and continuous data) | Strict; UCLN (mean ~ U(0, 100), coefficient of variation ~ Exp(mean 1/3)); or a 3-epoch clock with boundaries 509 and 497 Ma (heights 23 and 11 above 486 Ma) |
| Meristic / continuous | 2 meristic traits as a discrete partition with adjacent-state-only rate matrix (ordered); 6 continuous traits under Brownian diffusion with Wishart priors |
| MCMC | 10⁸ generations; stepping-stone marginal likelihoods |

The paper's result is rate constancy across the Cambrian. The strict clock fits better than the epoch clock, and UCLN rates show no early peak.

**BEAST2 (discrete data only, `BEAST2_executables_DiscreteDataOnly`)**: 106 discrete characters (104 unordered + 2 ordered), strict clock. One file uses the FBD prior with sampled ancestors (`SABD`); the other approximates birth–death serial sampling without ancestors (`BDSS`: `SABirthDeathModel` with removal probability 0.999999999999999). Diversification ~ LogNormal(−3.5, 1.5); turnover and sampling proportion ~ U(0, 1); origin ~ U(0, 1000); Γ shape ~ U(0, 10); ingroup monophyletic.

## Oddities

- **Spelling.** The matrix and the BEAST XMLs spell *Neocobboldia chilinica*; `Dataset_1` and `Dataset_3` spell *N. chinlinica*. The matrix spelling is kept as the key; `name_in_spreadsheet` preserves the other.
- **No age intervals.** The source gives one point age per species (42 distinct values, 486–519 Ma), derived from biozone correlation to GTS 2012. The `gts2012_stage` and `biozone` columns would allow stage-level intervals; none were derived here.
- **BEAST2 character count.** The BEAST2 executables carry 106 discrete characters against 107 in the NEXUS and the BEAST 1.8 files; which character was dropped was not established.
- **BEAST2 ρ.** In the BEAST2 XMLs the `rho` input of `SABirthDeathModel` is tied to `samplingProportion` (an XML comment says so deliberately), although there are no extant taxa.
- **Inapplicable codings.** The matrix uses `-` extensively (1835 cells). Mk software treats `-` as missing; the mean missing-or-inapplicable fraction per taxon is 36%.
- **Structure.** All tips are fossils, so the extant-sampling part of the model is uninformative. Meristic and continuous traits fall outside Mk.
