# Empirical tip-dated datasets

Curated inputs for the CET empirical comparison (phase 4): a 2 × 2 (β, γ) model family plus a UCLN comparator. The three flagship datasets were chosen in the dataset-tiers memo (penguins as the pilot, Hymenoptera as the main total-evidence test, trilobites as the Cambrian test).

| Dataset | Source study | Data DOI | Taxa | Fossils | Extant | Morph. chars analysed | Molecular sites | Tip ages (Ma) | Age format |
|---|---|---|---|---|---|---|---|---|---|
| [`penguins/`](penguins/README.md) | Gavryushkina et al. 2017 *Syst. Biol.* | 10.5061/dryad.44pf8 | 55 | 36 | 19 | 202 | 8145 (5 loci) | 0; fossils 2.5–61.6 | intervals |
| [`hymenoptera/`](hymenoptera/README.md) | Ronquist et al. 2012 *Syst. Biol.*; intervals from Zhang et al. 2016 | 10.5061/dryad.j2r64 | 113 | 45 | 68 (8 outgroups) | 343 of 353 | 5096 (7 markers) | 0; fossils 80–242 | intervals (and 2012 point ages) |
| [`trilobites/`](trilobites/README.md) | Paterson, Edgecombe & Lee 2019 *PNAS* | 10.5061/dryad.v7q827k | 107 | 107 | 0 | 107 (+2 meristic, 6 continuous) | none | 486–519 | points |

All three Dryad packages are CC0-1.0.

## Layout

Each dataset directory contains:

- `raw/`: byte-exact downloads, with `raw/MANIFEST.csv` giving bytes, MD5 (as reported by Dryad), SHA-256, licence and URL. A nested `.gitattributes` keeps git from normalizing their line endings.
- `curated/`: values transcribed from the papers, where the deposited files do not carry them.
- `<name>_morph.nex`: morphological matrix (NEXUS); `<name>_dna.nex`: molecular partitions, kept separate.
- `<name>_tips.csv`: one row per terminal. Columns `taxon`, `status` (`extant`/`fossil`), `min_age`, `max_age` (Ma) and `age_source`, plus dataset-specific columns.
- `<name>_constraints.csv`: topological constraints and outgroups used by the source study.
- `<name>_characters.csv`: per-character metadata.
- `README.md`: provenance, transformations, source analysis settings, licence and oddities.

`summary.csv` holds the counts in the table above as written by the build script.

## Rebuilding

```sh
Rscript data-raw/empirical/build.R
```

Run from the package root. The script needs `xml2`, `readxl` and `TreeTools`. It also needs `curl` and `jsonlite` if raw files must be downloaded, in which case Dryad requires `DRYAD_CLIENT_ID` and `DRYAD_CLIENT_SECRET` in the environment. It:

1. downloads any committed raw file that is missing and verifies every MD5 against `raw/MANIFEST.csv`;
2. rewrites all processed files;
3. re-reads every processed NEXUS file with `TreeTools::ReadCharacters()` and stops unless its taxa equal those of the tip table.

`data-raw/` is excluded from the package build by `.Rbuildignore`.
