# Penguins: Gavryushkina et al. 2017

Pilot dataset for the CET empirical tests.

- **Study:** Gavryushkina A., Heath T.A., Ksepka D.T., Stadler T., Welch D., Drummond A.J. 2017. Bayesian total-evidence dating reveals the recent crown radiation of penguins. *Systematic Biology* 66:57–73. doi:10.1093/sysbio/syw060
- **Data:** Dryad doi:10.5061/dryad.44pf8 (version 1). Morphology derives from the Ksepka et al. 2012 matrix.
- **Licence:** CC0-1.0 (Dryad). The processed files here are released under the same terms.

## Files

| File | Content |
|---|---|
| `penguins_morph.nex` | 55 taxa × 202 morphological characters (the characters analysed), with character and state labels; `CHARSET statesK` gives the Mkv partition (K assumed states) of each character |
| `penguins_characters.csv` | Per character: column in the XML matrix, states assumed by the model, states observed among the 55 penguins, label |
| `penguins_dna.nex` | 55 taxa × 8145 sites; `CHARSET`s per locus as named in the XML (fossils are all-missing) |
| `penguins_tips.csv` | Taxon, extant/fossil, `min_age`/`max_age` (Ma), source of the age, note |
| `penguins_constraints.csv` | Empty: the source analysis imposes no topological constraints and no outgroup |
| `curated/si-table1.csv` | SI Table 1 age intervals, transcribed by hand from `raw/SI.pdf` p. 7; used only as a cross-check |
| `raw/` | Byte-exact Dryad downloads; checksums in `raw/MANIFEST.csv` |

## Raw downloads

| File | Bytes | MD5 | SHA-256 | Committed | Source |
|---|---|---|---|---|---|
| `SI.pdf` | 317120 | `0497a483b8da555c4f70b43566ff6efd` | `89538902d13c60d177cca3272681ce2193bf51fb5e5172ef164d3fc50e58bb76` | yes | doi:10.5061/dryad.44pf8 |
| `xmls.zip` | 372935 | `651d83e5b74f7ca09234fb78da0accf8` | `c3038de2b7848394b30e8f9a8a4cda1f50bc54d85bc1055d186de2e696a3a44a` | yes | doi:10.5061/dryad.44pf8 |

## Transformations (all in `../build.R`)

1. Unzip `xmls.zip`; take the main total-evidence analysis `xmls/8+dna_R.xml` (analysis "8+dna R", SI Table 3; Fig. 4 of the paper).
2. Morphology: read `data[@id='penguins_morph']`. The XML has 252 columns: 245 matrix characters, then 7 dummy constant columns (states 0–6) that BEAST2 uses for the Mkv ascertainment correction. Drop the dummies (checked to be exactly 0–6 in every taxon).
3. Keep the 202 columns that appear in the six `FilteredAlignment` morphology partitions (`penguins2`–`penguins7`, i.e. 2–7 states). The other 43 are checked to be constant among the 55 penguins, matching the paper (outgroups removed, then constant characters dropped). Column order is preserved; `penguins_characters.csv` maps back to XML columns.
4. Polymorphic cells are written compactly (`{01}`); no recoding.
5. DNA: read `data[@id='penguins']`; `CHARSET`s from the `penguins_<locus>` filters.
6. Tip ages: fossils are the 36 taxa with a `SamplingDate` element; `min_age`/`max_age` are its `lower`/`upper`. Extant taxa are 0. The script checks that every other XML in the archive (except the crown-only analysis, which uses a 6-fossil subset) carries identical intervals.
7. Cross-check against SI Table 1; differences are recorded in the `note` column.

## Source analysis settings (8+dna_R.xml)

| Component | Setting |
|---|---|
| Tree prior | Fossilized birth–death with sampled ancestors (`SABirthDeathModel`), conditioned on ρ-sampling; ρ = 1 (all 19 extant species sampled), removal probability r = 0 |
| FBD parameterization | Net diversification d ~ LogNormal(M = −3.5, S = 1.5); turnover ν ~ U(0, 1); sampling proportion s ~ U(0, 1); origin ~ U(0, 160) Ma |
| Fossil ages | Uniform on the interval in `penguins_tips.csv` (sampled by BEAST2) |
| Morphology | Lewis Mkv, partitioned by number of states (2–7, numbers of states taken from the outgroup-inclusive Ksepka et al. 2012 matrix), Γ4 rate variation shared across partitions, shape ~ U(0, 10); all characters unordered (34 ordered in Ksepka et al. 2012 were run unordered) |
| Morphological clock | UCLN; mean ~ LogNormal(−5.5, 2); SD ~ Gamma(α = 0.5396, β = 0.3819) |
| Molecular model | GTR+Γ4 per locus (5 loci), unlinked |
| Molecular clock | UCLN shared across loci with relative rates; mean ~ LogNormal(−3.5, 1); SD ~ Gamma(0.5396, 0.3819). Analysis "8+dna S" uses a strict molecular clock instead |
| Topological constraints | None |
| MCMC | chainLength = 10⁸ |

Other analyses in `xmls.zip`: `1`–`8` (morphology only, the model comparison of Table 1), `1_LN281`/`1_U011`/`1_U1000` (diversification-rate prior sensitivity), `8+dna_S` (strict molecular clock), `8BD+dna_R` (FBD without sampled ancestors), `8+dna_R_crown_only` (19 extant + 6 crown fossils, Mk not Mkv).

## Oddities

- **Locus labels.** The XML partition lengths are RAG-1 2682, `cyt_b` 1143, `coi` 1551, `12srna` 1105, `16srna` 1664 sites. The SI table lists 12S 1143, 16S 1551, COI 1105, cyt b 1664. The XML lengths fit the genes (cyt b is 1143 bp), so the XML names are kept and the SI table appears to be mislabelled.
- **Spheniscus urbinai.** The XML and the SI text use 5.7–9.63 Ma; SI Table 1 prints 5.6–9.63 Ma. The XML value is used.
- **Perudyptes devriesi.** The SI text gives 38–46 Ma (≈42 Ma ± 4 Ma), but SI Table 1 and the XML use 38–42 Ma. The XML value is used.
- **Initial dates.** The XML's initial tip dates lie outside the sampling interval for *Archaeospheniscus lopdelli* (25.24 Ma), *Marplesornis novaezealandiae* (5.10 Ma) and *Palaeospheniscus patagonicus* (14.58 Ma). These are MCMC starting values only and are not used here.
- **Spelling.** *Palaeeudyptes klekowskii* is spelled "klekowski" in the SI text; the XML spelling is kept. Burnside "Palaeeudyptes" is an unnamed species, labelled `Burnside_Palaeudyptes` in the XML.
- **Structure.** The crown is extant-only and the stem is fossil-only. Fossils have no DNA and are on average 84% missing in the morphological matrix; extant taxa are 5% missing.
