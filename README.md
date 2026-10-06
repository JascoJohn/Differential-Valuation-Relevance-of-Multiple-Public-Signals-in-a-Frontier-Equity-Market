# Replication package: Differential Valuation Relevance of Public Signals in DSE-Listed Commercial Banks

**Authors:** Jasco Francis John (Ardhi University); Minja Ladislaus
**Package version:** 1.0 (October 2026)
**Manuscript:** *[title, journal and DOI to be added]*
**Dataset:** Tanzania_DSE_AnalyticalDataset_v1.3.1_LOCKED (included)

This package reproduces every table, test, diagnostic and robustness check reported in the manuscript from the locked analytical dataset, using R. One script (`R/master.R`) rebuilds all results from scratch and verifies them against the reference outputs supplied in `reference/`.

---

## 1. Contents

```
replication_package/
├── README.md                       this file
├── MANIFEST.csv                    path, size and SHA-256 of every file in the package
├── <project>.Rproj                 RStudio project file (sets the working directory)
├── .Rprofile, renv/activate.R,     renv project files (restore the pinned package library)
│   renv/settings.json
├── renv.lock                       exact versions of R and all packages
├── data/raw/
│   └── Tanzania_DSE_AnalyticalDataset_v1.3.1_LOCKED.xlsx
├── R/
│   ├── R00_setup.R                 original environment setup (documentation only; see Section 3)
│   ├── R01_import_verify.R         import and verify the locked dataset
│   ├── R02_construct.R             analysis variables and core sample
│   ├── R03_descriptives_diagnostics.R
│   ├── R04_H1.R                    H1: valuation relevance (Eqs. 2a–2c, 3)
│   ├── R05_H2.R                    H2: differential valuation relevance (Eqs. 4–6)
│   ├── R06_robustness.R            robustness, dominant-bank exclusions, PCSE, CR2
│   ├── R07_power.R                 Monte Carlo size and power
│   ├── R09_design_simulation.R     exploratory design simulation (Section 8; not run by master.R)
│   ├── wcb_functions.R             restricted wild cluster bootstrap (used by R04–R07)
│   └── master.R                    runs R01–R07 and verifies all outputs
├── logs/                           logs from the official run (R01–R07, master) and the exploratory R09 run
├── reference/                      expected output tables (T1–T9) from the official run
└── output/
    ├── OUTPUT_HASHES.csv           SHA-256 of the derived data and output tables
    └── exploratory/
        └── design_simulation.csv   results of the exploratory design simulation
```

Derived data (`data/derived/`) and output tables (`output/tables/`) are not included; `master.R` regenerates them.

## 2. Data

| Item | Value |
|---|---|
| File | `Tanzania_DSE_AnalyticalDataset_v1.3.1_LOCKED.xlsx` |
| Size | 134,968 bytes |
| File SHA-256 | `1fffff6c33475ce40dc3a8528ef4846fe4e6ae53352666d93a3f5537ed9eb289` |
| Panel-content SHA-256 (recorded in the workbook; `compute_panel_checksum.py`) | `92fc00fdf74a8333a6e823429f44e57af2c6943eb726d543d590df6e50681043` |
| Panel | 5 banks × 10 years (2016–2025), 22 variables; sheet `ANALYTICAL_PANEL` |
| Panel locked | 11 September 2026 (v1.3); documentation revision v1.3.1, 5 October 2026 |

v1.3.1 is a documentation-only revision of v1.3-LOCKED: no panel value changed, and the panel-content checksum is unchanged (see the workbook's `VERSION_LOG` and `DECISION_LOG` D-20). The workbook's own sheets document variable definitions, provenance of every cell, construction decisions and the version history. Its REPRODUCIBILITY sheet describes the analysis environment planned at the time of the lock; the analysis reported here is implemented entirely in R as documented in this README. The companion Research Archive, construction methodology and build scripts referred to in the workbook are deposited with the dataset record *[dataset DOI to be added]*.

## 3. Software environment

| Component | Version |
|---|---|
| R | 4.6.1 (Windows 11 x64 in the official run) |
| renv | 1.3.0 |
| fixest | 0.14.2 |
| clubSandwich | 0.7.0 |
| plm | 2.6-7 |
| sandwich | 3.1-3 |
| readxl | 1.5.0.1 |
| digest | 0.6.39 |

`renv.lock` pins every package and the R version. Restore the environment with `renv::restore()` rather than running `R00_setup.R`, which installs the newest available versions and is included only to document how the original environment was created. The environment of the official run is recorded in renv.lock; no log of the original R00 run is included.

`fwildclusterboot` (0.14.3) is recorded in `renv.lock` but is not called by the analysis scripts; the wild cluster bootstrap is implemented directly in `wcb_functions.R` (Section 5).

## 4. How to reproduce the results

1. Install R 4.6.1 (and, optionally, RStudio).
2. Open the project: double-click the `.Rproj` file, or start R with the package folder as the working directory.
3. Restore the packages (first time only):
   ```r
   renv::restore()
   ```
4. Run the full analysis:
   ```r
   source("R/master.R")
   ```

`master.R` deletes any derived data and output tables, runs R01–R07 in order (each in a clean environment, each writing its own log to `logs/`), records the SHA-256 of all outputs in `output/OUTPUT_HASHES.csv`, and compares every regenerated table with `reference/`. A successful run ends with:

```
ALL OUTPUTS REPRODUCED.
```

Run time is about 10–15 minutes on a standard laptop, mostly the confidence-interval inversion in R04–R05 and the simulations in R07. Results are deterministic: all random draws use R's default generator (Mersenne-Twister) with fixed seeds, so a rerun reproduces every value exactly on the same platform. On a different operating system or processor, coefficients reproduce exactly and bootstrap p-values may differ only at the level of floating-point rounding.

## 5. Specification summary

| Element | Specification |
|---|---|
| Sample | Core estimation sample of 40 bank-years (5 banks, 2017–2025), common to all principal models |
| Timing | All signals at fiscal year t−1; ln(MTB) at the end of year t |
| Dependent variable | ln(market capitalization ÷ total shareholders' equity) |
| Estimator | Two-way fixed effects (bank, year) |
| H1 models | Eqs. 2a–2c and 3: lagged CAR, NPL, ROA/ROE/EPS and FULLREP in percent or native units |
| H2 models | Eqs. 4–6: signals scaled by a common within-bank SD; NPL entered as −NPL; composites (banking strength, financial performance) built and rescaled on the same basis |
| Principal inference | Restricted wild cluster bootstrap, Webb six-point weights, 9,999 replications, clustered by bank; one hypothesis per call with `set.seed(20261004)` before each test; joint tests by bootstrap Wald statistic; confidence intervals by test inversion (reported as unbounded where the confidence set does not close within ±10 clustered SEs) |
| Testing protocol | Pairwise equality tests are interpreted only after a joint rejection; otherwise reported as exploratory |
| Secondary inference | Clustered SEs with t(G−1); CR2 with Satterthwaite / HTZ tests; panel-corrected SEs (Beck–Katz, cross-bank covariance from years in which all banks are observed) |
| Diagnostics | Variance shares by bank, VIF, modified Wald (groupwise heteroskedasticity), Pesaran CD, Wooldridge serial-correlation test, Hausman (Amemiya variance components; Swamy–Arora is not estimable with five banks) |
| Robustness | Size added; LDR added; residualized performance; winsorized signals (5th/95th percentiles); stale-price observations excluded; Eq. 5 with banking-risk signals; leave-one-bank-out; dominant-bank exclusions; PCSE; CR2 |
| Power | Monte Carlo size and power (500 samples per cell; 399 bootstrap draws per sample) for NPL, H2a, H2b and H2c |

## 6. Outputs

| File | Content |
|---|---|
| `T1_descriptives.csv` | Descriptive statistics, core sample |
| `T1b_between_within.csv` | Overall, between-bank and within-bank standard deviations |
| `T1c_within_correlations.csv` | Within-bank correlations |
| `T1d_variance_shares.csv` | Share of each signal's within-bank variation by bank |
| `T2_diagnostics.csv` | Residual diagnostics and Hausman test |
| `T3_H1_baseline.csv` | H1 models: coefficients, clustered and bootstrap inference, bootstrap CIs |
| `T4a_H2_coefficients.csv` | H2 models: coefficients and bootstrap inference |
| `T4b_H2_tests.csv` | H2 equality tests (joint and pairwise) |
| `T5_robustness_bootstrap.csv` | Robustness specifications and leave-one-bank-out |
| `T6_dominant_bank_exclusions.csv` | Dominant-bank exclusions and NPL bank-by-bank |
| `T7a_pcse_eq3.csv`, `T7b_pcse_H2.csv` | Panel-corrected SE results |
| `T8a_cr2_eq3.csv`, `T8b_cr2_H2.csv` | CR2 small-sample results |
| `T9_power.csv` | Monte Carlo size and power |

*[Mapping of these files to the manuscript's table numbers to be added.]*

## 7. Checksums

| File | SHA-256 |
|---|---|
| `data/raw/Tanzania_DSE_AnalyticalDataset_v1.3.1_LOCKED.xlsx` | `1fffff6c33475ce40dc3a8528ef4846fe4e6ae53352666d93a3f5537ed9eb289` |
| `data/derived/panel_imported.csv` (regenerated) | `0e7db0e18ce5c906f7af7eef5e71bccc3b310bae5fc3382170df961bbf35775a` |
| `data/derived/panel_analysis.csv` (regenerated) | `de5b2268e868832710b01c711ae081dc3248d933feab3fef31f9c224f40ec1f6` |
| `T1_descriptives.csv` | `c85155ce1b49f6372ce308f87bd769d015fa2e942f1884aacd16d5c142356b8f` |
| `T1b_between_within.csv` | `486d647bf9e012b644ca148a99ef2ff201152ed001b47ffb3f33d3ffa1a03006` |
| `T1c_within_correlations.csv` | `9d8f8fcf726ac08e2d0d3044d15f9c4ef0336fd8aa4258bfb10ce48229058907` |
| `T1d_variance_shares.csv` | `aa2bbd283a007e1400971b69a9f9aa9a1a92c1947655ba18f18c6a844dbb56ad` |
| `T2_diagnostics.csv` | `900bd822b0d3cb12db4a2ac145c8a440059fca0a9d63672676d6ba7d8a765178` |
| `T3_H1_baseline.csv` | `4db9af63eecff6f770940f50e328c743f668b921ae877cfe28cdab6424f7e3c2` |
| `T4a_H2_coefficients.csv` | `6bf8e23af4cc2cbd957a94402d30aaf0c7d1be6399a58f81f16a971caeefbd50` |
| `T4b_H2_tests.csv` | `d7c6ffb82b3442f60aa1b5c0ec25f2efdbb842089824063c3a219f0ad6c49f4e` |
| `T5_robustness_bootstrap.csv` | `447cfe2e4bb62475768d9672193a2617952fbd239682a8b7fa88b1a58ab5cea9` |
| `T6_dominant_bank_exclusions.csv` | `12bfd55bd18039e947193ccfe0fe3d99f6fc5c1a2d902b24ff4e8b080da77a88` |
| `T7a_pcse_eq3.csv` | `f52b58ff8f15ec9dba558b576624be48207d85f5366e14cddbb055c1c41b2827` |
| `T7b_pcse_H2.csv` | `827a861534afe3af945eee2a863e2b9ca302cf5c9f714c49b8838571efcabd0d` |
| `T8a_cr2_eq3.csv` | `757df56a42bb9f6021aef29c3cec2f299ec2e6af064fdc4cb428fad37084840e` |
| `T8b_cr2_H2.csv` | `32166cd64ed7b823de9f60c0f0bf072af6ca6a24d27579e9cdf9cee124e6e864` |
| `T9_power.csv` | `222cc7e6f35bfaa0ee145ba26d630fb18fb514de5dfd3f66012467ab83e94b31` |

`MANIFEST.csv` lists the size and SHA-256 of every file in the package. The output-table hashes above hold for the official platform (Windows, R 4.6.1); the numerical comparison in `master.R` is the platform-independent check.

## 8. Notes

- The H2 tests use a common within-bank standard deviation so that every bank shares one slope and the comparison is in units of typical deviations from each bank's own level.
- Within-bank variation in CAR comes almost entirely from one bank (T1d), which is why its bootstrap confidence intervals are unbounded.
- The panel-corrected SEs rely on the four years in which all five banks are observed (2021, 2023, 2024, 2025).
- With four clusters, leave-one-bank-out bootstrap p-values are indicative only.
- The panel was locked on 11 September 2026. Some earlier VERSION_LOG entries in the workbook carry the date 7 September 2026; this inconsistency is recorded in the workbook rather than altered retrospectively.
- `R09_design_simulation.R` is an exploratory analysis of how test power would change with 5–25 banks. It is not part of the verified pipeline: `master.R` does not run it, and its output (`output/exploratory/design_simulation.csv`) is not compared with `reference/`. It runs in about 15 minutes with `source("R/R09_design_simulation.R")` after `master.R`.
- AI-assisted tooling (Claude, Anthropic) was used during dataset construction and in writing the analysis code; all data, specifications and results were reviewed and decided by the authors.

## 9. Licence and citation

- Code: MIT License
- Data: Creative Commons Attribution 4.0 International (CC BY 4.0)

Please cite the article *[reference]* and this package *[Zenodo DOI]*.
