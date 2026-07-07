# Diabetes-care coverage under Chile’s Explicit Health Guarantees, 2011 to 2023: A longitudinal municipal analysis of equity and adequacy

This repository contains the code to reproduce the latent class mixed model (LCMM) analysis of diabetic retinopathy screening coverage (DRSC) and effective glycaemic control coverage (EGCC) trajectories across Chilean municipalities, 2011–2023.

The analysis uses [`targets`](https://docs.ropensci.org/targets/) to manage the data pipeline and ensure reproducibility.

## Variable naming note

The internal codebase uses `dgcc` for the effective glycaemic control coverage variable (proportion of HbA1c < 7% among monitored diabetes patients). In the manuscript and supplementary tables this variable is reported as **EGCC**. The two refer to the same quantity; the renaming is purely presentational.

| Variable in code | Variable in manuscript | Description |
|---|---|---|
| `dgcc` | EGCC | Effective glycaemic control coverage |
| `drsc` | DRSC | Diabetic retinopathy screening coverage |

## Requirements

### R version

R ≥ 4.2.0

### R packages

```r
install.packages(c(
  # Pipeline framework
  "targets", "tarchetypes",
  # Data manipulation
  "dplyr", "tidyr", "readr", "stringr", "stringi", "tibble", "rlang",
  "lubridate", "forcats", "purrr", "janitor", "readxl",
  # Trajectory modelling
  "lcmm",
  # Statistical analysis
  "nnet", "broom", "lme4", "lmerTest", "MASS", "psych",
  # Visualisation
  "ggplot2", "patchwork", "cowplot", "ggsci", "scales", "grid", "glue"
))
```

`LCTMtools` is not on CRAN. Install from GitHub:

```r
remotes::install_github("hlennon/LCTMtools")
```

## Input data

Five input files are required. None are included in this repository due to file size or access restrictions.

### 1. Annual primary care activity records — `SerieP20XX.txt` (2011–2023)

Thirteen semicolon-delimited files, one per year, containing monthly primary care activity records from Chilean public health facilities.

- **Source:** [DEIS — Ministerio de Salud de Chile](https://deis.minsal.cl/) → Estadísticas → Serie de Prestaciones.
- Some years may require a formal data request to DEIS if not available for direct download.
- Place all 13 files in `data/raw/` with the exact names `SerieP2011.txt` through `SerieP2023.txt`.

### 2. `datosf_2017.xlsx` — 2017 supplementary patch

The 2017 series file from DEIS is missing December records for a subset of establishments. This file provides the missing data.

- **How to obtain:** Contact the corresponding author.
- Place in `data/raw/datosf_2017.xlsx`.

### 3. `DPA2018.xls` — Municipality lookup table

The División Político-Administrativa (DPA) 2018 lookup table, used to harmonise municipality codes across the study period (including the creation of the Ñuble region in 2018).

- **Source:** INE — Instituto Nacional de Estadísticas, or contact the corresponding author.
- Place in `data/raw/DPA2018.xls`.

### 4. `indice_desarrollo_comunal_pdf_excel.xlsx` — Index of Communal Development

Used as a deprivation covariate in downstream regression analyses.

- **How to obtain:** Contact the corresponding author.
- Place in `data/raw/indice_desarrollo_comunal_pdf_excel.xlsx`.

### 5. `comunas_pndr.xlsx` — National Rural Development Policy classification

Used to assign each municipality to urban, mixed, or rural categories.

- **How to obtain:** Contact the corresponding author.
- Place in `data/raw/comunas_pndr.xlsx`.

## How to reproduce the analysis

### Step 1 — Clone the repository

```bash
git clone https://github.com/Rolo-Silva/diabetes-trajectories-chile.git
```

Open `diabetes-trajectories-chile.Rproj` in RStudio.

### Step 2 — Place input files

```
data/
└── raw/
    ├── SerieP2011.txt
    ├── SerieP2012.txt
    ├── ...
    ├── SerieP2023.txt
    ├── datosf_2017.xlsx
    ├── DPA2018.xls
    ├── indice_desarrollo_comunal_pdf_excel.xlsx
    └── comunas_pndr.xlsx
```

### Step 3 — Set the data path

Open `_targets.R` and update `base_dir` to the absolute path of your `data/raw/` folder:

```r
base_dir <- "/your/path/to/data/raw"
```

### Step 4 — Run the upstream data pipeline

This builds all targets from raw data through to `coverage_data`. Takes a few minutes.

```r
targets::tar_make()
```

You can inspect the pipeline before running:

```r
targets::tar_visnetwork()
```

### Step 5 — Fit the trajectory models

> ⚠️ This step takes approximately **4 days on a multi-core machine**. If you have been provided with `all_named_models.rds` by the corresponding author, place it in the project root and skip to Step 6.

```r
source("scripts/fit_models.R")
```

This script loads `coverage_data` from the `targets` cache, fits 140 `hlme` models in parallel (10 covariance structures × 2 outcomes × 7 class solutions), and saves the results as `all_named_models.rds` in the project root. Individual models are saved to `data/models/` as they complete — if the script is interrupted, it can be resumed without refitting completed models.

### Step 6 — Build the full pipeline

Once `all_named_models.rds` is in the project root, re-run the pipeline:

```r
targets::tar_make()
```

`targets` will build all downstream analyses automatically:

- Model adequacy table and selection diagnostics (Figures 2–6)
- Trajectory plots (Figure 1B, Figure 2B, SF3)
- Sensitivity validation: subsampling (Figure 7) and observation perturbation (Figure 8)
- Enrichment with municipal covariates (IDC, PNDR)
- Cross-classification of EGCC × DRSC (SF9)
- Multinomial regression of trajectory class membership (ST7)
- Mixed-effects regression with year-by-covariate interactions (ST8)

### Step 7 — Adjust sensitivity validation seeds (optional)

By default, the sensitivity analyses run with `sensitivity_seeds = 1L` (one random seed) for fast pipeline validation. For the production analysis reported in the manuscript, change the corresponding line in `_targets.R`:

```r
tar_target(sensitivity_seeds, 1:20)
```

This re-runs Stage 1 (subsampling) and Stage 2 (observation perturbation) across 20 seeds (approximately **21 hours total**).

## Repository structure

```
.
├── _targets.R                      # Pipeline definition
├── R/
│   ├── 00_data.R                   # Data assembly and coverage construction
│   ├── 01_models.R                 # Model adequacy and diagnostic figures
│   ├── 02_trajectories.R           # Manuscript trajectory figures
│   ├── 03_sensitivity.R            # Sensitivity validation (Stage 1 + Stage 2)
│   ├── 04_enrichment.R             # IDC and PNDR integration, class membership
│   ├── 05_multinom.R               # Multinomial regression with Rubin's rules
│   ├── 06_cross_class.R            # EGCC × DRSC cross-classification
│   └── 07_mixed_effects.R          # Mixed-effects linear analysis
├── scripts/
│   └── fit_models.R                # Standalone model fitting script
├── data/
│   ├── raw/                        # Input files (not tracked by git)
│   └── models/                     # Individual model .rds files (not tracked by git)
├── outputs/
│   └── figures/                    # Manuscript figures (generated by pipeline)
├── all_named_models.rds            # Fitted LCMM objects (not tracked by git)
└── README.md
```

## Key analytical outputs

After a complete pipeline run, the following targets are available via `targets::tar_load()`:

| Target | Description |
|---|---|
| `coverage_data` | 299 municipalities × annual observations, ready for LCMM |
| `coverage_demographics_class` | Enriched analytical data with covariates and class assignments |
| `model_adequacy_table` | LCMM diagnostics for all 140 candidate models |
| `subsampling_validation` | Stage 1 sensitivity: trajectory stability under municipality subsampling |
| `perturbation_validation` | Stage 2 sensitivity: trajectory stability under observation removal |
| `multinom_validation` | 8 multinomial tables (4 EGCC reference classes + 4 DRSC reference classes) |
| `cross_classification` | EGCC × DRSC contingency table, chi-square, Cramér's V, heatmap |
| `mixed_effects_analysis` | Main and interaction lmer models, LRT comparison, Supplementary Table 7 |

## Notes on reproducibility

- All analysis code uses `dplyr` and the native R pipe `|>`. Do not load `plyr` in the same session as it masks key `dplyr` functions.
- LCMM fitting uses `gridsearch(rep = 20)` with `set.seed(1234)` set both globally and within each parallel worker.
- DRSC values exceeding 1 (a small fraction of observations) are capped at 1 prior to modelling.
- Mixed-effects models require `lmerTest::lmer` (not `lme4::lmer`) to obtain Satterthwaite p-values.
- Cross-classification chi-square uses 10,000 Monte Carlo replicates and is stochastic; expected variation in `p-value` of approximately ±0.01 across runs.

## Citation

If you use this code or build upon this analysis, please cite:

> Silva-Jorquera R, Zett C, Haghparast-Bidgoli H, Nderitu P, Olvera-Barrios A, Warwick AN, Allel K. Trajectories of glycaemic control and diabetic retinopathy screening coverage across Chilean municipalities, 2011–2023: a latent trajectory analysis. *Journal of Public Health (Oxford)*. Manuscript submitted for publication.

Archived analysis pipeline: [Zenodo DOI: 10.5281/zenodo.20488260](https://doi.org/10.5281/zenodo.20488260).

## Contact

For access to `datosf_2017.xlsx`, `DPA2018.xls`, `indice_desarrollo_comunal_pdf_excel.xlsx`, `comunas_pndr.xlsx`, or `all_named_models.rds`, or for questions about the code, please open a [GitHub Issue](https://github.com/Rolo-Silva/diabetes-trajectories-chile/issues) or contact the corresponding author at `rolando.silva@ucl.ac.uk`.