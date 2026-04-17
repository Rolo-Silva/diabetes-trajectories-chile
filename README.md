# Diabetes Coverage Trajectories in Chilean Primary Care (2011–2023)

This repository contains the code to reproduce the latent class mixed model
(LCMM) analysis of diabetic retinopathy screening (DRS) and glycaemic control
(HbA1c < 7%) coverage trajectories across Chilean municipalities between 2011
and 2023.

The analysis uses [`targets`](https://books.ropensci.org/targets/) to manage
the data pipeline and ensure reproducibility.

---

## Requirements

### R version
R ≥ 4.2.0

### R packages

```r
install.packages(c(
  "targets", "tarchetypes",
  "tidyverse", "janitor", "readxl",
  "lcmm", "LCTMtoolkit",
  "future", "future.apply"
))
```

---

## Input data

Three types of input files are required. **None of them are included in this
repository** due to file size or access restrictions.

### 1. Annual primary care activity records — `SerieP20XX.txt` (2011–2023)

Thirteen semicolon-delimited files, one per year, containing monthly primary
care activity records from Chilean public health facilities.

- **Source:** [DEIS — Ministerio de Salud de Chile](https://deis.minsal.cl/)
  → *Estadísticas* → *Serie de Prestaciones*
- Some years may require a formal data request to DEIS if not available for
  direct download.
- Place all 13 files in `data/raw/` with the exact names `SerieP2011.txt`
  through `SerieP2023.txt`.

### 2. `datosf_2017.xlsx` — 2017 supplementary patch

The 2017 series file from DEIS is missing December records for a subset of
establishments. This file provides the missing data.

- **How to obtain:** Contact the corresponding author.
- Place in `data/raw/datosf_2017.xlsx`.

### 3. `DPA2018.xls` — Municipality lookup table

The División Político-Administrativa (DPA) 2018 lookup table, used to
harmonise municipality codes across the study period, including the creation
of Ñuble region in 2018.

- **Source:** [INE — Instituto Nacional de Estadísticas](https://www.ine.gob.cl/),
  or contact the corresponding author.
- Place in `data/raw/DPA2018.xls`.

---

## How to reproduce the analysis

### Step 1 — Clone the repository and open the R project

```bash
git clone https://github.com/YOUR_USERNAME/diabetes-trajectories-chile.git
```

Open `diabetes-trajectories-chile.Rproj` in RStudio.

### Step 2 — Place input files

Create the required folders and place all input files as described above:

```
data/
└── raw/
    ├── SerieP2011.txt
    ├── SerieP2012.txt
    ├── ...
    ├── SerieP2023.txt
    ├── datosf_2017.xlsx
    └── DPA2018.xls
```

### Step 3 — Set the data path

Open `_targets.R` and update `base_dir` to the absolute path of your `data/raw/`
folder:

```r
base_dir <- "/your/path/to/data/raw"
```

### Step 4 — Run the data pipeline

This builds all targets from raw data through to `coverage_data`. It takes a
few minutes.

```r
targets::tar_make()
```

You can inspect the pipeline before running:

```r
targets::tar_visnetwork()  # visual dependency graph
```

### Step 5 — Fit the trajectory models

> ⚠️ This step takes approximately **4 days** on a multi-core machine. If you
> have been provided with `all_named_models.rds` by the corresponding author,
> place it in the project root and skip to Step 6.

```r
source("scripts/fit_models.R")
```

This script loads `coverage_data` from the targets cache, fits 140 `hlme`
models in parallel (10 covariance structures × 2 outcomes × 7 class solutions),
and saves the results as `all_named_models.rds` in the project root. Individual
models are saved to `data/models/` as they complete — if the script is
interrupted, it can be resumed without refitting completed models.

### Step 6 — Build the model adequacy table

Once `all_named_models.rds` is in the project root, re-run the pipeline:

```r
targets::tar_make()
```

targets detects the new file and builds the `model_adequacy_table` target
automatically.

---

## Repository structure

```
├── _targets.R               # Pipeline definition
├── R/
│   └── Utils.R              # All pipeline functions
├── scripts/
│   └── fit_models.R         # Standalone model fitting script
├── data/
│   ├── raw/                 # Input files (not tracked by git)
│   └── models/              # Individual model .rds files (not tracked by git)
└── README.md
```

---

## Key variables

| Variable | Description |
|---|---|
| `drsc` | Diabetic retinopathy screening coverage (`dm_fo / dm`) |
| `dgcc` | Glycaemic control coverage — HbA1c < 7% (`dm_hg_menor7 / dm`) |
| `dm` | Diabetes patients enrolled in primary care per municipality per year |
| `year` | Years since 2011 (0 = 2011, 12 = 2023) |
| `id` | Municipality numeric identifier |

Municipalities in the lowest quintile of diabetes caseload (Q1, pooled across
all years) are excluded from modelling to avoid unstable coverage estimates in
very small practices.

---

## Notes

- All analysis code uses `dplyr` exclusively. Do not load `plyr` in the same
  session as it masks key `dplyr` functions.
- Model fitting uses `gridsearch(rep = 20)` for all multi-class solutions with
  `set.seed(1234)` set both globally and within each parallel worker.
- DRS coverage values greater than 1 are capped at 1 prior to modelling.

---

## Contact

For access to `datosf_2017.xlsx`, `DPA2018.xls`, or `all_named_models.rds`,
or for any questions about the code, please open a
[GitHub Issue](../../issues) or contact the corresponding author.
