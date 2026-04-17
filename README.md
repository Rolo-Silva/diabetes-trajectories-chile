# Diabetes Coverage Trajectories — R Project

Latent class mixed model (LCMM) analysis of diabetic retinopathy screening
(DRS) and glycaemic control (HbA1c < 7%) coverage trajectories across Chilean
municipalities, 2011–2023.

Built with [`targets`](https://books.ropensci.org/targets/) for reproducibility.

---

## Project structure

```
trajectories/
├── _targets.R                  # Pipeline definition (run with tar_make())
├── R/
│   └── Utils.R                 # All pipeline functions
├── scripts/
│   └── fit_models.R            # Standalone model fitting (~4 days, run once)
├── data/
│   ├── raw/
│   │   ├── SerieP2011.txt      # Annual series files — see section below
│   │   ├── SerieP2012.txt
│   │   ├── ...
│   │   ├── SerieP2023.txt
│   │   ├── datosf_2017.xlsx    # 2017 patch file (provided separately)
│   │   └── DPA2018.xls         # Municipality lookup table (provided separately)
│   └── models/                 # Individual model .rds files (created by fit_models.R)
├── all_named_models.rds        # Combined model list (created by fit_models.R)
└── _targets/                   # targets cache (auto-generated, do not edit)
```

---

## Input files you need before running

### 1. Annual series files — `SerieP20XX.txt` (2011–2023)

These are the primary care activity records from DEIS (Chilean Ministry of Health).

- **Source:** https://deis.minsal.cl/ → *Estadísticas* → *Serie de Prestaciones*
- **Format:** semicolon-delimited `.txt`, one file per year
- **Size:** several hundred MB per year — large files
- **Note:** Some years may require a formal data request to DEIS if not available
  for direct download from the portal.

Place all 13 files in `data/raw/` with the exact names: `SerieP2011.txt` through
`SerieP2023.txt`.

> ⚠️ **2011–2013 only:** These files have an `id_establecimiento` column with
> hyphens that requires special handling. This is already managed in the pipeline.

### 2. `datosf_2017.xlsx` — 2017 patch file

The 2017 series file is missing certain diabetes records. This Excel file
provides the missing December 2017 data for the affected establishments.

- **How to get it:** Contact the project author — this file can be shared directly.
- **Place it in:** `data/raw/datosf_2017.xlsx`

### 3. `DPA2018.xls` — Municipality lookup table

The División Político-Administrativa (DPA) lookup table, version 2018. Used to
harmonise municipality codes (`id_comuna`) across the study period, including
the creation of Ñuble region in 2018.

- **How to get it:** Contact the project author — this file can be shared directly.
  It is also available from the [INE (Instituto Nacional de Estadísticas)](https://www.ine.gob.cl/).
- **Place it in:** `data/raw/DPA2018.xls`

---

## How to run

### Step 0 — Update file paths

Open `_targets.R` and set `base_dir` to wherever you placed the `data/raw/`
folder on your machine:

```r
base_dir <- "/your/path/to/data/raw"
```

### Step 1 — Install R packages

```r
install.packages(c(
  "targets", "tarchetypes",
  "tidyverse", "janitor", "readxl",
  "lcmm", "LCTMtoolkit",
  "future", "future.apply"
))
```

### Step 2 — Build the data pipeline

This runs everything up to and including `coverage_data`. Fast (~minutes).

```r
targets::tar_make()
```

You can inspect what will run first:

```r
targets::tar_visnetwork()   # visual dependency graph
targets::tar_manifest()     # text summary
```

### Step 3 — Fit the models (manual, ~4 days)

> Skip this step if you have been provided with `all_named_models.rds` directly.

```r
source("scripts/fit_models.R")
```

This script:
- Loads `coverage_data` from the targets cache (`tar_load(coverage_data)`)
- Fits 140 `hlme` models in parallel (10 structures × 2 outcomes × 7 class solutions)
- Saves each model individually to `data/models/` as it completes (resumable if interrupted)
- Writes the combined `all_named_models.rds` to the project root

### Step 4 — Run the pipeline again to pick up the models

Once `all_named_models.rds` exists, re-run the pipeline to build the adequacy table:

```r
targets::tar_make()
```

targets tracks `all_named_models.rds` as a file target. It will automatically
detect if the file has changed and rerun downstream targets.

---

## Sharing results / running in a new environment

If you want a colleague to reproduce the analysis without re-fitting all models
(the most common scenario), share these files:

| File | Size | Notes |
|---|---|---|
| `all_named_models.rds` | ~1–2 GB | Pre-fitted models — skip `fit_models.R` |
| `datosf_2017.xlsx` | Small | 2017 patch — contact project author |
| `DPA2018.xls` | Small | Municipality lookup — contact project author |
| `SerieP20XX.txt` (×13) | Very large | Download from DEIS or request from project author |

Your colleague then:
1. Clones or copies the project folder
2. Places all input files in `data/raw/`
3. Places `all_named_models.rds` in the project root
4. Updates `base_dir` in `_targets.R`
5. Runs `targets::tar_make()` — targets will build everything except the models,
   then load the pre-fitted `all_named_models.rds` directly

---

## Pipeline overview

```
SerieP2011 ─┐
SerieP2012  │
...         ├──► all_series ──► all_series_diabetes ──► coverage_data
SerieP2023  │                        │
datosf_2017 │                    dpa2018
SerieP2017 ─┘

                              [fit_models.R — run manually]
                                        │
                              all_named_models.rds  ◄── targets tracks this file
                                        │
                              all_named_models ──► model_adequacy_table
```

---

## Key variables

| Variable | Description |
|---|---|
| `drsc` | Diabetic retinopathy screening coverage (`dm_fo / dm`) |
| `dgcc` | Glycaemic control coverage (`dm_hg_menor7 / dm`) |
| `dm` | Number of diabetes patients enrolled in primary care |
| `year` | Years since 2011 (0–12) |
| `id` | Municipality numeric ID (`comuna2`) |

Municipalities in the lowest quintile of diabetes caseload (Q1) are excluded
from modelling to avoid unstable coverage estimates in very small practices.

---

## Notes

- `plyr` is **not** used in this project. The 2017 patch is combined with
  `dplyr::bind_rows()`. If you load `plyr` elsewhere in your session, load it
  *before* `dplyr` to avoid masking issues.
- The quintile cut for Q1 exclusion is computed on the full pooled 2011–2023
  dataset. A municipality is excluded entirely if its diabetes caseload falls
  in Q1 in the pooled distribution (not year by year).
- Model fitting uses `gridsearch(rep = 20)` for all multi-class solutions.
  `set.seed(1234)` is set both globally and inside each parallel worker.

---

## Contact

For access to `datosf_2017.xlsx`, `DPA2018.xls`, or `all_named_models.rds`,
contact the project author directly.
