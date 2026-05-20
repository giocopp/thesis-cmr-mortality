# Analysis pipeline

End-to-end R pipeline for the CMR mortality thesis.

## Run

```r
source("analysis/R/run_all.R")
```

This restores the `renv` lockfile and runs every script under
`02_clean/` → `03_build/` → `04_descriptive/` → `05_analysis/` → `06_robustness/`
in alphanumeric order within each section.

## Layout

| Folder            | Purpose                                                           |
|-------------------|-------------------------------------------------------------------|
| `01_download/`    | External-data fetch scripts (ERA5). Other sources are documented. |
| `02_clean/`       | Raw → `data/processed/` cleaning.                                 |
| `03_build/`       | Daily and zone panels built into `analysis/data/`.                |
| `04_descriptive/` | Maps, time-series, distribution plots → `output/figures/`.        |
| `05_analysis/`    | Primary inference cited in the manuscript.                        |
| `06_robustness/`  | Robustness, sensitivity, source-comparison checks.                |

Shared utilities live in `_helpers.R` (sources `_constants.R`).

## Output

Every script writes to `output/figures/{section}/` or `output/tables/{section}/`.
The Quarto manuscript at `paper/thesis.qmd` references these paths directly.

## Raw data

See `01_download/README.md` for acquisition steps and the fixed snapshot dates
used to reproduce the manuscript.
