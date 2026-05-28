# When Rough Seas Become Deadlier: Migrant Mortality and the 2017 Italy–Libya Memorandum

Master's thesis, MDS in Public Policy, Hertie School, Class of 2026.  

Author: Giorgio Coppola. Supervisor: Prof. Dr. Asya Magazinnik.

<strong>Read the thesis</strong>
<a href="paper/thesis.pdf" target="_blank" rel="noopener">in <strong>PDF</strong></a>

## Abstract

The Central Mediterranean is the world's deadliest migration route, yet quantitative evidence on how the post-2017 EU migration policy regime changed the risk faced by migrants attempting the journey from North Africa to Europe remains limited. This thesis asks whether the same sea-state conditions are associated with more recorded migrant deaths on the Central Mediterranean after the 2 February 2017 Italy–Libya Memorandum of Understanding (MoU) than before, using a Central Mediterranean route-level panel from January 2014 to May 2023. In the primary negative-binomial model with month–year fixed effects, a one-metre rise in the lagged five-day mean of significant wave height (SWH) corresponds to a roughly 94% decrease in expected recorded deaths before the MoU but more than doubles them after (pre-MoU slope −2.829, post-MoU shift +3.629); the shift survives controls for crossing volume, boat composition, and conflicts in Libya and Tunisia. Mechanism checks point to search-and-rescue availability: in weeks with higher recorded SAR activity rough seas predict fewer additional deaths, and the post-MoU gap returns toward zero during the Lamorgese partial rollback. The findings reframe assessment of post-2017 sea-border policy: after the regime change the same sea conditions are associated with substantially more recorded deaths, and SAR availability buffers that association.

## Repository structure

```
.
├── analysis/
│   └── R/
│       ├── 01_download/      # raw-data fetch scripts (ERA5; others documented)
│       ├── 02_clean/         # raw → data/processed/
│       ├── 03_build/         # daily and zone panels into analysis/data/
│       ├── 04_descriptive/   # maps, time-series, distributions
│       ├── 05_analysis/      # primary count models, rolling β, mechanism, periods
│       ├── 06_robustness/    # FE, event placebos, push-factor, boat composition
│       ├── _constants.R      # shared constants
│       ├── _helpers.R        # shared utilities
│       ├── README.md         # pipeline notes
│       └── run_all.R         # end-to-end entry point
├── data/
│   ├── raw/                  # ACLED, IMO SAR, UNHCR tracked;
│   │                         # ERA5, UNITED, IOM, Frontex available on request
│   └── processed/            # intermediate RDS files used by analysis scripts
├── output/
│   ├── figures/              # PNGs produced by 04_descriptive and 05_analysis
│   └── tables/               # .tex tables included by paper/thesis.qmd
├── paper/
│   ├── thesis.qmd            # main Quarto document
│   ├── thesis.pdf            # rendered output
│   ├── references.bib        # bibliography
│   ├── apa.csl               # citation style
│   ├── assets/               # title-page logo and similar
│   └── wordcount_filter.lua  # word-count Lua filter
├── renv.lock                 # pinned R-package versions
└── README.md
```

## Data availability

`data/processed/` is tracked and contains all intermediate files the analysis scripts read. Inside `data/raw/`, the ACLED, IMO GISIS SAR boundaries, and UNHCR daily arrivals inputs are tracked; ERA5 sea-state files are not tracked because of their size, and UNITED, IOM MMP, and Frontex PAD-194 are not redistributed in line with each provider's terms. These four sources are available from the author on request: coppola.giorgio99@gmail.com.

## Reproducing the results

Requirements: R 4.5 or newer, the `renv` package, and a working Quarto + LuaLaTeX install for rendering the PDF.

1. Clone the repository and open `Thesis-MDS.Rproj` (or set the repo root as the working directory).
2. Request the ERA5, UNITED, IOM, and Frontex files from the author and place them under `data/raw/era5/`, `data/raw/united/`, `data/raw/iom/`, and `data/raw/frontex/` respectively. The other raw sources (ACLED, IMO SAR boundaries, UNHCR) ship with the repository.
3. Restore the R environment and run the full pipeline:
   ```r
   source("analysis/R/run_all.R")
   ```
   This calls `renv::restore(prompt = FALSE)` and then sources every script under `02_clean/` → `03_build/` → `04_descriptive/` → `05_analysis/` → `06_robustness/` in alphanumeric order. All figures and tables are written to `output/figures/` and `output/tables/`.
4. Render the paper:
   ```bash
   quarto render paper/thesis.qmd --to pdf
   ```
   The Quarto document reads tables and figures directly from `output/`, so the render reproduces the manuscript exactly as cited.
