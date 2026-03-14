# Thesis: Mortality on the Central Mediterranean Route

## Project overview

MDS thesis (Hertie School) studying whether the 2017 EU-Libya MoU made
the Central Mediterranean migration route more dangerous. Primary
analysis uses an event-level weather x post-MoU interaction design on
IOM Missing Migrants Project incident data matched to ERA5 weather.

## Research design

- **Primary (Model A):** NegBin triple-diff at event level:
  `deaths_i ~ NegBin(mu_i)`, `log(mu_i) = grid_FE + month_FE + beta_3(SWH x Post x CMR) + X'gamma`
  Uses EMR/WMR/WAAR as control routes.
- **Robustness:** DML (partially linear), Lee bounds (selection), causal forest (heterogeneity)
- **Complementary (monthly):** death counts ~ weather x post, crossings as outcome (diagnostic), rate regression (exploratory, collider caveat)

## Repository structure

- `analysis/` -- new event-level analysis (R scripts, python data downloads)
- `data/` -- all data: `raw/` (immutable), `processed/` (cleaned/merged)
- `output/` -- figures, tables, saved models
- `paper/` -- thesis draft
- `documents/` -- reasoning documents, DAG analysis, presentations
- `replication/` -- completed replication work (archival)
- `literature/` -- paper summaries

## Key conventions

- R is the primary language; Python for ERA5/Copernicus API downloads
- Numbered R scripts (01_, 02_, ...) run in order
- `_targets.R` for pipeline orchestration where used
- renv for package management (single lockfile at root)
- Raw data never modified; processed data is reproducible from code

## Data sources

- **IOM MMP:** incident-level records (geocoded, 2014-2025)
- **ERA5:** SWH, wind, SST, currents (daily, 0.25 degree)
- **UNHCR:** daily arrivals to Italy
- **IOM monthly:** crossings by route (arrivals + interceptions + deaths)
- **UNITED:** death list (supplementary)
- **Frontex JORA:** boat-type data (pending request)
