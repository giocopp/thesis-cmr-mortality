# Thesis: Mortality on the Central Mediterranean Route

## Project overview

MDS thesis (Hertie School) studying weather and mortality on the Central
Mediterranean migration route. Original design tested whether the 2017
EU-Libya MoU (treatment date: July 2017) changed the weather-mortality gradient (Deiana et al. 2024
design). Diagnostics revealed the gradient is **not stable** — it
oscillates substantially over the full period (2014-2025), and the MoU
date does not stand out in placebo tests. The thesis documents this
instability and discusses its implications.

## Current status (March 2026)

- Primary analysis complete: daily NegBin with month-year FE and NW(28) SEs
- Key finding: weather-mortality gradient oscillates over time
  (negative 2015-16, positive 2018+, with year-to-year variation)
- Placebo tests show MoU date is not special — many dates give
  similar or larger coefficients
- Falsification test passes: future weather (next-7d) has no effect
- SE diagnostic complete: HC SEs are anti-conservative; NW(28) is primary

## Analysis pipeline

### Current scripts (in `analysis/R/`)
- `00_define_sea_zones.R` — define sea zone boundaries
- `00b_build_daily_panel.R` — build daily panel
- `01_weather_danger_analysis.R` — weather danger analysis
- `02_fatality_rate_timeseries.R` — fatality rate time series
- `03_crossing_components.R` — crossing components
- `04_swh_vs_fatality_rate.R` — SWH vs fatality rate
- `05_reduced_form_model.R` — **primary model**: NegBin, SWH×PostMoU, NW(28) SEs
- `05b_reduced_form_plots.R` — reduced-form plots
- `05c_migrant_files_diagnostic.R` — Migrant Files data diagnostic
- `05d_fe_structure_diagnostic.R` — FE structure diagnostic
- `05e_weekly_panel_model.R` — weekly panel model
- `05f_se_diagnostic.R` — SE diagnostic (ACF, Breusch-Godfrey, NW vs HC comparison)

## Primary specification (05_reduced_form_model.R)

- **Design:** Daily reduced-form (not conditioning on crossings)
- **Outcome:** `n_dead_missing` (daily count of dead + missing, drowning + mixed/unknown)
- **Weather:** `swh_prevweek_z` — standardized previous-week average SWH (spatial mean from ERA5 daily panel)
- **Core corridor:** lon [10.0, 15.1] × lat [32.4, 37.8]
- **Interaction:** `swh_prevweek_z × post_mou` (β₃), treatment = 2017-07-01
- **FE:** Month-by-year (`month_year`)
- **Distribution:** Negative binomial (`fixest::fenegbin`)
- **SEs:** Newey-West(28) primary; NW(14) and HC as robustness
- **Timing:** Previous-week average SWH is primary; prev-3d as robustness
- **Falsification:** Next-7d future SWH (should be null — confirmed)
- **Sample periods:** 2014-2021 and 2014-2024 (both reported)
- **Relation to Deiana et al. (2024):** Shares the logic of testing whether
  a policy changed the weather-outcome elasticity. Key differences: they study
  crossings (not deaths), use a triple interaction with boat-type fraction,
  Poisson QMLE, and treat SAR expansion (not MoU) as the policy event.

## Key results (NW(28) SEs)

### Primary interaction (β₃ = swh_prevweek_z × post_mou)
- **2014-2021:** β₃ = +1.241 (SE=0.540, IRR=3.46, p=0.022)
- **2014-2024:** β₃ = +0.996 (SE=0.446, IRR=2.71, p=0.026)

### Robustness
- All CMR (no corridor restriction): p=0.011 (2014-2021), p=0.023 (2014-2024)
- No outliers (>100 deaths removed): p=0.024 (2014-2021), p=0.021 (2014-2024)
- Prev-3d SWH: weaker, not significant (p≈0.19 both periods)
- Next-7d falsification: null as expected (p=0.61 and p=0.89)

### SE sensitivity (2014-2024 primary spec)
- HC: SE=0.250, p<0.001 (anti-conservative)
- NW(14): SE=0.400, p=0.013
- NW(28): SE=0.446, p=0.026 (primary)
- Cluster(week): SE=0.337, p=0.003

### Year-by-year gradient (2014-2024, NW(28))
- Negative: 2015 (-1.75), 2016 (-0.56), 2017 (-0.71), 2022 (-0.88*)
- Positive: 2018 (+0.75), 2021 (+0.76), 2023 (+0.67)
- Wide CIs for most years; gradient oscillates, not a clean break at MoU

### Diagnostics
- **Placebo test:** MoU date is not uniquely special — many dates give
  similar or larger coefficients
- **Residual ACF:** individually small (lag-1 = -0.012) but Breusch-Godfrey
  rejects no serial correlation at all orders → NW SEs are appropriate

## Repository structure

- `analysis/` — R scripts (numbered pipeline), Python for ERA5 downloads
- `data/` — `raw/` (immutable), `processed/` (cleaned/merged)
- `output/` — figures, tables, saved models
- `paper/` — thesis draft
- `documents/` — reasoning documents, DAG analysis, presentations
- `replication/` — completed replication work (archival)
- `literature/` — paper summaries

## Key conventions

- R is the primary language; Python for ERA5/Copernicus API downloads
- Numbered R scripts (01_, 02_, ...) run in order
- renv for package management (single lockfile at root)
- Raw data never modified; processed data is reproducible from code

## Data sources

- **IOM MMP:** incident-level records (geocoded, 2014-2025)
- **ERA5:** SWH, wind (u10/v10), mean wave period, wind gust (0.25°/0.5°, daily)
- **IOM monthly:** crossings by route (arrivals + interceptions + deaths)
- **Replication data:** Rodriguez-Sanchez extension (monthly, 2011-2021)
