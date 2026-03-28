# Thesis: Mortality on the Central Mediterranean Route

## Project overview

MDS thesis (Hertie School) studying weather and mortality on the Central
Mediterranean migration route. Original design tested whether the 2017
EU-Libya MoU (treatment date: July 2017) changed the weather-mortality gradient (Deiana et al. 2024
design). Diagnostics revealed the gradient is **not stable** — it
oscillates substantially over the full period (2014-2025), and the MoU
date does not stand out in placebo tests. The thesis documents this
instability and discusses its implications.

## Current status (March 2025)

- Primary analysis complete: daily panel NegBin with week-year FE
- Key finding: weather-mortality gradient oscillates over time
  (negative 2015-16, positive 2017-19, negative again 2020-21)
- Placebo tests show MoU date is not special — many dates give
  similar or larger coefficients
- Thesis framing in progress — writing phase, one month to submission

## Analysis pipeline (run in order)

1. `01b_core_corridor_dataset.R` — build clean incident-level dataset for core corridor
2. `02_build_event_data.R` — build event-level dataset from IOM MMP + ERA5
3. `02b_daily_weather_panel.R` — build daily panel with spatial-mean weather
4. `03a_weather_overlap_diagnostic.R` — pre/post weather distribution overlap
5. `03b_placebo_stability.R` — placebo dates (day-0 spec, archival comparison)
6. `03b2_gradient_evolution.R` — rolling-window gradient + expanding-window β₃
7. `03b3_placebo_lag1.R` — placebo dates (lag-1 spec, primary diagnostic)
8. `03c_daily_panel_model.R` — **primary model**: core SWH lag-1 × Post, NegBin, week-year FE
9. `03d_event_model_revised.R` — complementary event-level analysis

Archived (in `analysis/R/archive/`): `03_negbin_model_a.R` (original event model),
`03e_fe_sensitivity.R`, `04a_dml_robustness.R`, `04b_sensitivity_analysis.R`.

Exploratory (not in pipeline): `autocorrelation_diagnostics.R`,
`explore_pca_sea_danger.R`, `monthly_count_negbin.R`,
`monthly_weather_mortality_regression.R`.

## Primary specification

- **Design:** Daily panel, Deiana et al. (2024) style
- **Outcome:** `n_dead_missing` (daily count of dead + missing, drowning only)
- **Weather:** `swh_core_lag1` — spatial mean SWH over core corridor [10.5,15.5]×[32.3,36.2], lag-1
- **Interaction:** `swh_core_lag1 × post_mou` (β₃)
- **FE:** Week-by-year (primary) and month-by-year (robustness)
- **Distribution:** Negative binomial (extreme overdispersion)
- **SEs:** Heteroskedasticity-robust (clustered SEs reported as robustness)
- **Timing:** Lag-1 is primary (IOM date = reporting date; lag-1 captures transit weather)
- **No wind control:** SWH already incorporates wind effects (Deiana et al. do not include wind separately)

## Key results

- Primary β₃ = +1.355 (IRR = 3.88, p = 0.011) with week-year FE
- Month-year FE: β₃ = +1.394 (p < 0.001, N = 4,143)
- Lag-2 confirms: β₃ = +1.380 (p = 0.012); lag-3/lag-7 null (timing falsification passes)
- Clustered SEs: p rises to 0.105
- Event-level model: same direction (+0.56 at lag-1, +0.88 at lag-2)
- **Placebo test: MoU date is not uniquely special** — among the strongest
  dates under week-year FE, but matched by 2015-04 and 2021-01
- **Gradient evolution:** rolling 2-year window shows gradient is always
  negative (deterrence) but oscillates; less negative around MoU period

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
