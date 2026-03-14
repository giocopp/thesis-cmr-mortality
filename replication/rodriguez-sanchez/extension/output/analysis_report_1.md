# Mortality Counterfactual Analysis: Extension of Rodriguez Sanchez et al. (2023)

## Analysis Report — Extension 1: BSTS Mortality Model

---

## 1. Research Question

This analysis extends Rodriguez Sanchez et al. (2023), "Search-and-rescue in the Central Mediterranean Route does not induce migration" (Nature Scientific Reports), by asking a different question using the same methodological framework:

**Did the EU-Libya Memorandum of Understanding (MoU) of February 2017 make the Central Mediterranean route more dangerous?**

While the original paper tests whether search-and-rescue (SAR) operations acted as a "pull factor" for migration (outcome: crossing attempts), this extension tests whether the policy shift toward Libyan Coast Guard cooperation increased mortality among those who attempted the crossing (outcome: mortality rate).

We additionally test whether the two preceding policy periods — Mare Nostrum (Oct 2013) and the NGO SAR era (Nov 2014) — had any effect on mortality, asking whether active rescue operations reduced the death rate.

---

## 2. Methodology

### 2.1 Causal Inference Framework

We use the same Bayesian Structural Time-Series (BSTS) approach as the original paper, implemented via Google's `CausalImpact` R package (Brodersen et al., 2015). The method:

1. Fits a state-space model to the **pre-intervention** time series, using a local linear trend component and a regression component with covariates selected via **spike-and-slab priors**
2. Generates a **counterfactual prediction** of what the outcome would have been in the post-intervention period, absent the intervention
3. Compares the observed outcome to this counterfactual to estimate the **causal effect** of the intervention

The spike-and-slab prior is a Bayesian variable selection mechanism: each candidate predictor has a "spike" probability of being included in the model and a "slab" prior on its coefficient if included. This allows automatic selection from a large set of candidate covariates without overfitting.

### 2.2 Key Assumptions

The CausalImpact method requires:

1. **No anticipation**: The intervention was not anticipated by the outcome variable
2. **Exogeneity of covariates**: The covariates used to build the counterfactual must not themselves be affected by the intervention
3. **Stable relationships**: The covariate-outcome relationships estimated in the pre-period must remain valid in the post-period (absent the intervention)

### 2.3 Model Parameters

Identical to the original paper for comparability:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `niter` | 10,000 | MCMC iterations for posterior sampling |
| `dynamic.regression` | FALSE | Static regression coefficients (not time-varying) |
| `standardize.data` | TRUE | Z-score standardization of all predictors |
| `set.seed` | 270488 | Reproducibility (same seed as original paper) |
| `alpha` | 0.05 | Significance level for credible intervals |

For the **curated specification** (primary), we set `max.flips = -1` (unlimited), allowing the MCMC sampler to consider flipping every variable's inclusion indicator at each iteration. For the **full specification** (robustness), we use `max.flips = 100` as in the original paper.

### 2.4 Outcome Variables

- **Primary**: Log-transformed mortality rate per 100 crossing attempts: `log(deaths/crossings * 100 + 0.01)`
- **Robustness**: Log-transformed death count: `log(deaths + 1)`

The small constant (0.01 for rates, 1 for counts) handles zero-value months under the log transformation.

### 2.5 Intervention Periods

Three interventions are tested, with the same pre/post definitions as the original paper:

| Model | Intervention | Pre-period | Post-period |
|-------|-------------|------------|-------------|
| A | Mare Nostrum (state-led SAR) | Feb 2011 – Sep 2013 | Oct 2013 – Sep 2021 |
| B | NGO SAR operations | Feb 2011 – Oct 2014 | Nov 2014 – Sep 2021 |
| C | EU-Libya / MoU cooperation | Feb 2011 – Jan 2017 | Feb 2017 – Sep 2021 |

---

## 3. Data

### 3.1 Original Dataset

The analysis starts from the dataset constructed by Rodriguez Sanchez et al. (`df.RDS`), which contains 156 monthly observations (January 2009 – December 2021) across 14,840 columns. Variables include:

- **Conflict indicators**: UCDP/PRIO Armed Conflict data for ~100 countries (battles, explosions/remote violence, protests, riots, strategic developments, violence against civilians) — with lags 01–24
- **Disaster counts**: EM-DAT disaster counts by country — with lags 01–24
- **Airport passenger flows**: Eurostat international air passenger data for ~80 countries — with lags 1–6
- **Exchange rates**: 38 African/MENA currencies to EUR — with lags 01–24
- **Commodity prices**: 86 IMF commodity price indices — with lags 01–24
- **Unemployment**: 38 EU/OECD countries — with lags 01–24
- **Weather**: Temperature, precipitation, and storm days for Italy and Malta — with lags 01–12
- **Google Trends**: Job/employment search indices for 5 North African countries
- **Migration data**: Frontex/IOM arrivals by route, IOM dead and missing by route, LCG/TCG pushback counts, geographic dispersion indices

### 3.2 New ERA5 Sea Condition Variables

We downloaded ERA5 reanalysis monthly mean data from the Copernicus Climate Data Store for the study period (Jan 2009 – Oct 2021). Three separate requests were made:

**Request 1: Atmospheric variables — Central Mediterranean (0.25° grid)**
- Bounding box: 31°N–38°N, 10°E–20°E (Central Med crossing zone)
- Variables: 10m wind u/v-components, sea surface temperature, 2m temperature, total cloud cover, low cloud cover, 2m dewpoint temperature

**Request 2: Wave variables — Central Mediterranean (0.5° grid)**
- Same bounding box
- Variables: significant wave height, mean wave period, mean wave direction

**Request 3: Atmospheric variables — North African departure coast (0.25° grid)**
- Bounding box: 30°N–34°N, 8°E–25°E (Libya/Tunisia coastal strip)
- Variables: 10m wind u/v-components, 2m temperature, total precipitation, total cloud cover

From these raw data, we computed spatial means over each grid and derived the following 9 base variables (plus lags 01–06 for each = 63 total new columns):

| Variable | Source | Description |
|----------|--------|-------------|
| `wave_height_central_med` | ERA5 waves | Significant height of combined wind waves and swell (m) |
| `wave_period_central_med` | ERA5 waves | Mean wave period (s); shorter periods = more dangerous |
| `wave_direction_central_med` | ERA5 waves | Mean wave direction (degrees) |
| `wind_speed_central_med` | ERA5 atmos | 10m wind speed magnitude over crossing zone (m/s) |
| `sst_central_med` | ERA5 atmos | Sea surface temperature (°C); affects hypothermia survival |
| `cloud_cover_central_med` | ERA5 atmos | Total cloud cover (fraction 0–1); visibility proxy |
| `dewpoint_depression_central_med` | ERA5 atmos | T2m minus dewpoint (°C); near 0 = fog risk |
| `wind_speed_departure_coast` | ERA5 coast | 10m wind speed at departure (Libya/Tunisia coast) |
| `cloud_cover_departure_coast` | ERA5 coast | Cloud cover at departure; affects departure decisions |

Additionally, temperature and precipitation from the ERA5 coast request were integrated as `temperature_departure_coast`, `precipitation_departure_coast`, and `temperature_central_med` — adding ~30 more columns (with lags). The total extension is **93 new columns**, bringing the dataset to 14,933 columns.

### 3.3 Mortality Outcome Construction

Following the original paper's convention:

```
crossings_CMR = arrivals_CMR + LCG_pushbacks_count + TCG_pushbacks_count
                + dead_and_missing_Central_Mediterranean

mortality_rate_100 = (dead_and_missing_Central_Mediterranean / crossings_CMR) * 100
```

Descriptive statistics for the analysis window (Feb 2011 – Sep 2021, 128 months):

| Statistic | Mortality rate (per 100) | Monthly deaths | Monthly crossings |
|-----------|:---:|:---:|:---:|
| Mean | 2.66 | 183.5 | 8,116 |
| Median | 1.74 | — | — |
| Std. dev. | 3.25 | — | — |
| Min | 0.00 | — | — |
| Max | 20.70 | — | — |
| Zero months | 8 | — | — |

**Pre/post MoU comparison** (raw means, not causal estimates):

| Period | Mortality rate | Monthly deaths | Monthly crossings |
|--------|:---:|:---:|:---:|
| Pre-MoU (Feb 2011 – Jan 2017) | 2.47 | 218.3 | 9,260 |
| Post-MoU (Feb 2017 – Sep 2021) | 2.90 | 131.1 | 6,433 |

The raw comparison shows a modest increase in the mortality *rate* (2.47 to 2.90) but a substantial decrease in absolute deaths (218 to 131) driven by the large drop in crossings (9,260 to 6,433).

---

## 4. Variable Selection: Three Attempts

### 4.1 Attempt 1 — Full Predictor Set (~4,819 predictors)

Our first attempt used the same variable filtering as the original paper: remove high-order lags (keep only 01–06), remove endogenous variables (crossing volumes, pushbacks, SAR indicators, deaths on other routes), and pass all remaining ~4,819 predictors to CausalImpact's spike-and-slab.

**This approach failed.** The spike-and-slab selected only **2 out of 4,819 predictors** with inclusion probability above 5%: `num_expvio_Algeria_lag_03` (99.9%) and `PSORG_lag_06` (99.8%) — both certainly spurious correlations with no causal mechanism for Central Med mortality. All 63 ERA5 sea condition variables had exactly 0.000 inclusion probability. Confidence intervals were extremely wide.

### 4.2 Attempt 2 — Curated Predictor Set (~808 predictors)

We reduced the predictor set to ~808 more theoretically motivated variables: ERA5 sea conditions, weather, conflicts for key Central Med origin/transit countries, exchange rates for key currencies, destination unemployment, commodity prices, and Google Trends. We also set `max.flips = -1` for exhaustive MCMC exploration.

**This produced statistically significant results** — but the predictor selection raises serious concerns (see Section 6).

### 4.3 Variable Exclusions (Exogeneity)

Both specifications exclude the following variables as endogenous to the intervention:

- **Crossing volumes** (arrivals by route): directly affected by MoU deterrence
- **Pushback counts** (LCG, TCG): created by MoU policy
- **SAR vessel/operation indicators**: affected by MoU-era restrictions
- **Deaths on other routes**: could reflect route substitution
- **Geographic dispersion indices**: reflect route changes caused by policy
- **Asylum application data**: downstream of crossing decisions

---

## 5. Results

### 5.1 Curated Specification — Mortality Rate

| Model | Intervention | Cum. effect | Rel. effect | 95% CI | p-value | Significant |
|:---:|-------------|:---:|:---:|:---:|:---:|:---:|
| A | Mare Nostrum | +601.7 | — | [-582, +1292] | 0.110 | No |
| B | NGO SAR | -75.6 | — | [-675, +747] | 0.387 | No |
| **C** | **EU-Libya/MoU** | **-165.7** | **-70.8%** | **[-608, +9.4]** | **0.029** | **Yes** |

### 5.2 Curated Specification — Death Count (Robustness)

| Model | Intervention | Cum. effect | Rel. effect | 95% CI | p-value | Significant |
|:---:|-------------|:---:|:---:|:---:|:---:|:---:|
| **C** | **EU-Libya/MoU** | **-326.3** | **-56.6%** | **[-785, -33]** | **0.012** | **Yes** |

### 5.3 Full Specification (Robustness) — 4,819 Predictors

| Model | Outcome | Cum. effect | Rel. effect | 95% CI | p-value | Significant |
|:---:|---------|:---:|:---:|:---:|:---:|:---:|
| C | Mortality rate | -23.9 | -40.8% | [-68, +18] | 0.094 | No (borderline) |
| C | Death count | -253.3 | -46.9% | [-686, +28] | 0.058 | No (borderline) |

### 5.4 Summary of Findings

1. **Model C (EU-Libya/MoU) shows a statistically significant effect in the curated specification** (p = 0.029 for mortality rate, p = 0.012 for death count), consistent at borderline significance in the full specification (p = 0.094 and p = 0.058).

2. **The direction is opposite to the hypothesis**: mortality was *lower* than the counterfactual prediction, not higher. The model estimates that the MoU period is associated with a 70.8% reduction in cumulative mortality rate relative to the counterfactual.

3. **Models A and B are not significant**: Mare Nostrum (p = 0.110) and NGO SAR (p = 0.387) show no statistically detectable effect on mortality at the 5% level.

---

## 6. Critical Evaluation: Why These Results Are Unreliable

### 6.1 The Predictor Selection Problem

The curated Model C selected the following top predictors:

| Rank | Variable | Inc. prob. |
|:---:|----------|:---:|
| 1 | `unem_MALTA` | 21.5% |
| 2 | `unem_euro_area_all` | 15.4% |
| 3 | `unem_euro_area_all_lag_01` | 14.2% |
| 4 | `num_violciv_Nigeria` | 13.8% |
| 5 | `unem_eu_27Acountries_lag_01` | 10.9% |
| 6 | `unem_SPAIN_lag_01` | 10.9% |
| 7 | `disas_count_Libya` | 10.8% |
| 8 | `num_violciv_South.Sudan_lag_01` | 10.0% |

**European unemployment dominates the predictor set.** This is a major red flag. There is no plausible mechanism by which unemployment in Malta or the EU would directly cause deaths at sea. Unemployment is almost certainly being selected because it correlates with the general time trend of mortality — both decline over parts of the period. The spike-and-slab is using unemployment as a **spurious trend proxy**, not as a genuine predictor of per-crossing danger.

This means the "significant" result (p = 0.029) may be an artifact: the model extrapolates a trend captured by unemployment, and when the post-MoU period diverges from this extrapolated trend, it registers as a "causal effect." But the trend was never causally driven by unemployment — it was driven by unobserved factors (rescue capacity, route choices, boat quality) that happen to correlate with economic trends.

### 6.2 ERA5 Variables Still Marginal

Even in the curated set, ERA5 sea condition variables received very low inclusion probabilities. The highest was `wind_speed_departure_coast_lag_06` at just 4.4%. Wave height, SST, and other variables were below 1%. This confirms that monthly-aggregated sea conditions are weak predictors of monthly mortality — mortality events depend on conditions during specific crossing days, not monthly averages.

### 6.3 The Fundamental Methodological Tension

The CausalImpact approach requires covariates that are (a) exogenous to the intervention AND (b) good predictors of the outcome. For mortality, these two requirements are in direct conflict:

- **Variables that genuinely predict per-crossing mortality** are endogenous to the MoU:
  - Rescue vessel proximity / response time (reduced by MoU)
  - Route choices and distances (changed by LCG enforcement)
  - Boat quality and overcrowding (affected by smuggler adaptation to enforcement)
  - Pushback behavior (created by MoU)

- **Variables that are exogenous to the MoU** do not genuinely predict mortality:
  - Macro-economic indicators (conflicts, unemployment, exchange rates) predict who *decides* to cross, not who *dies* while crossing
  - Sea conditions at monthly aggregation are too coarse to capture crossing-day risk
  - Time components capture seasonality but cannot separate it from the intervention effect

When the treatment works through *all* the plausible predictors of mortality, there is nothing left to build a meaningful counterfactual with. The model defaults to learning spurious time-trend correlations, producing unreliable causal estimates.

### 6.4 What This Means for the Results

The results from this extension should be treated as **exploratory and preliminary**, not as definitive causal estimates. The significant p-values (0.029, 0.012) likely overstate confidence because:

1. The counterfactual is built on spurious predictors (unemployment-as-trend-proxy)
2. The pre-treatment fit, while visually acceptable, may reflect overfitting to trends rather than genuine covariate relationships
3. The direction of the effect (lower mortality) is consistent across specifications, which suggests a real pattern in the data, but the magnitude and significance are sensitive to predictor choice

---

## 7. Paths Forward

The failure of the standard BSTS/CausalImpact approach for mortality motivates three alternative strategies:

### Path A: Minimal Exogenous Model (Extension 2)

Strip the predictor set to *only* variables that are both truly exogenous to the MoU AND have a direct theoretical link to per-crossing mortality risk:

- **ERA5 sea conditions**: wave height, wave period, wind speed, SST, cloud cover, fog index, departure coast conditions (base + lags 01-06)
- **Weather**: temperature, precipitation, storm days for Italy, Malta, Central Med, departure coast (base + lags 01-06)
- **Oil prices**: POILBRE (proxy for fuel costs, which affect boat range and smuggler economics)
- **Seasonality**: month, quarter, semester

This yields ~170 predictors. If these do not predict mortality well (as initial evidence suggests), the model reduces to a Bayesian interrupted time series — a simpler but *honest* test of whether mortality changed from its pre-intervention trend. We report this transparently as a limitation of monthly-aggregated data.

### Path B: Mechanism-Specific Models (Future Work)

Instead of testing "did mortality change?", test whether the *mechanisms* that drive mortality changed:

- **Geographic dispersion of deaths**: Did deaths shift further from the Libyan coast? (using `sd_lat`, `sd_lon`, `frac_index` as *outcomes*, not covariates)
- **Seasonal concentration**: Did crossings shift to more dangerous months?
- **Deaths per incident**: Did individual shipwrecks become deadlier? (requires incident-level IOM MMP data)

### Path C: Mediation Analysis (Extension 3)

Decompose the causal chain: MoU -> fewer crossings -> fewer deaths, vs. MoU -> more dangerous crossings -> higher per-crossing mortality. Include crossing volume as a mediator rather than excluding it, allowing decomposition of total vs. direct effects.

---

## 8. Files and Reproducibility

### Code

| File | Description |
|------|-------------|
| `code/download_era5.py` | Downloads ERA5 data from Copernicus CDS (3 requests) |
| `code/build_mortality_dataset.R` | Processes ERA5 NetCDF -> merges with df.RDS -> df_extended.RDS |
| `code/run_mortality_model.R` | Runs all CausalImpact models (curated + full) and generates figures |

### Data

| File | Description |
|------|-------------|
| `data/era5/*.nc` | Raw ERA5 NetCDF files (4 files, ~4.6 MB total) |
| `data/df_extended.RDS` | Extended dataset (156 rows x 14,933 columns) |

### Results

| File | Description |
|------|-------------|
| `results/models/model_curated_mortality_*.RDS` | Curated CausalImpact model objects (A, B, C) |
| `results/models/model_curated_deaths_sarlibya.RDS` | Curated death count model (C only) |
| `results/models/model_full_mortality_sarlibya.RDS` | Full-specification mortality model (C only) |
| `results/models/model_full_deaths_sarlibya.RDS` | Full-specification death count model (C only) |
| `results/data/mortality_model_results.RDS` | Extracted time series results for all models |
| `results/figures/png/Figure1_mortality.png` | Descriptive time series (deaths + mortality rate) |
| `results/figures/png/Figure2_mortality_counterfactual.png` | Counterfactual analysis (curated, 3 panels) |
| `results/figures/png/Figure3_mortality_pointwise.png` | Pointwise effects (curated, 3 panels) |

Seed: `set.seed(270488)` — identical to the original paper.
