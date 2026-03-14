---
output:
  html_document: default
  word_document: default
---
# Replication Guide: "Search-and-Rescue in the Central Mediterranean Route Does Not Induce Migration"

**Paper:** Rodriguez Sanchez, A., Wucherpfennig, J., Rischke, R., & Iacus, S.M. (2023). *Scientific Reports*, 13, 11014, https://doi.org/10.1038/s41598-023-38119-4  
**Replication by:** Giorgio Coppola

---

## 1. Motivation and Policy Context

### 1.1 The "Pull Factor" Claim

The Central Mediterranean Route (CMR) -- connecting North African departure points (primarily Libya) to Italy and Malta -- is the most dangerous irregular migration route into Europe. Thousands of migrants have drowned attempting the crossing.

In response to the humanitarian crisis, both state-led and private search-and-rescue (SAR) operations have been deployed in the Mediterranean since 2013. These operations have become the subject of a politically charged debate centered on the **"pull factor" hypothesis**: the claim that the presence of SAR boats *encourages* more migrants to attempt the crossing, thereby *increasing* migration flows. This claim has been invoked by EU policymakers to justify the scaling-back of rescue operations and the criminalization of NGO vessels.

The pull factor claim rests on three sub-arguments:
1. SAR boats encourage more crossing attempts by reducing the perceived danger.
2. SAR operations unintentionally help smugglers' businesses by reducing their operational costs (e.g., smugglers can use cheaper, more dangerous boats if they expect rescues).
3. SAR has the unintended consequence of making the trip more dangerous (e.g., through shifts to flimsier vessels).

### 1.2 Why Existing Evidence Was Insufficient

Prior to this paper, the empirical evidence for the pull factor claim was "both scant and methodologically compromised" (Supplementary Materials, Table S7). Previous studies suffered from:

- **No causal identification strategy.** Most studies simply compared arrival counts during high-SAR vs. low-SAR periods, without accounting for confounders.
- **Ignoring time-series properties.** Migration data exhibit serial autocorrelation, seasonality, and non-stationarity -- features that standard regression methods (OLS, Poisson) cannot properly handle.
- **Endogeneity.** SAR operations generally came *after* increases in migration flow (i.e., SAR is a response to higher flows, not a cause), making naive comparisons misleading.

### 1.3 What This Paper Contributes

Rodriguez Sanchez et al. provide the first **causal inference assessment** of the pull factor claim using a method specifically designed for time-series intervention analysis: Bayesian Structural Time-Series (BSTS) models with the `CausalImpact` framework. Their approach:

- Accounts for autocorrelation, trends, seasonality, and non-stationarity simultaneously.
- Builds a synthetic counterfactual from hundreds of exogenous predictors ("push-and-pull factors").
- Tests three distinct intervention periods, each corresponding to a major policy change in Mediterranean SAR.
- Uses spike-and-slab priors for automatic variable selection.

---

## 2. Research Design

### 2.1 The Three Interventions

The paper identifies three distinct intervention periods based on key changes in the politics of SAR:

**Model A -- Mare Nostrum (October 2013 -- October 2014)**
The largest state-led SAR operation in European history, launched by Italy after the Lampedusa shipwreck of October 3, 2013 (366 deaths). Mare Nostrum was not purely humanitarian: it also included anti-smuggler operations that disrupted smuggling networks and destroyed vessels. Pre-intervention period: February 2011 -- September 2013.

**Model B -- NGO Search-and-Rescue (November 2014 -- February 2017)**
After Mare Nostrum ended, private NGOs filled the gap, starting with MOAS (Migrant Offshore Aid Station) on August 26, 2014. Over time, up to 10+ NGO vessels operated simultaneously (see Figure S1). The intervention start is set at November 2014, immediately after Mare Nostrum's end. Pre-intervention period: February 2011 -- October 2014.

**Model C -- EU-Libya Cooperation (February 2017 onward)**
Beginning in early 2017, the EU cooperated with the Libyan Coast Guard (LCG) to intercept and return ("push back") migrants before they could reach European waters. This period also saw the extension of the Libyan SAR zone and the criminalization of NGO vessels. Pre-intervention period: February 2011 -- January 2017.

### 2.2 Causal Logic

The identification strategy is counterfactual:

1. **Pre-intervention period (t <= tau):** Fit a predictive model using exogenous covariates to learn the relationship between "push-and-pull factors" and crossing attempts.
2. **Post-intervention period (t > tau):** Use the fitted model to predict *what crossing attempts would have been* in the absence of the intervention (the counterfactual).
3. **Causal effect:** Compare the observed post-intervention crossings with the counterfactual prediction. A significant difference implies a causal effect of the intervention.

The key identifying assumption is that the control series (covariates) are **exogenous** to the interventions -- i.e., SAR operations did not cause changes in weather, commodity prices, African conflicts, etc. The paper tests this assumption via Granger causality tests (Tables S4-S6 in the supplementary materials).

---

## 3. Data

### 3.1 Dependent Variable

The outcome variable is **monthly attempted crossings** through the CMR, defined as:

```
Y_t = A_t + P_t + D_t
```

Where:  
- **A_t** = Arrivals to Italy and Malta (FRONTEX illegal border crossings + IOM arrival data)  
- **P_t** = Pushbacks by the Libyan and Tunisian Coast Guards (IOM Missing Migrants Project)  
- **D_t** = Dead and missing migrants in the Central Mediterranean (compiled from UNITED, IOM-MMP, and The Migrant Files)  

This composite measure captures *all* crossing attempts, not just successful arrivals. The variable is **log-transformed** before modeling (`log(Y_t)`) to stabilize variance, as the series contains periods of both increasing and decreasing variability.

In the code, this is constructed as:
```r
crossings_CMR = arrivals_CMR + LCG_pushbacks_count + TCG_pushbacks_count +
                dead_and_missing_Central_Mediterranean
```

NA values in pushbacks and deaths are imputed as 0 (pre-2016 for pushbacks; pre-2014 for IOM deaths). The mortality rate (`D_t / Y_t * 1000`) is also computed for Figure 1B but is not used in the models.

### 3.2 Covariates (Push-and-Pull Factors)

The predictor set includes ~4,726 variables at the original lag structure, drawn from 11 data source categories (Table S1 in supplementary materials):

| Category | Variables | Source | Period |
|---|---|---|---|
| Airport flows | Monthly passenger counts from African/MENA countries to Europe | Sabre Travel Data | 2011-2021 |
| Currency exchange rates | African and MENA currencies to EUR | ECB Statistical Data Warehouse | 2010-2020 |
| Job search indicators | Google Trends for "job", "work", "employment" in Arabic | Google Trends | 2008-2020 |
| Conflict indicators | Battles, explosions, violence, protests, riots in African countries | ACLED | 2008-2020 |
| Syrian conflict intensity | Google Trends for "Syrian war" | Google Trends | 2008-2020 |
| Commodity prices | Energy, agriculture, fertilizers, metals (30+ commodities) | IMF | 2000-2020 |
| Environmental disasters | Natural and technological disaster counts by country | EM-DAT | 1998-2020 |
| Weather | Temperature, precipitation, storm days in Italy and Malta | ECAD | 2010-2020 |
| EU unemployment | Unemployment rates in EU and specific countries | Eurostat | 2008-2022 |
| African unemployment | Unemployment rates in African countries | ILO/World Bank | various |
| SAR operation dates | Binary indicators for 22 SAR vessels + EU operations | Own elaboration | 2009-2021 |

Each time-varying covariate is lagged up to 6 months (`lag_01` through `lag_06`), reflecting that push-and-pull factors may operate with delay (e.g., a conflict shock in Libya may take months to translate into migration flows).

### 3.3 Variable Filtering

Before modeling, the predictor set is reduced:

1. **Time window:** Restricted to February 2011 -- September 2021 (128 months).  
2. **High-order lags removed:** Lags 07-24 are dropped, keeping only lags 01-06. This reduces the predictor count to 4,726 variables.
3. **Outcome-related variables removed:** All arrival counts for other routes (BSR, CRAG, EBR, EMR, OR, WAR, WBR, WMR), death counts from other Mediterranean routes, pushback counts, mortality rate, and geographic dispersion indices (`sd_lat_*`, `sd_lon_*`, `frac_index_*`).  
4. **Problematic variables removed:** Palestinian Territories airflow data (missing/problematic) and asylum application data.
5. **Time components added:** `month`, `semester`, and `quarter` indicators are appended.  
6. **Missing data:** Rows with any remaining NA are dropped via `na.omit()`.  

The resulting analysis dataset (`df_min_A`) contains ~113 monthly observations and ~600-700 predictors.

### 3.4 Data Construction Pipeline

The original codebase contains 22 R scripts in `sar_migration_aleja/code/` for building `df.RDS`:

- `arrivals_*.R` (4 scripts): Process FRONTEX, Eurostat, IOM, and JRC airport data.  
- `deaths_*.R` (3 scripts): Compile mortality data from The Migrant Files, IOM-MMP, and UNITED.  
- `covar_*.R` (9 scripts): Process weather, conflicts (ACLED, UCDP), exchange rates, unemployment, disasters, commodity prices, and Google Trends.  
- `dates_search_and_rescue_EU_NGOS_dates.R`: Defines daily binary indicators for 22 SAR vessels and EU operations.  
- `create_data_frame.R`: Merges all sources via `left_join()` operations into the final `df.RDS`.  

---

## 4. Methodology

### 4.1 Bayesian Structural Time-Series (BSTS)

BSTS is a state-space model that decomposes a time series into structural components (trend, seasonality) and regression effects from covariates. It consists of two equations:

**Observation equation:**
```
Y_t = Z_t^T * alpha_t + epsilon_t,    epsilon_t ~ N(0, sigma_t)
```

**State equation:**
```
alpha_{t+1} = T_t * alpha_t + R_t * eta_t,    eta_t ~ N(0, Q_t)
```

The *state equation* represents the hidden, evolving condition of the system at time $t$ that is not directly observed, but that produces what you do observe.

Where:  
- `Y_t` is the observed outcome (log crossings).  
- `alpha_t` is the latent state vector (captures trend, seasonality, etc.).   
- `Z_t` links the state to the observation (can include covariates: `Z_t = beta_t * X_t` when `alpha_t = 1`).  
- `T_t` is the state transition matrix.  
- `R_t * eta_t` is the state error structure.  

In the context of this paper, the paper fit a BSTS with a latent structural component (at least a local level) and static regression coefficients, estimated via MCMC with spike-and-slab variable selection. The model simultaneously estimates:  
- A **local level** trend component (random walk for the level of the series).  
- **Regression coefficients** for the covariates, with spike-and-slab priors for variable selection.  

The key advantage over standard regression is that the BSTS model handles autocorrelation, non-stationarity, and changing variance through the state-space structure, while the spike-and-slab prior handles the high dimensionality of the predictor set.

### 4.2 Spike-and-Slab Variable Selection

The BSTS model uses a **spike-and-slab prior**:

- The **spike** component places a point mass at zero for each regression coefficient, with some prior probability. This means each covariate has a prior probability of being *excluded* from the model.  
- The **slab** component is a diffuse (e.g., normal) prior for the coefficient value *when included*.  

At each MCMC iteration, the algorithm:  
1. Updates the inclusion vector and is restricted to at most `max.flips = 100` inclusion toggles per iteration
2. For included covariates, draws a coefficient value from the posterior.  
3. The **inclusion probability** for each covariate is the fraction of MCMC iterations in which it was included.  

This produces a form of Bayesian model averaging over an enormous model space (2^700 possible models), automatically selecting the most predictive covariates. The parameter `standardize.data = TRUE` ensures that the spike-and-slab prior treats all covariates on the same scale.  

### 4.3 CausalImpact Framework

The `CausalImpact` R package (Google, Brodersen et al. 2015) wraps the BSTS model for causal inference:  

1. The user specifies pre- and post-intervention periods.   
2. The post-intervention outcome values are set to `NA` in the BSTS model.  
3. The BSTS model is fitted on the pre-intervention data only.  
4. The fitted model generates posterior predictive draws for the post-intervention period -- this is the **counterfactual**.  
5. The package computes:  
   - **Pointwise effects:** `observed_t - predicted_t` for each post-intervention month.  
   - **Cumulative effects:** Running sum of pointwise effects.  
   - **Average effects** (absolute and relative) with 95% credible intervals.  
   - A **posterior tail-area probability** (Bayesian p-value) testing whether the cumulative effect is significantly different from zero.  

### 4.4 Model Configuration

All three models use identical MCMC settings:

```r
model.args = list(
  dynamic.regression = FALSE,   # Static regression coefficients (not time-varying)
  standardize.data = TRUE,      # Standardize predictors before fitting
  max.flips = 100,              # Max covariate inclusion changes per MCMC iteration
  niter = 10000                 # 10,000 MCMC iterations
)
alpha = 0.05                    # 95% credible intervals
set.seed(270488)                # Reproducibility seed
```

- `dynamic.regression = FALSE` means the regression coefficients `beta` are constant over time (not varying with each time step). This is a deliberate choice: the paper wants to test whether the relationship between push-pull factors and crossings changed at the intervention, not to model a smoothly evolving relationship.  
- `niter = 10000` provides 10,000 posterior draws for inference.
- `max.flips = 100` limits the spike-and-slab sampler to 100 inclusion/exclusion proposals per iteration, which is necessary for computational tractability with ~600-700 predictors.  

---

## 5. Results

### 5.1 Descriptive Findings (Figure 1)

**Figure 1A** shows attempted crossings over time with shaded intervention periods. Key patterns:  
- Low and rising crossings from 2009 to mid-2011 (Syrian civil war onset).  
- Dramatic peak around 2014-2016 (height of the crisis).  
- Sharp decline after mid-2017 (EU-Libya cooperation / pushbacks).  
- Strong seasonality (summer peaks, winter troughs).  

**Figure 1B** shows the mortality rate per 100 attempted crossings. During the SAR period (Mare Nostrum + NGOs, shaded purple), the mortality rate was *lower* and *less volatile* than before or after. This descriptive finding already challenges the pull factor claim: SAR is associated with lower, not higher, mortality rates.

### 5.2 Counterfactual Analysis (Figure 2)

The three panels of Figure 2 show observed log crossings (black line) against the counterfactual prediction (dashed blue line) with 95% credible intervals (light blue ribbon):

**Panel A (Mare Nostrum):** The observed series remains within the confidence intervals of the counterfactual throughout the entire post-intervention period (October 2013 onward). There is a slight bump where observed crossings exceed the prediction around Spring-Summer 2014, but this difference is not statistically significant. **Result: No significant effect.**

**Panel B (NGO SAR):** The observed and predicted counterfactual series follow each other closely throughout the post-intervention period (November 2014 onward), including during the peak NGO activity in 2015-2016. The model predicts slightly *higher* crossings than were observed in later years, suggesting that, if anything, crossings were lower than expected during NGO SAR. **Result: No significant effect.**

**Panel C (EU-Libya Cooperation):** Here, the counterfactual diverges from the observed series: the model predicts substantially *higher* crossings than were actually observed after February 2017. The observed series drops below the lower confidence bound, indicating that the EU-Libya cooperation (pushbacks, LCG SAR zone extension) **significantly reduced** crossing attempts. **Result: Significant negative effect.**

### 5.3 Pointwise Effects (Figure 3)

Figure 3 shows the difference between observed and predicted crossings (pointwise effects) over time:

- **Panels A and B:** Effects fluctuate around zero with no systematic pattern, confirming that Mare Nostrum and NGO SAR did not shift crossing attempts away from what the push-pull factors would predict.  
- **Panel C:** Effects are consistently negative after February 2017, with the credible interval mostly below zero. This confirms a statistically significant reduction in crossings due to EU-Libya cooperation.  

### 5.4 Summary of Causal Effects

| Intervention | Period | Absolute Effect | Significant? |
|---|---|---|---|
| Mare Nostrum | Oct 2013 -- Oct 2014 | Near zero | No (p > 0.05) |
| NGO SAR | Nov 2014 -- Feb 2017 | Near zero | No (p > 0.05) |
| EU-Libya cooperation | Feb 2017 onward | Negative | Yes (p < 0.05) |

### 5.5 Interpretation

The paper concludes that **search-and-rescue does not function as a pull factor for migration**. The observed flow of migrants can be fully explained by exogenous push-and-pull factors (conflicts, economics, weather, etc.) without any additional effect attributable to the presence of SAR boats. In other words, migrants would have crossed regardless of whether rescue operations were active.

The EU-Libya cooperation period, by contrast, *did* significantly reduce crossings -- but the authors emphasize that this came "at the expense of significant deterioration of the human rights situation of migrants in Libya" (detention centers, pushbacks, violence). The reduction reflects deterrence through increased risk and interception, not a resolution of the underlying migration drivers.

---

## 6. Model Validation and Robustness (Supplementary Materials)

### 6.1 Time-Series Decomposition (Figure S2)

An additive decomposition of `log(crossings_CMR)` confirms strong seasonal and trend components. The trend peaks around 2014-2016 and declines thereafter. The seasonal component shows regular summer highs and winter lows, changing amplitude over time (non-stationary seasonality). This motivated the choice of BSTS over methods that assume stationary seasonality.

The decomposition is performed in the code as:
```r
y <- ts(log(df$crossings_CMR), start = c(2009, 1), frequency = 12)
decomp <- decompose(y, type = "additive")
```

### 6.2 Structural Breakpoints (Figures S4-S5)

Using the Bai-Perron test for multiple structural breaks (`strucchange::breakpoints()`), the paper identifies four breakpoints in the log-crossings series:  
1. ~2011: Start of the Syrian civil war and associated migration increase.  
2. ~2014: Shortly after Mare Nostrum, coinciding with the European migration crisis peak.  
3. ~2017: Start of EU-Libya cooperation.  
4. ~2020: COVID-19 pandemic effects.  

These breakpoints roughly align with the three intervention periods but do not coincide exactly, which the authors note is consistent with the breakpoints reflecting broader geopolitical shifts rather than SAR-specific effects.

For deaths, the breakpoints occur at different times, reflecting the distinct dynamics of mortality vs. crossing attempts.

Code (for crossings):
```r
bp_crossings <- breakpoints(log(y_crossings) ~ 1, h = 12)
```
The `h = 12` parameter sets a minimum segment length of 12 months.

### 6.3 Alternative State-Space Specifications (Figures S6-S7, Table S3)

The paper tests 8 different state-space model specifications to assess robustness:

| Model | State Components |
|---|---|
| LLT | Local Linear Trend |
| SLLT | Semi-Local Linear Trend |
| LL | Local Level |
| LLT + AR | Local Linear Trend + AutoAR(3) |
| SLLT + AR | Semi-Local Linear Trend + AutoAR(3) |
| LL + AR | Local Level + AutoAR(3) |
| LL + SLLT + AR | Local Level + Semi-Local Linear Trend + AutoAR(3) |
| Default | CausalImpact default (Local Level) |

The mathematical definitions of these components are:

- **Local Level:** `alpha_{t+1} = alpha_t + epsilon_t` -- a random walk. Captures a stochastic level that drifts over time.  
- **Local Linear Trend:** `mu_{t+1} = mu_t + delta_t + epsilon_t`, `delta_{t+1} = delta_t + eta_t` -- both level and slope follow random walks. Captures evolving trends.  
- **Semi-Local Linear Trend:** Same as LLT but the slope follows a stationary AR(1) process centered on D: `delta_{t+1} = D + phi*(delta_t - D) + eta_t`. This prevents the trend from drifting indefinitely, making it better for long-horizon predictions.  
- **AutoAR(3):** An autoregressive component with up to 3 lags, allowing the model to capture short-term serial correlation.  

The authors deliberately excluded an explicit seasonality component (`AddSeasonal`) from the final model because the seasonality in the data is non-stationary (changing over time), and adding it worsened fit. Instead, seasonality is captured indirectly through the weather covariates and the `month`/`quarter`/`semester` indicators.

**Cross-validation results** (Table S3): 5-fold time-series cross-validation (expanding window) shows that simpler models (LLT, SLLT, LL) generally outperform more complex ones on RMSE, MAPE, and MASE metrics. The Default specification (Local Level) is competitive but not always best. The horizon is 3 months for Mare Nostrum (short pre-period) and 6 months for the other two interventions.

**Effect size comparison** (Figure S7): All 8 model specifications produce similar point estimates for the causal effect across all three interventions. The key difference is in the width of credible intervals (more complex models produce wider CIs). The qualitative conclusion -- no effect for SAR, significant negative effect for EU-Libya -- is robust across all specifications.

**Cumulative prediction errors** (Figure S6): The `CompareBstsModels()` function compares the one-step-ahead prediction error accumulated over the sample period. The Default model shows competitive performance, with errors of similar magnitude to the more complex specifications.

### 6.4 Residual Diagnostics (Figures S8-S9)

**ACF plots** (Figure S8): The autocorrelation function of the posterior residuals shows the characteristic decay toward zero, indicating that the model adequately captures the serial correlation structure. The boxplots show the distribution of ACF values across MCMC draws.

**Q-Q plots** (Figure S9): Quantile-quantile plots of the posterior residuals against the standard normal distribution show reasonable fit, with some deviation in the tails (likely due to the heavy-tailed nature of migration data and variable omission). The authors conclude the model provides "a good fit in general and not a strong deviation from the assumed normal distribution."

### 6.5 Covariate Inclusion Probabilities (Figure S10)

The spike-and-slab prior produces **inclusion probabilities** for each covariate -- the fraction of MCMC iterations in which the covariate was included in the regression. Figure S10 shows the top covariates by inclusion probability for each intervention:

- **Mare Nostrum (Panel A):** Airport flows (Lithuania, Luxembourg), weather variables, commodity prices, conflict indicators from various African countries.  
- **NGO SAR (Panel B):** Exchange rates (APK, LBL), airport flows, weather indicators, and some commodity prices dominate.  
- **EU-Libya (Panel C):** Commodity prices (PWOTAN, PNFUEL), airport flows, unemployment rates (Greece, Bulgaria), and exchange rates.  

Notably, no single covariate dominates the model -- inclusion probabilities are generally moderate (10-50%), reflecting the Bayesian model averaging across many possible predictor combinations.

### 6.6 Mediation Analysis (Sensitivity Test)

The most important robustness check addresses the **exogeneity assumption**: could any control series be affected by the interventions? The authors focus on labor market indicators (European unemployment rates, African job searches) as the most plausible mediators. If SAR is a pull factor, more migrants arrive, changing labor supply and unemployment.

Two tests are performed:

1. **Exclusion of labor market indicators** (Figure S11): Re-running the BSTS models without unemployment and Google Trends job search variables produces "very similar effects" with "only small variations in the predicted counterfactual."  

2. **Granger causality tests** (Tables S4-S6): For each control series in the model, a bivariate Granger test checks whether the intervention Granger-causes the control series (which would violate exogeneity). With Bonferroni correction across hundreds of tests, only 2-8 control series per intervention show significant associations -- and none pass the reverse causality test (the control series does not Granger-cause the intervention). The authors interpret these as "spurious correlations."  

---

## 7. Code Architecture

### 7.1 Original Code (`sar_migration_aleja/code/`)

The original repository contains 22 R scripts organized as follows:

**Data construction layer** (22 scripts):  
- 4 arrival scripts, 3 death scripts, 9 covariate scripts assemble raw data sources.  
- `dates_search_and_rescue_EU_NGOS_dates.R` defines daily binary operation indicators for 22 SAR vessels (SeaWatch 1/2/3/4, Lifeline, Sea-Eye, Open Arms, Mare Liberum, Mediterranea, SMH, Louise Michel, Refugee Rescue, MOAS, Jugend Rettet, MSF/SOS, Save the Children, MSF Bourbon Argos, MSF Dignity1, MSF Vos Prudence, Resqship, Lifeboat) and 8 EU operations (FRONTEX POSEIDON, HERA, TRITON, THEMIS, MINERVA; EUNAVFOR SOPHIA, IRINI; Mare Nostrum; Libyan Coast Guard).  
- `create_data_frame.R` merges everything into `df.RDS` (~4,726 columns).  

**Modeling layer** (4 scripts):  
- `run_bsts_models.R`: The main analysis script. Fits the 3 CausalImpact models with 10,000 MCMC iterations, produces plots, and saves results. Also includes exploratory analyses (spike-and-slab variable pre-screening, deaths models, custom BSTS specifications) that are not part of the final paper results.  
- `model_validation_marenostrum.R`, `model_validation_sarngos.R`, `model_validation_lcg_eu.R`: Fit 8 alternative state-space specifications per intervention with 2,500 MCMC iterations, run 5-fold cross-validation, produce cumulative error comparisons, effect size plots, and diagnostic figures.  

### 7.2 Replication Code (Root Directory)

The replication consolidates the original's scattered analysis into two self-contained scripts:

**`replicate_main.R`** -- replicates Figures 1-3 and the core results:  
1. Loads `df.RDS` and constructs `crossings_CMR` and `mortality_rate`.  
2. Applies the variable filtering (removes high-order lags, outcome-related variables, etc.).  
3. Creates Figure 1 (descriptive time series with intervention shading).  
4. Defines the three intervention periods (identical to the original).  
5. Fits three CausalImpact models with `niter=10000, max.flips=100, seed=270488`.  
6. Extracts results and creates Figures 2-3 (counterfactual and pointwise effects).  
7. Prints a summary of absolute effects, relative effects, p-values, and significance.  

**`replicate_supplementary.R`** -- replicates Figures S1-S10:  
1. Figure S1: NGO operations timeline (daily resolution, 22 vessels).  
2. Figure S2: Additive decomposition of log crossings.  
3. Figure S3: Death ratio over time.  
4. Figures S4-S5: Structural breakpoints (deaths and crossings).  
5. Figures S6-S7: 8-model comparison (cumulative errors and effect sizes).  
6. Figures S8-S9: ACF and Q-Q diagnostic plots.  
7. Figure S10: Covariate inclusion probabilities.  

The supplementary script requires pre-fitted models from `replicate_main.R` (loaded from `replicated_results/models/`). For the 8-model comparison (S6-S7), it runs 7 alternative BSTS specifications with `niter=2500` (reduced for computational tractability) plus the default from `replicate_main.R`.  

### 7.3 Output Structure

```
replicated_results/
  models/                          # Fitted CausalImpact model objects (.RDS)
    model_marenostrum_full.RDS
    model_sarngos_full.RDS
    model_sarlibya_full.RDS
  figures/
    png/                           # Main figures at 300 DPI
      Figure1_full_model.png
      Figure2_full_model.png
      Figure3_full_model.png
    pdf/                           # Main figures in vector format
    supplementary/
      png/                         # Figures S1-S10
      pdf/
  data/                            # Extracted results as data frames
    model_results_all_full.RDS
  validation/                      # Model comparison outputs
    effs_marenostrum.RDS           # Effect sizes from 8 model specifications
    effs_sarngos.RDS
    effs_sarlibya.RDS
    compare_mod_long_*.RDS         # Cumulative error data
```

---

## 8. Reproducing the Analysis

### 8.1 Prerequisites

**R packages required:**. 
- Core: `CausalImpact`, `bsts`, `tidyverse`, `lubridate`, `scales`, `gridExtra`. 
- Supplementary: `strucchange`, `forecast`. 
- Data construction (only if rebuilding `df.RDS`): `janitor`, `countrycode`, `imputeTS`, `readxl`, `readstata13`, `Hmisc`, `priceR`, `gtrendsR`. 

### 8.2 Running the Replication

```bash
# Step 1: Fit the three main BSTS models and generate Figures 1-3
# WARNING: This takes ~8-15 hours
Rscript replicate_main.R

# Step 2: Generate supplementary figures S1-S10
# Requires fitted models from Step 1
Rscript replicate_supplementary.R
```

### 8.3 Key Parameters to Understand

| Parameter | Value | Purpose |
|---|---|---|
| `set.seed(270488)` | Fixed | Ensures MCMC reproducibility |
| `niter = 10000` | Main models | Number of MCMC iterations (10,000 for main, 2,500 for validation) |
| `max.flips = 100` | All models | Max spike-and-slab inclusion changes per MCMC iteration |
| `dynamic.regression = FALSE` | All models | Static (not time-varying) regression coefficients |
| `standardize.data = TRUE` | All models | Standardize predictors before fitting |
| `alpha = 0.05` | All models | 95% credible intervals |

---

## 9. Implications and Limitations

### 9.1 Policy Implications

The paper provides strong evidence against the pull factor claim:  
- SAR operations (both state-led and NGO) did not increase crossing attempts beyond what exogenous factors would predict.  
- Migration flows are driven by structural push-and-pull factors (conflicts, economics, weather) that operate independently of SAR.  
- Deterrence policies (EU-Libya cooperation / pushbacks) did reduce crossings, but at a severe humanitarian cost.  

### 9.2 Methodological Contribution

The paper demonstrates that BSTS + CausalImpact is a powerful tool for causal inference in time-series settings where:  
- The outcome exhibits strong autocorrelation, trends, and seasonality.  
- The predictor space is very high-dimensional.  
- Traditional approaches (DiD, synthetic control) are infeasible because there is no suitable "control unit."  

### 9.3 Limitations Acknowledged by the Authors

1. **Post-intervention contamination.** The post-intervention periods for earlier interventions overlap with later interventions. For example, the post-period for Mare Nostrum includes the NGO SAR and EU-Libya periods. This means the models cannot isolate the effect of one intervention from the delayed effects of another.  

2. **Confidence interval widening.** The counterfactual predictions become increasingly unreliable as we move further from the intervention start. Confidence intervals widen substantially after ~2 years.  

3. **Aggregate data.** The analysis uses monthly aggregate counts, precluding individual-level causal statements. The ecological fallacy applies: aggregate non-effects do not rule out micro-level behavioral responses.  

4. **Death underreporting.** The true number of deaths in the CMR is unknown. Underreporting may affect the crossing attempts measure, though deaths constitute ~4% of total crossings, limiting the impact.  

5. **Endogeneity of some covariates.** While Granger causality tests show minimal evidence of feedback effects, the exogeneity assumption cannot be definitively proven. Labor market indicators are the most plausible channel for endogeneity (addressed by the sensitivity analysis in Figure S11).  

6. **No explicit seasonality component.** The non-stationary seasonality in the data made an explicit seasonal state component counterproductive. Seasonality is instead captured indirectly through weather covariates and calendar indicators, which is somewhat unconventional.  

---

## 10. Key References

- Brodersen, K.H., Gallusser, F., Koehler, J., Remy, N., & Scott, S.L. (2015). Inferring causal impact using Bayesian structural time-series models. *Annals of Applied Statistics*, 9(1), 247-274. (The CausalImpact methodology.). 
- Scott, S.L. & Varian, H.R. (2014). Predicting the present with Bayesian structural time series. *International Journal of Mathematical Modelling and Numerical Optimisation*, 5(1-2), 4-23. (The BSTS framework.). 
- George, E.I. & McCulloch, R.E. (1993). Variable selection via Gibbs sampling. *Journal of the American Statistical Association*, 88(423), 881-889. (Spike-and-slab priors.). 
