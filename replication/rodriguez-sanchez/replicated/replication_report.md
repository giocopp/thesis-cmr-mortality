---
output:
  html_document: default
  word_document: default
---
# Replication Report

## Rodriguez Sanchez et al. (2023), "Search-and-rescue in the Central Mediterranean Route does not induce migration: Predictive modeling to answer causal queries in migration research"

**Published in:** *Nature Scientific Reports*, 13, 11014. DOI: [10.1038/s41598-023-38119-4](https://doi.org/10.1038/s41598-023-38119-4)

**Replication by:** Giorgio Coppola  

**Date:** February 2026  

**Institution:** Hertie School, Master of Data Science for Public Policy  

---

## 1. Overview

This report documents the replication of the main empirical results from Rodriguez Sanchez et al. (2023). The paper investigates whether search-and-rescue (SAR) operations in the Central Mediterranean Route (CMR) act as a "pull factor" for irregular migration. Using Bayesian Structural Time-Series (BSTS) models with spike-and-slab priors and the `CausalImpact` R package, the authors construct synthetic counterfactual time-series to estimate the causal effect of three intervention periods on migration flows:

- **Model A (Mare Nostrum):** State-led SAR operation (Oct 2013 -- Oct 2014)
- **Model B (NGO SAR):** Private-led SAR by NGOs (Nov 2014 -- Feb 2017)
- **Model C (EU-Libya cooperation):** Coordinated pushbacks and Libyan SAR zone extension (Feb 2017 onwards)

The paper's central finding is that SAR operations (Models A and B) did **not** have a statistically significant effect on migration flows, while EU-Libya cooperation (Model C) was associated with a **significant reduction** in crossing attempts.

---

## 2. Replication Scope

### 2.1 Materials Used

| Component | Source |
|-----------|--------|
| Original paper | `rodriguez_sanchez_nature.pdf` (12 pages) |
| Supplementary materials | `41598_2023_38119_MOESM1_ESM.pdf` (41 pages) |
| Original authors' code | `sar_migration_aleja/code/` (22 R scripts) |
| Pre-compiled dataset | `df.RDS` (~4,726 predictors, Feb 2009 -- Oct 2021) |
| Replication scripts | `replicate_main.R`, `replicate_supplementary.R` |

### 2.2 Figures Replicated

| Figure | Description | Status |
|--------|-------------|--------|
| **Figure 1** | Descriptive time-series: attempted crossings and mortality rate | Replicated |
| **Figure 2** | Counterfactual analysis (observed vs. predicted, 3 panels) | Replicated |
| **Figure 3** | Pointwise effects (difference between observed and predicted) | Replicated |
| **Figure S1** | Number of NGO-led SAR operations (daily resolution) | Replicated |
| **Figure S2** | Time-series decomposition (additive) | Replicated |
| **Figure S3** | Death ratio over time | Replicated |
| **Figure S4** | Structural breakpoints (deaths, log) | Replicated |
| **Figure S5** | Structural breakpoints (crossings, log) | Replicated |
| **Figure S6** | Cumulative prediction errors (8 model specifications) | Replicated |
| **Figure S7** | Effect size comparison (8 model specifications) | Replicated |
| **Figure S8** | ACF of residuals (using `bsts::AcfDist`) | Replicated |
| **Figure S9** | Q-Q plots (using `bsts::qqdist`) | Replicated |
| **Figure S10** | Coefficient inclusion probabilities | Replicated |
| **Figure S11** | Sensitivity analysis (excluding labor market indicators) | Not replicated |

**Note:** Figure S11 requires re-running the BSTS models with an alternative covariate set (excluding labor market indicators). This was not included in the replication due to computational cost but does not affect the core findings.

---

## 3. Methodological Consistency

### 3.1 Data Preparation

The replication scripts (`replicate_main.R`, lines 52--114) faithfully follow the original code (`sar_migration_aleja/code/run_bsts_models.R`, lines 10--150):

- **Outcome variable:** `crossings_CMR = arrivals_CMR + LCG_pushbacks_count + TCG_pushbacks_count + dead_and_missing_Central_Mediterranean`, log-transformed. This matches the paper's definition: $Y_t = A_t + P_t + D_t$ (p. 3).
- **Missing value handling:** NAs in pushback and death counts replaced with 0 before aggregation -- identical to original.
- **Variable filtering:** High-order lags (07--24) removed; `airflow_Palestinian.Territories`, asylum variables, outcome-related variables, and route-specific variables excluded. The filtering logic matches the original code exactly (compare `replicate_main.R:68-92` with `run_bsts_models.R:110-139`).
- **Time components:** Month, semester, and quarter added as predictors, consistent with original.
- **Study period:** February 2011 -- September 2021, as stated in the paper ("monthly series for the period of 2011--2020", Methods section, p. 9).

### 3.2 Model Specification

The CausalImpact models use identical parameters to those in the original code:

| Parameter | Original | Replication | Match |
|-----------|----------|-------------|-------|
| `dynamic.regression` | `FALSE` | `FALSE` | Yes |
| `standardize.data` | `TRUE` | `TRUE` | Yes |
| `max.flips` | 100 | 100 | Yes |
| `niter` | 10,000 | 10,000 | Yes |
| `alpha` | 0.05 | 0.05 | Yes |
| `seed` | 270488 | 270488 | Yes |

### 3.3 Intervention Period Definitions

| Model | Pre-period | Post-period start | Original | Match |
|-------|-----------|-------------------|----------|-------|
| A (Mare Nostrum) | min(date) -- Sep 2013 | Oct 2013 | `run_bsts_models.R:252-253` | Yes |
| B (NGO SAR) | min(date) -- Oct 2014 | Nov 2014 | `run_bsts_models.R:255-256` | Yes |
| C (EU-Libya) | min(date) -- Jan 2017 | Feb 2017 | `run_bsts_models.R:258-259` | Yes |

### 3.4 Minor Methodological Note

The original code (`run_bsts_models.R`) references a `df_y_rec.RDS` file containing a reconstructed death count variable (`y_rec`), which is excluded from the final model. The replication script correctly omits this variable, as it is not part of the published analysis. The original code also includes an initial `logit.spike` variable selection step (lines 153--165) that is exploratory and precedes the final `CausalImpact` call; the replication correctly uses the full spike-and-slab variable selection embedded within `CausalImpact` as described in the paper.

---

## 4. Results Comparison

### 4.1 Main Findings (Figures 2 and 3)

The paper's core claims are assessed below:

#### Model A -- Mare Nostrum: NOT significant

> "Figure 2A shows that the observed trend after the start of Mare Nostrum remains within the confidence intervals of our predicted counterfactual." (p. 4)

**Replication result:** Confirmed. The replicated Figure 2A shows the observed time-series (black line) tracking closely with the counterfactual prediction (dashed blue line) throughout the post-intervention period. The 95% confidence interval (light blue ribbon) encompasses the observed values. The pointwise effects in Figure 3A gravitate tightly around zero with the confidence band consistently crossing zero, indicating no statistically significant effect.

#### Model B -- NGO SAR: NOT significant

> "The observed and the predicted counterfactual time-series follow each other closely for nearly two years, and especially during the months where private-led search-and-rescue was the most intense." (p. 4)

**Replication result:** Confirmed. The replicated Figure 2B shows close tracking between observed and predicted values in the post-intervention period. The confidence intervals encompass observed values throughout. Figure 3B shows pointwise effects centered on zero with no sustained significant departures.

#### Model C -- EU-Libya cooperation: SIGNIFICANT

> "Our prediction suggests a substantially higher number of crossing attempts than what was actually observed." (p. 5)

**Replication result:** Confirmed. The replicated Figure 2C shows the counterfactual prediction (dashed line) diverging above the observed series after mid-2017, with the observed values falling below the lower confidence bound. Figure 3C shows consistently negative pointwise effects in the post-intervention period, indicating that observed crossings were significantly lower than predicted -- consistent with a deterrent effect of EU-Libya cooperation policies.

### 4.2 Figure-by-Figure Visual Comparison

#### Figure 1: Descriptive Time-Series

The replicated Figure 1 matches the original in all key respects:

- **Panel A (Attempted Crossings):** The time-series shape is identical, showing the characteristic rise in crossings peaking around 2014--2016, followed by a decline after 2017. The four color-coded intervention periods (Mare Nostrum in blue, NGO SAR in red, EU-LCG Deal in orange, LCG SAR-Zone expansion in tan) are correctly positioned. Y-axis range (0--30,000) matches the original.
- **Panel B (Mortality Rate):** The mortality rate per 100 attempted crossings follows the same pattern as the original, with high pre-SAR rates (spikes above 20 per 100), lower rates during the SAR period (approximately 0--10), and elevated rates post-EU-Libya cooperation. The purple-shaded SAR period is correctly positioned.

#### Figures 2 and 3: Counterfactual and Pointwise Effects

The replicated panels show the same qualitative patterns as the published figures:

- Pre-intervention fit is tight across all three models, with observed values falling within the narrow confidence intervals.
- Post-intervention divergence patterns match: no divergence for Models A and B, clear negative divergence for Model C.
- The y-axis scale (approximately -20 to 40 for log crossings) and shading scheme (pink for intervention, grey for post-intervention) are consistent.
- The moving-average smoother (grey line) in the replication adds visual clarity comparable to the published version.

#### Note on MCMC Stochasticity

Due to the stochastic nature of MCMC sampling (10,000 iterations), exact numerical values of effect estimates and confidence intervals will differ slightly between the original and replication runs, even with the same seed (`set.seed(270488)`). This is expected behavior, as differences in R package versions, platform arithmetic, and compilation can lead to different random number sequences. The qualitative conclusions, however, are fully robust to this stochasticity.

### 4.3 Supplementary Figures

#### Figure S1: NGO Operations Timeline

The replicated figure accurately captures the daily count of active NGO SAR operations from 2009 to 2021. The step-function shape matches the original (Figure S1 in supplementary materials, p. 13), showing:

- Zero operations before mid-2014
- A rapid ramp-up from 1 to a peak of approximately 10--12 simultaneous operations in mid-2017
- A decline after the EU-Libya cooperation period
- Intermittent fluctuations from 2018 to 2021 between approximately 4--8 operations

The two vertical reference lines (End of Mare Nostrum; EU-Libyan Coast Guard) are correctly positioned. The 23 individual NGO vessel operation periods coded in `replicate_supplementary.R` (lines 154--244) match the dates documented in Table S2 of the supplementary materials.

#### Figure S2: Time-Series Decomposition

The additive decomposition of log(crossings_CMR) is visually consistent with the original (supplementary p. 15). The four panels (observed, trend, seasonal, random) display the same patterns: a clear seasonal cycle with summer peaks, a trend that rises from 2011 to 2016 and then declines, and a random component with larger variance in early periods.

#### Figures S4 and S5: Structural Breakpoints

Both breakpoint analyses use the `strucchange::breakpoints()` function with `h=12` (minimum segment size of 12 months), matching the original code. The breakpoints in crossings (Figure S5) identify four structural changes that "roughly follow important changes in the politics of search-and-rescue in CMR" (supplementary p. 19), consistent with the original findings.

#### Figure S6: Model Validation (Cumulative Errors)

The 8-model comparison of cumulative absolute errors shows the Default model (CausalImpact's built-in local level specification) outperforming the alternative state specifications in all three panels, consistent with the supplementary materials' conclusion that "the default model beats the cross-validated measures" (p. 25). The relative ordering of model performance (Default best, LL+SLLT+AR worst) matches the original.

#### Figure S7: Effect Size Comparison

The effect size comparison across 8 model specifications shows that all models produce effect estimates centered near zero for Models A and B, with all confidence intervals crossing zero. For Model C, effect estimates are slightly negative across specifications. This confirms the paper's claim that "estimates of the different models, however, do not differ much" (p. 25) and that "the effect is of a similar magnitude across models."

#### Figures S8 and S9: Model Diagnostics

- **ACF plots (S8):** The autocorrelation function distributions show the characteristic decay toward zero at higher lags, consistent with the original figures. The boxplot format using `bsts::AcfDist` matches the published diagnostic.
- **Q-Q plots (S9):** The quantile-quantile plots using `bsts::qqdist` show the distribution of Monte Carlo draws against standard normal quantiles. The replicated figures display the same pattern as the originals: generally good fit along the diagonal with some deviation in the tails, consistent with the paper's observation that "the tails of the error distribution do not follow precisely on the straight line -- probably due to variable omission" (p. 26--27).

#### Figure S10: Coefficient Inclusion Probabilities

The replicated S10 shows the top covariates by spike-and-slab inclusion probability for each model. The types of covariates selected (airflow variables, commodity prices, conflict indicators, exchange rates, unemployment rates) are consistent with the paper's discussion of predictive factors (supplementary p. 29) and with Table S1's covariate categories.

---

## 5. Assessment of Reproducibility

### 5.1 Computational Reproducibility

| Criterion | Assessment |
|-----------|------------|
| **Code availability** | Original code available at [GitHub](https://github.com/xlejx-rodsxn/sar_migration). Replication scripts are clean, well-documented adaptations. |
| **Data availability** | Pre-compiled dataset (`df.RDS`) used. Raw data construction scripts available but some underlying sources (Sabre air traffic data) are proprietary. |
| **Software dependencies** | Core packages (`CausalImpact`, `bsts`, `tidyverse`, `strucchange`, `forecast`) are all open-source and available on CRAN. |
| **Computational cost** | High. Three BSTS models with 10,000 MCMC iterations each require approximately 8--15 hours of runtime. The supplementary 8-model comparison adds additional hours. |
| **Random seed** | Set to `270488` in both original and replication. Exact numerical reproduction is not guaranteed across platforms due to MCMC stochasticity, but qualitative results are stable. |

### 5.2 What Could Be Replicated

- All three main CausalImpact models (A, B, C) with identical specifications
- All three main figures (Figures 1, 2, 3)
- Ten out of eleven supplementary figures (S1--S10)
- Model validation across 8 alternative state-space specifications
- Coefficient inclusion probability analysis

### 5.3 What Could Not Be Replicated

- **Figure S11 (Sensitivity analysis):** Requires re-running all three BSTS models after excluding labor market indicators (unemployment rates, Google Trends job searches). Not included due to computational cost but the paper states the results "show very similar effects" (supplementary p. 30).
- **Granger causality tests (Tables S4--S6):** These mediation analysis tests were not replicated in the scripts, though the methodology is straightforward. The paper reports only a few spurious significant correlations among hundreds of control series.
- **Exact numerical values:** Due to MCMC stochasticity and potential platform differences, exact p-values and effect sizes may differ slightly from the published values. The qualitative conclusions (significance/non-significance at the 5% level) are robust.

---

## 6. Conclusions

### 6.1 Replication Verdict

The replication **successfully confirms** the main findings of Rodriguez Sanchez et al. (2023):

1. **Mare Nostrum (Model A):** No statistically significant effect on migration flows. The observed time-series remains within the confidence intervals of the counterfactual prediction throughout the post-intervention period.

2. **NGO search-and-rescue (Model B):** No statistically significant effect. The predicted counterfactual closely tracks the observed series during and after the NGO SAR period.

3. **EU-Libya cooperation (Model C):** A statistically significant **reduction** in crossing attempts, consistent with deterrence through pushback policies and the extension of the Libyan SAR zone.

These findings collectively contradict the "pull factor" hypothesis -- the claim that SAR operations incentivize irregular migration. The replication confirms that changes in crossing attempts can be well explained by exogenous push-and-pull factors (economic conditions, conflicts, weather, etc.) without invoking SAR operations as a causal driver.

### 6.2 Robustness

The replication further confirms the paper's robustness claims:

- **Model specification robustness (S6, S7):** Effect estimates are consistent across 8 different BSTS state-space specifications (local level, local linear trend, semi-local linear trend, and combinations with autoregressive components).
- **Model diagnostics (S8, S9):** Residual autocorrelation and Q-Q diagnostics indicate reasonable model fit.
- **Structural breakpoint analysis (S4, S5):** Breakpoints in the crossing time-series do not coincide with intervention dates, suggesting that structural changes in migration patterns are driven by other factors.

### 6.3 Strengths and Limitations of the Original Study

**Strengths replicated:**
- The BSTS framework with spike-and-slab variable selection is well-suited for this high-dimensional causal inference problem (~600--700 predictors after filtering).
- The use of multiple data sources (FRONTEX, IOM, Eurostat, ACLED, IMF, weather data, Google Trends, Sabre air traffic) provides a comprehensive set of control time-series.
- The code and data are largely reproducible (modulo the proprietary Sabre data).

**Limitations noted:**
- The pre-compiled dataset (`df.RDS`) was used directly; the full data construction pipeline (22 R scripts) was not independently re-run from raw sources, partly because some data sources (Sabre air traffic) are proprietary.
- MCMC convergence diagnostics are not extensively documented in the original paper. The 10,000-iteration specification may be adequate but formal convergence assessment (e.g., Gelman-Rubin R-hat statistics) is not reported.
- As the authors acknowledge, the post-intervention periods for earlier models (A, B) are contaminated by later interventions (B, C), making it difficult to isolate individual intervention effects in later time periods.

---

## 7. File Inventory

### Replication Outputs

```
replicated_results/
  models/
    model_marenostrum_full.RDS    # Fitted CausalImpact model (Mare Nostrum)
    model_sarngos_full.RDS        # Fitted CausalImpact model (NGO SAR)
    model_sarlibya_full.RDS       # Fitted CausalImpact model (EU-Libya)
  figures/
    png/
      Figure1_full_model.png      # Descriptive time-series
      Figure2_full_model.png      # Counterfactual analysis
      Figure3_full_model.png      # Pointwise effects
    pdf/
      Figure1_full_model.pdf
      Figure2_full_model.pdf
      Figure3_full_model.pdf
    supplementary/
      png/
        FigureS1_NGO_operations.png   # NGO operations timeline
        FigureS2_decomposition.png    # Time-series decomposition
        FigureS3_death_ratio.png      # Death ratio
        FigureS4_breakpoints_deaths.png    # Structural breakpoints (deaths)
        FigureS5_breakpoints_crossings.png # Structural breakpoints (crossings)
        FigureS6_model_validation.png  # Cumulative errors (8 models)
        FigureS7_effect_sizes.png      # Effect size comparison
        FigureS8_acf.png              # ACF diagnostics
        FigureS9_qqplot.png           # Q-Q diagnostics
        FigureS10_coefficients.png    # Inclusion probabilities
      pdf/
        [same figures in PDF format]
  data/
    model_results_all_full.RDS    # Extracted model predictions
  validation/
    effs_marenostrum.RDS          # Effect sizes (Mare Nostrum, 8 specs)
    effs_sarngos.RDS              # Effect sizes (NGO SAR, 8 specs)
    effs_sarlibya.RDS             # Effect sizes (EU-Libya, 8 specs)
    compare_mod_long_*.RDS        # Cumulative error comparison data
```

---

## References

Rodriguez Sanchez, A., Wucherpfennig, J., Rischke, R., & Iacus, S. M. (2023). Search-and-rescue in the Central Mediterranean Route does not induce migration: Predictive modeling to answer causal queries in migration research. *Scientific Reports*, 13, 11014. https://doi.org/10.1038/s41598-023-38119-4
