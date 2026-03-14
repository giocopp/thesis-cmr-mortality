# Session Notes — 26 February 2026

## Full Session: Six-Model Expansion, Volume Channel Analysis, Moon Removal, Covariate Correlations

---

## 1. Presentation Rendering Fix

The Quarto presentation (`progress_presentation.qmd`) failed to render
because absolute image paths were invalid under LuaTeX. Converted all
image paths to relative paths. Rendering verified.

---

## 2. Death Count Models Added for Models A and B

Previously, death count models existed only for Model C. We added
`A_deaths` and `B_deaths` to produce a full 6-model comparison.

### Changes made

- **`R/02_data_prepare.R`** — added `A_deaths` and `B_deaths`
  intervention specs to `define_interventions()`
- **`_targets.R`** — added 4 new targets: `model_a_deaths`,
  `model_b_deaths`, `placebos_a_deaths`, `placebos_b_deaths`;
  updated figure and snapshot calls to pass all 6 models
- **`R/05_figures.R`** — expanded counterfactual, pointwise, and
  placebo figures from 1-panel or 3-panel to full 6-panel layouts
- **`data/raw/df.RDS`** — fixed broken symlink (folder name had
  hyphens vs spaces)

### 6-model results

| Model | Outcome | p-value | Cum. effect | RMSE/SD |
|---|---|---:|---:|---:|
| A | rate | 0.085 | +54.7 | 0.75 |
| A | deaths | **0.0004** | **+127.1** | **0.62** |
| B | rate | 0.002 | +100.4 | 0.86 |
| B | deaths | **0.0001** | **+127.4** | **0.77** |
| C | rate | **0.008** | +45.7 | 0.97 |
| C | deaths | 0.176 | +18.1 | 0.92 |

**Key finding — reversed pattern:** A/B deaths are highly significant
but C deaths fails; C rate is significant but A/B rates fail. This
reveals a volume channel: A/B death models detect the MN/NGO SAR
crossing surge, not a per-crossing danger change. C rate detects
per-crossing danger because the rate removes volume.

---

## 3. Extension 1 (Naive Approach) — Complete Results

Ran curated-specification (809 predictors) death count models for
A and B to complete the Extension 1 comparison.

| Model | Outcome | p-value | Rel. effect | Significant? |
|---|---|---:|---:|---|
| A | rate | 0.110 | −102% | No |
| A | deaths | 0.436 | −193% | No |
| B | rate | 0.387 | −96% | No |
| B | deaths | 0.096 | −49% | Marginal |
| C | rate | **0.029** | **−71%** | **Yes** |
| C | deaths | **0.012** | **−57%** | **Yes** |

Top predictors: Malta unemployment (21%), EU unemployment (15%),
Nigeria violence (14%). ERA5 sea variables max 4.4%.

**Problem:** these covariates predict crossing *volume* (unemployment,
conflict), not per-crossing danger. Model C deaths is significant
because volume-predicting covariates forecast continued high crossings
while actual crossings dropped — a volume artifact, not a danger
signal.

Models saved to `Extension-1-BSTS-mortality/results/models/`.

---

## 4. Investigating the Reversed Pattern

### Why A/B deaths work

Pre-period is the low-volume era (~3,600 crossings/mo). After MN/NGO
SAR, crossings jump to ~14,000 → deaths jump mechanically (+1.86 log
units). The model captures this volume shift well.

### Why C deaths fails

The C death-count pre-period (Feb 2011–Jan 2017) **contains** the
MN/NGO SAR volume surge as a structural break:

| Sub-period | Log deaths | Log rate |
|---|---:|---:|
| Pre-MN (Feb 11–Sep 13) | 2.98 | −0.55 |
| Mare Nostrum (Oct 13–Oct 14) | 4.10 (+1.12) | −0.63 (flat) |
| NGO SAR (Nov 14–Jan 17) | 4.84 (+0.74) | +0.12 (gradual) |
| Post-MoU (Feb 17+) | 4.34 | +0.55 |

Death counts jump +1.86 log units within the pre-period. The actual
MoU effect is only +0.32 — dwarfed by the pre-period break. Placebos
detect this break at any fake date before the volume jump → 89% FPR.

### Why C rate works

The rate shows no sharp break at MN (−0.55 → −0.63, flat). With 72
months of stationary pre-period, the model calibrates well. Post-MoU
rate shift (2.22% → 2.90%) is genuine relative to the trajectory.

---

## 5. Conditioning on Crossing Attempts — Discussion

### The idea

Including log(crossings) as a covariate absorbs the volume channel.
The model learns deaths = f(crossings, weather, ...) in the
pre-period, then checks whether post-MoU deaths are higher than
predicted at the *observed* (reduced) crossing level.

This resembles the **controlled direct effect** (CDE):
CDE(m) = E[Y(1,m)] − E[Y(0,m)], where CausalImpact estimates
E[Y(0,m)] from the pre-period relationship at observed post-MoU
crossings m.

### Benefits

1. Absorbs the volume channel
2. Fixes the C deaths placebo problem (pre-period structural break
   explained by crossings → pre-period becomes stationary)
3. Avoids Kronmal's (1993) ratio problem

### Problems (from literature review)

1. **Unmeasured mediator–outcome confounders** (VanderWeele 2015):
   CDE requires no unmeasured common causes of crossings and deaths.
   Likely violated: smuggler behavior, seasonal composition, Libyan
   conflict dynamics jointly affect volume and per-crossing danger.
2. **CDE requires a fixed mediator value** (Pearl 2001; Robins &
   Greenland 1992): CausalImpact conditions on time-varying observed
   crossings — yielding a weighted average of CDEs.
3. **CausalImpact not designed for this**: Brodersen et al. (2015)
   require covariates unaffected by the intervention.

### Honest framing

What we can say: conditioning on crossings absorbs the volume channel
and tests for a residual per-crossing danger change. Suggestive but
not formally identified as CDE.

**Added to presentation:** 4 slides covering the idea, mechanics,
problems, and honest framing.

---

## 6. Volume Channel Evidence (Script: `R/07_volume_channel_evidence.R`)

Empirical tests of two claims previously stated as assertions.

### Claim 1: Moon predicts crossing volume — NOT SUPPORTED

- cor(moon_lag1, log_crossings) = −0.124, p = 0.163 — not significant
- Pre-MoU subsample correlation is stronger (−0.313) but likely
  reflects seasonal confounding
- After controlling for crossings, moon non-significant (p = 0.34)
- High BSTS inclusion (86–99%) likely captures seasonal/volume
  patterns jointly with trend, not a clean channel

**Decision:** remove moon illumination from the model entirely.

### Claim 2: Death structural break driven by crossings — SUPPORTED

| | Unconditional | Conditional on crossings |
|---|---|---|
| NGO SAR dummy | coef = +1.87, p = 0.0008 | coef = −0.26, p = 0.64 |
| log(crossings) | — | coef = +1.19, p < 0.0001 |
| R² | 0.15 | 0.46 |

Structural break **disappears entirely** when conditioning on
crossings. Rate has no structural break even unconditionally
(F-test p = 0.39).

---

## 7. Moon Illumination Removed from Model

**File:** `R/02_data_prepare.R` — removed `"moon_illumination_frac"`
from `select_predictors()`. Reduces predictors from 40 to 38
exogenous (+ 11 month dummies = 49 total).

### Pipeline re-run results (completed in 3 min 13 sec)

| Model | Before (with moon) | After (no moon) |
|---|---|---|
| C rate | p = 0.0082 | p = 0.0079 |
| C deaths | p = 0.1759 | p = 0.1965 |
| A deaths | p = 0.0004 | p = 0.0003 |
| B deaths | p = 0.0001 | p = 0.0009 |

Removing moon changed nothing substantively — confirms it wasn't
doing real work.

---

## 8. Covariate Correlation Analysis (Script: `R/08_covariate_correlations.R`)

### Finding 1: Wrong-sign correlations with deaths

All sea-danger variables correlate **negatively** with deaths:

| Variable | r(deaths) | p | Expected | Actual |
|---|---:|---:|---|---|
| Wave SD (lag 1) | −0.278 | 0.002 | + | **−** |
| Wave max (lag 1) | −0.277 | 0.002 | + | **−** |
| Current speed (lag 1) | −0.266 | 0.003 | + | **−** |
| Wave days >2m (lag 1) | −0.261 | 0.003 | + | **−** |
| Temperature (departure) | +0.237 | 0.007 | ? | + |

Bad weather → fewer crossings → fewer deaths. The volume channel
dominates.

### Finding 2: Same variables are 2–3× stronger crossing predictors

| Variable | r(crossings) | r(deaths) | Ratio |
|---|---:|---:|---:|
| Wave height | −0.519*** | −0.227* | 2.3× |
| Wave days >2m | −0.498*** | −0.215* | 2.3× |
| Wave max | −0.488*** | −0.210* | 2.3× |
| Wave SD | −0.478*** | −0.224* | 2.1× |
| Temperature | +0.447*** | +0.237** | 1.9× |

### Finding 3: Near-zero signal for mortality rate

| Outcome | # sig. covariates (of 38) | Max |r| |
|---|---:|---:|
| Crossings | 30 | 0.519 |
| Deaths (and missing) | 22 | 0.278 |
| Mortality rate | 1 | 0.198 |

Covariates predict *when people cross*, not *how dangerous each
crossing is*.

---

## 9. Why No Model Produces a Defensible Causal Estimate

### Death counts: covariates predict volume, not danger

The correlation analysis proves that all sea-danger covariates
correlate with deaths in the wrong direction (rougher seas → fewer
deaths). The signal flows entirely through crossing volume. Death
count models (A, B) detect the MN/NGO SAR volume surge — a change
in how many people crossed, not in how dangerous each crossing was.
Model C deaths sees two effects cancel (higher per-crossing risk ×
fewer crossings ≈ flat deaths) and fails placebos.

### Mortality rate: three compounding problems

1. **Kronmal (1993) ratio bias.** Mortality rate = deaths/crossings.
   We have shown that 30/38 covariates significantly predict the
   denominator (crossings). When covariates predict the denominator
   of a ratio outcome, the regression coefficients are biased by
   spurious arithmetic correlation — the covariates appear related
   to the rate through the denominator, not through any genuine
   relationship with per-crossing danger. Even the single covariate
   significant with the rate (dewpoint depression, r = −0.198)
   may be a Kronmal artifact.

2. **Monthly aggregation destroys the danger signal.** The
   correlation asymmetry is the empirical proof: monthly weather
   averages correlate strongly with crossings (r up to 0.52)
   because the decision to cross integrates over the month — if
   weather is bad most of the month, fewer people depart. But
   mortality depends on conditions on the *specific departure day*:
   the one storm, the one overcrowded boat, whether a distress call
   was answered that night. A month with 29 calm days and one deadly
   storm looks average in the monthly mean but produces a mass
   casualty event. Monthly aggregation preserves the volume signal
   but destroys the danger signal — exactly what the 30/38 vs 1/38
   asymmetry shows.

3. **No covariate signal → trend-driven counterfactual.** With
   covariates contributing nothing (−2% marginal in the forecasting
   diagnostic, 1/38 significant with the rate), CausalImpact's
   counterfactual is identified by the local linear trend alone.
   The p-value (0.008 for Model C rate) tells us a structural break
   occurred, not what caused it. The result is sensitive to
   pre-period length and trend specification rather than anchored
   to observable exogenous variation.

### The fundamental tension

CausalImpact requires covariates that are both **exogenous** (not
affected by the intervention) and **predictive** (correlated with
the outcome). For mortality at monthly frequency:

- The variables that are **exogenous** (weather, sea state) predict
  **crossing volume**, not per-crossing danger
- The variables that would **predict** per-crossing danger (SAR
  capacity, vessel type, smuggler behavior, rescue response time)
  are **endogenous** to the policy
- No covariates satisfying both requirements exist at this
  aggregation level

None of the six models — across two outcomes and three interventions
— produce a defensible causal estimate.

---

## 10. Presentation Updates — Complete List

### Slides added (7 new)

1. **"Extension: Curated Specification (809 Predictors)"** — Extension
   1 results with 6 models
2. **"Curated Specification — Implications"** — why volume-channel
   predictors bias results
3. **"The Reversed Pattern: Volume Channel"** — period-by-period
   crossings/deaths/rate table
4. **"The Reversed Pattern: Mortality Rate"** — why rate works for C
5. **"Why C Deaths Fails: Pre-Period Non-Stationarity"** — structural
   break table
6. **"Why C Deaths Fails: Placebo Mechanics"** (+ continuation slide)
   — step-by-step placebo explanation with concrete example
7. **"A Possible Fix: Conditioning on Crossing Attempts"** —
   mechanics, CDE interpretation, problems, honest framing (4 slides)
8. **"Evidence: Structural Break Driven by Crossings"** — regression
   table proving break is volume-driven
9. **"Why Covariates Fail: Wrong Sign"** — wrong-sign correlations
10. **"Why Covariates Fail: They Predict Crossings"** — crossing vs
    death correlations
11. **"The Ratio Problem (Kronmal 1993)"** — Kronmal critique

### Slides modified

- Updated results table to 6 models
- Updated placebo table to 6 models
- Updated death count figures from 1-panel to 3-panel
- Updated covariate table (removed moon, updated counts)
- Removed moon references throughout
- Updated truncation test values
- Updated outline

---

## 10. Complete File Inventory

### Files created

| File | Description |
|---|---|
| `R/07_volume_channel_evidence.R` | Volume channel evidence script |
| `R/08_covariate_correlations.R` | Covariate correlation analysis |
| `output/session_notes/2026-02-26_*.md` | This document |
| `Extension-1/results/models/model_curated_deaths_A.rds` | Ext. 1 Model A deaths |
| `Extension-1/results/models/model_curated_deaths_B.rds` | Ext. 1 Model B deaths |

### Files modified

| File | Changes |
|---|---|
| `R/02_data_prepare.R` | Added A/B deaths interventions; removed moon |
| `_targets.R` | Added 4 targets (A/B death models + placebos); updated figure/snapshot calls |
| `R/05_figures.R` | Expanded all figures to 6-model layout |
| `data/raw/df.RDS` | Fixed broken symlink |
| `Presentation 2/progress_presentation.qmd` | ~20 edits, 10+ new slides |

### Pipeline runs

1. First `tar_make()` — failed (broken symlink), fixed, re-run
2. Second `tar_make()` — 22 targets, 3 min 14 sec (with moon)
3. Third `tar_make()` — 22 targets, 3 min 13 sec (without moon)
