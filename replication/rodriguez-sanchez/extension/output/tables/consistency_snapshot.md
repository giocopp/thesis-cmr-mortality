# Extension-2 Consistency Snapshot

Generated: 2026-02-26 16:17:38

Pipeline: `targets::tar_make()`

## Model Results

Note: AbsEffect and CI are on the model's transformed outcome scale (log units).

| Model | Treatment | Outcome | p-value | AbsEffect (log units) | 95% CI lower | 95% CI upper | RMSE/SD (pre) | Burn used (draws) |
|---|---|---|---:|---:|---:|---:|---:|---:|
| A_mortality | Mare Nostrum | mortality_rate | 0.0881 | 53.0243 | -22.3584 | 133.1337 | 0.7514 | 123 |
| B_mortality | NGO SAR | mortality_rate | 0.0046 | 81.6618 | 21.0311 | 140.3617 | 0.9222 | 1 |
| C_mortality | EU-Libya MoU | mortality_rate | 0.0079 | 45.2984 | 9.0555 | 82.2641 | 0.9769 | 27 |
| A_deaths | Mare Nostrum | death_count | 3e-04 | 109.6973 | 45.0372 | 178.5749 | 0.7246 | 45 |
| B_deaths | NGO SAR | death_count | 9e-04 | 96.0513 | 39.3731 | 152.8314 | 0.8739 | 3 |
| C_deaths | EU-Libya MoU | death_count | 0.1965 | 16.6144 | -21.4419 | 55.1192 | 0.9244 | 1 |

## Placebo Diagnostics

| Model | # Placebos | # Significant | FPR | Actual p | Placebos p<=actual | Share p<=actual |
|---|---:|---:|---:|---:|---:|---:|
| A_mortality | 7 | 7 | 1.000 | 0.0881 | 7/7 | 1.000 |
| B_mortality | 9 | 6 | 0.667 | 0.0046 | 5/9 | 0.556 |
| C_mortality | 18 | 8 | 0.444 | 0.0079 | 3/18 | 0.167 |
| A_deaths | 7 | 7 | 1.000 | 3e-04 | 2/7 | 0.286 |
| B_deaths | 9 | 5 | 0.556 | 9e-04 | 1/9 | 0.111 |
| C_deaths | 18 | 16 | 0.889 | 0.1965 | 17/18 | 0.944 |

## Model A Truncation

| Label | Post months | p-value | Significant |
|---|---:|---:|---|
| Mare Nostrum window only (13 mo) | 13 | 0.3888 | NO |
| +1yr after MN (25 mo) | 25 | 0.3610 | NO |
| Pre-MoU window (37 mo) | 37 | 0.3654 | NO |
| Full post-period (96 mo) | 96 | 0.0877 | NO |

## Forecast Diagnostic (weighted averages)

| Model | wt.RMSE | wt.MAE | RMSE/SD | vs M0a | vs M0c |
|---|---:|---:|---:|---:|---:|
| M0a (month mean) | 1.870 | 1.581 | 0.893 | +0.0% | -0.5% |
| M0c (state-only) | 1.861 | 1.584 | 0.888 | +0.5% | +0.0% |
| M1-A1 (LL+dummies) | 1.595 | 1.379 | 0.761 | +14.7% | +14.3% |
| M1-A2 (LLT+dummies) | 1.858 | 1.576 | 0.887 | +0.6% | +0.2% |
| M1-B (LL+Seasonal) | 1.863 | 1.588 | 0.889 | +0.4% | -0.1% |
