# run_weather_mortality_regressions.R
# ====================================
# Quick diagnostic: does weather predict mortality after controlling
# for seasonality and volume? Tests whether bivariate r ≈ 0 is because
# (a) weather truly doesn't matter, or (b) confounders hide the signal.

library(dplyr)
library(lubridate)

# --- Load data from targets ---
message("Loading data from targets store...")
df_model <- targets::tar_read(df_model, store = "_targets")
df_full  <- df_model$df_full

# --- Build analysis dataset ---
PRE_START <- as.Date("2011-02-01")
END_DATE  <- as.Date("2021-09-01")
MOU_DATE  <- as.Date("2017-02-01")

df <- df_full %>%
  filter(date >= PRE_START & date <= END_DATE) %>%
  mutate(
    log_rate     = log(mortality_rate_100 + 0.01),
    log_deaths   = log(dead_and_missing_Central_Mediterranean + 1),
    log_cross    = log(crossings_CMR + 1),
    post_mou     = as.integer(date >= MOU_DATE),
    month_factor = factor(month(date))
  )

# Weather variables (use key ones, not all 19)
weather_vars <- c(
  "wave_height_central_med",
  "wind_speed_central_med",
  "sst_central_med",
  "current_speed_central_med",
  "wave_max_central_med",
  "wave_days_above_2m"
)

# Check availability
available <- weather_vars[weather_vars %in% names(df)]
message("Weather variables available: ", length(available), "/", length(weather_vars))

# =========================================================================
# TEST 1: Bivariate (confirmation of weak correlations)
# =========================================================================
message("\n========================================")
message("TEST 1: Bivariate correlations (no controls)")
message("========================================\n")

for (v in available) {
  r_rate   <- cor(df[[v]], df$log_rate, use = "complete.obs")
  r_deaths <- cor(df[[v]], df$log_deaths, use = "complete.obs")
  message(sprintf("  %-30s  r(rate)=%+.3f  r(deaths)=%+.3f", v, r_rate, r_deaths))
}

# =========================================================================
# TEST 2: Partial correlation controlling for month dummies
# =========================================================================
message("\n========================================")
message("TEST 2: Weather → mortality, controlling for SEASONALITY")
message("========================================")
message("(month dummies absorb the seasonal confound)\n")

for (v in available) {
  fml_rate   <- as.formula(paste0("log_rate ~ ", v, " + month_factor"))
  fml_deaths <- as.formula(paste0("log_deaths ~ ", v, " + month_factor"))

  m_rate   <- lm(fml_rate, data = df)
  m_deaths <- lm(fml_deaths, data = df)

  coef_rate   <- summary(m_rate)$coefficients[v, ]
  coef_deaths <- summary(m_deaths)$coefficients[v, ]

  message(sprintf("  %-30s", v))
  message(sprintf("    Rate:   β=%+.4f  SE=%.4f  t=%.2f  p=%.4f %s",
    coef_rate[1], coef_rate[2], coef_rate[3], coef_rate[4],
    ifelse(coef_rate[4] < 0.05, " *", "")))
  message(sprintf("    Deaths: β=%+.4f  SE=%.4f  t=%.2f  p=%.4f %s",
    coef_deaths[1], coef_deaths[2], coef_deaths[3], coef_deaths[4],
    ifelse(coef_deaths[4] < 0.05, " *", "")))
}

# =========================================================================
# TEST 3: Weather → deaths, controlling for seasonality + volume
# =========================================================================
message("\n========================================")
message("TEST 3: Weather → deaths, controlling for SEASONALITY + VOLUME")
message("========================================")
message("(month dummies + log(crossings) — isolates per-crossing danger)\n")

for (v in available) {
  fml <- as.formula(paste0("log_deaths ~ ", v, " + log_cross + month_factor"))
  m   <- lm(fml, data = df)
  coef_v <- summary(m)$coefficients[v, ]

  message(sprintf("  %-30s  β=%+.4f  SE=%.4f  t=%.2f  p=%.4f %s",
    v, coef_v[1], coef_v[2], coef_v[3], coef_v[4],
    ifelse(coef_v[4] < 0.05, " *", "")))
}

# =========================================================================
# TEST 4: Weather × post-MoU interaction (the key test)
# =========================================================================
message("\n========================================")
message("TEST 4: Weather × Post-MoU INTERACTION")
message("========================================")
message("(does weather's effect on mortality CHANGE after the MoU?)\n")

for (v in available) {
  # Rate model
  fml_rate <- as.formula(paste0("log_rate ~ ", v, " * post_mou + month_factor"))
  m_rate   <- lm(fml_rate, data = df)
  s_rate   <- summary(m_rate)$coefficients
  interact_name <- paste0(v, ":post_mou")

  if (interact_name %in% rownames(s_rate)) {
    ci <- s_rate[interact_name, ]
    main <- s_rate[v, ]
    message(sprintf("  %-30s", v))
    message(sprintf("    Rate main:     β=%+.4f  p=%.4f", main[1], main[4]))
    message(sprintf("    Rate interact: β=%+.4f  p=%.4f %s", ci[1], ci[4],
      ifelse(ci[4] < 0.05, " *", "")))
  }

  # Deaths model (with volume control)
  fml_deaths <- as.formula(paste0("log_deaths ~ ", v, " * post_mou + log_cross + month_factor"))
  m_deaths   <- lm(fml_deaths, data = df)
  s_deaths   <- summary(m_deaths)$coefficients

  if (interact_name %in% rownames(s_deaths)) {
    ci <- s_deaths[interact_name, ]
    main <- s_deaths[v, ]
    message(sprintf("    Deaths main:     β=%+.4f  p=%.4f", main[1], main[4]))
    message(sprintf("    Deaths interact: β=%+.4f  p=%.4f %s", ci[1], ci[4],
      ifelse(ci[4] < 0.05, " *", "")))
  }
  message("")
}

# =========================================================================
# TEST 5: Joint F-test — do ALL weather vars together explain mortality?
# =========================================================================
message("\n========================================")
message("TEST 5: Joint F-test — all weather vars together")
message("========================================\n")

weather_formula_rate <- as.formula(
  paste0("log_rate ~ ", paste(available, collapse = " + "), " + month_factor"))
null_formula_rate <- log_rate ~ month_factor

m_full_rate <- lm(weather_formula_rate, data = df)
m_null_rate <- lm(null_formula_rate, data = df)
f_rate <- anova(m_null_rate, m_full_rate)

message(sprintf("  Rate: null R²=%.4f, full R²=%.4f, ΔR²=%.4f",
  summary(m_null_rate)$r.squared, summary(m_full_rate)$r.squared,
  summary(m_full_rate)$r.squared - summary(m_null_rate)$r.squared))
message(sprintf("  F-test: F=%.3f, p=%.4f %s",
  f_rate$F[2], f_rate$`Pr(>F)`[2],
  ifelse(f_rate$`Pr(>F)`[2] < 0.05, " *", "")))

weather_formula_deaths <- as.formula(
  paste0("log_deaths ~ ", paste(available, collapse = " + "), " + log_cross + month_factor"))
null_formula_deaths <- log_deaths ~ log_cross + month_factor

m_full_deaths <- lm(weather_formula_deaths, data = df)
m_null_deaths <- lm(null_formula_deaths, data = df)
f_deaths <- anova(m_null_deaths, m_full_deaths)

message(sprintf("\n  Deaths (controlling for volume): null R²=%.4f, full R²=%.4f, ΔR²=%.4f",
  summary(m_null_deaths)$r.squared, summary(m_full_deaths)$r.squared,
  summary(m_full_deaths)$r.squared - summary(m_null_deaths)$r.squared))
message(sprintf("  F-test: F=%.3f, p=%.4f %s",
  f_deaths$F[2], f_deaths$`Pr(>F)`[2],
  ifelse(f_deaths$`Pr(>F)`[2] < 0.05, " *", "")))
