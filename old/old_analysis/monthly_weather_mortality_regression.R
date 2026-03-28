# monthly_weather_mortality_regression.R
# ======================================
# Monthly death rate regression with SWH + wind speed as the two
# weather variables, interacted with post-MoU indicator.
#
# Based on the PCA exploration: at monthly frequency with spatial means,
# sea danger is essentially two-dimensional (SWH + wind). Wave period
# is redundant with SWH (r = 0.98), wave extremes are mechanically
# derived from SWH. PCA deferred to event-level analysis.

library(dplyr)
library(lubridate)

BASE_DIR <- file.path("replication", "rodriguez-sanchez", "extension")
PROJECT_DIR <- here::here()

# ============================================================
# 1. Load and prepare data
# ============================================================
df <- readRDS(file.path(BASE_DIR, "data", "df_extended.RDS"))

PRE_START <- as.Date("2011-02-01")
END_DATE  <- as.Date("2021-09-01")
MOU_DATE  <- as.Date("2017-07-01")

df <- df %>%
  filter(date >= PRE_START & date <= END_DATE) %>%
  mutate(
    log_rate     = log(mortality_rate_100 + 0.01),
    log_rate_1   = log(mortality_rate_100 + 1),
    log_deaths   = log(dead_and_missing_Central_Mediterranean + 1),
    log_cross    = log(crossings_CMR + 1),
    post_mou     = as.integer(date >= MOU_DATE),
    month_fac    = factor(month(date)),
    # Standardize weather for comparable coefficients
    swh_z        = scale(wave_height_central_med)[, 1],
    wind_z       = scale(wind_speed_central_med)[, 1],
    year_fac     = factor(year(date))
  )

# ============================================================
# 1b. Merge corridor/core weather from daily ERA5 panel
# ============================================================
daily_panel <- readRDS(file.path(PROJECT_DIR, "data", "processed",
                                  "cmr_daily_weather_panel.RDS"))

monthly_wx <- daily_panel %>%
  mutate(ym_date = lubridate::floor_date(date, "month")) %>%
  group_by(ym_date) %>%
  summarise(
    swh_corridor  = mean(swh_mean, na.rm = TRUE),
    swh_core      = mean(swh_core, na.rm = TRUE),
    wind_corridor = mean(wind_mean, na.rm = TRUE),
    wind_core     = mean(wind_core, na.rm = TRUE),
    .groups = "drop"
  )

df <- df %>% left_join(monthly_wx, by = c("date" = "ym_date"))

cat("Analysis sample:", nrow(df), "months\n")
cat("  Pre-MoU:", sum(df$post_mou == 0), "months\n")
cat("  Post-MoU:", sum(df$post_mou == 1), "months\n\n")

# ============================================================
# 2. Newey-West HAC standard errors
# ============================================================
newey_west_vcov <- function(model, max_lag = 4) {
  X <- model.matrix(model)
  e <- residuals(model)
  n <- length(e)
  k <- ncol(X)
  meat <- matrix(0, k, k)
  for (j in 0:max_lag) {
    w <- if (j == 0) 1 else 1 - j / (max_lag + 1)
    for (t in (j + 1):n) {
      if (j == 0) {
        meat <- meat + e[t]^2 * tcrossprod(X[t, ])
      } else {
        meat <- meat + w * e[t] * e[t - j] *
          (tcrossprod(X[t, ], X[t - j, ]) +
             tcrossprod(X[t - j, ], X[t, ]))
      }
    }
  }
  bread <- solve(crossprod(X))
  n / (n - k) * bread %*% meat %*% bread
}

hac_summary <- function(model, coef_patterns, max_lag = 4) {
  vcov_hac <- newey_west_vcov(model, max_lag)
  hac_se <- sqrt(diag(vcov_hac))
  betas <- coef(model)
  all_names <- names(betas)

  # Match patterns to actual coefficient names
  coef_names <- character(0)
  for (pat in coef_patterns) {
    matched <- all_names[all_names == pat]
    if (length(matched) == 0) {
      # Try matching with reversed interaction order
      parts <- strsplit(pat, ":")[[1]]
      if (length(parts) == 2) {
        alt <- paste(parts[2], parts[1], sep = ":")
        matched <- all_names[all_names == alt]
      }
    }
    if (length(matched) > 0) coef_names <- c(coef_names, matched[1])
  }

  results <- data.frame(
    beta = betas[coef_names],
    ols_se = summary(model)$coefficients[coef_names, 2],
    ols_p  = summary(model)$coefficients[coef_names, 4],
    hac_se = hac_se[coef_names],
    hac_t  = betas[coef_names] / hac_se[coef_names],
    hac_p  = 2 * pt(-abs(betas[coef_names] / hac_se[coef_names]),
                     df = model$df.residual),
    row.names = coef_names
  )
  results
}

# ============================================================
# 3. Main specification: log(rate) ~ SWH * post + Wind * post + month FE
# ============================================================
cat("============================================================\n")
cat("MODEL 1: JOINT SWH + WIND INTERACTION (raw scale)\n")
cat("log(rate + 0.01) ~ SWH * PostMoU + Wind * PostMoU + month FE\n")
cat("============================================================\n\n")

m1 <- lm(log_rate ~ wave_height_central_med * post_mou +
           wind_speed_central_med * post_mou + month_fac,
         data = df)

coefs_of_interest <- c(
  "wave_height_central_med",
  "wind_speed_central_med",
  "post_mou",
  "wave_height_central_med:post_mou",
  "wind_speed_central_med:post_mou"
)

res1 <- hac_summary(m1, coefs_of_interest)
cat(sprintf("%-38s %10s %8s %8s %8s\n",
    "Coefficient", "Beta", "OLS p", "HAC SE", "HAC p"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in seq_len(nrow(res1))) {
  cat(sprintf("%-38s %+10.4f %8.4f %8.4f %8.4f%s\n",
      rownames(res1)[i], res1$beta[i], res1$ols_p[i],
      res1$hac_se[i], res1$hac_p[i],
      ifelse(res1$hac_p[i] < 0.01, " ***",
             ifelse(res1$hac_p[i] < 0.05, " **",
                    ifelse(res1$hac_p[i] < 0.1, " *", "")))))
}
cat(sprintf("\nR-squared: %.4f, Adj R-squared: %.4f, df.residual: %d\n",
    summary(m1)$r.squared, summary(m1)$adj.r.squared, m1$df.residual))

# Joint F-test: both interactions = 0
m1_no_int <- lm(log_rate ~ wave_height_central_med + wind_speed_central_med +
                  post_mou + month_fac, data = df)
f1 <- anova(m1_no_int, m1)
cat(sprintf("Joint F-test (both interactions = 0): F = %.3f, p = %.4f%s\n",
    f1$F[2], f1$`Pr(>F)`[2],
    ifelse(f1$`Pr(>F)`[2] < 0.05, " **", "")))

# ============================================================
# 4. Standardized version (for comparing effect sizes)
# ============================================================
cat("\n============================================================\n")
cat("MODEL 2: STANDARDIZED (z-scored weather)\n")
cat("============================================================\n\n")

m2 <- lm(log_rate ~ swh_z * post_mou + wind_z * post_mou + month_fac,
         data = df)

coefs_z <- c("swh_z", "wind_z", "post_mou",
             "swh_z:post_mou", "wind_z:post_mou")
res2 <- hac_summary(m2, coefs_z)

cat(sprintf("%-28s %10s %8s %8s\n",
    "Coefficient", "Beta", "HAC SE", "HAC p"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in seq_len(nrow(res2))) {
  cat(sprintf("%-28s %+10.4f %8.4f %8.4f%s\n",
      rownames(res2)[i], res2$beta[i], res2$hac_se[i], res2$hac_p[i],
      ifelse(res2$hac_p[i] < 0.01, " ***",
             ifelse(res2$hac_p[i] < 0.05, " **",
                    ifelse(res2$hac_p[i] < 0.1, " *", "")))))
}

# ============================================================
# 5. Robustness: alternative log(rate + 1) constant
# ============================================================
cat("\n============================================================\n")
cat("MODEL 3: ROBUSTNESS -- log(rate + 1) instead of log(rate + 0.01)\n")
cat("============================================================\n\n")

m3 <- lm(log_rate_1 ~ wave_height_central_med * post_mou +
           wind_speed_central_med * post_mou + month_fac,
         data = df)

res3 <- hac_summary(m3, coefs_of_interest)
cat(sprintf("%-38s %10s %8s\n", "Coefficient", "Beta", "HAC p"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in seq_len(nrow(res3))) {
  cat(sprintf("%-38s %+10.4f %8.4f%s\n",
      rownames(res3)[i], res3$beta[i], res3$hac_p[i],
      ifelse(res3$hac_p[i] < 0.05, " **",
             ifelse(res3$hac_p[i] < 0.1, " *", ""))))
}

# ============================================================
# 6. Robustness: quarter FE instead of month FE
# ============================================================
cat("\n============================================================\n")
cat("MODEL 4: ROBUSTNESS -- quarter FE (conserve df)\n")
cat("============================================================\n\n")

df$quarter_fac <- factor(quarter(df$date))
m4 <- lm(log_rate ~ wave_height_central_med * post_mou +
           wind_speed_central_med * post_mou + quarter_fac,
         data = df)

res4 <- hac_summary(m4, coefs_of_interest)
cat(sprintf("%-38s %10s %8s\n", "Coefficient", "Beta", "HAC p"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in seq_len(nrow(res4))) {
  cat(sprintf("%-38s %+10.4f %8.4f%s\n",
      rownames(res4)[i], res4$beta[i], res4$hac_p[i],
      ifelse(res4$hac_p[i] < 0.05, " **",
             ifelse(res4$hac_p[i] < 0.1, " *", ""))))
}
cat(sprintf("df.residual: %d (vs %d with month FE)\n",
    m4$df.residual, m1$df.residual))

# ============================================================
# 6b. CORRIDOR / CORE WEATHER (ERA5 subsample, 2014+)
# ============================================================
cat("\n============================================================\n")
cat("CORRIDOR / CORE WEATHER RATE REGRESSIONS (ERA5 subsample)\n")
cat("============================================================\n\n")

df_era5 <- df %>% filter(!is.na(swh_corridor))
cat(sprintf("ERA5 subsample: %d months\n\n", nrow(df_era5)))

# Core SWH + Wind joint (month FE)
m_core <- lm(log_rate ~ swh_core * post_mou + wind_core * post_mou + month_fac,
              data = df_era5)
coefs_core <- c("swh_core", "wind_core", "post_mou",
                "swh_core:post_mou", "wind_core:post_mou")
res_core <- hac_summary(m_core, coefs_core)
cat("--- Core SWH + Wind x Post (month FE) ---\n")
cat(sprintf("%-30s %10s %8s %8s\n", "Coefficient", "Beta", "HAC SE", "HAC p"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in seq_len(nrow(res_core))) {
  cat(sprintf("%-30s %+10.4f %8.4f %8.4f%s\n",
      rownames(res_core)[i], res_core$beta[i], res_core$hac_se[i], res_core$hac_p[i],
      ifelse(res_core$hac_p[i] < 0.05, " **", ifelse(res_core$hac_p[i] < 0.1, " *", ""))))
}

# Core SWH + Wind joint (year + month FE)
m_core_yr <- lm(log_rate ~ swh_core * post_mou + wind_core * post_mou +
                   year_fac + month_fac, data = df_era5)
coefs_core_yr <- c("swh_core:post_mou", "wind_core:post_mou")
res_core_yr <- hac_summary(m_core_yr, coefs_core_yr)
cat("\n--- Core SWH + Wind x Post (year + month FE) ---\n")
for (i in seq_len(nrow(res_core_yr))) {
  cat(sprintf("  %-28s %+10.4f  HAC p=%8.4f%s\n",
      rownames(res_core_yr)[i], res_core_yr$beta[i], res_core_yr$hac_p[i],
      ifelse(res_core_yr$hac_p[i] < 0.05, " **", ifelse(res_core_yr$hac_p[i] < 0.1, " *", ""))))
}

# Broad SWH on ERA5 subsample (for comparison)
m_broad_era5 <- lm(log_rate ~ wave_height_central_med * post_mou +
                      wind_speed_central_med * post_mou + month_fac,
                    data = df_era5)
res_broad_era5 <- hac_summary(m_broad_era5, coefs_of_interest)

# Geography comparison table
# Use grep to find interaction rows (R may reverse interaction ordering)
find_row <- function(res, pattern) {
  idx <- grep(pattern, rownames(res))
  if (length(idx) > 0) idx[1] else NA
}

cat("\n--- Geography comparison (same ERA5 sample, month FE) ---\n")
cat(sprintf("  %-12s  %10s  %8s  |  %10s  %8s\n",
    "Geography", "SWH×Post", "HAC p", "Wind×Post", "HAC p"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (geo_spec in list(
  list("Broad", res_broad_era5, "wave_height.*post|post.*wave_height",
       "wind_speed.*post|post.*wind_speed"),
  list("Core", res_core, "swh_core.*post|post.*swh_core",
       "wind_core.*post|post.*wind_core")
)) {
  r <- geo_spec[[2]]
  swh_i <- find_row(r, geo_spec[[3]])
  wnd_i <- find_row(r, geo_spec[[4]])
  cat(sprintf("  %-12s  %+10.4f  %8.4f  |  %+10.4f  %8.4f\n",
      geo_spec[[1]],
      ifelse(!is.na(swh_i), r$beta[swh_i], NA),
      ifelse(!is.na(swh_i), r$hac_p[swh_i], NA),
      ifelse(!is.na(wnd_i), r$beta[wnd_i], NA),
      ifelse(!is.na(wnd_i), r$hac_p[wnd_i], NA)))
}

# ============================================================
# 7. Individual variable models (for comparison)
# ============================================================
cat("\n============================================================\n")
cat("COMPARISON: INDIVIDUAL vs JOINT models\n")
cat("============================================================\n\n")

# SWH only
m_swh <- lm(log_rate ~ wave_height_central_med * post_mou + month_fac,
            data = df)
s_swh <- hac_summary(m_swh, c("wave_height_central_med:post_mou"))

# Wind only
m_wind <- lm(log_rate ~ wind_speed_central_med * post_mou + month_fac,
             data = df)
s_wind <- hac_summary(m_wind, c("wind_speed_central_med:post_mou"))

# Joint (from m1) — use grep since R may reverse interaction ordering
swh_j <- find_row(res1, "wave_height.*post|post.*wave_height")
wnd_j <- find_row(res1, "wind_speed.*post|post.*wind_speed")

cat(sprintf("%-20s | %12s %8s | %12s %8s\n",
    "", "Individual", "HAC p", "Joint", "HAC p"))
cat(paste(rep("-", 70), collapse = ""), "\n")
cat(sprintf("%-20s | %+12.4f %8.4f | %+12.4f %8.4f\n",
    "SWH x PostMoU",
    s_swh$beta, s_swh$hac_p,
    res1$beta[swh_j], res1$hac_p[swh_j]))
cat(sprintf("%-20s | %+12.4f %8.4f | %+12.4f %8.4f\n",
    "Wind x PostMoU",
    s_wind$beta, s_wind$hac_p,
    res1$beta[wnd_j], res1$hac_p[wnd_j]))

# ============================================================
# 8. Deaths model (with volume control) -- note mediator caveat
# ============================================================
cat("\n============================================================\n")
cat("MODEL 5: DEATHS (with volume control)\n")
cat("log(deaths+1) ~ SWH*post + Wind*post + log(cross) + month FE\n")
cat("NOTE: log(crossings) is a post-treatment mediator. Interpret with caution.\n")
cat("============================================================\n\n")

m5 <- lm(log_deaths ~ wave_height_central_med * post_mou +
           wind_speed_central_med * post_mou + log_cross + month_fac,
         data = df)

coefs_d <- c("wave_height_central_med", "wind_speed_central_med",
             "log_cross", "post_mou",
             "wave_height_central_med:post_mou",
             "wind_speed_central_med:post_mou")
res5 <- hac_summary(m5, coefs_d)

cat(sprintf("%-38s %10s %8s %8s\n",
    "Coefficient", "Beta", "HAC SE", "HAC p"))
cat(paste(rep("-", 70), collapse = ""), "\n")
for (i in seq_len(nrow(res5))) {
  cat(sprintf("%-38s %+10.4f %8.4f %8.4f%s\n",
      rownames(res5)[i], res5$beta[i], res5$hac_se[i], res5$hac_p[i],
      ifelse(res5$hac_p[i] < 0.05, " **",
             ifelse(res5$hac_p[i] < 0.1, " *", ""))))
}

# ============================================================
# 9. Autocorrelation diagnostic
# ============================================================
cat("\n============================================================\n")
cat("AUTOCORRELATION DIAGNOSTIC\n")
cat("============================================================\n\n")

dw1 <- sum(diff(residuals(m1))^2) / sum(residuals(m1)^2)
cat(sprintf("Durbin-Watson (joint rate model): %.3f %s\n",
    dw1, ifelse(dw1 < 1.5, "(positive autocorrelation)", "(OK)")))

# Residual ACF at lags 1-3
resid1 <- residuals(m1)
n <- length(resid1)
for (lag in 1:3) {
  r <- cor(resid1[(lag + 1):n], resid1[1:(n - lag)])
  cat(sprintf("  ACF lag %d: %.3f\n", lag, r))
}

# ============================================================
# 10. Summary table for export
# ============================================================
cat("\n============================================================\n")
cat("SUMMARY TABLE\n")
cat("============================================================\n\n")

cat(sprintf("%-25s | %10s %8s | %10s %8s | %10s %8s\n",
    "", "Rate", "p", "Rate(c=1)", "p", "Rate(qtr)", "p"))
cat(paste(rep("-", 90), collapse = ""), "\n")

# SWH x PostMoU (use grep to handle reversed interaction names)
swh_pat <- "wave_height.*post|post.*wave_height"
wnd_pat <- "wind_speed.*post|post.*wind_speed"

get_val <- function(res, pattern, col) {
  idx <- find_row(res, pattern)
  if (!is.na(idx)) res[[col]][idx] else NA
}

cat(sprintf("%-25s | %+10.4f %8.4f | %+10.4f %8.4f | %+10.4f %8.4f\n",
    "SWH x PostMoU",
    get_val(res1, swh_pat, "beta"), get_val(res1, swh_pat, "hac_p"),
    get_val(res3, swh_pat, "beta"), get_val(res3, swh_pat, "hac_p"),
    get_val(res4, swh_pat, "beta"), get_val(res4, swh_pat, "hac_p")))

# Wind x PostMoU
cat(sprintf("%-25s | %+10.4f %8.4f | %+10.4f %8.4f | %+10.4f %8.4f\n",
    "Wind x PostMoU",
    get_val(res1, wnd_pat, "beta"), get_val(res1, wnd_pat, "hac_p"),
    get_val(res3, wnd_pat, "beta"), get_val(res3, wnd_pat, "hac_p"),
    get_val(res4, wnd_pat, "beta"), get_val(res4, wnd_pat, "hac_p")))

cat(sprintf("\n%-25s | %10.4f %8s | %10.4f %8s | %10.4f\n",
    "R-squared",
    summary(m1)$r.squared, "",
    summary(m3)$r.squared, "",
    summary(m4)$r.squared))
cat(sprintf("%-25s | %10d %8s | %10d %8s | %10d\n",
    "df.residual",
    m1$df.residual, "", m3$df.residual, "", m4$df.residual))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
