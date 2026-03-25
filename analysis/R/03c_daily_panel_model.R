# 03c_daily_panel_model.R
# ========================
# Deiana-style daily panel model: deaths ~ weather x post + week-by-year FE.
#
# Design logic (following Deiana, Maheshri & Mastrobuoni 2024):
# - Daily time series of death counts on the CMR
# - Weather = daily spatial mean over CMR sea area (from ERA5)
# - Week-by-year FE absorb ALL slow-moving confounders (seasonality, trends,
#   policy shifts WITHIN a week are negligible)
# - Identification from day-to-day weather variation WITHIN each week
# - The weather x post interaction captures whether the same weather shock
#   kills more people after the MoU than before
# - No gradient stability assumption across years needed
#
# Why this works when event-level week-year FE failed:
# - Event-level: 1,381 incidents / 628 cells = median 2 per cell → singular
# - Daily panel: 4,383 days / 628 cells = 7 per cell → identified
#
# Input:  data/processed/cmr_daily_weather_panel.RDS
# Output: output/tables/daily_panel_results.csv
#         output/figures/daily_panel_coefplot.pdf
#         printed diagnostics

library(fixest)
library(data.table)
library(ggplot2)

BASE_DIR <- here::here()
d <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                      "cmr_daily_weather_panel.RDS")))

MOU_DATE <- as.Date("2017-02-01")

cat("============================================================\n")
cat("DAILY PANEL MODEL: Deaths ~ Weather x Post + Week-Year FE\n")
cat("============================================================\n\n")

cat(sprintf("Panel: %d days (%s to %s)\n", nrow(d), min(d$date), max(d$date)))
cat(sprintf("Week-year cells: %d (mean %.1f days/cell)\n",
    uniqueN(d$week_year), nrow(d) / uniqueN(d$week_year)))
cat(sprintf("Outcome: n_dead_missing (mean=%.2f, max=%d, zeros=%d [%.1f%%])\n",
    mean(d$n_dead_missing), max(d$n_dead_missing),
    sum(d$n_dead_missing == 0), 100 * mean(d$n_dead_missing == 0)))
cat(sprintf("Variance/mean: %.1f (confirms overdispersion)\n\n",
    var(d$n_dead_missing) / mean(d$n_dead_missing)))

# ============================================================
# 1. Primary specification: Week-year FE
# ============================================================
cat("============================================================\n")
cat("1. PRIMARY SPEC: Week-Year FE (Deiana-style)\n")
cat("============================================================\n\n")

# Spec A: SWH x Post | week-year FE
mA <- fenegbin(n_dead_missing ~ swh_mean * post_mou | week_year_fac,
               data = d, vcov = "hetero")

# Spec B: Wind x Post | week-year FE
mB <- fenegbin(n_dead_missing ~ wind_mean * post_mou | week_year_fac,
               data = d, vcov = "hetero")

# Spec C: Gust x Post | week-year FE
mC <- fenegbin(n_dead_missing ~ gust_mean * post_mou | week_year_fac,
               data = d, vcov = "hetero")

# Spec D: SWH x Post + Wind | week-year FE
mD <- fenegbin(n_dead_missing ~ swh_mean * post_mou + wind_mean | week_year_fac,
               data = d, vcov = "hetero")

# Spec E: Wind x Post + SWH | week-year FE
mE <- fenegbin(n_dead_missing ~ wind_mean * post_mou + swh_mean | week_year_fac,
               data = d, vcov = "hetero")

cat("--- Results (hetero-robust SEs) ---\n\n")

report_model <- function(m, label) {
  ct <- summary(m, vcov = "hetero")$coeftable
  cat(sprintf("  %s:\n", label))
  # Show non-FE coefficients
  for (i in seq_len(nrow(ct))) {
    stars <- ifelse(ct[i, 4] < 0.01, "***",
                    ifelse(ct[i, 4] < 0.05, "**",
                           ifelse(ct[i, 4] < 0.1, "*", "")))
    cat(sprintf("    %-35s %+8.4f (SE=%6.4f) p=%6.4f %s  IRR=%.4f\n",
        rownames(ct)[i], ct[i, 1], ct[i, 2], ct[i, 4], stars, exp(ct[i, 1])))
  }
  cat(sprintf("    N=%d, LogLik=%.1f\n\n", nobs(m), logLik(m)))
}

report_model(mA, "Spec A: SWH x Post")
report_model(mB, "Spec B: Wind x Post")
report_model(mC, "Spec C: Gust x Post")
report_model(mD, "Spec D: SWH x Post + Wind")
report_model(mE, "Spec E: Wind x Post + SWH")


# ============================================================
# 1b. LAG STRUCTURE: Prior-day weather (Camarena-style)
# ============================================================
cat("============================================================\n")
cat("1b. LAG STRUCTURE: Prior-day weather (corridor)\n")
cat("============================================================\n\n")

# SWH lag1 (yesterday)
mL1 <- fenegbin(n_dead_missing ~ swh_mean_lag1 * post_mou + wind_mean_lag1 | week_year_fac,
                data = d[!is.na(swh_mean_lag1)], vcov = "hetero")

# SWH prev3d (prior 3-day mean)
mL2 <- fenegbin(n_dead_missing ~ swh_mean_prev3d * post_mou + wind_mean_prev3d | week_year_fac,
                data = d[!is.na(swh_mean_prev3d)], vcov = "hetero")

# SWH prev7d (prior 7-day mean — Camarena's main spec)
mL3 <- fenegbin(n_dead_missing ~ swh_mean_prev7d * post_mou + wind_mean_prev7d | week_year_fac,
                data = d[!is.na(swh_mean_prev7d)], vcov = "hetero")

cat("--- Corridor lags ---\n\n")
report_model(mL1, "Lag: SWH lag1 x Post + Wind lag1")
report_model(mL2, "Lag: SWH prev3d x Post + Wind prev3d")
report_model(mL3, "Lag: SWH prev7d x Post + Wind prev7d")


# ============================================================
# 1c. CORE WEATHER: Tighter geography [11,15]x[32,36]
# ============================================================
cat("============================================================\n")
cat("1c. CORE WEATHER + LAGS [11,15]x[32,36]\n")
cat("============================================================\n\n")

# Core day-0
mK0 <- fenegbin(n_dead_missing ~ swh_core * post_mou + wind_core | week_year_fac,
                data = d, vcov = "hetero")

# Core lags
mK1 <- fenegbin(n_dead_missing ~ swh_core_lag1 * post_mou + wind_core_lag1 | week_year_fac,
                data = d[!is.na(swh_core_lag1)], vcov = "hetero")
mK2 <- fenegbin(n_dead_missing ~ swh_core_prev3d * post_mou + wind_core_prev3d | week_year_fac,
                data = d[!is.na(swh_core_prev3d)], vcov = "hetero")
mK3 <- fenegbin(n_dead_missing ~ swh_core_prev7d * post_mou + wind_core_prev7d | week_year_fac,
                data = d[!is.na(swh_core_prev7d)], vcov = "hetero")

cat("--- Core geography ---\n\n")
report_model(mK0, "Core: SWH day0 x Post + Wind")
report_model(mK1, "Core: SWH lag1 x Post + Wind lag1")
report_model(mK2, "Core: SWH prev3d x Post + Wind prev3d")
report_model(mK3, "Core: SWH prev7d x Post + Wind prev7d")


# ============================================================
# 2. Robustness: alternative FE structures
# ============================================================
cat("============================================================\n")
cat("2. ROBUSTNESS: Alternative FE structures\n")
cat("============================================================\n\n")

# Quarter-year FE (coarser)
mR1 <- fenegbin(n_dead_missing ~ swh_mean * post_mou + wind_mean | quarter_fac,
                data = d, vcov = "hetero")
report_model(mR1, "SWH x Post + Wind | quarter-year FE")

# Year + month-of-year FE (Rodriguez-Sanchez style)
mR2 <- fenegbin(n_dead_missing ~ swh_mean * post_mou + wind_mean | year_fac + month_fac,
                data = d, vcov = "hetero")
report_model(mR2, "SWH x Post + Wind | year + month FE")

# Week-year + DOW FE (adds within-week day pattern)
d[, dow_fac := factor(dow)]
mR3 <- tryCatch(
  fenegbin(n_dead_missing ~ swh_mean * post_mou + wind_mean | week_year_fac + dow_fac,
           data = d, vcov = "hetero"),
  error = function(e) { cat("  week-year + DOW failed:", e$message, "\n"); NULL }
)
if (!is.null(mR3)) {
  report_model(mR3, "SWH x Post + Wind | week-year + DOW FE")
}

# Wind as primary with week-year FE
mR4 <- fenegbin(n_dead_missing ~ wind_mean * post_mou + swh_mean | quarter_fac,
                data = d, vcov = "hetero")
report_model(mR4, "Wind x Post + SWH | quarter-year FE")


# ============================================================
# 3. Poisson QMLE comparison
# ============================================================
cat("============================================================\n")
cat("3. POISSON QMLE (robustness to distributional assumption)\n")
cat("============================================================\n\n")

mP1 <- fepois(n_dead_missing ~ swh_mean * post_mou + wind_mean | week_year_fac,
              data = d, vcov = "hetero")
report_model(mP1, "Poisson QMLE: SWH x Post + Wind | week-year FE")

mP2 <- fepois(n_dead_missing ~ wind_mean * post_mou + swh_mean | week_year_fac,
              data = d, vcov = "hetero")
report_model(mP2, "Poisson QMLE: Wind x Post + SWH | week-year FE")


# ============================================================
# 4. Restrict sample period (2014-2021, matching replication data)
# ============================================================
cat("============================================================\n")
cat("4. RESTRICTED SAMPLE: 2014-2021\n")
cat("============================================================\n\n")

d_restr <- d[year <= 2021]
cat(sprintf("Restricted sample: %d days, %d week-year cells\n",
    nrow(d_restr), uniqueN(d_restr$week_year)))

mS1 <- fenegbin(n_dead_missing ~ swh_mean * post_mou + wind_mean | week_year_fac,
                data = d_restr, vcov = "hetero")
report_model(mS1, "SWH x Post + Wind | week-year FE (2014-2021)")

mS2 <- fenegbin(n_dead_missing ~ wind_mean * post_mou + swh_mean | week_year_fac,
                data = d_restr, vcov = "hetero")
report_model(mS2, "Wind x Post + SWH | week-year FE (2014-2021)")


# ============================================================
# 5. Crossings diagnostic (if available — daily deaths only for now)
# ============================================================
cat("============================================================\n")
cat("5. INCIDENTS AS OUTCOME (volume channel diagnostic)\n")
cat("============================================================\n\n")

# If weather x post is null for incident COUNTS, the death effect is DANGER
mV1 <- fenegbin(n_incidents ~ swh_mean * post_mou + wind_mean | week_year_fac,
                data = d[n_incidents > 0 | n_dead_missing > 0 | TRUE],
                vcov = "hetero")
report_model(mV1, "Incidents ~ SWH x Post + Wind | week-year FE")

# Fatal incidents only
mV2 <- fenegbin(n_fatal ~ swh_mean * post_mou + wind_mean | week_year_fac,
                data = d, vcov = "hetero")
report_model(mV2, "Fatal incidents ~ SWH x Post + Wind | week-year FE")


# ============================================================
# 6. Extract and save all results
# ============================================================
cat("============================================================\n")
cat("6. SUMMARY TABLE\n")
cat("============================================================\n\n")

extract_results <- function(model, spec_label, fe_label, sample_label = "Full") {
  ct <- summary(model, vcov = "hetero")$coeftable
  # Find interaction row
  int_rows <- grep(":post_mou|post_mou:", rownames(ct))
  if (length(int_rows) == 0) return(NULL)

  rbindlist(lapply(int_rows, function(r) {
    varname <- sub(":post_mou|post_mou:", "", rownames(ct)[r])
    data.table(
      spec = spec_label, fe = fe_label, sample = sample_label,
      weather_var = varname,
      beta = ct[r, 1], se = ct[r, 2], p = ct[r, 4],
      irr = exp(ct[r, 1]),
      ci_lo = ct[r, 1] - 1.96 * ct[r, 2],
      ci_hi = ct[r, 1] + 1.96 * ct[r, 2],
      n_obs = nobs(model),
      loglik = as.numeric(logLik(model))
    )
  }))
}

results <- rbindlist(list(
  # Day-0, corridor
  extract_results(mA, "SWH only", "week-year"),
  extract_results(mB, "Wind only", "week-year"),
  extract_results(mC, "Gust only", "week-year"),
  extract_results(mD, "SWH + Wind ctrl", "week-year"),
  extract_results(mE, "Wind + SWH ctrl", "week-year"),
  # Corridor lags
  extract_results(mL1, "SWH lag1 + Wind", "week-year"),
  extract_results(mL2, "SWH prev3d + Wind", "week-year"),
  extract_results(mL3, "SWH prev7d + Wind", "week-year"),
  # Core day-0 + lags
  extract_results(mK0, "SWH core day0", "week-year"),
  extract_results(mK1, "SWH core lag1", "week-year"),
  extract_results(mK2, "SWH core prev3d", "week-year"),
  extract_results(mK3, "SWH core prev7d", "week-year"),
  # Robustness FE
  extract_results(mR1, "SWH + Wind ctrl", "quarter-year"),
  extract_results(mR2, "SWH + Wind ctrl", "year+month"),
  if (!is.null(mR3)) extract_results(mR3, "SWH + Wind ctrl", "week-year+DOW"),
  extract_results(mR4, "Wind + SWH ctrl", "quarter-year"),
  # Poisson
  extract_results(mP1, "SWH + Wind (Poisson)", "week-year"),
  extract_results(mP2, "Wind + SWH (Poisson)", "week-year"),
  # Restricted sample
  extract_results(mS1, "SWH + Wind ctrl", "week-year", "2014-2021"),
  extract_results(mS2, "Wind + SWH ctrl", "week-year", "2014-2021")
), fill = TRUE)

# Print summary
cat(sprintf("%-35s %-15s %-10s %8s %8s %8s\n",
    "Specification", "FE", "Sample", "Beta", "IRR", "p"))
cat(paste(rep("-", 92), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i]
  stars <- ifelse(r$p < 0.01, "***",
                  ifelse(r$p < 0.05, "**",
                         ifelse(r$p < 0.1, "*", "")))
  cat(sprintf("%-35s %-15s %-10s %+8.4f %8.4f %8.4f %s\n",
      paste0(r$spec, " (", r$weather_var, ")"),
      r$fe, r$sample, r$beta, r$irr, r$p, stars))
}

# Save
fwrite(results, file.path(BASE_DIR, "output", "tables", "daily_panel_results.csv"))
cat("\nSaved: output/tables/daily_panel_results.csv\n")

# ============================================================
# 7. Coefficient plot
# ============================================================
cat("\n============================================================\n")
cat("7. COEFFICIENT PLOT\n")
cat("============================================================\n\n")

# Plot primary specs (NegBin, full sample)
plot_dt <- results[sample == "Full" & !grepl("Poisson", spec)]
plot_dt[, label := paste0(spec, "\n(", fe, ")")]
plot_dt[, label := factor(label, levels = rev(unique(label)))]

p <- ggplot(plot_dt, aes(x = beta, y = label, colour = weather_var)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2,
                 position = position_dodge(width = 0.4)) +
  labs(
    title = "Daily panel: Weather x Post interaction coefficients",
    subtitle = "NegBin, hetero-robust 95% CI. Deiana-style week-year FE and alternatives.",
    x = expression(hat(beta) ~ "(Weather × Post interaction)"),
    y = NULL, colour = "Weather variable"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(BASE_DIR, "output", "figures", "daily_panel_coefplot.pdf"),
       p, width = 12, height = 8)
ggsave(file.path(BASE_DIR, "output", "figures", "daily_panel_coefplot.png"),
       p, width = 12, height = 8, dpi = 200)
cat("Saved: output/figures/daily_panel_coefplot.pdf + .png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
