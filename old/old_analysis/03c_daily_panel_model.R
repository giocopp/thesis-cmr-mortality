# 03c_daily_panel_model.R
# ========================
# Primary analysis: daily panel NegBin model (adapted from Deiana et al. 2024).
#
# Design:
#   n_dead_missing_t ~ NegBin(mu_t)
#   log(mu_t) = week-year_FE + beta_1 * SWH_t + beta_3 * (SWH_t x post_t)
#
# Outcome filtering (Zambiasi & Albarosa 2025 precedent):
#   Primary: drowning incidents only (SWH does not cause sickness/violence)
#   Robustness 1: drowning + dead/missing > 1 (drops body-found records)
#   Robustness 2: all incidents (unfiltered, for comparability)
#
# Identification:
#   Week-by-year FE absorb all slow-moving confounders (seasonality,
#   trends, policy shifts within a week are negligible). beta_3 is
#   identified from day-to-day weather variation WITHIN each week.
#
# Geography:
#   Weather = daily spatial mean over core CMR channel [11,15]x[32,36],
#   where ~80% of incidents occur.  This is the weather actually
#   experienced by boats on the Libya-to-Lampedusa corridor.
#
# Timing:
#   Day-0 (contemporaneous) is primary: weather on the day the incident
#   is recorded.  Lag-1 (yesterday) is a robustness variant, reflecting
#   that boats depart the prior day.
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

MOU_DATE <- as.Date("2017-07-01")

# ============================================================
# 1. Sample description
# ============================================================
cat("============================================================\n")
cat("DAILY PANEL MODEL: Deaths ~ Weather x Post + Week-Year FE\n")
cat("============================================================\n\n")

cat(sprintf("Panel: %d days (%s to %s)\n", nrow(d), min(d$date), max(d$date)))
cat(sprintf("Week-year cells: %d (mean %.1f days/cell)\n",
    uniqueN(d$week_year), nrow(d) / uniqueN(d$week_year)))
cat(sprintf("Outcome: n_dead_missing (mean=%.2f, max=%d, zeros=%d [%.1f%%])\n",
    mean(d$n_dead_missing), max(d$n_dead_missing),
    sum(d$n_dead_missing == 0), 100 * mean(d$n_dead_missing == 0)))
cat(sprintf("Variance/mean: %.1f (overdispersion → NegBin)\n\n",
    var(d$n_dead_missing) / mean(d$n_dead_missing)))


# ============================================================
# 2. Primary specification: SWH lag-1
# ============================================================
cat("============================================================\n")
cat("2. PRIMARY SPEC: SWH lag-1 x Post | week-year FE\n")
cat("============================================================\n\n")

# IOM MMP incident_date is the REPORTING date, not the crossing date.
# Lag-1 (yesterday's weather) captures conditions during departure/transit,
# which is 1-2 days before the incident is documented.
m_primary <- fenegbin(
  n_dead_missing ~ swh_core_lag1 + swh_core_lag1:post_mou | week_year_fac,
  data = d[!is.na(swh_core_lag1)], vcov = "hetero"
)

# --- Reporting helper ---
report <- function(m, label, vcov_type = "hetero") {
  ct <- summary(m, vcov = vcov_type)$coeftable
  theta_str <- if (!is.null(m$theta)) sprintf(", theta=%.3f", m$theta) else ""
  cat(sprintf("  %s (N=%d%s):\n", label, nobs(m), theta_str))
  for (i in seq_len(nrow(ct))) {
    stars <- ifelse(ct[i, 4] < 0.01, "***",
             ifelse(ct[i, 4] < 0.05, "**",
             ifelse(ct[i, 4] < 0.1, "*", "")))
    cat(sprintf("    %-30s %+8.4f (SE=%6.4f) p=%6.4f %s  IRR=%.4f\n",
        rownames(ct)[i], ct[i, 1], ct[i, 2], ct[i, 4], stars, exp(ct[i, 1])))
  }
  cat("\n")
}

# --- Extraction helper ---
extract_int <- function(model, spec_label, vcov_type = "hetero") {
  ct <- summary(model, vcov = vcov_type)$coeftable
  int_rows <- grep(":post_mou|post_mou:", rownames(ct))
  if (length(int_rows) == 0) return(NULL)
  rbindlist(lapply(int_rows, function(r) {
    data.table(
      spec = spec_label,
      coef = rownames(ct)[r],
      beta = ct[r, 1], se = ct[r, 2], p = ct[r, 4],
      irr = exp(ct[r, 1]),
      ci_lo = ct[r, 1] - 1.96 * ct[r, 2],
      ci_hi = ct[r, 1] + 1.96 * ct[r, 2],
      n_obs = nobs(model),
      theta = if (!is.null(model$theta)) model$theta else NA_real_
    )
  }))
}

report(m_primary, "PRIMARY: SWH lag-1 x Post | week-year FE")


# ============================================================
# 3. Monthly FE (month-year)
# ============================================================
cat("============================================================\n")
cat("3. MONTHLY FE: SWH lag-1 x Post | month-year FE\n")
cat("============================================================\n\n")

# Monthly FE (~144 cells) are less aggressive than weekly FE (~628 cells).
# More data survives singleton removal — important given 84% zeros.
# Identification still from day-to-day weather variation, but within
# a month rather than within a week.
m_monthly <- fenegbin(
  n_dead_missing ~ swh_core_lag1 + swh_core_lag1:post_mou | month_year_fac,
  data = d[!is.na(swh_core_lag1)], vcov = "hetero"
)
report(m_monthly, "SWH lag-1 x Post | month-year FE")


# ============================================================
# 4. Timing robustness: lag-2, lag-3 (each alone)
# ============================================================
cat("============================================================\n")
cat("3. TIMING ROBUSTNESS: Lag-2, Lag-3\n")
cat("============================================================\n\n")

m_lag2_wk <- fenegbin(
  n_dead_missing ~ swh_core_lag2 + swh_core_lag2:post_mou | week_year_fac,
  data = d[!is.na(swh_core_lag2)], vcov = "hetero"
)
report(m_lag2_wk, "SWH lag-2 x Post | week-year FE")

m_lag2_mo <- fenegbin(
  n_dead_missing ~ swh_core_lag2 + swh_core_lag2:post_mou | month_year_fac,
  data = d[!is.na(swh_core_lag2)], vcov = "hetero"
)
report(m_lag2_mo, "SWH lag-2 x Post | month-year FE")

m_lag3_wk <- fenegbin(
  n_dead_missing ~ swh_core_lag3 + swh_core_lag3:post_mou | week_year_fac,
  data = d[!is.na(swh_core_lag3)], vcov = "hetero"
)
report(m_lag3_wk, "SWH lag-3 x Post | week-year FE")

m_lag3_mo <- fenegbin(
  n_dead_missing ~ swh_core_lag3 + swh_core_lag3:post_mou | month_year_fac,
  data = d[!is.na(swh_core_lag3)], vcov = "hetero"
)
report(m_lag3_mo, "SWH lag-3 x Post | month-year FE")

# Lag-7: diagnostic. No plausible mechanism — crossing takes 1-2 days.
# If significant, suggests FE structure is too loose or autocorrelation
# in weather creates spurious significance.
m_lag7_wk <- fenegbin(
  n_dead_missing ~ swh_core_lag7 + swh_core_lag7:post_mou | week_year_fac,
  data = d[!is.na(swh_core_lag7)], vcov = "hetero"
)
report(m_lag7_wk, "SWH lag-7 x Post | week-year FE")

m_lag7_mo <- fenegbin(
  n_dead_missing ~ swh_core_lag7 + swh_core_lag7:post_mou | month_year_fac,
  data = d[!is.na(swh_core_lag7)], vcov = "hetero"
)
report(m_lag7_mo, "SWH lag-7 x Post | month-year FE")


# ============================================================
# 5. Alternative weather variable: wind
# ============================================================
cat("============================================================\n")
cat("6. ALTERNATIVE WEATHER: Wind\n")
cat("============================================================\n\n")

m_wind <- fenegbin(
  n_dead_missing ~ wind_core + wind_core:post_mou | week_year_fac,
  data = d, vcov = "hetero"
)
report(m_wind, "Wind_core x Post | week-year FE")


# ============================================================
# 6. Clustered SEs
# ============================================================
cat("============================================================\n")
cat("7. CLUSTERED STANDARD ERRORS (primary spec)\n")
cat("============================================================\n\n")

report(m_primary, "Cluster by week-year", vcov_type = ~week_year_fac)


# ============================================================
# 7. Restricted sample: 2014-2021
# ============================================================
cat("============================================================\n")
cat("8. RESTRICTED SAMPLE: 2014-2021\n")
cat("============================================================\n\n")

d_restr <- d[year <= 2021]
cat(sprintf("Restricted: %d days, %d week-year cells\n\n",
    nrow(d_restr), uniqueN(d_restr$week_year)))

m_restr <- fenegbin(
  n_dead_missing ~ swh_core_lag1 + swh_core_lag1:post_mou | week_year_fac,
  data = d_restr[!is.na(d_restr$swh_core_lag1)], vcov = "hetero"
)
report(m_restr, "SWH lag-1 x Post | week-year FE [2014-2021]")


# ============================================================
# 8. All incidents (unfiltered, for comparability)
# ============================================================
cat("============================================================\n")
cat("9. ALL INCIDENTS (unfiltered)\n")
cat("============================================================\n\n")

cat(sprintf("All incidents outcome: var/mean = %.1f\n\n",
    var(d$all_n_dead_missing) / mean(d$all_n_dead_missing)))

m_all <- fenegbin(
  all_n_dead_missing ~ swh_core_lag1 + swh_core_lag1:post_mou | week_year_fac,
  data = d[!is.na(swh_core_lag1)], vcov = "hetero"
)
report(m_all, "All incidents (unfiltered) | week-year FE")


# ============================================================
# 9. Summary table
# ============================================================
cat("============================================================\n")
cat("10. SUMMARY TABLE\n")
cat("============================================================\n\n")

results <- rbindlist(list(
  extract_int(m_primary,  "SWH lag-1 | week-year FE"),
  extract_int(m_monthly,  "SWH lag-1 | month-year FE"),
  extract_int(m_lag2_wk,  "SWH lag-2 | week-year FE"),
  extract_int(m_lag2_mo,  "SWH lag-2 | month-year FE"),
  extract_int(m_lag3_wk,  "SWH lag-3 | week-year FE"),
  extract_int(m_lag3_mo,  "SWH lag-3 | month-year FE"),
  extract_int(m_lag7_wk,  "SWH lag-7 | week-year FE"),
  extract_int(m_lag7_mo,  "SWH lag-7 | month-year FE"),
  extract_int(m_wind,     "Wind lag-0 | week-year FE"),
  extract_int(m_primary,  "SWH lag-1 | week-year (clustered SE)", vcov_type = ~week_year_fac),
  extract_int(m_restr,    "SWH lag-1 | week-year [2014-2021]"),
  extract_int(m_all,      "All incidents (unfiltered)")
), fill = TRUE)

cat(sprintf("%-35s %+8s %7s %7s %8s %5s\n",
    "Specification", "Beta", "SE", "p", "IRR", "N"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i]
  stars <- ifelse(r$p < 0.01, "***",
           ifelse(r$p < 0.05, "**",
           ifelse(r$p < 0.1, "*", "")))
  cat(sprintf("%-35s %+8.4f %7.4f %7.4f %8.4f %5d %s\n",
      r$spec, r$beta, r$se, r$p, r$irr, r$n_obs, stars))
}

fwrite(results, file.path(BASE_DIR, "output", "tables", "daily_panel_results.csv"))
cat("\nSaved: output/tables/daily_panel_results.csv\n")


# ============================================================
# 10. Coefficient plot
# ============================================================
cat("\n============================================================\n")
cat("10. COEFFICIENT PLOT\n")
cat("============================================================\n\n")

# Plot the main weather x post interaction specs
plot_dt <- results[!grepl("clustered|Outcome", spec)]
plot_dt[, label := factor(spec, levels = rev(unique(spec)))]

p <- ggplot(plot_dt, aes(x = beta, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.5, colour = "#2166AC") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2,
                 colour = "#2166AC") +
  labs(
    title = "Daily panel: Weather x Post interaction (NegBin, week-year FE)",
    subtitle = "Core geography [11-15, 32-36]. Hetero-robust 95% CI.",
    x = expression(hat(beta)[3] ~ "(Weather × Post interaction)"),
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(BASE_DIR, "output", "figures", "daily_panel_coefplot.pdf"),
       p, width = 10, height = 6)
ggsave(file.path(BASE_DIR, "output", "figures", "daily_panel_coefplot.png"),
       p, width = 10, height = 6, dpi = 200)
cat("Saved: output/figures/daily_panel_coefplot.pdf + .png\n")


# ============================================================
# 11. Save model objects
# ============================================================
model_dir <- file.path(BASE_DIR, "output", "models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

models <- list(
  primary   = m_primary,
  monthly   = m_monthly,
  lag2_wk   = m_lag2_wk,
  lag2_mo   = m_lag2_mo,
  lag3_wk   = m_lag3_wk,
  lag3_mo   = m_lag3_mo,
  lag7_wk   = m_lag7_wk,
  lag7_mo   = m_lag7_mo,
  wind      = m_wind,
  restr     = m_restr,
  all_inc   = m_all
)
saveRDS(models, file.path(model_dir, "daily_panel_models.RDS"))
cat("Saved: output/models/daily_panel_models.RDS\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
