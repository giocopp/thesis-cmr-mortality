# 04b_sensitivity_analysis.R
# =========================
# Cinelli & Hazlett (2020) omitted variable sensitivity analysis.
#
# Purpose: Quantify how strong an unobserved confounder would need to be
#          to explain away the SWH x Post interaction.
#
# Method: sensemakr computes:
#   - Robustness Value (RV): minimum strength of confounding (in terms of
#     partial R^2 with treatment and outcome) needed to reduce the estimate
#     to zero or non-significance
#   - Benchmarks: compare required confounding strength to observed covariates
#
# Limitation: sensemakr works with OLS. Our primary model is NegBin.
# We run the OLS version (linear model of counts) as a conservative proxy.
# If the OLS result is robust to omitted variables, the NegBin result
# (which fits the DGP better) is likely robust too.
#
# Input:  data/processed/cmr_daily_weather_panel.RDS
#         data/processed/cmr_events_with_weather.RDS
# Output: output/tables/sensitivity_analysis.csv
#         output/figures/sensitivity_contour.pdf
#         printed diagnostics

library(sensemakr)
library(fixest)
library(data.table)
library(ggplot2)

BASE_DIR <- here::here()

# ============================================================
# 1. Load data
# ============================================================
cat("============================================================\n")
cat("CINELLI & HAZLETT SENSITIVITY ANALYSIS\n")
cat("============================================================\n\n")

d <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                      "cmr_daily_weather_panel.RDS")))

cat(sprintf("Panel: %d days\n\n", nrow(d)))


# ============================================================
# 2. OLS version of primary spec (sensemakr needs lm object)
# ============================================================
# sensemakr requires a base R lm() object. We include week-year as
# factor dummies. With 628 levels this is feasible for lm().
#
# The OLS model: n_dead_missing ~ swh_core * post_mou + wind_core + week_year_fac
# Treatment of interest: the swh_core:post_mou interaction

cat("============================================================\n")
cat("2. OLS MODEL FOR SENSEMAKR\n")
cat("============================================================\n\n")

m_ols <- lm(n_dead_missing ~ swh_core * post_mou + wind_core + week_year_fac,
             data = d)

# Check the interaction coefficient
s <- summary(m_ols)$coefficients
int_name <- grep("swh_core:post_mou|post_mou:swh_core", rownames(s), value = TRUE)
cat(sprintf("OLS interaction: %s\n", int_name))
cat(sprintf("  beta=%+.4f, SE=%.4f, t=%.3f, p=%.4f\n",
    s[int_name, 1], s[int_name, 2], s[int_name, 3], s[int_name, 4]))


# ============================================================
# 3. sensemakr: sensitivity to omitted variables
# ============================================================
cat("\n============================================================\n")
cat("3. SENSITIVITY ANALYSIS\n")
cat("============================================================\n\n")

# Run sensemakr
# treatment: the interaction term
# benchmark_covariates: observed covariates to benchmark against
sens <- sensemakr(
  model = m_ols,
  treatment = int_name,
  benchmark_covariates = "wind_core",
  kd = c(1, 2, 3, 5)  # multiples of benchmark strength
)

cat("--- sensemakr summary ---\n\n")
print(summary(sens))

# Extract key quantities
cat("\n--- Key robustness quantities ---\n\n")
cat(sprintf("Robustness Value (RV, q=1):    %.4f\n", sens$sensitivity_stats$rv_q))
cat(sprintf("  An unobserved confounder would need partial R^2 of %.2f%%\n",
    100 * sens$sensitivity_stats$rv_q))
cat(sprintf("  with BOTH treatment and outcome to bring beta to zero.\n\n"))

cat(sprintf("Robustness Value (RV, q=1, alpha=0.05): %.4f\n",
    sens$sensitivity_stats$rv_qa))
cat(sprintf("  Confounder strength needed to make result non-significant.\n\n"))

# Benchmark interpretation
cat("--- Benchmark interpretation ---\n")
cat("How strong would the confounder need to be relative to wind_core?\n\n")
if (!is.null(sens$bounds)) {
  bounds <- sens$bounds
  print(bounds)
}


# ============================================================
# 4. sensemakr contour plot
# ============================================================
cat("\n--- Generating contour plot ---\n")

pdf(file.path(BASE_DIR, "output", "figures", "sensitivity_contour.pdf"),
    width = 10, height = 7)
plot(sens, sensitivity.of = "estimate")
dev.off()

pdf(file.path(BASE_DIR, "output", "figures", "sensitivity_contour_tval.pdf"),
    width = 10, height = 7)
plot(sens, sensitivity.of = "t-value")
dev.off()

cat("Saved: output/figures/sensitivity_contour.pdf\n")
cat("Saved: output/figures/sensitivity_contour_tval.pdf\n")


# ============================================================
# 5. Repeat with corridor SWH (broader geography)
# ============================================================
cat("\n============================================================\n")
cat("5. SENSITIVITY: CORRIDOR SWH (broader geography)\n")
cat("============================================================\n\n")

m_ols_corr <- lm(n_dead_missing ~ swh_mean * post_mou + wind_mean + week_year_fac,
                  data = d)
int_corr <- grep("swh_mean:post_mou|post_mou:swh_mean",
                  rownames(summary(m_ols_corr)$coefficients), value = TRUE)
s_corr <- summary(m_ols_corr)$coefficients
cat(sprintf("Corridor OLS: %s = %+.4f, p=%.4f\n",
    int_corr, s_corr[int_corr, 1], s_corr[int_corr, 4]))

sens_corr <- sensemakr(
  model = m_ols_corr,
  treatment = int_corr,
  benchmark_covariates = "wind_mean",
  kd = c(1, 2, 3)
)
cat(sprintf("Corridor RV (q=1): %.4f\n", sens_corr$sensitivity_stats$rv_q))
cat(sprintf("Corridor RV (q=1, alpha=0.05): %.4f\n",
    sens_corr$sensitivity_stats$rv_qa))


# ============================================================
# 6. Event-level sensitivity
# ============================================================
cat("\n============================================================\n")
cat("6. EVENT-LEVEL SENSITIVITY\n")
cat("============================================================\n\n")

df_ev <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                          "cmr_events_with_weather.RDS")))
df_ev[, grid_1deg := paste0(sprintf("%.0f", round(grid_lat)), "_",
                              sprintf("%.0f", round(grid_lon)))]
df_ev[, month_fac := factor(month(date))]

# OLS at event level
m_ev <- lm(dead_missing ~ swh_day0 * post_mou + wind_day0 +
             grid_1deg + year_fac + month_fac,
           data = df_ev[!is.na(swh_day0)])

int_ev <- grep("swh_day0:post_mou|post_mou:swh_day0",
                rownames(summary(m_ev)$coefficients), value = TRUE)
s_ev <- summary(m_ev)$coefficients
cat(sprintf("Event-level OLS: %s = %+.4f, p=%.4f\n",
    int_ev, s_ev[int_ev, 1], s_ev[int_ev, 4]))

sens_ev <- sensemakr(
  model = m_ev,
  treatment = int_ev,
  benchmark_covariates = "wind_day0",
  kd = c(1, 2, 3)
)
cat(sprintf("Event-level RV (q=1): %.4f\n", sens_ev$sensitivity_stats$rv_q))
cat(sprintf("Event-level RV (q=1, alpha=0.05): %.4f\n",
    sens_ev$sensitivity_stats$rv_qa))
print(summary(sens_ev))


# ============================================================
# 7. Summary table
# ============================================================
cat("\n============================================================\n")
cat("7. SUMMARY\n")
cat("============================================================\n\n")

results <- data.table(
  model = c("Daily panel, core SWH (primary)",
            "Daily panel, corridor SWH",
            "Event-level, incident SWH"),
  beta = c(s[int_name, 1], s_corr[int_corr, 1], s_ev[int_ev, 1]),
  se = c(s[int_name, 2], s_corr[int_corr, 2], s_ev[int_ev, 2]),
  p = c(s[int_name, 4], s_corr[int_corr, 4], s_ev[int_ev, 4]),
  rv_q1 = c(sens$sensitivity_stats$rv_q,
            sens_corr$sensitivity_stats$rv_q,
            sens_ev$sensitivity_stats$rv_q),
  rv_q1_a05 = c(sens$sensitivity_stats$rv_qa,
                sens_corr$sensitivity_stats$rv_qa,
                sens_ev$sensitivity_stats$rv_qa)
)

cat(sprintf("%-40s %8s %8s %8s %8s\n",
    "Model", "beta", "p", "RV(q=1)", "RV(sig)"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i]
  cat(sprintf("%-40s %+8.4f %8.4f %8.4f %8.4f\n",
      r$model, r$beta, r$p, r$rv_q1, r$rv_q1_a05))
}

cat("\nRV(q=1): partial R^2 of confounder with both D and Y needed to explain away beta\n")
cat("RV(sig): partial R^2 needed to make result non-significant at alpha=0.05\n")

fwrite(results, file.path(BASE_DIR, "output", "tables", "sensitivity_analysis.csv"))
cat("\nSaved: output/tables/sensitivity_analysis.csv\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
