# 05f_se_diagnostic.R
# ===================
# Diagnose serial correlation in residuals and compare SE approaches.
#
# 1. Residual ACF from the primary NegBin model
# 2. Breusch-Godfrey test for serial correlation
# 3. SE comparison: HC vs NW(7) vs NW(14) vs NW(28) vs cluster(iso_week)
#
# Input:  analysis/data/daily_panel.RDS, data/processed/iom_mmp_incidents.RDS
# Output: output/figures/se_diagnostic_acf.png
#         output/tables/se_diagnostic_comparison.txt

library(tidyverse)
library(lubridate)
library(fixest)
library(lmtest)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
YEAR_START <- 2014
YEAR_END <- 2024
SEA_CAUSES <- c("Drowning", "Mixed or unknown")
CORE <- list(lon_min = 10.0, lon_max = 15.1, lat_min = 32.4, lat_max = 37.8)

cat("============================================================\n")
cat("SE DIAGNOSTIC: autocorrelation + SE comparison\n")
cat("============================================================\n\n")

# ── 1. Build data (same as 05_reduced_form_model.R) ────────
cat("--- 1. Data preparation ---\n")

iom <- readRDS(file.path(BASE_DIR, "data", "processed",
                           "iom_mmp_incidents.RDS")) %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) == "incident",
         `Cause of death (category)` %in% SEA_CAUSES,
         Longitude >= CORE$lon_min, Longitude <= CORE$lon_max,
         Latitude  >= CORE$lat_min, Latitude  <= CORE$lat_max) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date))

daily_deaths <- iom %>%
  group_by(date) %>%
  summarise(n_dead_missing = sum(dead_missing), .groups = "drop")

daily_full <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS")) %>%
  select(date, swh, swh_prevweek, iso_week) %>%
  left_join(daily_deaths, by = "date") %>%
  replace_na(list(n_dead_missing = 0)) %>%
  arrange(date) %>%
  mutate(
    post_mou   = as.integer(date >= MOU_DATE),
    year       = year(date),
    month_year = factor(format(date, "%Y-%m")),
    unit       = 1L
  ) %>%
  filter(!is.na(swh_prevweek), year >= YEAR_START, year <= YEAR_END)

d <- daily_full %>%
  mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)))

cat(sprintf("  N = %d days | Pre: %d | Post: %d\n",
    nrow(d), sum(d$post_mou == 0), sum(d$post_mou == 1)))

# ── 2. Estimate primary model ─────────────────────────────
cat("\n--- 2. Primary model (HC SEs) ---\n")

m <- fenegbin(n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year,
              data = d, vcov = "hetero", panel.id = ~unit + date)

cat("  Primary estimates (HC):\n")
print(summary(m))

# ── 3. Residual autocorrelation ───────────────────────────
cat("\n--- 3. Residual autocorrelation ---\n")

# Pearson residuals — model drops singletons, so filter to estimation sample
obs_removed <- m$obs_selection$obsRemoved
d_est <- d[setdiff(seq_len(nrow(d)), abs(obs_removed)), ]
d_est$resid <- as.numeric(residuals(m, type = "pearson"))

# ACF
acf_vals <- acf(d_est$resid, lag.max = 35, plot = FALSE)

png(file.path(BASE_DIR, "output", "figures", "se_diagnostic_acf.png"),
    width = 10, height = 6, units = "in", res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# ACF plot
plot(acf_vals, main = "ACF of Pearson residuals (NegBin primary model)",
     xlab = "Lag (days)", ylab = "Autocorrelation")

# PACF plot
pacf_vals <- pacf(d_est$resid, lag.max = 35, plot = FALSE)
plot(pacf_vals, main = "PACF of Pearson residuals",
     xlab = "Lag (days)", ylab = "Partial autocorrelation")

dev.off()
cat("Saved: output/figures/se_diagnostic_acf.png\n")

# Report first few autocorrelations
cat("\n  Autocorrelations (lags 1-7):\n")
for (k in 1:7) {
  cat(sprintf("    Lag %d: %.4f\n", k, acf_vals$acf[k + 1]))
}

# Breusch-Godfrey via a linear proxy (OLS on same spec)
cat("\n  Breusch-Godfrey test (linear OLS proxy):\n")
m_ols <- lm(n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou + month_year,
            data = d_est)
for (order in c(7, 14, 28)) {
  bg <- bgtest(m_ols, order = order)
  cat(sprintf("    Order %2d: LM = %.2f, p = %.6f\n", order, bg$statistic, bg$p.value))
}

# ── 4. SE comparison ──────────────────────────────────────
cat("\n--- 4. SE comparison ---\n")

vcov_list <- list(
  "HC (current)"  = "hetero",
  "NW(7)"         = NW(7),
  "NW(14)"        = NW(14),
  "NW(28)"        = NW(28),
  "Cluster(week)" = ~iso_week
)

cat("\n  Interaction coefficient: swh_prevweek_z x post_mou\n")
cat("  ─────────────────────────────────────────────────────\n")
cat(sprintf("  %-18s %8s %8s %8s %8s\n", "VCOV", "beta", "SE", "z", "p"))
cat("  ─────────────────────────────────────────────────────\n")

comparison_rows <- list()

for (nm in names(vcov_list)) {
  s <- summary(m, vcov = vcov_list[[nm]])
  ct <- coeftable(s)
  # interaction is the second row
  b  <- ct[2, 1]
  se <- ct[2, 2]
  z  <- b / se
  p  <- 2 * pnorm(-abs(z))
  cat(sprintf("  %-18s %+8.4f %8.4f %+8.3f %8.4f\n", nm, b, se, z, p))
  comparison_rows[[nm]] <- tibble(vcov = nm, beta = b, se = se, z = z, p = p)
}

comparison_df <- bind_rows(comparison_rows)

cat("\n  SE inflation relative to HC:\n")
hc_se <- comparison_df$se[1]
for (i in seq_len(nrow(comparison_df))) {
  cat(sprintf("    %-18s  SE ratio = %.2f\n",
      comparison_df$vcov[i], comparison_df$se[i] / hc_se))
}

# ── 5. Save results ───────────────────────────────────────
cat("\n--- 5. Saving results ---\n")

sink(file.path(BASE_DIR, "output", "tables", "se_diagnostic_comparison.txt"))
cat("SE DIAGNOSTIC: primary NegBin model\n")
cat("deaths ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year\n")
cat(sprintf("Sample: %d-%d, N = %d\n\n", YEAR_START, YEAR_END, nrow(d)))

cat("Interaction coefficient (swh_prevweek_z x post_mou):\n")
cat(sprintf("%-18s %8s %8s %8s %8s\n", "VCOV", "beta", "SE", "z", "p"))
cat(strrep("-", 56), "\n")
for (i in seq_len(nrow(comparison_df))) {
  r <- comparison_df[i, ]
  cat(sprintf("%-18s %+8.4f %8.4f %+8.3f %8.4f\n", r$vcov, r$beta, r$se, r$z, r$p))
}
cat("\nSE inflation relative to HC:\n")
for (i in seq_len(nrow(comparison_df))) {
  cat(sprintf("  %-18s  %.2fx\n", comparison_df$vcov[i], comparison_df$se[i] / hc_se))
}

cat("\nResidual autocorrelations (lags 1-7):\n")
for (k in 1:7) {
  cat(sprintf("  Lag %d: %.4f\n", k, acf_vals$acf[k + 1]))
}

cat("\nBreusch-Godfrey (OLS proxy):\n")
for (order in c(7, 14, 28)) {
  bg <- bgtest(m_ols, order = order)
  cat(sprintf("  Order %2d: LM = %.2f, p = %.6f\n", order, bg$statistic, bg$p.value))
}
sink()
cat("Saved: output/tables/se_diagnostic_comparison.txt\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
