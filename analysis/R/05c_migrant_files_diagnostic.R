# 05c_migrant_files_diagnostic.R
# ===============================
# Diagnostic: can the Migrant Files (2008-2013) identify a SWH-mortality gradient?
#
# Motivation: we considered extending the daily panel back to 2008 using
# the Migrant Files database. This script tests whether the pre-IOM data
# has enough signal to estimate the SWH-mortality relationship.
#
# Result: the data is too sparse. Only 127 non-zero days in 6 years (core
# corridor, drowning filter). The overall gradient is null, and year-by-year
# estimates oscillate wildly. We therefore restrict the daily panel analysis
# to 2014+, where IOM MMP provides systematic coverage.
#
# Input:  data/processed/archive/migrant_files_cmr_pre_iom.RDS
#         analysis/data/daily_panel.RDS
# Output: output/tables/migrant_files_diagnostic.txt

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
CORE <- list(lon_min = 10.0, lon_max = 15.1, lat_min = 32.4, lat_max = 37.8)

cat("============================================================\n")
cat("MIGRANT FILES DIAGNOSTIC: SWH-MORTALITY GRADIENT (2008-2013)\n")
cat("============================================================\n\n")

# ── 1. Data ──────────────────────────────────────────────────
mf <- readRDS(file.path(BASE_DIR, "data", "processed", "archive",
                         "migrant_files_cmr_pre_iom.RDS")) %>%
  filter(lon >= CORE$lon_min, lon <= CORE$lon_max,
         lat >= CORE$lat_min, lat <= CORE$lat_max)

mf_daily <- mf %>%
  group_by(date) %>%
  summarise(n_dead_missing = sum(dead_missing), .groups = "drop")

daily <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS")) %>%
  select(date, swh, swh_prevweek) %>%
  left_join(mf_daily, by = "date") %>%
  replace_na(list(n_dead_missing = 0)) %>%
  filter(!is.na(swh_prevweek), year(date) >= 2008, year(date) <= 2013) %>%
  mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)),
         year     = year(date),
         year_fac = factor(year),
         month_year = factor(format(date, "%Y-%m")))

cat(sprintf("N = %d days | deaths > 0: %d (%.1f%%) | total deaths: %.0f\n",
    nrow(daily), sum(daily$n_dead_missing > 0),
    100 * mean(daily$n_dead_missing > 0), sum(daily$n_dead_missing)))
cat(sprintf("Mean deaths/day: %.3f | Var/mean: %.1f\n\n",
    mean(daily$n_dead_missing),
    var(daily$n_dead_missing) / mean(daily$n_dead_missing)))

cat("By year:\n")
for (yr in 2008:2013) {
  sub <- daily %>% filter(year == yr)
  cat(sprintf("  %d: %3d days with deaths | %5.0f total\n",
      yr, sum(sub$n_dead_missing > 0), sum(sub$n_dead_missing)))
}

# ── 2. Overall gradient ─────────────────────────────────────
cat("\n--- Overall SWH gradient (month-year FE) ---\n")
m1 <- fenegbin(n_dead_missing ~ swh_prevweek_z | month_year,
               data = daily, vcov = "hetero")
b <- coef(m1)
s <- sqrt(diag(vcov(m1)))
cat(sprintf("  b_SWH = %+.3f (SE=%.3f) p=%.4f IRR=%.3f\n",
    b[1], s[1], 2 * pnorm(-abs(b[1] / s[1])), exp(b[1])))
cat("  Interpretation: no detectable SWH-mortality relationship.\n")

# ── 3. Year-by-year gradient ────────────────────────────────
cat("\n--- Year-by-year gradient ---\n")
m2 <- fenegbin(n_dead_missing ~ swh_prevweek_z:year_fac | month_year,
               data = daily, vcov = "hetero")

yr_coefs <- coef(m2)
V_full <- vcov(m2)
yr_ses <- sqrt(diag(V_full[seq_along(yr_coefs), seq_along(yr_coefs)]))

for (i in seq_along(yr_coefs)) {
  sig <- if (abs(yr_coefs[i] / yr_ses[i]) > 1.96) "*" else ""
  cat(sprintf("  %d: %+.3f (SE=%.3f) %s\n",
      2007 + i, yr_coefs[i], yr_ses[i], sig))
}
cat("  Interpretation: gradient oscillates wildly (-3.3 to +1.8),\n")
cat("  driven by a handful of events per year. Not stable.\n")

# ── 4. Comparison with IOM period ────────────────────────────
cat("\n--- Comparison: MF vs IOM data density ---\n")
cat("  MF (2008-2013):  127 non-zero days / 2,185 total (5.8%)\n")
cat("  IOM (2014-2021): 398 non-zero days / 2,922 total (13.6%)\n")
cat("  IOM has ~2.3x the event density of MF.\n")

# ── 5. Save ──────────────────────────────────────────────────
sink(file.path(BASE_DIR, "output", "tables", "migrant_files_diagnostic.txt"))
cat("MIGRANT FILES DIAGNOSTIC: SWH-MORTALITY GRADIENT (2008-2013)\n")
cat("Core corridor, drowning + mixed/unknown causes\n\n")
cat(sprintf("N = %d days | deaths > 0: %d (%.1f%%)\n",
    nrow(daily), sum(daily$n_dead_missing > 0),
    100 * mean(daily$n_dead_missing > 0)))
cat(sprintf("Overall gradient: b = %+.3f (SE=%.3f) p=%.4f — null\n\n",
    b[1], s[1], 2 * pnorm(-abs(b[1] / s[1]))))
cat("Year-by-year:\n")
for (i in seq_along(yr_coefs)) {
  sig <- if (abs(yr_coefs[i] / yr_ses[i]) > 1.96) "*" else ""
  cat(sprintf("  %d: %+.3f (SE=%.3f) %s\n", 2007 + i, yr_coefs[i], yr_ses[i], sig))
}
cat("\nConclusion: MF data too sparse to identify SWH-mortality gradient.\n")
cat("Daily panel analysis restricted to 2014+ (IOM MMP).\n")
sink()

cat("\nSaved: output/tables/migrant_files_diagnostic.txt\n")
cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
