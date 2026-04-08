# 07b_temporal_disaggregation.R
# =============================
# Temporal disaggregation of monthly LCG/TCG interceptions to daily frequency.
#
# Method: Proportional Denton disaggregation (Denton, 1971).
#   Monthly interceptions are distributed to daily values in proportion to a
#   high-frequency indicator series, preserving monthly totals exactly.
#
#   Formula:
#     daily_d = monthly_m * (indicator_d / sum(indicator in month m))
#
#   When the indicator is zero for an entire month but interceptions are positive,
#   the script falls back to uniform distribution (monthly / days_in_month).
#
# Indicator choice:
#   LCG interceptions -> Frontex daily persons departing from Libya
#   TCG interceptions -> Frontex daily persons departing from Tunisia
#   Rationale: coast guard interceptions track the temporal pattern of departure
#   attempts from the respective country.
#
# Why proportional Denton over Chow-Lin or Litterman:
#   (a) The tempdisagg R package (Sax & Steiner 2013) does not natively support
#       monthly-to-daily conversion (only annual->quarterly, annual->monthly,
#       quarterly->monthly).
#   (b) The indicator-target relationship is conceptually direct (departures
#       predict interceptions), so GLS estimation adds complexity without
#       clear benefit.
#   (c) The formula is transparent and auditable for academic work.
#   (d) No new package dependencies are required.
#
# References:
#   Denton, F.T. (1971). "Adjustment of monthly or quarterly series to annual
#     totals." JASA, 66(333), 99-102.
#   Sax, C. & Steiner, P. (2013). "Temporal disaggregation of time series."
#     The R Journal, 5(2), 80-87.
#   Chow, G. & Lin, A. (1971). "Best linear unbiased interpolation, distribution,
#     and extrapolation of time series by related series." Review of Economics
#     and Statistics, 53, 372-375.
#
# Input:
#   data/processed/iom_med_crossings_monthly.RDS  (monthly LCG/TCG interceptions)
#   data/processed/frontex_incidents.RDS           (incident-level, for daily indicator)
#
# Output:
#   analysis/data/interceptions_daily_disagg.RDS
#
# Panel scope: 2014-01-01 to Frontex end date (~2023-06-11).
#   Months with NA interceptions are set to 0.
#   Post-Frontex months are excluded (no daily indicator available).

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()
CMR_DEPARTURES <- c("Libya", "Tunisia", "Algeria")

cat("============================================================\n")
cat("TEMPORAL DISAGGREGATION: MONTHLY INTERCEPTIONS -> DAILY\n")
cat("============================================================\n\n")

# ── 1. Load monthly interceptions ────────────────────────────
cat("--- 1. Loading monthly interceptions ---\n")

monthly_raw <- readRDS(file.path(BASE_DIR, "data", "processed",
                                  "iom_med_crossings_monthly.RDS"))

monthly <- monthly_raw %>%
  transmute(
    ym          = as.Date(date),
    lcg_monthly = replace_na(as.numeric(interceptions_by_libyan_coast_guard), 0),
    tcg_monthly = replace_na(as.numeric(interceptions_by_tunisian_coast_guard), 0)
  )

cat(sprintf("  %d months (%s to %s)\n", nrow(monthly),
    min(monthly$ym), max(monthly$ym)))
cat(sprintf("  LCG total: %.0f | TCG total: %.0f\n",
    sum(monthly$lcg_monthly), sum(monthly$tcg_monthly)))

# ── 2. Load Frontex and build daily indicators ───────────────
cat("\n--- 2. Building daily indicators from Frontex ---\n")

frx_raw <- readRDS(file.path(BASE_DIR, "data", "processed",
                              "frontex_incidents.RDS"))

frx <- frx_raw %>%
  filter(country_of_departure %in% CMR_DEPARTURES)

FRX_END <- max(frx$date)
cat(sprintf("  Frontex CMR incidents: %d (%s to %s)\n",
    nrow(frx), min(frx$date), FRX_END))

# Aggregate to daily persons by departure country
frx_daily <- frx %>%
  group_by(date) %>%
  summarise(
    indicator_lcg = sum(num_persons[country_of_departure == "Libya"], na.rm = TRUE),
    indicator_tcg = sum(num_persons[country_of_departure == "Tunisia"], na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("  Daily indicator: %d days with Frontex activity\n", nrow(frx_daily)))
cat(sprintf("  Libya persons total:   %.0f\n", sum(frx_daily$indicator_lcg)))
cat(sprintf("  Tunisia persons total: %.0f\n", sum(frx_daily$indicator_tcg)))

# ── 3. Build date spine and attach indicators ────────────────
cat("\n--- 3. Building date spine ---\n")

# Truncate monthly data to Frontex coverage
# Use the last full month within Frontex coverage
last_full_month <- floor_date(FRX_END, "month")
# If Frontex ends mid-month, include that partial month only if >14 days
days_in_last <- as.integer(FRX_END - last_full_month) + 1
if (days_in_last < 15) {
  panel_end_month <- last_full_month - months(1)
} else {
  panel_end_month <- last_full_month
}
panel_end_date <- min(FRX_END, ceiling_date(panel_end_month, "month") - days(1))

cat(sprintf("  Frontex ends: %s\n", FRX_END))
cat(sprintf("  Last month included: %s (%d days of Frontex data)\n",
    panel_end_month, days_in_last))

spine <- tibble(date = seq(as.Date("2014-01-01"), panel_end_date, by = "day")) %>%
  mutate(ym = floor_date(date, "month")) %>%
  left_join(frx_daily, by = "date") %>%
  replace_na(list(indicator_lcg = 0, indicator_tcg = 0))

cat(sprintf("  Date spine: %d days (%s to %s)\n",
    nrow(spine), min(spine$date), max(spine$date)))

# Attach monthly interceptions to spine
spine <- spine %>%
  left_join(monthly, by = "ym") %>%
  replace_na(list(lcg_monthly = 0, tcg_monthly = 0))

# ── 4. Proportional Denton disaggregation ────────────────────
cat("\n--- 4. Applying proportional Denton disaggregation ---\n")

#' Proportional Denton disaggregation
#'
#' For each month m:
#'   S_m = sum(indicator over days in month m)
#'   If S_m > 0:  daily = monthly * (indicator_day / S_m)
#'   If S_m == 0 and monthly > 0:  daily = monthly / days_in_month  (uniform fallback)
#'   If monthly == 0:  daily = 0

disagg <- spine %>%
  group_by(ym) %>%
  mutate(
    # Monthly sums of indicator within each month
    S_lcg = sum(indicator_lcg),
    S_tcg = sum(indicator_tcg),
    n_days = n(),

    # LCG disaggregation
    lcg_daily = case_when(
      lcg_monthly == 0              ~ 0,
      S_lcg > 0                     ~ lcg_monthly * (indicator_lcg / S_lcg),
      TRUE                          ~ lcg_monthly / n_days  # uniform fallback
    ),
    disagg_method_lcg = case_when(
      lcg_monthly == 0              ~ "zero",
      S_lcg > 0                     ~ "denton",
      TRUE                          ~ "uniform"
    ),

    # TCG disaggregation
    tcg_daily = case_when(
      tcg_monthly == 0              ~ 0,
      S_tcg > 0                     ~ tcg_monthly * (indicator_tcg / S_tcg),
      TRUE                          ~ tcg_monthly / n_days  # uniform fallback
    ),
    disagg_method_tcg = case_when(
      tcg_monthly == 0              ~ "zero",
      S_tcg > 0                     ~ "denton",
      TRUE                          ~ "uniform"
    )
  ) %>%
  ungroup()

# Combine interceptions
disagg <- disagg %>%
  mutate(interceptions_daily = lcg_daily + tcg_daily)

# Report method usage
cat("\n  LCG disaggregation method by month:\n")
lcg_methods <- disagg %>%
  distinct(ym, disagg_method_lcg) %>%
  count(disagg_method_lcg)
print(lcg_methods)

cat("\n  TCG disaggregation method by month:\n")
tcg_methods <- disagg %>%
  distinct(ym, disagg_method_tcg) %>%
  count(disagg_method_tcg)
print(tcg_methods)

# ── 5. Validation tests ─────────────────────────────────────
cat("\n--- 5. Validation tests ---\n")

# Test 1: Temporal consistency — daily sums must match monthly totals
cat("\n  Test 1: Temporal consistency (daily sums == monthly totals)\n")
monthly_check <- disagg %>%
  group_by(ym) %>%
  summarise(
    lcg_daily_sum   = sum(lcg_daily),
    lcg_monthly_val = first(lcg_monthly),
    tcg_daily_sum   = sum(tcg_daily),
    tcg_monthly_val = first(tcg_monthly),
    .groups = "drop"
  ) %>%
  mutate(
    lcg_diff = abs(lcg_daily_sum - lcg_monthly_val),
    tcg_diff = abs(tcg_daily_sum - tcg_monthly_val)
  )

stopifnot(
  "LCG daily sums do not match monthly totals" =
    all(monthly_check$lcg_diff < 0.01),
  "TCG daily sums do not match monthly totals" =
    all(monthly_check$tcg_diff < 0.01)
)
cat("    PASS: All monthly sums match within tolerance (0.01)\n")

# Test 2: Non-negativity
cat("  Test 2: Non-negativity\n")
stopifnot(
  "Negative LCG daily values found" = all(disagg$lcg_daily >= 0),
  "Negative TCG daily values found" = all(disagg$tcg_daily >= 0)
)
cat("    PASS: All daily values >= 0\n")

# Test 3: No NAs in output
cat("  Test 3: No NAs in daily estimates\n")
stopifnot(
  "NA values in lcg_daily" = !any(is.na(disagg$lcg_daily)),
  "NA values in tcg_daily" = !any(is.na(disagg$tcg_daily)),
  "NA values in interceptions_daily" = !any(is.na(disagg$interceptions_daily))
)
cat("    PASS: No NAs in daily interception estimates\n")

# Test 4: Row integrity
cat("  Test 4: Row integrity\n")
expected_days <- as.integer(panel_end_date - as.Date("2014-01-01")) + 1
stopifnot(
  "Row count mismatch" = nrow(disagg) == expected_days,
  "Duplicate dates found" = !any(duplicated(disagg$date))
)
cat(sprintf("    PASS: %d rows, no duplicates\n", nrow(disagg)))

# Test 5: Correlation diagnostic (informational)
cat("  Test 5: Correlation diagnostic\n")
lcg_nonzero <- disagg %>% filter(lcg_daily > 0 & indicator_lcg > 0)
tcg_nonzero <- disagg %>% filter(tcg_daily > 0 & indicator_tcg > 0)
if (nrow(lcg_nonzero) > 10) {
  cat(sprintf("    LCG cor(daily, indicator) = %.3f (N=%d non-zero days)\n",
      cor(lcg_nonzero$lcg_daily, lcg_nonzero$indicator_lcg), nrow(lcg_nonzero)))
}
if (nrow(tcg_nonzero) > 10) {
  cat(sprintf("    TCG cor(daily, indicator) = %.3f (N=%d non-zero days)\n",
      cor(tcg_nonzero$tcg_daily, tcg_nonzero$indicator_tcg), nrow(tcg_nonzero)))
}

# Test 6: Uniform fallback count
cat("  Test 6: Uniform fallback summary\n")
n_uniform_lcg <- sum(lcg_methods$n[lcg_methods$disagg_method_lcg == "uniform"])
n_uniform_tcg <- sum(tcg_methods$n[tcg_methods$disagg_method_tcg == "uniform"])
cat(sprintf("    LCG months using uniform fallback: %d\n",
    ifelse(length(n_uniform_lcg) == 0, 0L, n_uniform_lcg)))
cat(sprintf("    TCG months using uniform fallback: %d\n",
    ifelse(length(n_uniform_tcg) == 0, 0L, n_uniform_tcg)))

# ── 6. Summary statistics ───────────────────────────────────
cat("\n--- 6. Summary ---\n")

annual <- disagg %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    lcg_daily_total = sum(lcg_daily),
    tcg_daily_total = sum(tcg_daily),
    total_interceptions = sum(interceptions_daily),
    days_with_interceptions = sum(interceptions_daily > 0),
    .groups = "drop"
  )
cat("\n  Annual totals (disaggregated daily):\n")
print(annual, n = 15)

# ── 7. Save ─────────────────────────────────────────────────
cat("\n--- 7. Saving ---\n")

output <- disagg %>%
  select(date, lcg_daily, tcg_daily, interceptions_daily,
         disagg_method_lcg, disagg_method_tcg)

out_dir <- file.path(BASE_DIR, "analysis", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(output, file.path(out_dir, "interceptions_daily_disagg.RDS"))
cat(sprintf("Saved: analysis/data/interceptions_daily_disagg.RDS\n"))
cat(sprintf("  %d rows x %d columns\n", nrow(output), ncol(output)))
cat(sprintf("  Date range: %s to %s\n", min(output$date), max(output$date)))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
