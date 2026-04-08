# 08_build_frontex_iom_panel.R
# ============================
# Build integrated daily panel merging Frontex Themis with IOM MMP.
#
# Crossing attempts follow the literature convention (Deiana et al. 2024,
# Battiston 2022, Rodriguez-Sanchez et al. 2023), extended to four components:
#
#   crossing_attempts = frontex_persons          (detected during operations)
#                     + interceptions_daily      (LCG + TCG, Denton-disaggregated)
#                     + iom_deaths               (dead + missing, IOM MMP)
#
#   Undetected arrivals (UNHCR - Frontex) are excluded from the daily formula
#   because daily-level subtraction inflates the gap 3-27x due to timing
#   mismatches. The annual gap is 1-12% of arrivals. crossing_attempts is
#   therefore a lower bound. See Section 6 comments for details.
#
#   fatality_rate     = iom_deaths / crossing_attempts
#
# LCG/TCG interceptions are disaggregated from monthly to daily by
# 07b_temporal_disaggregation.R using the proportional Denton method
# (Denton 1971) with Frontex daily departures as indicators.
#
# Input:
#   data/processed/frontex_incidents.RDS
#   data/processed/iom_mmp_incidents.RDS
#   data/processed/era5_swh_daily.RDS
#   data/processed/acled_daily.RDS
#   data/processed/unhcr_daily_arrivals.RDS
#   data/processed/archive/sar_ngo_ops_daily_RS.RDS
#   analysis/data/interceptions_daily_disagg.RDS   (from 07b)
#
# Output:
#   analysis/data/daily_panel_complete.RDS

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
SEA_CAUSES <- c("Drowning", "Mixed or unknown")
CMR_INCIDENT_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia")
CMR_DEPARTURES <- c("Libya", "Tunisia", "Algeria")

cat("============================================================\n")
cat("BUILD INTEGRATED FRONTEX + IOM DAILY PANEL\n")
cat("============================================================\n\n")

# ── 1. Load and filter Frontex Themis ─────────────────────
cat("--- 1. Loading Frontex Themis ---\n")

frx_all <- readRDS(file.path(BASE_DIR, "data", "processed", "frontex_incidents.RDS"))

cat(sprintf("  Total Frontex: %d incidents (%s to %s)\n",
    nrow(frx_all), min(frx_all$date), max(frx_all$date)))

# Filter to CMR departures
frx <- frx_all %>%
  filter(country_of_departure %in% CMR_DEPARTURES)

FRX_END <- max(frx$date)
cat(sprintf("  After CMR filter (Libya/Tunisia/Algeria): %d incidents\n", nrow(frx)))
cat(sprintf("  Excluded: %d non-CMR incidents\n", nrow(frx_all) - nrow(frx)))
cat(sprintf("  Frontex end date: %s\n", FRX_END))

cat("  Boat type distribution:\n")
print(table(frx$boat_category))

# ── 2. Aggregate Frontex to daily ─────────────────────────
cat("\n--- 2. Aggregating Frontex to daily ---\n")

frx_daily <- frx %>%
  group_by(date) %>%
  summarise(
    frx_incidents             = n(),
    frx_persons               = sum(num_persons, na.rm = TRUE),
    frx_n_sar                 = sum(sar_flag, na.rm = TRUE),
    # SAR by actor type
    frx_n_sar_ngo             = sum(event_type == "SAR: NGO"),
    frx_n_sar_ita             = sum(event_type == "SAR: Italian authorities"),
    frx_n_sar_eu              = sum(event_type == "SAR: EU operations (IRINI)"),
    frx_n_sar_commercial      = sum(event_type == "SAR: Commercial vessels"),
    frx_n_sar_other           = sum(event_type == "SAR: Other"),
    # Detection method (detected_by can have compound values like "CPB;CPV")
    frx_n_det_fwa             = sum(grepl("FWA", detected_by), na.rm = TRUE),
    frx_n_det_helo            = sum(grepl("HELO", detected_by), na.rm = TRUE),
    frx_n_det_ngo             = sum(grepl("NGO vessel", detected_by), na.rm = TRUE),
    frx_n_det_call            = sum(grepl("Call-Migrant", detected_by), na.rm = TRUE),
    frx_n_det_patrol          = sum(grepl("CPV|CPB|Land Patrol", detected_by), na.rm = TRUE),
    # Non-SAR categories
    frx_n_notsar_cg           = sum(event_type == "Not SAR: Coast Guard"),
    frx_n_notsar_land         = sum(event_type == "Not SAR: Land patrol"),
    frx_n_selfarrived         = sum(event_type == "Not SAR: Self-arrived"),
    frx_n_notsar_other        = sum(event_type == "Not SAR: Other"),
    # Boat type
    frx_n_inflatable          = sum(boat_category == "Inflatable"),
    frx_persons_inflatable    = sum(num_persons[boat_category == "Inflatable"], na.rm = TRUE),
    frx_n_wooden              = sum(boat_category == "Wooden"),
    frx_persons_wooden        = sum(num_persons[boat_category == "Wooden"], na.rm = TRUE),
    frx_dep_libya             = sum(country_of_departure == "Libya"),
    frx_dep_tunisia           = sum(country_of_departure == "Tunisia"),
    frx_dep_algeria           = sum(country_of_departure == "Algeria"),
    frx_n_in_oparea           = sum(in_op_area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    frx_sar_share             = frx_n_sar / frx_incidents,
    frx_sar_ngo_share         = frx_n_sar_ngo / frx_incidents,
    frx_sar_ita_share         = frx_n_sar_ita / frx_incidents,
    frx_sar_eu_share          = frx_n_sar_eu / frx_incidents,
    frx_sar_commercial_share  = frx_n_sar_commercial / frx_incidents,
    frx_det_fwa_share         = frx_n_det_fwa / frx_incidents,
    frx_det_helo_share        = frx_n_det_helo / frx_incidents,
    frx_det_ngo_share         = frx_n_det_ngo / frx_incidents,
    frx_det_call_share        = frx_n_det_call / frx_incidents,
    frx_det_patrol_share      = frx_n_det_patrol / frx_incidents,
    frx_inflatable_share      = frx_n_inflatable / frx_incidents,
    frx_frac_inflatable_persons = ifelse(frx_persons > 0,
                                          frx_persons_inflatable / frx_persons, NA_real_),
    frx_wooden_share          = frx_n_wooden / frx_incidents,
    frx_libya_share           = frx_dep_libya / frx_incidents,
    frx_in_oparea_share       = frx_n_in_oparea / frx_incidents
  )

cat(sprintf("  Frontex daily: %d days with activity\n", nrow(frx_daily)))

# ── 3. Load and filter IOM MMP ────────────────────────────
cat("\n--- 3. Loading IOM MMP ---\n")

iom_raw <- readRDS(file.path(BASE_DIR, "data", "processed", "iom_mmp_incidents.RDS"))

# CMR sea deaths (drowning + mixed/unknown, filtered by country of incident)
iom_core <- iom_raw %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) %in% c("incident", "split incident"),
         `Cause of death (category)` %in% SEA_CAUSES,
         `Country of Incident` %in% CMR_INCIDENT_COUNTRIES) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date))

cat(sprintf("  IOM CMR deaths: %d incidents\n", nrow(iom_core)))

# ── 4. Aggregate IOM to daily ─────────────────────────────
cat("\n--- 4. Aggregating IOM to daily ---\n")

iom_daily <- iom_core %>%
  group_by(date) %>%
  summarise(iom_deaths    = sum(dead_missing),
            iom_incidents = n(), .groups = "drop")

cat(sprintf("  IOM daily: %d days with deaths\n", nrow(iom_daily)))

# ── 5. Build date spine and merge ─────────────────────────
cat("\n--- 5. Building integrated panel ---\n")

# Date spine: 2014-01-01 to Frontex end date
spine <- tibble(date = seq(as.Date("2014-01-01"), FRX_END, by = "day"))

# Weather from ERA5 clean data
weather <- readRDS(file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS")) %>%
  select(date, swh, swh_prev3days, swh_prevweek)

# ACLED weekly conflict
acled <- readRDS(file.path(BASE_DIR, "data", "processed", "acled_daily.RDS")) %>%
  select(date, week_date,
         libya_conflict, libya_conflict_fatalities,
         libya_battles, libya_expvio, libya_violciv,
         tunisia_conflict, tunisia_conflict_fatalities,
         tunisia_battles, tunisia_expvio, tunisia_violciv)
cat(sprintf("  ACLED: %d days loaded\n", nrow(acled)))

# Daily LCG/TCG interceptions (from 07b_temporal_disaggregation.R, Denton method)
disagg_path <- file.path(BASE_DIR, "analysis", "data", "interceptions_daily_disagg.RDS")
if (file.exists(disagg_path)) {
  interceptions_daily <- readRDS(disagg_path) %>%
    select(date, lcg_daily, tcg_daily, interceptions_daily)
  cat(sprintf("  Interceptions (Denton disagg): %d days loaded\n", nrow(interceptions_daily)))
} else {
  interceptions_daily <- tibble(date = as.Date(character()),
                                 lcg_daily = numeric(), tcg_daily = numeric(),
                                 interceptions_daily = numeric())
  cat("  WARNING: interceptions_daily_disagg.RDS not found — run 07b first\n")
}

# UNHCR daily arrivals (for arrivals not detected during Frontex operations)
unhcr <- readRDS(file.path(BASE_DIR, "data", "processed",
                            "unhcr_daily_arrivals.RDS"))
cat(sprintf("  UNHCR arrivals: %d days (%s to %s)\n",
    nrow(unhcr), min(unhcr$date), max(unhcr$date)))

# Merge everything
panel <- spine %>%
  mutate(iso_week = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))) %>%
  left_join(frx_daily, by = "date") %>%
  left_join(iom_daily, by = "date") %>%
  left_join(weather, by = "date") %>%
  left_join(acled, by = "date") %>%
  left_join(interceptions_daily, by = "date") %>%
  left_join(unhcr, by = "date")

# Replace NA counts with 0 (Frontex and IOM counts)
count_cols <- c("frx_incidents", "frx_persons", "frx_n_sar",
                "frx_n_sar_ngo", "frx_n_sar_ita", "frx_n_sar_eu",
                "frx_n_sar_commercial", "frx_n_sar_other",
                "frx_n_det_fwa", "frx_n_det_helo", "frx_n_det_ngo",
                "frx_n_det_call", "frx_n_det_patrol",
                "frx_n_notsar_cg", "frx_n_notsar_land", "frx_n_selfarrived",
                "frx_n_notsar_other",
                "frx_n_inflatable", "frx_persons_inflatable",
                "frx_n_wooden", "frx_persons_wooden",
                "frx_dep_libya", "frx_dep_tunisia", "frx_dep_algeria",
                "frx_n_in_oparea",
                "iom_deaths", "iom_incidents",
                "lcg_daily", "tcg_daily", "interceptions_daily")
panel <- panel %>%
  mutate(across(all_of(count_cols), ~ replace_na(.x, 0L)))

cat(sprintf("  Panel: %d days (%s to %s)\n",
    nrow(panel), min(panel$date), max(panel$date)))

# ── 6. Derive variables ───────────────────────────────────
cat("\n--- 6. Deriving variables ---\n")

# crossing_attempts components (daily):
#   1. frx_persons            — persons detected during Frontex operations
#   2. interceptions_daily    — LCG + TCG interceptions (Denton-disaggregated from monthly)
#   3. iom_deaths             — deaths and missing (IOM MMP, CMR countries)
#
# NOTE on undetected arrivals:
#   UNHCR daily arrivals to Italy exceed Frontex daily detections by 1-12%
#   annually, reflecting arrivals not recorded during Frontex operations.
#   However, computing undetected_arrivals = max(UNHCR - Frontex, 0) at the
#   DAILY level inflates this gap 3-27x because of timing mismatches:
#   UNHCR records arrivals on landing day, Frontex records detections on
#   detection day. A rescue on day t may produce a Frontex count on day t
#   and a UNHCR count on day t+1. The daily pmax captures all "UNHCR-high"
#   days but discards offsetting "Frontex-high" days.
#
#   Example (2017): annual gap = 3,040 persons, but daily pmax sum = 83,435.
#
#   We therefore exclude undetected arrivals from the daily crossing formula.
#   The gap is small at the monthly/annual level (see validation below) and
#   crossing_attempts should be understood as a LOWER BOUND.
#   The UNHCR daily arrivals column is retained in the panel for reference.
panel <- panel %>%
  mutate(
    crossing_attempts = frx_persons + interceptions_daily + iom_deaths,
    fatality_rate     = ifelse(crossing_attempts > 0,
                                iom_deaths / crossing_attempts, NA_real_),
    post_mou          = as.integer(date >= MOU_DATE),
    year              = year(date),
    month_year        = factor(format(date, "%Y-%m"))
  )

# ── 7. Validation diagnostics ─────────────────────────────
cat("\n--- 7. Validation ---\n")

cat("\n  Annual totals:\n")
annual <- panel %>%
  group_by(year) %>%
  summarise(
    frx_persons       = sum(frx_persons),
    iom_deaths        = sum(iom_deaths),
    crossing_attempts = sum(crossing_attempts),
    frx_incidents     = sum(frx_incidents),
    days_with_frx     = sum(frx_incidents > 0),
    .groups = "drop"
  )
print(annual, n = 12)

# Undetected arrivals diagnostic (monthly level, where timing washes out)
cat("\n  Undetected arrivals (monthly UNHCR - Frontex, for reference):\n")
monthly_gap <- panel %>%
  filter(!is.na(arrivals)) %>%
  mutate(ym = floor_date(date, "month")) %>%
  group_by(ym) %>%
  summarise(unhcr = sum(arrivals), frx = sum(frx_persons), .groups = "drop") %>%
  mutate(gap = pmax(unhcr - frx, 0))
cat(sprintf("    Monthly gap total: %.0f persons (%.1f%% of UNHCR)\n",
    sum(monthly_gap$gap),
    sum(monthly_gap$gap) / sum(monthly_gap$unhcr) * 100))

excess_days <- panel %>%
  filter(iom_deaths > frx_persons, iom_deaths > 0)
cat(sprintf("\n  Days where IOM deaths > Frontex persons: %d\n", nrow(excess_days)))

rate_days <- panel %>% filter(!is.na(fatality_rate), crossing_attempts > 0)
cat(sprintf("\n  Fatality rate (days with crossings, N=%d):\n", nrow(rate_days)))
cat(sprintf("    Mean:   %.4f (%.2f%%)\n", mean(rate_days$fatality_rate),
    mean(rate_days$fatality_rate) * 100))
cat(sprintf("    Pre-MoU mean:  %.4f (%.2f%%)\n",
    mean(rate_days$fatality_rate[rate_days$post_mou == 0]),
    mean(rate_days$fatality_rate[rate_days$post_mou == 0]) * 100))
cat(sprintf("    Post-MoU mean: %.4f (%.2f%%)\n",
    mean(rate_days$fatality_rate[rate_days$post_mou == 1]),
    mean(rate_days$fatality_rate[rate_days$post_mou == 1]) * 100))

# ── 8. Save ───────────────────────────────────────────────
cat("\n--- 8. Saving ---\n")

out_dir <- file.path(BASE_DIR, "analysis", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(panel, file.path(out_dir, "daily_panel_complete.RDS"))
cat(sprintf("Saved: analysis/data/daily_panel_complete.RDS\n"))
cat(sprintf("  %d rows x %d columns\n", nrow(panel), ncol(panel)))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
