# 02_build_daily_panel.R
# ======================
# Build integrated daily panel merging Frontex Themis with IOM MMP.
#
# Crossing attempts follow the literature convention (Deiana et al. 2024,
# Battiston 2022, Rodriguez-Sanchez et al. 2023), extended to three components:
#
#   crossing_attempts = frx_persons              (events engaged by Frontex)
#                     + lcg_tcg_pushbacks        (LCG + TCG pushbacks to Africa,
#                                                 Denton-disaggregated from monthly,
#                                                 rounded to integers in 01)
#                     + n_dead_missing           (dead + missing, IOM MMP)
#
#   Undetected arrivals (UNHCR - Frontex) are excluded from the daily formula
#   because daily-level subtraction inflates the gap 3-27x due to timing
#   mismatches. The annual gap is 1-12% of arrivals. crossing_attempts is
#   therefore a lower bound. See Section 6 comments for details.
#
#   fatality_rate     = n_dead_missing / crossing_attempts
#
# IOM death count — single inclusive series:
#   The panel aggregates all CMR dead+missing to daily frequency using the
#   broadest reasonable filter so the panel serves as the descriptive
#   reference (matches 04_descriptive_statistics.R exactly) and the base
#   for downstream analyses that narrow the sample further:
#     - Route == "Central Mediterranean"
#     - Incident Type in {incident, split incident}  (cumulative / sub-
#       incident excluded because they often double-count the main incident)
#     - Country of Incident in {Algeria, Italy, Libya, Malta, Tunisia}
#     - No cause-of-death restriction
#     - No spatial restriction to the core corridor polygon
#   The single column is `n_dead_missing` (deaths plus missing).
#   Analytical scripts (05_reduced_form_primary.R, 052–056) rebuild a
#   narrower series from raw IOM: incident-only (for robust dates) + optional
#   geographic and cause filters.
#
# Terminology note:
#   - "interception/rescue" (Frontex frx_n_intcp_* columns) refers to
#      who physically engaged with the boat on the EUROPEAN side — these
#      people end up arriving in Europe.
#   - "pushbacks" (lcg_pushbacks / tcg_pushbacks / lcg_tcg_pushbacks) refers
#      to coast guard pullbacks that RETURN people to Libya/Tunisia.
#   These are opposite outcomes and must not be conflated.
#
# LCG/TCG pushbacks are disaggregated from monthly to daily by
# 01_temporal_disaggregation.R using the proportional Denton method
# (Denton 1971) with Frontex daily departures as indicators, then rounded
# to integers per month via largest-remainder (monthly sums preserved exactly).
#
# Input:
#   data/processed/frontex_incidents.RDS
#   data/processed/iom_mmp_incidents.RDS
#   data/processed/era5_swh_daily.RDS
#   data/processed/acled_daily.RDS
#   data/processed/unhcr_daily_arrivals.RDS
#   analysis/data/interceptions_daily_disagg.RDS   (from 01, integer-valued)
#
# Output:
#   analysis/data/daily_panel_complete.RDS

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
CMR_DEPARTURES <- c("Libya", "Tunisia", "Algeria")
CMR_INCIDENT_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia")

# Shared helper for IOM filtering — see analysis/R/_helpers.R
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

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

# Interception/rescue aggregation:
#   frx_n_intcp_* and frx_persons_intcp_* count Frontex events and persons
#   classified by who physically intercepted/rescued the boat on the European
#   side (Italian authorities, NGO, EU, Commercial, Coast Guard, Land patrol,
#   etc.) combined with whether it was a SAR event. These are events where
#   migrants ended up arriving in Europe. Do NOT confuse with
#   lcg_pushbacks / tcg_pushbacks, which are LCG/TCG pushbacks to Africa.
#   The 9 categories come from the Frontex `event_type` field.
frx_daily <- frx %>%
  group_by(date) %>%
  summarise(
    frx_incidents             = n(),
    frx_persons               = sum(num_persons, na.rm = TRUE),
    frx_n_sar                 = sum(sar_flag, na.rm = TRUE),
    # Interception/rescue event counts (SAR + Not SAR x interceptor)
    frx_n_intcp_sar_ngo           = sum(event_type == "SAR: NGO"),
    frx_n_intcp_sar_ita           = sum(event_type == "SAR: Italian authorities"),
    frx_n_intcp_sar_eu            = sum(event_type == "SAR: EU operations (IRINI)"),
    frx_n_intcp_sar_commercial    = sum(event_type == "SAR: Commercial vessels"),
    frx_n_intcp_sar_other         = sum(event_type == "SAR: Other"),
    frx_n_intcp_notsar_cg         = sum(event_type == "Not SAR: Coast Guard"),
    frx_n_intcp_notsar_land       = sum(event_type == "Not SAR: Land patrol"),
    frx_n_intcp_notsar_self       = sum(event_type == "Not SAR: Self-arrived"),
    frx_n_intcp_notsar_other      = sum(event_type == "Not SAR: Other"),
    # Interception/rescue persons (same categories, persons instead of events)
    frx_persons_intcp_sar_ngo         = sum(num_persons[event_type == "SAR: NGO"], na.rm = TRUE),
    frx_persons_intcp_sar_ita         = sum(num_persons[event_type == "SAR: Italian authorities"], na.rm = TRUE),
    frx_persons_intcp_sar_eu          = sum(num_persons[event_type == "SAR: EU operations (IRINI)"], na.rm = TRUE),
    frx_persons_intcp_sar_commercial  = sum(num_persons[event_type == "SAR: Commercial vessels"], na.rm = TRUE),
    frx_persons_intcp_sar_other       = sum(num_persons[event_type == "SAR: Other"], na.rm = TRUE),
    frx_persons_intcp_notsar_cg       = sum(num_persons[event_type == "Not SAR: Coast Guard"], na.rm = TRUE),
    frx_persons_intcp_notsar_land     = sum(num_persons[event_type == "Not SAR: Land patrol"], na.rm = TRUE),
    frx_persons_intcp_notsar_self     = sum(num_persons[event_type == "Not SAR: Self-arrived"], na.rm = TRUE),
    frx_persons_intcp_notsar_other    = sum(num_persons[event_type == "Not SAR: Other"], na.rm = TRUE),
    # Detection method (detected_by can have compound values like "CPB;CPV")
    frx_n_det_fwa             = sum(grepl("FWA", detected_by), na.rm = TRUE),
    frx_n_det_helo            = sum(grepl("HELO", detected_by), na.rm = TRUE),
    frx_n_det_ngo             = sum(grepl("NGO vessel", detected_by), na.rm = TRUE),
    frx_n_det_call            = sum(grepl("Call-Migrant", detected_by), na.rm = TRUE),
    frx_n_det_patrol          = sum(grepl("CPV|CPB|Land Patrol", detected_by), na.rm = TRUE),
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
    frx_sar_share                 = frx_n_sar / frx_incidents,
    frx_intcp_sar_ngo_share        = frx_n_intcp_sar_ngo / frx_incidents,
    frx_intcp_sar_ita_share        = frx_n_intcp_sar_ita / frx_incidents,
    frx_intcp_sar_eu_share         = frx_n_intcp_sar_eu / frx_incidents,
    frx_intcp_sar_commercial_share = frx_n_intcp_sar_commercial / frx_incidents,
    frx_det_fwa_share             = frx_n_det_fwa / frx_incidents,
    frx_det_helo_share            = frx_n_det_helo / frx_incidents,
    frx_det_ngo_share             = frx_n_det_ngo / frx_incidents,
    frx_det_call_share            = frx_n_det_call / frx_incidents,
    frx_det_patrol_share          = frx_n_det_patrol / frx_incidents,
    frx_inflatable_share          = frx_n_inflatable / frx_incidents,
    frx_frac_inflatable_persons   = ifelse(frx_persons > 0,
                                            frx_persons_inflatable / frx_persons, NA_real_),
    frx_wooden_share              = frx_n_wooden / frx_incidents,
    frx_libya_share               = frx_dep_libya / frx_incidents,
    frx_in_oparea_share           = frx_n_in_oparea / frx_incidents
  )

cat(sprintf("  Frontex daily: %d days with activity\n", nrow(frx_daily)))

# ── 3. Build daily IOM aggregate ──────────────────────────
cat("\n--- 3. Building daily IOM aggregate (broad/descriptive filter) ---\n")

# Inclusive filter (see header): CMR route + incident/split + CMR countries,
# no cause or spatial restriction. Matches 04_descriptive_statistics.R.
# Filter is centralised in analysis/R/_helpers.R::build_iom_daily(); change
# parameters there or in the call below to swap variants for sensitivity.
iom_daily <- build_iom_daily(
  incident_types = c("incident", "split incident"),
  spatial        = "all_cmr",
  causes         = "all",
  countries      = CMR_INCIDENT_COUNTRIES,
  base_dir       = BASE_DIR
)

cat(sprintf("  IOM daily: %d days, %.0f dead+missing\n",
            nrow(iom_daily), sum(iom_daily$n_dead_missing)))

# ── 5. Build date spine and merge ─────────────────────────
cat("\n--- 5. Building integrated panel ---\n")

# Date spine: 2014-01-01 to the last date with COMPLETE coverage of all sources.
# That is bounded by interceptions_daily_disagg.RDS (built by 01), which
# truncates at the last full month inside the Frontex window — currently
# 2023-05-31. Reading the cap from the disagg file rather than hard-coding it
# means the panel adapts automatically if 01's logic changes.
disagg_path <- file.path(BASE_DIR, "analysis", "data", "interceptions_daily_disagg.RDS")
if (!file.exists(disagg_path)) {
  stop("interceptions_daily_disagg.RDS not found — run 01 first.")
}
PANEL_END <- max(readRDS(disagg_path)$date)
cat(sprintf("  Panel end (from disagg file): %s\n", PANEL_END))

spine <- tibble(date = seq(as.Date("2014-01-01"), PANEL_END, by = "day"))

# Weather from ERA5 clean data
weather <- readRDS(file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS")) %>%
  select(date, swh, swh_prev3days, swh_prevweek,
         swh_w, swh_w_prev3days, swh_w_prevweek)

# ACLED weekly conflict
acled <- readRDS(file.path(BASE_DIR, "data", "processed", "acled_daily.RDS")) %>%
  select(date, week_date,
         libya_conflict, libya_conflict_fatalities,
         libya_battles, libya_expvio, libya_violciv,
         tunisia_conflict, tunisia_conflict_fatalities,
         tunisia_battles, tunisia_expvio, tunisia_violciv)
cat(sprintf("  ACLED: %d days loaded\n", nrow(acled)))

# Daily LCG/TCG pushbacks (from 01_temporal_disaggregation.R, Denton method).
# disagg_path/PANEL_END were already verified and read above (Section 5).
lcg_tcg_daily <- readRDS(disagg_path) %>%
  select(date, lcg_pushbacks, tcg_pushbacks, lcg_tcg_pushbacks)
cat(sprintf("  LCG/TCG pushbacks (Denton disagg): %d days loaded\n", nrow(lcg_tcg_daily)))

# UNHCR daily arrivals (for arrivals not detected during Frontex operations)
unhcr <- readRDS(file.path(BASE_DIR, "data", "processed",
                            "unhcr_daily_arrivals.RDS"))
cat(sprintf("  UNHCR arrivals: %d days (%s to %s)\n",
    nrow(unhcr), min(unhcr$date), max(unhcr$date)))

# Merge everything
panel <- spine %>%
  mutate(iso_week = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))) %>%
  left_join(frx_daily,     by = "date") %>%
  left_join(iom_daily,     by = "date") %>%
  left_join(weather,       by = "date") %>%
  left_join(acled,         by = "date") %>%
  left_join(lcg_tcg_daily, by = "date") %>%
  left_join(unhcr,         by = "date")

# Replace NA counts with 0 (Frontex and IOM counts)
count_cols <- c("frx_incidents", "frx_persons", "frx_n_sar",
                "frx_n_intcp_sar_ngo", "frx_n_intcp_sar_ita", "frx_n_intcp_sar_eu",
                "frx_n_intcp_sar_commercial", "frx_n_intcp_sar_other",
                "frx_n_intcp_notsar_cg", "frx_n_intcp_notsar_land",
                "frx_n_intcp_notsar_self", "frx_n_intcp_notsar_other",
                "frx_persons_intcp_sar_ngo", "frx_persons_intcp_sar_ita",
                "frx_persons_intcp_sar_eu", "frx_persons_intcp_sar_commercial",
                "frx_persons_intcp_sar_other", "frx_persons_intcp_notsar_cg",
                "frx_persons_intcp_notsar_land", "frx_persons_intcp_notsar_self",
                "frx_persons_intcp_notsar_other",
                "frx_n_det_fwa", "frx_n_det_helo", "frx_n_det_ngo",
                "frx_n_det_call", "frx_n_det_patrol",
                "frx_n_inflatable", "frx_persons_inflatable",
                "frx_n_wooden", "frx_persons_wooden",
                "frx_dep_libya", "frx_dep_tunisia", "frx_dep_algeria",
                "frx_n_in_oparea",
                "n_dead_missing",
                "lcg_pushbacks", "tcg_pushbacks", "lcg_tcg_pushbacks")
panel <- panel %>%
  mutate(across(all_of(count_cols), ~ replace_na(.x, 0L)))

cat(sprintf("  Panel: %d days (%s to %s)\n",
    nrow(panel), min(panel$date), max(panel$date)))

# ── 6. Derive variables ───────────────────────────────────
cat("\n--- 6. Deriving variables ---\n")

# crossing_attempts components (daily):
#   1. frx_persons         — persons in Frontex events (arrivals to Europe)
#   2. lcg_tcg_pushbacks   — LCG + TCG pushbacks to Africa (Denton-disaggregated)
#   3. n_dead_missing      — deaths and missing (IOM MMP, CMR countries)
#
# NOTE on undetected arrivals:
#   UNHCR daily arrivals to Italy exceed Frontex daily events by 1-12%
#   annually, reflecting arrivals not recorded during Frontex operations.
#   However, computing undetected_arrivals = max(UNHCR - Frontex, 0) at the
#   DAILY level inflates this gap 3-27x because of timing mismatches:
#   UNHCR records arrivals on landing day, Frontex records events on
#   event day. A rescue on day t may produce a Frontex count on day t
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
    crossing_attempts = frx_persons + lcg_tcg_pushbacks + n_dead_missing,
    fatality_rate     = ifelse(crossing_attempts > 0,
                                n_dead_missing / crossing_attempts, NA_real_),
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
    n_dead_missing    = sum(n_dead_missing),
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
  filter(n_dead_missing > frx_persons, n_dead_missing > 0)
cat(sprintf("\n  Days where dead+missing > Frontex persons: %d\n", nrow(excess_days)))

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
