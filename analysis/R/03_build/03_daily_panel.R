# ── Daily panel: Frontex + IOM + UNHCR + LCG/TCG + ERA5 + ACLED ────────────
# crossing_attempts = frx_persons + lcg_tcg_pushbacks + n_dead_missing.
# Daily spine 2014-01-01 to last day of interceptions_daily_disagg.RDS.
# Broad IOM filter here (incident + split, no cause/spatial cut) is the
# volume denominator; analytical scripts narrow via build_iom_daily().

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()
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
frx <- frx_all |>
  filter(country_of_departure %in% CMR_DEPARTURES)

FRX_END <- max(frx$date)
cat(sprintf("  After CMR filter (Libya/Tunisia/Algeria): %d incidents\n", nrow(frx)))
cat(sprintf("  Excluded: %d non-CMR incidents\n", nrow(frx_all) - nrow(frx)))
cat(sprintf("  Frontex end date: %s\n", FRX_END))

cat("  Boat type distribution:\n")
print(table(frx$boat_category))

# ── 2. Aggregate Frontex to daily ─────────────────────────
cat("\n--- 2. Aggregating Frontex to daily ---\n")

# Interception aggregation — full matrix:
#   frx_n_{sar|notsar}_{interceptor} and frx_persons_{sar|notsar}_{interceptor}
#   count Frontex events and persons by interceptor_type × SAR status.
#   SAR bucket: sar_ops == TRUE.
#   Not-SAR bucket: sar_ops == FALSE OR sar_ops == NA (raw SAR = "N/A" —
#     mostly pre-2020 land-patrol events; grouped with Not-SAR per design).
#   Interceptor types (9): NGO, EU_ops, Ita_ops, Commercial, EU_Coast_Guard,
#                          Land_patrol, No_intercept, Other, NA.
#   Do NOT confuse with lcg_pushbacks / tcg_pushbacks (LCG/TCG pushbacks
#   to Africa, separate data source).
frx_daily <- frx |>
  mutate(
    sar_bucket = ifelse(!is.na(sar_ops) & sar_ops, "sar", "notsar"),
    int_lab = dplyr::recode(interceptor_type,
      "NGO" = "ngo", "EU_ops" = "eu", "Ita_ops" = "ita",
      "Commercial" = "comm", "EU_Coast_Guard" = "cg",
      "Land_patrol" = "land", "No_intercept" = "noint",
      "Other" = "other", "NA" = "na")
  ) |>
  group_by(date) |>
  summarise(
    frx_incidents             = n(),
    frx_persons               = sum(num_persons, na.rm = TRUE),
    frx_n_sar                 = sum(sar_bucket == "sar"),

    # Counts: interceptor × SAR matrix (9 × 2 = 18 cells)
    frx_n_sar_ngo        = sum(sar_bucket == "sar"    & int_lab == "ngo"),
    frx_n_sar_eu         = sum(sar_bucket == "sar"    & int_lab == "eu"),
    frx_n_sar_ita        = sum(sar_bucket == "sar"    & int_lab == "ita"),
    frx_n_sar_comm       = sum(sar_bucket == "sar"    & int_lab == "comm"),
    frx_n_sar_cg         = sum(sar_bucket == "sar"    & int_lab == "cg"),
    frx_n_sar_land       = sum(sar_bucket == "sar"    & int_lab == "land"),
    frx_n_sar_noint      = sum(sar_bucket == "sar"    & int_lab == "noint"),
    frx_n_sar_other      = sum(sar_bucket == "sar"    & int_lab == "other"),
    frx_n_sar_na         = sum(sar_bucket == "sar"    & int_lab == "na"),
    frx_n_notsar_ngo     = sum(sar_bucket == "notsar" & int_lab == "ngo"),
    frx_n_notsar_eu      = sum(sar_bucket == "notsar" & int_lab == "eu"),
    frx_n_notsar_ita     = sum(sar_bucket == "notsar" & int_lab == "ita"),
    frx_n_notsar_comm    = sum(sar_bucket == "notsar" & int_lab == "comm"),
    frx_n_notsar_cg      = sum(sar_bucket == "notsar" & int_lab == "cg"),
    frx_n_notsar_land    = sum(sar_bucket == "notsar" & int_lab == "land"),
    frx_n_notsar_noint   = sum(sar_bucket == "notsar" & int_lab == "noint"),
    frx_n_notsar_other   = sum(sar_bucket == "notsar" & int_lab == "other"),
    frx_n_notsar_na      = sum(sar_bucket == "notsar" & int_lab == "na"),

    # Persons: same matrix
    frx_persons_sar_ngo      = sum(num_persons[sar_bucket == "sar"    & int_lab == "ngo"],   na.rm = TRUE),
    frx_persons_sar_eu       = sum(num_persons[sar_bucket == "sar"    & int_lab == "eu"],    na.rm = TRUE),
    frx_persons_sar_ita      = sum(num_persons[sar_bucket == "sar"    & int_lab == "ita"],   na.rm = TRUE),
    frx_persons_sar_comm     = sum(num_persons[sar_bucket == "sar"    & int_lab == "comm"],  na.rm = TRUE),
    frx_persons_sar_cg       = sum(num_persons[sar_bucket == "sar"    & int_lab == "cg"],    na.rm = TRUE),
    frx_persons_sar_land     = sum(num_persons[sar_bucket == "sar"    & int_lab == "land"],  na.rm = TRUE),
    frx_persons_sar_noint    = sum(num_persons[sar_bucket == "sar"    & int_lab == "noint"], na.rm = TRUE),
    frx_persons_sar_other    = sum(num_persons[sar_bucket == "sar"    & int_lab == "other"], na.rm = TRUE),
    frx_persons_sar_na       = sum(num_persons[sar_bucket == "sar"    & int_lab == "na"],    na.rm = TRUE),
    frx_persons_notsar_ngo   = sum(num_persons[sar_bucket == "notsar" & int_lab == "ngo"],   na.rm = TRUE),
    frx_persons_notsar_eu    = sum(num_persons[sar_bucket == "notsar" & int_lab == "eu"],    na.rm = TRUE),
    frx_persons_notsar_ita   = sum(num_persons[sar_bucket == "notsar" & int_lab == "ita"],   na.rm = TRUE),
    frx_persons_notsar_comm  = sum(num_persons[sar_bucket == "notsar" & int_lab == "comm"],  na.rm = TRUE),
    frx_persons_notsar_cg    = sum(num_persons[sar_bucket == "notsar" & int_lab == "cg"],    na.rm = TRUE),
    frx_persons_notsar_land  = sum(num_persons[sar_bucket == "notsar" & int_lab == "land"],  na.rm = TRUE),
    frx_persons_notsar_noint = sum(num_persons[sar_bucket == "notsar" & int_lab == "noint"], na.rm = TRUE),
    frx_persons_notsar_other = sum(num_persons[sar_bucket == "notsar" & int_lab == "other"], na.rm = TRUE),
    frx_persons_notsar_na    = sum(num_persons[sar_bucket == "notsar" & int_lab == "na"],    na.rm = TRUE),

    # Detection method (grepl on raw detected_by — preserves multi-method events)
    frx_n_det_fwa             = sum(grepl("FWA",  detected_by), na.rm = TRUE),
    frx_n_det_helo            = sum(grepl("HELO", detected_by), na.rm = TRUE),
    frx_n_det_rpas            = sum(grepl("RPAS", detected_by), na.rm = TRUE),
    frx_n_det_mas             = sum(grepl("MAS",  detected_by), na.rm = TRUE),
    frx_n_det_aerial          = sum(grepl("FWA|HELO|RPAS|MAS", detected_by), na.rm = TRUE),
    frx_n_det_ngo             = sum(grepl("NGO vessel", detected_by), na.rm = TRUE),
    frx_n_det_call            = sum(grepl("Call-", detected_by), na.rm = TRUE),
    frx_n_det_cg              = sum(grepl("CPV|CPB|OPV", detected_by), na.rm = TRUE),
    frx_n_det_land            = sum(grepl("Land Patrol", detected_by), na.rm = TRUE),

    # Boat type
    frx_n_inflatable          = sum(boat_category == "Inflatable"),
    frx_persons_inflatable    = sum(num_persons[boat_category == "Inflatable"], na.rm = TRUE),
    frx_n_wooden              = sum(boat_category == "Wooden"),
    frx_persons_wooden        = sum(num_persons[boat_category == "Wooden"], na.rm = TRUE),
    frx_n_fibreglass          = sum(boat_category == "Fibre glass"),
    frx_persons_fibreglass    = sum(num_persons[boat_category == "Fibre glass"], na.rm = TRUE),

    # Departure country
    frx_dep_libya             = sum(country_of_departure == "Libya"),
    frx_dep_tunisia           = sum(country_of_departure == "Tunisia"),
    frx_dep_algeria           = sum(country_of_departure == "Algeria"),

    # Operational area
    frx_n_in_oparea           = sum(in_op_area, na.rm = TRUE),

    # Multi-actor flag (any event with compound actors of multiple categories)
    frx_n_multi_actors        = sum(multi_actors_inv, na.rm = TRUE),

    .groups = "drop"
  ) |>
  mutate(
    frx_sar_share                 = frx_n_sar / frx_incidents,
    frx_det_aerial_share          = frx_n_det_aerial / frx_incidents,
    frx_det_ngo_share             = frx_n_det_ngo / frx_incidents,
    frx_det_call_share            = frx_n_det_call / frx_incidents,
    frx_det_cg_share              = frx_n_det_cg / frx_incidents,
    frx_det_land_share            = frx_n_det_land / frx_incidents,
    frx_inflatable_share          = frx_n_inflatable / frx_incidents,
    frx_wooden_share              = frx_n_wooden / frx_incidents,
    frx_fibreglass_share          = frx_n_fibreglass / frx_incidents,
    frx_frac_inflatable_persons   = ifelse(frx_persons > 0,
                                            frx_persons_inflatable / frx_persons, NA_real_),
    frx_libya_share               = frx_dep_libya / frx_incidents,
    frx_in_oparea_share           = frx_n_in_oparea / frx_incidents
  )

cat(sprintf("  Frontex daily: %d days with activity\n", nrow(frx_daily)))

# ── 3. Build daily IOM aggregate ──────────────────────────
cat("\n--- 3. Building daily IOM aggregate (broad/descriptive filter) ---\n")

# Inclusive filter (see header): CMR route + incident/split + CMR countries,
# no cause or spatial restriction. Matches 04_descriptive_statistics.R.
# NOTE: this build deliberately KEEPS split incidents (incident_types
# passed explicitly below) because n_dead_missing here feeds the broad
# crossing_attempts volume LOWER BOUND, not the IOM death OUTCOME. The
# analysis death series uses build_iom_daily()'s default, which now
# EXCLUDES split incidents (date-smearing flattens the SWH gradient and
# breaks IOM/UNITED comparability — see _helpers.R). These are different
# variables; the split choice is intentionally opposite for each.
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

# Weather from ERA5 clean data (all 16 SWH columns: 1d / 3d / 5d / 7d
# windows, mean + max, uniform + death-weighted). See
# data/scripts/01_clean_era5_swh.R for the column definitions.
weather <- readRDS(file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS"))

# ACLED weekly conflict
acled <- readRDS(file.path(BASE_DIR, "data", "processed", "acled_daily.RDS")) |>
  select(date, week_date,
         libya_conflict, libya_conflict_fatalities,
         libya_battles, libya_expvio, libya_violciv,
         tunisia_conflict, tunisia_conflict_fatalities,
         tunisia_battles, tunisia_expvio, tunisia_violciv)
cat(sprintf("  ACLED: %d days loaded\n", nrow(acled)))

# Daily LCG/TCG pushbacks (from 01_temporal_disaggregation.R, Denton method).
# disagg_path/PANEL_END were already verified and read above (Section 5).
lcg_tcg_daily <- readRDS(disagg_path) |>
  select(date, lcg_pushbacks, tcg_pushbacks, lcg_tcg_pushbacks)
cat(sprintf("  LCG/TCG pushbacks (Denton disagg): %d days loaded\n", nrow(lcg_tcg_daily)))

# UNHCR daily arrivals (for arrivals not detected during Frontex operations)
unhcr <- readRDS(file.path(BASE_DIR, "data", "processed",
                            "unhcr_daily_arrivals.RDS"))
cat(sprintf("  UNHCR arrivals: %d days (%s to %s)\n",
    nrow(unhcr), min(unhcr$date), max(unhcr$date)))

# Merge everything
panel <- spine |>
  mutate(iso_week = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))) |>
  left_join(frx_daily,     by = "date") |>
  left_join(iom_daily,     by = "date") |>
  left_join(weather,       by = "date") |>
  left_join(acled,         by = "date") |>
  left_join(lcg_tcg_daily, by = "date") |>
  left_join(unhcr,         by = "date")

# Replace NA counts with 0 (Frontex and IOM counts).
# Covers the interceptor × SAR matrix (18 n_* + 18 persons_*), detection
# method counts, boat-type counts, departure counts, op-area, multi-actor,
# IOM deaths, and LCG/TCG pushbacks.
int_types <- c("ngo","eu","ita","comm","cg","land","noint","other","na")
frx_matrix_cols <- c(
  paste0("frx_n_sar_",       int_types),
  paste0("frx_n_notsar_",    int_types),
  paste0("frx_persons_sar_",    int_types),
  paste0("frx_persons_notsar_", int_types)
)
count_cols <- c("frx_incidents", "frx_persons", "frx_n_sar",
                frx_matrix_cols,
                "frx_n_det_fwa", "frx_n_det_helo", "frx_n_det_rpas",
                "frx_n_det_mas", "frx_n_det_aerial", "frx_n_det_ngo",
                "frx_n_det_call", "frx_n_det_cg", "frx_n_det_land",
                "frx_n_inflatable", "frx_persons_inflatable",
                "frx_n_wooden", "frx_persons_wooden",
                "frx_n_fibreglass", "frx_persons_fibreglass",
                "frx_dep_libya", "frx_dep_tunisia", "frx_dep_algeria",
                "frx_n_in_oparea", "frx_n_multi_actors",
                "n_dead_missing",
                "lcg_pushbacks", "tcg_pushbacks", "lcg_tcg_pushbacks")
panel <- panel |>
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
panel <- panel |>
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
annual <- panel |>
  group_by(year) |>
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
monthly_gap <- panel |>
  filter(!is.na(arrivals)) |>
  mutate(ym = floor_date(date, "month")) |>
  group_by(ym) |>
  summarise(unhcr = sum(arrivals), frx = sum(frx_persons), .groups = "drop") |>
  mutate(gap = pmax(unhcr - frx, 0))
cat(sprintf("    Monthly gap total: %.0f persons (%.1f%% of UNHCR)\n",
    sum(monthly_gap$gap),
    sum(monthly_gap$gap) / sum(monthly_gap$unhcr) * 100))

excess_days <- panel |>
  filter(n_dead_missing > frx_persons, n_dead_missing > 0)
cat(sprintf("\n  Days where dead+missing > Frontex persons: %d\n", nrow(excess_days)))

rate_days <- panel |> filter(!is.na(fatality_rate), crossing_attempts > 0)
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
