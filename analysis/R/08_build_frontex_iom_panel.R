# 08_build_frontex_iom_panel.R
# ============================
# Build integrated daily panel merging Frontex Themis with IOM MMP.
#
# Frontex provides: crossing denominator (persons), boat type, SAR flag,
#                    departure country, operational area.
# IOM provides:     death numerator (comprehensive, captures unrescued drownings).
#
# The two sources are combined at the daily level following the literature convention
# (Deiana et al. 2024, Battiston 2022, Rodriguez-Sanchez et al. 2023):
#   crossing_attempts = frontex_persons + iom_deaths
#   fatality_rate     = iom_deaths / crossing_attempts
#
# Input:
#   data/processed/frontex_incidents.RDS
#   data/processed/iom_mmp_incidents.RDS
#   data/processed/era5_swh_daily.RDS
#   data/processed/acled_daily.RDS
#   data/processed/archive/sar_ngo_ops_daily_RS.RDS
#   data/processed/archive/conflicts_interceptions_monthly_RS.RDS
#
# Output:
#   analysis/data/daily_panel_complete.RDS

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
SEA_CAUSES <- c("Drowning", "Mixed or unknown")
CORE <- list(lon_min = 10.0, lon_max = 15.1, lat_min = 32.4, lat_max = 37.8)
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

# Core corridor deaths (model convention from 05_reduced_form_model.R)
iom_core <- iom_raw %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) == "incident",
         `Cause of death (category)` %in% SEA_CAUSES,
         Longitude >= CORE$lon_min, Longitude <= CORE$lon_max,
         Latitude  >= CORE$lat_min, Latitude  <= CORE$lat_max) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date))

# Broad CMR deaths (all except sub-incidents)
iom_broad <- iom_raw %>%
  filter(Route == "Central Mediterranean",
         !tolower(`Incident Type`) %in% c("sub-incident")) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date))

cat(sprintf("  IOM core corridor: %d incidents\n", nrow(iom_core)))
cat(sprintf("  IOM broad CMR:     %d incidents\n", nrow(iom_broad)))

# ── 4. Aggregate IOM to daily ─────────────────────────────
cat("\n--- 4. Aggregating IOM to daily ---\n")

iom_daily_core <- iom_core %>%
  group_by(date) %>%
  summarise(iom_deaths_core    = sum(dead_missing),
            iom_incidents_core = n(), .groups = "drop")

iom_daily_broad <- iom_broad %>%
  group_by(date) %>%
  summarise(iom_deaths_allcmr    = sum(dead_missing),
            iom_incidents_allcmr = n(), .groups = "drop")

cat(sprintf("  IOM daily core:  %d days with deaths\n", nrow(iom_daily_core)))
cat(sprintf("  IOM daily broad: %d days with deaths\n", nrow(iom_daily_broad)))

# ── 5. Build date spine and merge ─────────────────────────
cat("\n--- 5. Building integrated panel ---\n")

# Date spine: 2014-01-01 to Frontex end date
spine <- tibble(date = seq(as.Date("2014-01-01"), FRX_END, by = "day"))

# Weather from ERA5 clean data
weather <- readRDS(file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS")) %>%
  select(date, swh, swh_prev3days, swh_prevweek)

# SAR (NGO vessel counts from Rodriguez-Sanchez)
sar_path <- file.path(BASE_DIR, "data", "processed", "archive", "sar_ngo_ops_daily_RS.RDS")
if (file.exists(sar_path)) {
  sar <- readRDS(sar_path)
  cat(sprintf("  SAR: %d days loaded\n", nrow(sar)))
} else {
  sar <- tibble(date = as.Date(character()), n_ngo_vessels = integer(), n_gov_operations = integer())
  cat("  SAR data not found — skipping\n")
}

# ACLED weekly conflict
acled <- readRDS(file.path(BASE_DIR, "data", "processed", "acled_daily.RDS")) %>%
  select(date, week_date,
         libya_conflict, libya_conflict_fatalities,
         libya_battles, libya_expvio, libya_violciv,
         tunisia_conflict, tunisia_conflict_fatalities,
         tunisia_battles, tunisia_expvio, tunisia_violciv)
cat(sprintf("  ACLED: %d days loaded\n", nrow(acled)))

# Monthly LCG/TCG interceptions
ic_path <- file.path(BASE_DIR, "data", "processed", "archive",
                      "conflicts_interceptions_monthly_RS.RDS")
if (file.exists(ic_path)) {
  interceptions_monthly <- readRDS(ic_path) %>%
    select(date, lcg_interceptions, tcg_interceptions)

  interceptions_daily <- spine %>%
    mutate(month_date = floor_date(date, "month")) %>%
    left_join(interceptions_monthly, by = c("month_date" = "date")) %>%
    select(date, lcg_interceptions, tcg_interceptions)
  cat(sprintf("  LCG/TCG interceptions: %d months with data\n",
      sum(!is.na(interceptions_daily$lcg_interceptions))))
} else {
  interceptions_daily <- spine %>%
    mutate(lcg_interceptions = NA_real_, tcg_interceptions = NA_real_)
  cat("  LCG/TCG interceptions data not found — skipping\n")
}

# Merge everything
panel <- spine %>%
  mutate(iso_week = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))) %>%
  left_join(frx_daily, by = "date") %>%
  left_join(iom_daily_core, by = "date") %>%
  left_join(iom_daily_broad, by = "date") %>%
  left_join(weather, by = "date") %>%
  left_join(sar, by = "date") %>%
  left_join(acled, by = "date") %>%
  left_join(interceptions_daily, by = "date")

# Replace NA counts with 0 (Frontex and IOM counts)
count_cols <- c("frx_incidents", "frx_persons", "frx_n_sar",
                "frx_n_inflatable", "frx_persons_inflatable",
                "frx_n_wooden", "frx_persons_wooden",
                "frx_dep_libya", "frx_dep_tunisia", "frx_dep_algeria",
                "frx_n_in_oparea",
                "iom_deaths_core", "iom_incidents_core",
                "iom_deaths_allcmr", "iom_incidents_allcmr")
panel <- panel %>%
  mutate(across(all_of(count_cols), ~ replace_na(.x, 0L)))

cat(sprintf("  Panel: %d days (%s to %s)\n",
    nrow(panel), min(panel$date), max(panel$date)))

# ── 6. Derive variables ───────────────────────────────────
cat("\n--- 6. Deriving variables ---\n")

panel <- panel %>%
  mutate(
    crossing_attempts = frx_persons + iom_deaths_core,
    fatality_rate     = ifelse(crossing_attempts > 0,
                                iom_deaths_core / crossing_attempts, NA_real_),
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
    iom_deaths_core   = sum(iom_deaths_core),
    iom_deaths_allcmr = sum(iom_deaths_allcmr),
    crossing_attempts = sum(crossing_attempts),
    frx_incidents     = sum(frx_incidents),
    days_with_frx     = sum(frx_incidents > 0),
    .groups = "drop"
  )
print(annual, n = 12)

excess_days <- panel %>%
  filter(iom_deaths_core > frx_persons, iom_deaths_core > 0)
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
