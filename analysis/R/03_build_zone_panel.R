# 03_build_zone_panel.R
# =====================
# Build zone-level daily panel by STARTING from analysis/data/daily_panel_complete.RDS
# (the canonical daily-agg panel built by 02_build_daily_panel.R) and
# expanding it to (date x country) long format. Each zone row inherits the
# corridor-wide SWH from the base panel — ALL four zone rows on a given day
# share the same SWH, swh_prevweek, etc.
#
# Specification (option A2):
#   - SWH: single corridor-wide series from era5_swh_daily.RDS (already in
#     daily_panel_complete.RDS). Same value for all 4 zone rows on a given
#     day. Rationale: zone-specific SWH series are 77%-99% correlated with
#     the corridor-wide series (see diagnostic), so the marginal identifying
#     variation is small, and small-N zone means (Italy 22 cells, Tunisia 15
#     cells) add noise. The 4-country model still identifies zone-specific
#     responses because different zones have different death DGPs.
#
#   - Deaths: spatially joined to CORRIDOR-INTERSECTED SAR polygons
#     (sar_zones_in_corridor.RDS). Incidents outside the CMR core corridor
#     are dropped (~14% of CMR incidents). This matches the "we only care
#     about CMR" criterion and keeps Libya/Tunisia intact (99%+ retained).
#
#   - Date range: inherited from daily_panel_complete.RDS, i.e., 2014-01-01
#     to 2023-06-09 (bounded by Frontex). Zone-only data after 2023-06-09
#     is not used, accepting the trade-off for consistency with daily-agg.
#
#   - zone_area_km2: area of the CORRIDOR-INTERSECTED polygon (not the full
#     SAR), via lwgeom::st_area. Used for area-weighted 2-bloc collapses
#     (though area weighting is vacuous here since all zones share the same
#     SWH — kept for structural consistency and future extensibility).
#
# Column semantics:
#   - n_dead_missing      : deaths assigned to this SAR zone (spatial join to
#                           corridor-intersected polygon). Filter matches the
#                           primary analytical spec in 05_reduced_form_primary.R:
#                           incident-only, Cause = Drowning or Mixed/unknown.
#                           Zero if no matching deaths in this (date, country).
#   - n_dead_missing_nat  : national daily count inherited from the base panel
#                           built by 02 (broad descriptive filter: incident +
#                           split incident, all CMR countries, no cause filter).
#                           This is BROADER than n_dead_missing and is kept as a
#                           reference column. Same value on all 4 zone rows for
#                           a given date.
#   - swh, swh_prev3days,
#     swh_prevweek        : corridor-wide SWH inherited from base. Same on
#                           all 4 zone rows.
#   - zone_area_km2       : static area of each corridor-intersected polygon
#   - crossings, frx_*,
#     acled_*, etc.       : all inherited from base
#
# Input:
#   analysis/data/daily_panel_complete.RDS
#   data/processed/sar_zones_in_corridor.RDS
#   data/processed/iom_mmp_incidents.RDS
#
# Output: analysis/data/daily_panel_zone.RDS

library(dplyr)
library(tidyr)
library(sf)
library(lubridate)
# Note: lwgeom is not loaded here — zone areas are pre-computed and stored
# in sar_zones_in_corridor.RDS by analysis/R/00_define_sea_zones.R.

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")

sf_use_s2(FALSE)  # Italian SRR has self-intersecting edges

cat("============================================================\n")
cat("BUILD ZONE-LEVEL DAILY PANEL (from daily_panel_complete base)\n")
cat("Option A2: corridor SWH + corridor-intersected death assignment\n")
cat("============================================================\n\n")

# ── 1. Load the canonical daily base ────────────────────────
cat("--- 1. Loading daily_panel_complete.RDS as base ---\n")

base <- readRDS(file.path(BASE_DIR, "analysis", "data",
                           "daily_panel_complete.RDS")) %>%
  rename(n_dead_missing_nat = n_dead_missing)

cat(sprintf("  Base: %d rows x %d columns (%s to %s)\n",
    nrow(base), ncol(base), min(base$date), max(base$date)))

# ── 2. Load corridor-intersected SAR zones ──────────────────
cat("\n--- 2. Loading corridor-intersected SAR zones ---\n")

zones_in <- readRDS(file.path(BASE_DIR, "data", "processed",
                                "sar_zones_in_corridor.RDS"))
cat("  Intersected zone areas (km^2):\n")
print(st_drop_geometry(zones_in)[, c("country", "zone_area_km2")])

# ── 3. Load and filter IOM incidents ─────────────────────────
cat("\n--- 3. Loading and filtering IOM incidents ---\n")

# Filter matches 05_reduced_form_primary.R: incident + split incident,
# Cause of death in {Drowning, Mixed or unknown} (the cases most directly
# tied to the act of crossing the sea — see 05 header for the rationale).
# Geographic membership is resolved below by spatial join to the
# corridor-intersected SAR polygons.
iom <- readRDS(file.path(BASE_DIR, "data", "processed", "iom_mmp_incidents.RDS")) %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) %in% c("incident", "split incident"),
         `Cause of death (category)` %in% c("Drowning", "Mixed or unknown")) %>%
  transmute(date         = as.Date(incident_date_clean),
            dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE),
            lon          = as.numeric(Longitude),
            lat          = as.numeric(Latitude)) %>%
  drop_na(date, lon, lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

cat(sprintf("  IOM CMR incidents (incident + split, drowning + mixed): %d\n", nrow(iom)))

# ── 4. Spatial join to corridor-intersected zones ──────────
cat("\n--- 4. Spatial join to corridor-intersected zones ---\n")

iom_zoned <- iom %>%
  st_join(zones_in %>% select(country), left = TRUE) %>%
  st_drop_geometry()

n_in_corridor <- sum(!is.na(iom_zoned$country))
deaths_in_corridor <- sum(iom_zoned$dead_missing[!is.na(iom_zoned$country)])

cat(sprintf("  Inside a corridor-intersected SAR: %d incidents, %.0f deaths\n",
    n_in_corridor, deaths_in_corridor))
cat(sprintf("  Outside all intersected SARs:      %d incidents, %.0f deaths\n",
    sum(is.na(iom_zoned$country)),
    sum(iom_zoned$dead_missing[is.na(iom_zoned$country)])))

deaths_daily_zone <- iom_zoned %>%
  filter(!is.na(country)) %>%
  group_by(date, country) %>%
  summarise(n_dead_missing = sum(dead_missing),
            n_incidents    = n(),
            .groups = "drop")

cat(sprintf("  %d (date x zone) cells with deaths\n", nrow(deaths_daily_zone)))

# ── 5. Expand base to (date x country), merge deaths ────────
cat("\n--- 5. Building zone panel ---\n")

panel_zone <- base %>%
  tidyr::crossing(country = c("Italy", "Libya", "Malta", "Tunisia")) %>%
  left_join(deaths_daily_zone, by = c("date", "country")) %>%
  tidyr::replace_na(list(n_dead_missing = 0, n_incidents = 0)) %>%
  left_join(st_drop_geometry(zones_in)[, c("country", "zone_area_km2")],
             by = "country") %>%
  mutate(sar_bloc = if_else(country %in% c("Italy", "Malta"), "EU", "AFR"))

# Strip dim attribute from date
dim(panel_zone$date) <- NULL

cat(sprintf("  Panel: %d rows x %d cols (%s to %s)\n",
    nrow(panel_zone), ncol(panel_zone),
    min(panel_zone$date), max(panel_zone$date)))
cat(sprintf("  %d dates x 4 zones = %d (expected match)\n",
    length(unique(panel_zone$date)),
    length(unique(panel_zone$date)) * 4))

# ── 6. Sanity checks ────────────────────────────────────────
cat("\n--- 6. Sanity checks ---\n")

stopifnot(nrow(panel_zone) == length(unique(panel_zone$date)) * 4)
cat("  [OK] rows = days x 4\n")

req_cols <- c("swh", "swh_prevweek", "n_dead_missing", "n_dead_missing_nat",
               "zone_area_km2", "sar_bloc", "country", "post_mou",
               "month_year", "iso_week", "crossing_attempts")
missing_cols <- setdiff(req_cols, names(panel_zone))
stopifnot(length(missing_cols) == 0)
cat("  [OK] all required columns present\n")

# Verify SWH is the SAME across zones on a given day (defining feature of A2)
swh_check <- panel_zone %>%
  group_by(date) %>%
  summarise(swh_sd = sd(swh, na.rm = TRUE), .groups = "drop")
max_within_day_sd <- max(swh_check$swh_sd, na.rm = TRUE)
stopifnot(is.na(max_within_day_sd) || max_within_day_sd < 1e-10)
cat("  [OK] SWH is constant across zones within each day (A2 property)\n")

cat("  Deaths by bloc (spatial-join assignment, corridor-intersected):\n")
bloc_deaths <- panel_zone %>%
  group_by(sar_bloc) %>%
  summarise(total_deaths = sum(n_dead_missing),
            n_rows       = n(), .groups = "drop")
print(bloc_deaths)

cat("  Deaths by country:\n")
country_deaths <- panel_zone %>%
  group_by(country) %>%
  summarise(total_deaths = sum(n_dead_missing), .groups = "drop")
print(country_deaths)

# ── 7. Save ─────────────────────────────────────────────────
cat("\n--- 7. Saving ---\n")

out_path <- file.path(BASE_DIR, "analysis", "data", "daily_panel_zone.RDS")
saveRDS(panel_zone, out_path)
cat(sprintf("Saved: %s\n", out_path))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
