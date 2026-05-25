# Zone-level daily panel: expand daily_panel_complete to (date × SAR country)
# long format, with deaths assigned by spatial join to corridor-intersected zones.

library(dplyr)
library(tidyr)
library(sf)
library(lubridate)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

sf_use_s2(FALSE)

# ── 1. Base panel ───────────────────────────────────────────────────────────
base <- readRDS(file.path(BASE_DIR, "analysis", "data",
                           "daily_panel_complete.RDS")) |>
  rename(n_dead_missing_nat = n_dead_missing)

# ── 2. SAR zones intersected with corridor ──────────────────────────────────
zones_in <- readRDS(file.path(BASE_DIR, "data", "processed",
                                "sar_zones_in_corridor.RDS"))

# ── 3. IOM incidents (analytical filter, geometry kept for spatial join) ────
iom <- iom_incidents(spatial = "all_cmr") |>
  drop_na(lon, lat) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# ── 3b. UNITED incidents (analytical filter, geometry kept for spatial join) ─
united <- united_incidents(spatial = "all_cmr") |>
  tidyr::drop_na(longitude, latitude) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# ── 4. Assign deaths to corridor-intersected zones ──────────────────────────
iom_zoned <- iom |>
  st_join(zones_in |> select(country), left = TRUE) |>
  st_drop_geometry()

deaths_daily_zone <- iom_zoned |>
  filter(!is.na(country)) |>
  group_by(date, country) |>
  summarise(n_dead_missing = sum(dead_missing),
            n_incidents    = n(),
            .groups = "drop")

united_zoned <- united |>
  st_join(zones_in |> select(country), left = TRUE) |>
  st_drop_geometry()

deaths_daily_zone_united <- united_zoned |>
  filter(!is.na(country)) |>
  group_by(date, country) |>
  summarise(n_dead_united = sum(n_deaths, na.rm = TRUE),
            .groups = "drop")

# ── 5. Expand to (date × country) and merge ─────────────────────────────────
panel_zone <- base |>
  tidyr::crossing(country = c("Italy", "Libya", "Malta", "Tunisia")) |>
  left_join(deaths_daily_zone, by = c("date", "country")) |>
  left_join(deaths_daily_zone_united, by = c("date", "country")) |>
  tidyr::replace_na(list(n_dead_missing = 0, n_incidents = 0,
                          n_dead_united = 0)) |>
  left_join(st_drop_geometry(zones_in)[, c("country", "zone_area_km2")],
             by = "country") |>
  mutate(sar_bloc = if_else(country %in% c("Italy", "Malta"), "EU", "AFR"))

dim(panel_zone$date) <- NULL

# ── 6. Sanity checks ────────────────────────────────────────────────────────
stopifnot(nrow(panel_zone) == length(unique(panel_zone$date)) * 4)

req_cols <- c("swh", "swh_prev5days", "n_dead_missing", "n_dead_missing_nat",
               "n_dead_united", "zone_area_km2", "sar_bloc", "country",
               "post_mou", "month_year", "iso_week", "crossing_attempts")
stopifnot(length(setdiff(req_cols, names(panel_zone))) == 0)

swh_check <- panel_zone |>
  group_by(date) |>
  summarise(swh_sd = sd(swh, na.rm = TRUE), .groups = "drop")
max_within_day_sd <- max(swh_check$swh_sd, na.rm = TRUE)
stopifnot(is.na(max_within_day_sd) || max_within_day_sd < 1e-10)

# ── 7. Save ─────────────────────────────────────────────────────────────────
saveRDS(panel_zone,
        file.path(BASE_DIR, "analysis", "data", "daily_panel_zone.RDS"))
