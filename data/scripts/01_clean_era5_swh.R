# 01_clean_era5_swh.R
# ===================
# Build ERA5 daily SWH panel: spatial-mean significant wave height over
# the CMR core corridor sea zone.
#
# Sea zone definition (Core Corridor):
#   Follows Tunisian/Libyan coastlines, diagonal Cap Bon → W Sicily,
#   Sicily south coast, east edge at ~16-17E.
#   Used ONLY for wave height calculation — analysis scripts include
#   all CMR incidents regardless of location.
#
# Input:  data/raw/era5/era5_daily_cmr_wave_[YEAR].nc  (2008-2025)
# Output: data/processed/era5_swh_daily.RDS
#
# Columns: date, swh, swh_prev3days, swh_prevweek

library(dplyr)
library(sf)
library(ncdf4)
library(rnaturalearth)
library(zoo)
library(lubridate)
library(purrr)

BASE_DIR <- here::here()
ERA5_DIR <- file.path(BASE_DIR, "data", "raw", "era5")

cat("============================================================\n")
cat("CLEAN ERA5 DAILY SWH PANEL\n")
cat("============================================================\n\n")

# ── 1. Define core corridor sea zone ─────────────────────
cat("--- 1. Defining core corridor sea zone ---\n")

coords_core <- matrix(c(
    17, 30.0,
  15.1, 36.7,
  12.4, 37.8,
  11.0, 37.1,
   9.0, 34.0,
   9.0, 31.0,
    17, 30.0
), ncol = 2, byrow = TRUE)

outer_poly <- st_sfc(st_polygon(list(coords_core)), crs = 4326)

world <- ne_countries(scale = "medium", returnclass = "sf")
land  <- st_union(world)

sea_raw  <- st_difference(outer_poly, land)
sea_proj <- st_transform(sea_raw, 3857)
sea_buf  <- st_buffer(sea_proj, dist = 5000)
sea_back <- st_transform(sea_buf, 4326)
sea_vis  <- st_intersection(sea_back, outer_poly)

cat("  Sea zone polygon built\n")

# ── 2. Build spatial mask from ERA5 grid ─────────────────
cat("--- 2. Building spatial mask ---\n")

nc0 <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
wave_lon <- ncvar_get(nc0, "longitude")
wave_lat <- ncvar_get(nc0, "latitude")
nc_close(nc0)

grid_pts <- expand.grid(lon = wave_lon, lat = wave_lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
in_zone <- st_intersects(grid_pts, sea_vis, sparse = FALSE)[, 1]

# Also check for ocean (non-NA SWH) on first day
nc_chk <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
swh_day1 <- ncvar_get(nc_chk, "swh", start = c(1, 1, 1), count = c(-1, -1, 1))
nc_close(nc_chk)
mask <- matrix(!is.na(as.vector(swh_day1)) & in_zone, nrow = length(wave_lon))

cat(sprintf("  ERA5 grid: %d lon x %d lat\n", length(wave_lon), length(wave_lat)))
cat(sprintf("  Ocean cells in core zone: %d\n", sum(mask)))

# ── 3. Extract daily SWH from netCDF files ───────────────
cat("--- 3. Extracting daily SWH ---\n")

get_nc_dates <- function(nc) {
  tn <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  tv <- ncvar_get(nc, tn)
  tu <- ncatt_get(nc, tn, "units")$value
  if (grepl("seconds since 1970", tu)) {
    as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC"))
  } else if (grepl("hours since 1900", tu)) {
    as.Date(as.POSIXct("1900-01-01", tz = "UTC") + tv * 3600)
  } else {
    ref <- sub("(hours|seconds|days) since ", "", tu)
    mult <- ifelse(grepl("hours", tu), 3600, ifelse(grepl("days", tu), 86400, 1))
    as.Date(as.POSIXct(ref, tz = "UTC") + tv * mult)
  }
}

weather <- map_dfr(2008:2025, function(yr) {
  f <- file.path(ERA5_DIR, paste0("era5_daily_cmr_wave_", yr, ".nc"))
  if (!file.exists(f)) {
    cat(sprintf("  %d: file not found, skipping\n", yr))
    return(tibble())
  }
  nc <- nc_open(f)
  dates <- get_nc_dates(nc)
  swh_3d <- ncvar_get(nc, "swh")
  nc_close(nc)

  out <- tibble(
    date = dates,
    swh = map_dbl(seq_along(dates), ~ mean(swh_3d[, , .x][mask], na.rm = TRUE))
  )
  cat(sprintf("  %d: %d days\n", yr, nrow(out)))
  out
})

# ── 4. Compute rolling averages ──────────────────────────
cat("--- 4. Computing rolling averages ---\n")

weather <- weather %>%
  arrange(date) %>%
  mutate(
    swh_prev3days = zoo::rollmeanr(lag(swh, 1), k = 3, fill = NA),
    swh_prevweek  = zoo::rollmeanr(lag(swh, 1), k = 7, fill = NA)
  )

cat(sprintf("  Total: %d weather days (%s to %s)\n",
    nrow(weather), min(weather$date), max(weather$date)))
cat(sprintf("  SWH range: %.2f to %.2f\n",
    min(weather$swh, na.rm = TRUE), max(weather$swh, na.rm = TRUE)))

# ── 5. Save ──────────────────────────────────────────────
saveRDS(weather, file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS"))
cat("\nSaved: data/processed/era5_swh_daily.RDS\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
