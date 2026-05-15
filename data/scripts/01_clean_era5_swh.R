# 01_clean_era5_swh.R
# ===================
# Build ERA5 daily SWH panel over the CMR core corridor sea zone.
# Single series: `swh` = simple spatial mean over all ocean cells in
# the core-corridor polygon.
#
# Input:  data/raw/era5/era5_daily_cmr_wave_[YEAR].nc  (2008-2025)
#         data/processed/core_corridor.RDS
# Output: data/processed/era5_swh_daily.RDS
#
# Columns: date, swh,
#          swh_lag1, swh_prev3days, swh_prev5days, swh_prevweek
# Windows: 1-day (lag1), 1-3 day, 1-5 day, 1-7 day rolling means, each
# lagged by one day so day-t regressors exclude day-t weather.

library(dplyr)
library(sf)
library(ncdf4)
library(zoo)
library(lubridate)
library(purrr)

sf_use_s2(FALSE)

BASE_DIR <- here::here()
ERA5_DIR <- file.path(BASE_DIR, "data", "raw", "era5")

cat("============================================================\n")
cat("CLEAN ERA5 DAILY SWH PANEL\n")
cat("============================================================\n\n")

# ── 1. Load core corridor polygon ─────────────────────────
# Single source of truth: analysis/R/00_define_sea_zones.R builds the
# polygon and saves it to data/processed/core_corridor.RDS. Load it so
# the SWH spatial mean and the death filter always use the same bounds.
cat("--- 1. Loading core corridor polygon ---\n")

outer_poly <- readRDS(file.path(BASE_DIR, "data", "processed",
                                  "core_corridor.RDS"))
sea_vis    <- outer_poly  # kept for any downstream visualization hooks

cat("  Core Corridor polygon loaded from core_corridor.RDS\n")

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
cat("--- 3. Extracting daily SWH (corridor-wide mean) ---\n")

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

  swh_day <- vapply(seq_along(dates),
                     \(t) mean(swh_3d[, , t][mask], na.rm = TRUE),
                     numeric(1))

  cat(sprintf("  %d: %d days\n", yr, length(dates)))
  tibble(date = dates, swh = swh_day)
})

# ── 4. Compute lagged rolling means ───────────────────────
# Windows: 1-day (lag1), 1-3 day, 1-5 day, 1-7 day rolling means,
# each lagged one day so day-t regressors exclude day-t weather.
cat("--- 4. Computing lagged rolling means ---\n")

weather <- weather %>%
  arrange(date) %>%
  mutate(
    swh_lag1      = lag(swh, 1),
    swh_prev3days = zoo::rollmeanr(lag(swh, 1), k = 3, fill = NA),
    swh_prev5days = zoo::rollmeanr(lag(swh, 1), k = 5, fill = NA),
    swh_prevweek  = zoo::rollmeanr(lag(swh, 1), k = 7, fill = NA)
  )

cat(sprintf("  Total: %d weather days (%s to %s)\n",
    nrow(weather), min(weather$date), max(weather$date)))
cat(sprintf("  Columns: %d\n", ncol(weather)))
cat(sprintf("  swh range: %.2f to %.2f\n",
    min(weather$swh, na.rm = TRUE), max(weather$swh, na.rm = TRUE)))

cat("  Window summary (lagged rolling means):\n")
for (nm in c("swh_lag1", "swh_prev3days", "swh_prev5days", "swh_prevweek")) {
  v <- weather[[nm]]
  cat(sprintf("    %-14s mean=%.3f  range=[%.2f, %.2f]\n",
              nm, mean(v, na.rm = TRUE),
              min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
}

# Sanity check: wider window smooths variance (sd monotonically non-increasing)
sds <- sapply(c("swh_lag1", "swh_prev3days", "swh_prev5days", "swh_prevweek"),
              function(nm) sd(weather[[nm]], na.rm = TRUE))
cat(sprintf("  SDs (1d, 3d, 5d, 7d): %s\n",
    paste(sprintf("%.3f", sds), collapse = ", ")))
stopifnot(all(diff(sds) <= 0))

# ── 5. Save ──────────────────────────────────────────────
saveRDS(weather, file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS"))
cat("\nSaved: data/processed/era5_swh_daily.RDS\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
