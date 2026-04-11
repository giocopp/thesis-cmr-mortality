# 01_clean_era5_swh.R
# ===================
# Build ERA5 daily SWH panel over the CMR core corridor sea zone.
# Two parallel series:
#   swh   — simple spatial mean over all ocean cells in the polygon
#   swh_w — death-weighted spatial mean. Weights are the historical
#           count of CMR deaths per cell (all years, static), snapped
#           from IOM incidents to the nearest ERA5 cell. Cells with
#           zero deaths contribute nothing to swh_w.
#
# Both series are exogenous to any single day's outcome: swh weights
# cells uniformly, swh_w weights by pre-computed historical density.
#
# Input:  data/raw/era5/era5_daily_cmr_wave_[YEAR].nc  (2008-2025)
#         data/processed/core_corridor.RDS
#         data/processed/iom_mmp_incidents.RDS
# Output: data/processed/era5_swh_daily.RDS
#
# Columns: date,
#          swh,   swh_prev3days,   swh_prevweek
#          swh_w, swh_w_prev3days, swh_w_prevweek

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

# ── 2b. Build per-cell death-weight from IOM incidents ────
# Static weights: historical death density per ERA5 cell, pooled over
# the whole sample. Weights are time-invariant so they don't create
# feedback between day-t outcome and day-t regressor.
cat("--- 2b. Building per-cell death weights ---\n")

iom <- readRDS(file.path(BASE_DIR, "data", "processed",
                           "iom_mmp_incidents.RDS")) %>%
  dplyr::filter(
    Route == "Central Mediterranean",
    tolower(`Incident Type`) %in% c("incident", "split incident"),
    `Country of Incident` %in% c("Algeria","Italy","Libya","Malta","Tunisia"),
    `Cause of death (category)` %in% c("Drowning", "Mixed or unknown")
  ) %>%
  dplyr::mutate(
    lon  = as.numeric(Longitude),
    lat  = as.numeric(Latitude),
    dead = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)
  ) %>%
  tidyr::drop_na(lon, lat)

iom$lon_idx <- vapply(iom$lon, \(x) which.min(abs(wave_lon - x)), integer(1))
iom$lat_idx <- vapply(iom$lat, \(x) which.min(abs(wave_lat - x)), integer(1))

weight_mat <- matrix(0, nrow = length(wave_lon), ncol = length(wave_lat))
for (i in seq_len(nrow(iom))) {
  li <- iom$lon_idx[i]; la <- iom$lat_idx[i]
  if (mask[li, la]) {
    weight_mat[li, la] <- weight_mat[li, la] + iom$dead[i]
  }
}

cat(sprintf("  IOM incidents loaded (full 5-country + sea-cause filter): %d\n",
            nrow(iom)))
cat(sprintf("  Total deaths mapped to mask cells: %.0f\n", sum(weight_mat)))
cat(sprintf("  Mask cells with >=1 death: %d / %d\n",
            sum(weight_mat > 0 & mask), sum(mask)))

# ── 3. Extract daily SWH from netCDF files ───────────────
cat("--- 3. Extracting daily SWH (mean + death-weighted) ---\n")

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

# Weight vector flattened to match mask cell order
w_mask <- weight_mat[mask]  # length = sum(mask)
sum_w  <- sum(w_mask)

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

  swh_day  <- vapply(seq_along(dates),
                      \(t) mean(swh_3d[, , t][mask], na.rm = TRUE),
                      numeric(1))
  swh_w_day <- vapply(seq_along(dates), function(t) {
    v <- swh_3d[, , t][mask]
    ok <- !is.na(v) & w_mask > 0
    if (!any(ok)) return(NA_real_)
    sum(v[ok] * w_mask[ok]) / sum(w_mask[ok])
  }, numeric(1))

  cat(sprintf("  %d: %d days\n", yr, length(dates)))
  tibble(date = dates, swh = swh_day, swh_w = swh_w_day)
})

# ── 4. Compute rolling averages (both series) ─────────────
cat("--- 4. Computing rolling averages ---\n")

weather <- weather %>%
  arrange(date) %>%
  mutate(
    swh_prev3days   = zoo::rollmeanr(lag(swh,   1), k = 3, fill = NA),
    swh_prevweek    = zoo::rollmeanr(lag(swh,   1), k = 7, fill = NA),
    swh_w_prev3days = zoo::rollmeanr(lag(swh_w, 1), k = 3, fill = NA),
    swh_w_prevweek  = zoo::rollmeanr(lag(swh_w, 1), k = 7, fill = NA)
  )

cat(sprintf("  Total: %d weather days (%s to %s)\n",
    nrow(weather), min(weather$date), max(weather$date)))
cat(sprintf("  swh    range: %.2f to %.2f\n",
    min(weather$swh,   na.rm = TRUE), max(weather$swh,   na.rm = TRUE)))
cat(sprintf("  swh_w  range: %.2f to %.2f\n",
    min(weather$swh_w, na.rm = TRUE), max(weather$swh_w, na.rm = TRUE)))
cat(sprintf("  cor(swh, swh_w): %.4f\n",
    cor(weather$swh, weather$swh_w, use = "pairwise.complete.obs")))

# ── 5. Save ──────────────────────────────────────────────
saveRDS(weather, file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS"))
cat("\nSaved: data/processed/era5_swh_daily.RDS\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
