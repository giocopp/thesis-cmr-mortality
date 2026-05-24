# Build ERA5 daily SWH panel over the CMR core-corridor sea zone.

library(dplyr)
library(sf)
library(ncdf4)
library(zoo)
library(lubridate)
library(purrr)

sf_use_s2(FALSE)

BASE_DIR <- here::here()
ERA5_DIR <- file.path(BASE_DIR, "data", "raw", "era5")

# ── 1. Load core corridor polygon ───────────────────────────────────────────
outer_poly <- readRDS(file.path(BASE_DIR, "data", "processed",
                                  "core_corridor.RDS"))
sea_vis    <- outer_poly

# ── 2. Build spatial mask from ERA5 grid ────────────────────────────────────
nc0 <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
wave_lon <- ncvar_get(nc0, "longitude")
wave_lat <- ncvar_get(nc0, "latitude")
nc_close(nc0)

grid_pts <- expand.grid(lon = wave_lon, lat = wave_lat) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
in_zone <- st_intersects(grid_pts, sea_vis, sparse = FALSE)[, 1]

nc_chk <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
swh_day1 <- ncvar_get(nc_chk, "swh", start = c(1, 1, 1), count = c(-1, -1, 1))
nc_close(nc_chk)
mask <- matrix(!is.na(as.vector(swh_day1)) & in_zone, nrow = length(wave_lon))

# ── 3. Extract daily SWH from netCDF files ──────────────────────────────────
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
  if (!file.exists(f)) return(tibble())
  nc <- nc_open(f)
  dates <- get_nc_dates(nc)
  swh_3d <- ncvar_get(nc, "swh")
  nc_close(nc)

  swh_day <- vapply(seq_along(dates),
                     \(t) mean(swh_3d[, , t][mask], na.rm = TRUE),
                     numeric(1))

  tibble(date = dates, swh = swh_day)
})

# ── 4. Lagged rolling means (day-t regressors exclude day-t weather) ────────
weather <- weather |>
  arrange(date) |>
  mutate(
    swh_lag1      = lag(swh, 1),
    swh_prev3days = zoo::rollmeanr(lag(swh, 1), k = 3, fill = NA),
    swh_prev5days = zoo::rollmeanr(lag(swh, 1), k = 5, fill = NA),
    swh_prevweek  = zoo::rollmeanr(lag(swh, 1), k = 7, fill = NA)
  )

sds <- sapply(c("swh_lag1", "swh_prev3days", "swh_prev5days", "swh_prevweek"),
              function(nm) sd(weather[[nm]], na.rm = TRUE))
stopifnot(all(diff(sds) <= 0))

# ── 5. Save ─────────────────────────────────────────────────────────────────
saveRDS(weather, file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS"))
