# 02b_daily_weather_panel.R
# ========================
# Compute daily CMR spatial-mean weather from ERA5 grids.
# This creates a daily panel with weather for ALL days (not just incident days),
# which is needed for the daily count model.
#
# GEOGRAPHIC RESTRICTION:
#   Weather is averaged only over ocean cells within the CMR corridor,
#   defined from the incident data footprint (where crossings actually happen).
#
#   Corridor: lon [10, 18], lat [31, 38]
#   This captures ~90% of CMR incidents. The core Libya-Lampedusa channel
#   is lon [11, 15], lat [32, 36] (captures ~80%).
#   We also compute a "core" version for robustness.
#
# Input:  data/raw/era5/era5_daily_cmr_{atm_instant,atm_accum,atm_gust,wave}_YYYY.nc
# Output: data/processed/cmr_daily_weather_panel.RDS
#
# Variables computed (for corridor and core):
#   - swh_mean / swh_core:   spatial mean SWH (m) over CMR ocean cells
#   - wind_mean / wind_core: spatial mean wind speed (m/s)
#   - gust_mean / gust_core: spatial mean wind gust (m/s)
#   - mwp_mean / mwp_core:   spatial mean wave period (s)

library(ncdf4)
library(dplyr)
library(lubridate)

BASE_DIR <- here::here()
ERA5_DIR <- file.path(BASE_DIR, "data", "raw", "era5")
IOM_PATH <- file.path(BASE_DIR, "data", "processed", "cmr_events_with_weather.RDS")
OUTPUT_PATH <- file.path(BASE_DIR, "data", "processed", "cmr_daily_weather_panel.RDS")

YEARS <- 2014:2025

# ============================================================
# 0. Define CMR corridor bounds
# ============================================================
# Corridor: covers ~90% of CMR incidents
CORRIDOR <- list(lon_min = 10, lon_max = 18, lat_min = 31, lat_max = 38)
# Core: covers the Libya-Lampedusa channel (~80% of incidents)
CORE     <- list(lon_min = 11, lon_max = 15, lat_min = 32, lat_max = 36)

# ============================================================
# 1. Helpers
# ============================================================
get_nc_dates <- function(nc) {
  time_name <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  time_vals <- ncvar_get(nc, time_name)
  time_units <- ncatt_get(nc, time_name, "units")$value
  if (grepl("seconds since 1970", time_units)) {
    as.Date(as.POSIXct(time_vals, origin = "1970-01-01", tz = "UTC"))
  } else if (grepl("hours since 1900", time_units)) {
    as.Date(as.POSIXct("1900-01-01", tz = "UTC") + time_vals * 3600)
  } else {
    ref <- sub("(hours|seconds|days) since ", "", time_units)
    multiplier <- ifelse(grepl("hours", time_units), 3600,
                         ifelse(grepl("days", time_units), 86400, 1))
    as.Date(as.POSIXct(ref, tz = "UTC") + time_vals * multiplier)
  }
}

# Build a logical mask for grid cells within a bounding box
make_mask <- function(lon_vec, lat_vec, bounds) {
  lon_in <- lon_vec >= bounds$lon_min & lon_vec <= bounds$lon_max
  lat_in <- lat_vec >= bounds$lat_min & lat_vec <= bounds$lat_max
  # Outer product: [lon, lat] matrix
  outer(lon_in, lat_in, "&")
}

# Masked spatial mean: average only cells where mask=TRUE (and not NA)
masked_mean <- function(slice_2d, mask) {
  vals <- slice_2d[mask]
  mean(vals, na.rm = TRUE)
}

# ============================================================
# 2. Read grids and build masks
# ============================================================
cat("============================================================\n")
cat("Computing daily CMR corridor weather from ERA5\n")
cat("============================================================\n\n")

# Wave grid (0.5°)
nc_test_w <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
wave_lon <- ncvar_get(nc_test_w, "longitude")
wave_lat <- ncvar_get(nc_test_w, "latitude")
nc_close(nc_test_w)

wave_mask_corridor <- make_mask(wave_lon, wave_lat, CORRIDOR)
wave_mask_core     <- make_mask(wave_lon, wave_lat, CORE)

cat(sprintf("Wave grid: %d x %d (0.5°)\n", length(wave_lon), length(wave_lat)))
cat(sprintf("  Corridor mask: %d cells\n", sum(wave_mask_corridor)))
cat(sprintf("  Core mask:     %d cells\n", sum(wave_mask_core)))

# Atm grid (0.25°)
nc_test_a <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_atm_instant_2014.nc"))
atm_lon <- ncvar_get(nc_test_a, "longitude")
atm_lat <- ncvar_get(nc_test_a, "latitude")
nc_close(nc_test_a)

atm_mask_corridor <- make_mask(atm_lon, atm_lat, CORRIDOR)
atm_mask_core     <- make_mask(atm_lon, atm_lat, CORE)

cat(sprintf("Atm grid:  %d x %d (0.25°)\n", length(atm_lon), length(atm_lat)))
cat(sprintf("  Corridor mask: %d cells\n", sum(atm_mask_corridor)))
cat(sprintf("  Core mask:     %d cells\n", sum(atm_mask_core)))

# Check how many ocean cells are in each wave mask (NAs = land)
nc_chk <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
swh_chk <- ncvar_get(nc_chk, "swh", start = c(1, 1, 1), count = c(-1, -1, 1))
nc_close(nc_chk)

ocean_corridor <- sum(!is.na(swh_chk[wave_mask_corridor]))
ocean_core     <- sum(!is.na(swh_chk[wave_mask_core]))
cat(sprintf("\nOcean cells in corridor: %d (out of %d masked)\n",
    ocean_corridor, sum(wave_mask_corridor)))
cat(sprintf("Ocean cells in core:     %d (out of %d masked)\n",
    ocean_core, sum(wave_mask_core)))

# ============================================================
# 3. Compute daily spatial means from ERA5
# ============================================================
cat("\n--- Extracting daily weather ---\n")

all_days <- list()

for (yr in YEARS) {
  cat("  ", yr, ": ")

  f_wave <- file.path(ERA5_DIR, paste0("era5_daily_cmr_wave_", yr, ".nc"))
  f_atm  <- file.path(ERA5_DIR, paste0("era5_daily_cmr_atm_instant_", yr, ".nc"))
  f_gust <- file.path(ERA5_DIR, paste0("era5_daily_cmr_atm_gust_", yr, ".nc"))

  if (!file.exists(f_wave) || !file.exists(f_atm)) {
    cat("MISSING FILES\n"); next
  }

  # --- Wave variables (SWH, MWP) ---
  nc_w <- nc_open(f_wave)
  dates <- get_nc_dates(nc_w)
  n_days <- length(dates)
  swh_3d <- ncvar_get(nc_w, "swh")
  mwp_3d <- ncvar_get(nc_w, "mwp")
  nc_close(nc_w)

  swh_corr <- swh_core <- mwp_corr <- mwp_core <- numeric(n_days)
  for (t in seq_len(n_days)) {
    swh_corr[t] <- masked_mean(swh_3d[, , t], wave_mask_corridor)
    swh_core[t] <- masked_mean(swh_3d[, , t], wave_mask_core)
    mwp_corr[t] <- masked_mean(mwp_3d[, , t], wave_mask_corridor)
    mwp_core[t] <- masked_mean(mwp_3d[, , t], wave_mask_core)
  }

  # --- Atm variables (wind speed from u10, v10) ---
  nc_a <- nc_open(f_atm)
  u10_3d <- ncvar_get(nc_a, "u10")
  v10_3d <- ncvar_get(nc_a, "v10")
  nc_close(nc_a)

  wind_corr <- wind_core <- numeric(n_days)
  for (t in seq_len(n_days)) {
    ws <- sqrt(u10_3d[, , t]^2 + v10_3d[, , t]^2)
    wind_corr[t] <- masked_mean(ws, atm_mask_corridor)
    wind_core[t] <- masked_mean(ws, atm_mask_core)
  }

  # --- Gust (i10fg) ---
  gust_corr <- gust_core <- rep(NA_real_, n_days)
  if (file.exists(f_gust)) {
    nc_g <- nc_open(f_gust)
    i10fg_3d <- ncvar_get(nc_g, "i10fg")
    nc_close(nc_g)
    for (t in seq_len(n_days)) {
      gust_corr[t] <- masked_mean(i10fg_3d[, , t], atm_mask_corridor)
      gust_core[t] <- masked_mean(i10fg_3d[, , t], atm_mask_core)
    }
  }

  all_days[[as.character(yr)]] <- data.frame(
    date = dates,
    swh_mean  = swh_corr,  swh_core  = swh_core,
    wind_mean = wind_corr, wind_core = wind_core,
    gust_mean = gust_corr, gust_core = gust_core,
    mwp_mean  = mwp_corr,  mwp_core  = mwp_core
  )

  cat(n_days, "days\n")
}

daily <- bind_rows(all_days)
cat("\nTotal days:", nrow(daily), "\n")

# ============================================================
# 3b. Temporal lags (prior-day averages for transit time)
# ============================================================
# Boats take 1-7 days to cross (Camarena et al. 2020).
# We compute lagged weather to capture conditions during transit:
#   lag1:   yesterday (t-1)
#   prev3d: mean of t-1 to t-3
#   prev7d: mean of t-1 to t-7 (Camarena's main spec)

cat("\n--- Computing temporal lags ---\n")
daily <- daily %>% arrange(date)

add_lags <- function(df, var) {
  x <- df[[var]]
  df[[paste0(var, "_lag1")]]   <- dplyr::lag(x, 1)
  df[[paste0(var, "_prev3d")]] <- rowMeans(sapply(1:3, function(k) dplyr::lag(x, k)))
  df[[paste0(var, "_prev7d")]] <- rowMeans(sapply(1:7, function(k) dplyr::lag(x, k)))
  df
}

weather_vars <- c("swh_mean", "swh_core", "wind_mean", "wind_core",
                  "gust_mean", "gust_core", "mwp_mean", "mwp_core")
for (v in weather_vars) {
  daily <- add_lags(daily, v)
}

cat(sprintf("  Added lag1/prev3d/prev7d for %d variables (%d new columns)\n",
    length(weather_vars), length(weather_vars) * 3))
cat(sprintf("  NAs from lagging: lag1=%d, prev3d=%d, prev7d=%d\n",
    sum(is.na(daily$swh_mean_lag1)),
    sum(is.na(daily$swh_mean_prev3d)),
    sum(is.na(daily$swh_mean_prev7d))))

# ============================================================
# 4. Merge with incident counts
# ============================================================
cat("\n============================================================\n")
cat("Merging with IOM incident counts\n")
cat("============================================================\n\n")

df_incidents <- readRDS(IOM_PATH)

daily_counts <- df_incidents %>%
  group_by(date) %>%
  summarise(
    n_incidents = n(),
    n_dead = sum(dead, na.rm = TRUE),
    n_dead_missing = sum(dead_missing, na.rm = TRUE),
    n_fatal = sum(dead > 0 | dead_missing > 0, na.rm = TRUE),
    .groups = "drop"
  )

daily <- daily %>%
  left_join(daily_counts, by = "date") %>%
  mutate(
    n_incidents = ifelse(is.na(n_incidents), 0L, as.integer(n_incidents)),
    n_dead = ifelse(is.na(n_dead), 0, n_dead),
    n_dead_missing = ifelse(is.na(n_dead_missing), 0, n_dead_missing),
    n_fatal = ifelse(is.na(n_fatal), 0L, as.integer(n_fatal)),
    year = year(date),
    month = month(date),
    dow = wday(date),
    week = isoweek(date),
    post_mou = as.integer(date >= as.Date("2017-02-01")),
    month_fac = factor(month),
    year_fac = factor(year),
    week_year = paste0(year, "_w", sprintf("%02d", isoweek(date))),
    week_year_fac = factor(week_year)
  )

cat("Days with incidents:", sum(daily$n_incidents > 0), "\n")
cat("Days without incidents:", sum(daily$n_incidents == 0), "\n")
cat("Total incidents:", sum(daily$n_incidents), "\n")
cat("Total deaths:", sum(daily$n_dead), "\n")
cat("Total dead+missing:", sum(daily$n_dead_missing), "\n")
cat("Week-year cells:", length(unique(daily$week_year)), "\n")

# ============================================================
# 5. Time aggregation variables
# ============================================================
daily$half_year <- paste0(daily$year, ifelse(daily$month <= 6, "H1", "H2"))
daily$half_year_fac <- factor(daily$half_year)
daily$quarter <- paste0(daily$year, "Q", ceiling(daily$month / 3))
daily$quarter_fac <- factor(daily$quarter)

# ============================================================
# 6. Save
# ============================================================
cat("\n============================================================\n")
cat("Saving\n")
cat("============================================================\n\n")

saveRDS(daily, OUTPUT_PATH)
cat("Saved:", OUTPUT_PATH, "\n")
cat("Rows:", nrow(daily), ", Cols:", ncol(daily), "\n")

cat("\n--- Daily panel summary ---\n")
cat(sprintf("  Date range: %s to %s\n", min(daily$date), max(daily$date)))
cat(sprintf("  Corridor [%d,%d]x[%d,%d]:\n",
    CORRIDOR$lon_min, CORRIDOR$lon_max, CORRIDOR$lat_min, CORRIDOR$lat_max))
cat(sprintf("    SWH:  [%.2f, %.2f]\n", min(daily$swh_mean), max(daily$swh_mean)))
cat(sprintf("    Wind: [%.2f, %.2f]\n", min(daily$wind_mean), max(daily$wind_mean)))
cat(sprintf("    Gust: [%.2f, %.2f]\n", min(daily$gust_mean, na.rm=TRUE),
    max(daily$gust_mean, na.rm=TRUE)))
cat(sprintf("  Core [%d,%d]x[%d,%d]:\n",
    CORE$lon_min, CORE$lon_max, CORE$lat_min, CORE$lat_max))
cat(sprintf("    SWH:  [%.2f, %.2f]\n", min(daily$swh_core), max(daily$swh_core)))
cat(sprintf("    Wind: [%.2f, %.2f]\n", min(daily$wind_core), max(daily$wind_core)))
cat(sprintf("    Gust: [%.2f, %.2f]\n", min(daily$gust_core, na.rm=TRUE),
    max(daily$gust_core, na.rm=TRUE)))
cat(sprintf("  Correlation corridor vs core: SWH=%.3f, Wind=%.3f, Gust=%.3f\n",
    cor(daily$swh_mean, daily$swh_core),
    cor(daily$wind_mean, daily$wind_core),
    cor(daily$gust_mean, daily$gust_core, use = "complete.obs")))
cat(sprintf("  Dead+missing/day: mean=%.2f, max=%d, zeros=%d (%.1f%%)\n",
    mean(daily$n_dead_missing), max(daily$n_dead_missing),
    sum(daily$n_dead_missing == 0), 100 * mean(daily$n_dead_missing == 0)))
