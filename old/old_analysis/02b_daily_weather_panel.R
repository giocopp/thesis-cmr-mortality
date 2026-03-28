# 02b_daily_weather_panel.R
# ========================
# Compute daily CMR spatial-mean weather from ERA5 grids.
# This creates a daily panel with weather for ALL days (not just incident days),
# which is needed for the daily count model.
#
# GEOGRAPHIC RESTRICTION:
#   Weather is averaged only over ocean cells within the CMR core corridor,
#   defined as [10.5, 15.5] x [32.3, 36.2] — the Libya/Tunisia to
#   Lampedusa crossing channel.
#
# Input:  data/raw/era5/era5_daily_cmr_{atm_instant,atm_accum,atm_gust,wave}_YYYY.nc
# Output: data/processed/cmr_daily_weather_panel.RDS
#
# Variables computed:
#   - swh_core:  spatial mean SWH (m) over core corridor ocean cells
#   - wind_core: spatial mean wind speed (m/s)
#   - gust_core: spatial mean wind gust (m/s)
#   - mwp_core:  spatial mean wave period (s)

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
# Core CMR corridor: the full Libya/Tunisia to Lampedusa/Sicily crossing channel
CORE <- list(lon_min = 10.5, lon_max = 15.5, lat_min = 32.3, lat_max = 36.2)

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

wave_mask_core <- make_mask(wave_lon, wave_lat, CORE)

cat(sprintf("Wave grid: %d x %d (0.5°)\n", length(wave_lon), length(wave_lat)))
cat(sprintf("  Core mask: %d cells\n", sum(wave_mask_core)))

# Atm grid (0.25°)
nc_test_a <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_atm_instant_2014.nc"))
atm_lon <- ncvar_get(nc_test_a, "longitude")
atm_lat <- ncvar_get(nc_test_a, "latitude")
nc_close(nc_test_a)

atm_mask_core <- make_mask(atm_lon, atm_lat, CORE)

cat(sprintf("Atm grid:  %d x %d (0.25°)\n", length(atm_lon), length(atm_lat)))
cat(sprintf("  Core mask: %d cells\n", sum(atm_mask_core)))

# Check how many ocean cells are in the wave mask (NAs = land)
nc_chk <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
swh_chk <- ncvar_get(nc_chk, "swh", start = c(1, 1, 1), count = c(-1, -1, 1))
nc_close(nc_chk)

ocean_core <- sum(!is.na(swh_chk[wave_mask_core]))
cat(sprintf("\nOcean cells in core: %d (out of %d masked)\n",
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

  swh_core <- mwp_core <- numeric(n_days)
  for (t in seq_len(n_days)) {
    swh_core[t] <- masked_mean(swh_3d[, , t], wave_mask_core)
    mwp_core[t] <- masked_mean(mwp_3d[, , t], wave_mask_core)
  }

  # --- Atm variables (wind speed from u10, v10) ---
  nc_a <- nc_open(f_atm)
  u10_3d <- ncvar_get(nc_a, "u10")
  v10_3d <- ncvar_get(nc_a, "v10")
  nc_close(nc_a)

  wind_core <- numeric(n_days)
  for (t in seq_len(n_days)) {
    ws <- sqrt(u10_3d[, , t]^2 + v10_3d[, , t]^2)
    wind_core[t] <- masked_mean(ws, atm_mask_core)
  }

  # --- Gust (i10fg) ---
  gust_core <- rep(NA_real_, n_days)
  if (file.exists(f_gust)) {
    nc_g <- nc_open(f_gust)
    i10fg_3d <- ncvar_get(nc_g, "i10fg")
    nc_close(nc_g)
    for (t in seq_len(n_days)) {
      gust_core[t] <- masked_mean(i10fg_3d[, , t], atm_mask_core)
    }
  }

  all_days[[as.character(yr)]] <- data.frame(
    date = dates,
    swh_core  = swh_core,
    wind_core = wind_core,
    gust_core = gust_core,
    mwp_core  = mwp_core
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
  df[[paste0(var, "_lag2")]]   <- dplyr::lag(x, 2)
  df[[paste0(var, "_lag3")]]   <- dplyr::lag(x, 3)
  df[[paste0(var, "_lag7")]]   <- dplyr::lag(x, 7)
  df[[paste0(var, "_prev3d")]] <- rowMeans(sapply(1:3, function(k) dplyr::lag(x, k)))
  df[[paste0(var, "_prev7d")]] <- rowMeans(sapply(1:7, function(k) dplyr::lag(x, k)))
  df
}

weather_vars <- c("swh_core", "wind_core", "gust_core", "mwp_core")
for (v in weather_vars) {
  daily <- add_lags(daily, v)
}

cat(sprintf("  Added lag1-3/prev3d/prev7d for %d variables (%d new columns)\n",
    length(weather_vars), length(weather_vars) * 5))
cat(sprintf("  NAs from lagging: lag1=%d, lag2=%d, lag3=%d, prev7d=%d\n",
    sum(is.na(daily$swh_core_lag1)),
    sum(is.na(daily$swh_core_lag2)),
    sum(is.na(daily$swh_core_lag3)),
    sum(is.na(daily$swh_core_prev7d))))

# ============================================================
# 4. Merge with incident counts
# ============================================================
cat("\n============================================================\n")
cat("Merging with IOM incident counts\n")
cat("============================================================\n\n")

df_incidents <- readRDS(IOM_PATH)

# --- Geographic restriction to core corridor ---
# Weather is the spatial mean over [11,15]x[32,36]. Only incidents
# within this box are matched to the correct weather. Incidents off
# Egypt, near Greece/Italy arrival points, etc. experience different
# sea conditions than what SWH_core measures.
n_before_geo <- nrow(df_incidents)
df_incidents <- df_incidents %>%
  filter(lon >= CORE$lon_min & lon <= CORE$lon_max &
         lat >= CORE$lat_min & lat <= CORE$lat_max)

cat(sprintf("Geographic restriction to core [%.1f,%.1f]x[%.1f,%.1f]:\n",
    CORE$lon_min, CORE$lon_max, CORE$lat_min, CORE$lat_max))
cat(sprintf("  Before: %d incidents\n", n_before_geo))
cat(sprintf("  After:  %d incidents (dropped %d outside core)\n\n",
    nrow(df_incidents), n_before_geo - nrow(df_incidents)))

# --- Cause-of-death filtering ---
# Primary: drowning + suspected drowning (Drowning + Mixed/unknown).
#   On the CMR, "Mixed or unknown" are overwhelmingly bodies found at sea
#   or on boats — sea crossing deaths where exact cause was not determined.
#   Excludes clearly non-maritime causes: vehicle accidents, sickness, violence.
# Robustness: drowning/suspected + dead_missing > 1 (drops body-found records
#   where date = recovery date, not crossing date, so weather is noise).

sea_causes <- c("Drowning", "Mixed or unknown")
df_sea <- df_incidents %>%
  filter(Cause.of.death..category. %in% sea_causes)
df_sea_ship <- df_sea %>%
  filter(dead_missing > 1)

cat(sprintf("Cause-of-death filtering (within core corridor):\n"))
cat(sprintf("  All incidents in core:              %d\n", nrow(df_incidents)))
cat(sprintf("  Drowning + suspected drowning:      %d (dropped %d non-maritime)\n",
    nrow(df_sea), nrow(df_incidents) - nrow(df_sea)))
cat(sprintf("  Above + dead/missing > 1:           %d (dropped %d single-body records)\n",
    nrow(df_sea_ship), nrow(df_sea) - nrow(df_sea_ship)))

# Helper: aggregate incidents to daily counts
aggregate_daily <- function(inc_df, prefix = "") {
  counts <- inc_df %>%
    group_by(date) %>%
    summarise(
      n_incidents = n(),
      n_dead = sum(dead, na.rm = TRUE),
      n_dead_missing = sum(dead_missing, na.rm = TRUE),
      n_fatal = sum(dead > 0 | dead_missing > 0, na.rm = TRUE),
      .groups = "drop"
    )
  if (nchar(prefix) > 0) {
    names(counts)[-1] <- paste0(prefix, "_", names(counts)[-1])
  }
  counts
}

# Primary outcome: drowning + suspected drowning
daily_primary <- aggregate_daily(df_sea)

# Robustness outcome: above + dead/missing > 1
daily_ship <- aggregate_daily(df_sea_ship, prefix = "ship")
daily_all  <- aggregate_daily(df_incidents, prefix = "all")

daily <- daily %>%
  left_join(daily_primary, by = "date") %>%
  left_join(daily_ship, by = "date") %>%
  left_join(daily_all, by = "date") %>%
  mutate(
    # Primary: drowning only
    n_incidents = ifelse(is.na(n_incidents), 0L, as.integer(n_incidents)),
    n_dead = ifelse(is.na(n_dead), 0, n_dead),
    n_dead_missing = ifelse(is.na(n_dead_missing), 0, n_dead_missing),
    n_fatal = ifelse(is.na(n_fatal), 0L, as.integer(n_fatal)),
    # Robustness 1: drowning + dead/missing > 1
    ship_n_incidents = ifelse(is.na(ship_n_incidents), 0L, as.integer(ship_n_incidents)),
    ship_n_dead_missing = ifelse(is.na(ship_n_dead_missing), 0, ship_n_dead_missing),
    # Robustness 2: all incidents
    all_n_incidents = ifelse(is.na(all_n_incidents), 0L, as.integer(all_n_incidents)),
    all_n_dead_missing = ifelse(is.na(all_n_dead_missing), 0, all_n_dead_missing),
    year = year(date),
    month = month(date),
    dow = wday(date),
    week = isoweek(date),
    post_mou = as.integer(date >= as.Date("2017-07-01")),
    month_fac = factor(month),
    year_fac = factor(year),
    week_year = paste0(year, "_w", sprintf("%02d", isoweek(date))),
    week_year_fac = factor(week_year),
    month_year = paste0(year, "_m", sprintf("%02d", month)),
    month_year_fac = factor(month_year)
  )

cat("\n--- Primary outcome: drowning + suspected drowning ---\n")
cat("Days with incidents:", sum(daily$n_incidents > 0), "\n")
cat("Days without incidents:", sum(daily$n_incidents == 0), "\n")
cat("Total incidents:", sum(daily$n_incidents), "\n")
cat("Total dead+missing:", sum(daily$n_dead_missing), "\n")
cat(sprintf("Var/mean: %.1f\n", var(daily$n_dead_missing) / mean(daily$n_dead_missing)))

cat("\n--- Robustness: drowning/suspected + dead/missing > 1 ---\n")
cat("Days with incidents:", sum(daily$ship_n_incidents > 0), "\n")
cat("Total incidents:", sum(daily$ship_n_incidents), "\n")
cat("Total dead+missing:", sum(daily$ship_n_dead_missing), "\n")

cat("\n--- Robustness 2: all incidents ---\n")
cat("Days with incidents:", sum(daily$all_n_incidents > 0), "\n")
cat("Total incidents:", sum(daily$all_n_incidents), "\n")
cat("Total dead+missing:", sum(daily$all_n_dead_missing), "\n")

cat("\nWeek-year cells:", length(unique(daily$week_year)), "\n")

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
cat(sprintf("  Core corridor [%.1f,%.1f]x[%.1f,%.1f]:\n",
    CORE$lon_min, CORE$lon_max, CORE$lat_min, CORE$lat_max))
cat(sprintf("    SWH:  [%.2f, %.2f]\n", min(daily$swh_core), max(daily$swh_core)))
cat(sprintf("    Wind: [%.2f, %.2f]\n", min(daily$wind_core), max(daily$wind_core)))
cat(sprintf("    Gust: [%.2f, %.2f]\n", min(daily$gust_core, na.rm=TRUE),
    max(daily$gust_core, na.rm=TRUE)))
cat(sprintf("  Drowning/suspected dead+missing/day: mean=%.2f, max=%d, zeros=%d (%.1f%%)\n",
    mean(daily$n_dead_missing), max(daily$n_dead_missing),
    sum(daily$n_dead_missing == 0), 100 * mean(daily$n_dead_missing == 0)))
cat(sprintf("  All dead+missing/day:      mean=%.2f, max=%d, zeros=%d (%.1f%%)\n",
    mean(daily$all_n_dead_missing), max(daily$all_n_dead_missing),
    sum(daily$all_n_dead_missing == 0), 100 * mean(daily$all_n_dead_missing == 0)))
