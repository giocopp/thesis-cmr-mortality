# 02_build_event_data.R
# ====================
# Merge IOM MMP Central Mediterranean incidents with daily ERA5 weather.
#
# Input:
#   - data/processed/iom_mmp_incidents_2014_2025_reg.csv (cleaned IOM data)
#   - data/raw/era5/era5_daily_cmr_atm_instant_YYYY.nc (u10, v10, sst — 0.25°)
#   - data/raw/era5/era5_daily_cmr_atm_accum_YYYY.nc   (tp — 0.25°)
#   - data/raw/era5/era5_daily_cmr_wave_YYYY.nc         (swh, mwp — 0.5°)
#
# Output:
#   - data/processed/cmr_events_with_weather.RDS
#
# Steps:
#   1. Load and filter IOM data (CMR only, drop Atlantic outliers)
#   2. Load ERA5 grids (two resolutions: 0.25° atm, 0.5° wave)
#   3. For each incident: extract weather at nearest grid cell,
#      incident day + 7 preceding days
#   4. Compute derived variables (wind speed, grid ID, post-MoU)
#   5. Diagnostics

library(ncdf4)
library(dplyr)
library(lubridate)

BASE_DIR <- here::here()
IOM_PATH <- file.path(BASE_DIR, "data", "processed",
                       "iom_mmp_incidents_2014_2025_reg.csv")
ERA5_DIR <- file.path(BASE_DIR, "data", "raw", "era5")
OUTPUT_PATH <- file.path(BASE_DIR, "data", "processed",
                          "cmr_events_with_weather.RDS")

MOU_DATE <- as.Date("2017-02-01")

# ============================================================
# 1. Load and filter IOM data
# ============================================================
cat("============================================================\n")
cat("1. LOADING IOM MMP DATA\n")
cat("============================================================\n\n")

df <- read.csv(IOM_PATH, stringsAsFactors = FALSE)

# Filter to CMR only
df <- df %>% filter(Route == "Central Mediterranean")
cat("CMR incidents:", nrow(df), "\n")

# Drop Atlantic outliers (lon < 5)
n_atlantic <- sum(df$Longitude < 5, na.rm = TRUE)
df <- df %>% filter(Longitude >= 5)
cat("Dropped Atlantic outliers (lon < 5):", n_atlantic, "\n")

# Flag Egypt coordinates (lon > 25) for inspection
n_egypt <- sum(df$Longitude > 25, na.rm = TRUE)
df$flag_egypt <- df$Longitude > 25
cat("Egypt-area incidents flagged (lon > 25):", n_egypt, "\n")

# Parse dates
df$date <- as.Date(df$incident_date_clean)
df <- df %>% filter(!is.na(date))
cat("With valid dates:", nrow(df), "\n")

# Key columns
df <- df %>%
  mutate(
    dead = as.numeric(No..dead),
    dead = ifelse(is.na(dead), 0, dead),
    dead_missing = as.numeric(No..dead.missing),
    dead_missing = ifelse(is.na(dead_missing), 0, dead_missing),
    lat = as.numeric(Latitude),
    lon = as.numeric(Longitude),
    year = year(date),
    month = month(date),
    post_mou = as.integer(date >= MOU_DATE)
  ) %>%
  filter(!is.na(lat) & !is.na(lon))

cat("Final sample:", nrow(df), "incidents\n")
cat("  Pre-MoU:", sum(df$post_mou == 0), "\n")
cat("  Post-MoU:", sum(df$post_mou == 1), "\n\n")

# ============================================================
# 2. Load ERA5 grids
# ============================================================
cat("============================================================\n")
cat("2. LOADING ERA5 GRIDS\n")
cat("============================================================\n\n")

# Five file types per year, two grid resolutions
FILE_TYPES <- c("atm_instant", "atm_accum", "atm_gust", "wave", "wave_hmax")

# Which variables live in which file type
VAR_FILE_MAP <- list(
  u10   = "atm_instant", v10 = "atm_instant", sst = "atm_instant",
  tp    = "atm_accum",
  i10fg = "atm_gust",
  swh   = "wave",         mwp = "wave",
  hmax  = "wave_hmax"
)
WEATHER_VARS <- names(VAR_FILE_MAP)

# Check files exist for at least one year
test_year <- 2014
for (ft in FILE_TYPES) {
  f <- file.path(ERA5_DIR, paste0("era5_daily_cmr_", ft, "_", test_year, ".nc"))
  if (!file.exists(f))
    stop("Missing ERA5 file: ", f, "\n  Run analysis/python/download_era5_daily.py first.")
}

# Read the two grids: atm (0.25°) and wave (0.5°)
nc_atm <- nc_open(file.path(ERA5_DIR, paste0("era5_daily_cmr_atm_instant_", test_year, ".nc")))
atm_lon <- ncvar_get(nc_atm, "longitude")
atm_lat <- ncvar_get(nc_atm, "latitude")
nc_close(nc_atm)

nc_wav <- nc_open(file.path(ERA5_DIR, paste0("era5_daily_cmr_wave_", test_year, ".nc")))
wave_lon <- ncvar_get(nc_wav, "longitude")
wave_lat <- ncvar_get(nc_wav, "latitude")
nc_close(nc_wav)

cat(sprintf("Atm grid (0.25°): %d x %d, lon [%.1f, %.1f], lat [%.1f, %.1f]\n",
    length(atm_lon), length(atm_lat), min(atm_lon), max(atm_lon),
    min(atm_lat), max(atm_lat)))
cat(sprintf("Wave grid (0.5°): %d x %d, lon [%.1f, %.1f], lat [%.1f, %.1f]\n",
    length(wave_lon), length(wave_lat), min(wave_lon), max(wave_lon),
    min(wave_lat), max(wave_lat)))

# ============================================================
# 3. Helpers
# ============================================================

find_nearest <- function(target_lon, target_lat, grid_lon, grid_lat) {
  ix <- which.min(abs(grid_lon - target_lon))
  iy <- which.min(abs(grid_lat - target_lat))
  list(ix = ix, iy = iy,
       grid_lon = grid_lon[ix], grid_lat = grid_lat[iy])
}

extract_weather <- function(nc, ix, iy, time_idx, var_name) {
  tryCatch(
    ncvar_get(nc, var_name, start = c(ix, iy, time_idx), count = c(1, 1, 1)),
    error = function(e) NA_real_
  )
}

get_nc_time <- function(nc) {
  time_name <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  if (is.na(time_name)) stop("No time dimension found")
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

# NetCDF cache: keyed by "year_filetype" (e.g. "2017_wave")
nc_cache <- list()

open_nc <- function(year, file_type) {
  key <- paste0(year, "_", file_type)
  if (is.null(nc_cache[[key]])) {
    f <- file.path(ERA5_DIR, paste0("era5_daily_cmr_", file_type, "_", year, ".nc"))
    if (!file.exists(f)) return(NULL)
    nc_cache[[key]] <<- nc_open(f)
  }
  nc_cache[[key]]
}

# ============================================================
# 4. Main extraction loop
# ============================================================
cat("\n============================================================\n")
cat("3. EXTRACTING WEATHER FOR EACH INCIDENT\n")
cat("============================================================\n\n")

# Pre-compute nearest grid cells on BOTH grids for all incidents
df$atm_ix  <- NA_integer_; df$atm_iy  <- NA_integer_
df$wave_ix <- NA_integer_; df$wave_iy <- NA_integer_
df$grid_lon <- NA_real_;   df$grid_lat <- NA_real_

for (i in seq_len(nrow(df))) {
  nn_a <- find_nearest(df$lon[i], df$lat[i], atm_lon, atm_lat)
  nn_w <- find_nearest(df$lon[i], df$lat[i], wave_lon, wave_lat)
  df$atm_ix[i]  <- nn_a$ix;  df$atm_iy[i]  <- nn_a$iy
  df$wave_ix[i] <- nn_w$ix;  df$wave_iy[i] <- nn_w$iy
  # Primary grid ID from 0.25° atm grid
  df$grid_lon[i] <- nn_a$grid_lon
  df$grid_lat[i] <- nn_a$grid_lat
}

df$grid_id <- paste0(sprintf("%.2f", df$grid_lat), "_",
                      sprintf("%.2f", df$grid_lon))
cat("Unique grid cells (0.25°):", length(unique(df$grid_id)), "\n")

# Initialize weather columns (lag 0 through 7)
for (v in WEATHER_VARS) {
  for (lag in 0:7) df[[paste0(v, "_lag", lag)]] <- NA_real_
}

# Extract weather per year
years_in_data <- sort(unique(df$year))
cat("Extracting weather for years:", paste(years_in_data, collapse = ", "), "\n")

n_missing <- 0

# Pre-cache time vectors to avoid repeated parsing
time_cache <- list()
get_dates <- function(year, file_type) {
  key <- paste0(year, "_", file_type)
  if (is.null(time_cache[[key]])) {
    nc <- open_nc(year, file_type)
    if (is.null(nc)) return(NULL)
    time_cache[[key]] <<- get_nc_time(nc)
  }
  time_cache[[key]]
}

for (yr in years_in_data) {
  cat("  ", yr, ": ")

  # Open all 3 file types for this year (and previous year for lags)
  for (ft in FILE_TYPES) {
    open_nc(yr, ft)
    open_nc(yr - 1, ft)  # may be NULL for 2014
  }

  idx_yr <- which(df$year == yr)
  cat(length(idx_yr), "incidents... ")

  for (i in idx_yr) {
    incident_date <- df$date[i]

    for (lag in 0:7) {
      target_date <- incident_date - lag
      target_year <- year(target_date)

      # Check year availability
      if (target_year != yr && target_year != yr - 1) {
        n_missing <- n_missing + 1
        next
      }

      for (v in WEATHER_VARS) {
        ft <- VAR_FILE_MAP[[v]]
        nc_use <- open_nc(target_year, ft)
        if (is.null(nc_use)) { n_missing <- n_missing + 1; next }

        dates_use <- get_dates(target_year, ft)
        t_idx <- which(dates_use == target_date)
        if (length(t_idx) == 0) { n_missing <- n_missing + 1; next }

        # Use correct grid indices for this variable
        # wave and wave_hmax both use the 0.5-degree wave grid
        if (ft %in% c("wave", "wave_hmax")) {
          ix <- df$wave_ix[i]; iy <- df$wave_iy[i]
        } else {
          ix <- df$atm_ix[i]; iy <- df$atm_iy[i]
        }

        df[[paste0(v, "_lag", lag)]][i] <- extract_weather(nc_use, ix, iy, t_idx[1], v)
      }
    }
  }
  cat("done\n")
}

# Close all cached NetCDF files
for (nc in nc_cache) nc_close(nc)

cat("\nMissing extractions:", n_missing, "\n")

# ============================================================
# 6. Compute derived variables
# ============================================================
cat("\n============================================================\n")
cat("4. COMPUTING DERIVED VARIABLES\n")
cat("============================================================\n\n")

# Wind speed at each lag
for (lag in 0:7) {
  u_col <- paste0("u10_lag", lag)
  v_col <- paste0("v10_lag", lag)
  ws_col <- paste0("wind_speed_lag", lag)
  df[[ws_col]] <- sqrt(df[[u_col]]^2 + df[[v_col]]^2)
}

# Helper: safe row max/min that returns NA (not -Inf/Inf) when all values are NA
safe_row_max <- function(mat) {
  apply(mat, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else max(x)
  })
}

# SST: convert Kelvin to Celsius BEFORE computing summaries
for (lag in 0:7) {
  col <- paste0("sst_lag", lag)
  df[[col]] <- df[[col]] - 273.15
}

# --- Day-0 values (incident day) ---
df$swh_day0   <- df$swh_lag0
df$wind_day0  <- df$wind_speed_lag0
df$sst_day0   <- df$sst_lag0
df$i10fg_day0 <- df$i10fg_lag0
df$hmax_day0  <- df$hmax_lag0
df$mwp_day0   <- df$mwp_lag0
df$tp_day0    <- df$tp_lag0

# --- 3-day summaries (lag 0-2): transit/distress period ---
# SWH
df$swh_mean_3d <- rowMeans(df[, paste0("swh_lag", 0:2)], na.rm = TRUE)
df$swh_mean_3d[rowSums(!is.na(df[, paste0("swh_lag", 0:2)])) == 0] <- NA
df$swh_max_3d  <- safe_row_max(df[, paste0("swh_lag", 0:2)])

# Wind speed
df$wind_mean_3d <- rowMeans(df[, paste0("wind_speed_lag", 0:2)], na.rm = TRUE)
df$wind_mean_3d[rowSums(!is.na(df[, paste0("wind_speed_lag", 0:2)])) == 0] <- NA
df$wind_max_3d  <- safe_row_max(df[, paste0("wind_speed_lag", 0:2)])

# Wind gust
df$i10fg_mean_3d <- rowMeans(df[, paste0("i10fg_lag", 0:2)], na.rm = TRUE)
df$i10fg_mean_3d[rowSums(!is.na(df[, paste0("i10fg_lag", 0:2)])) == 0] <- NA
df$i10fg_max_3d  <- safe_row_max(df[, paste0("i10fg_lag", 0:2)])

# Mean wave period
df$mwp_mean_3d <- rowMeans(df[, paste0("mwp_lag", 0:2)], na.rm = TRUE)
df$mwp_mean_3d[rowSums(!is.na(df[, paste0("mwp_lag", 0:2)])) == 0] <- NA

# SST (hypothermia relevance — mean is appropriate)
df$sst_mean_3d <- rowMeans(df[, paste0("sst_lag", 0:2)], na.rm = TRUE)
df$sst_mean_3d[rowSums(!is.na(df[, paste0("sst_lag", 0:2)])) == 0] <- NA

# Precipitation (sum over 3 days)
df$tp_sum_3d <- rowSums(df[, paste0("tp_lag", 0:2)], na.rm = TRUE)
df$tp_sum_3d[rowSums(!is.na(df[, paste0("tp_lag", 0:2)])) == 0] <- NA

# Hmax
df$hmax_mean_3d <- rowMeans(df[, paste0("hmax_lag", 0:2)], na.rm = TRUE)
df$hmax_mean_3d[rowSums(!is.na(df[, paste0("hmax_lag", 0:2)])) == 0] <- NA
df$hmax_max_3d  <- safe_row_max(df[, paste0("hmax_lag", 0:2)])

# --- 7-day summaries (lag 0-7): full transit window ---
# SWH
df$swh_mean_7d <- rowMeans(df[, paste0("swh_lag", 0:7)], na.rm = TRUE)
df$swh_mean_7d[rowSums(!is.na(df[, paste0("swh_lag", 0:7)])) == 0] <- NA
df$swh_max_7d  <- safe_row_max(df[, paste0("swh_lag", 0:7)])

# Wind speed
df$wind_mean_7d <- rowMeans(df[, paste0("wind_speed_lag", 0:7)], na.rm = TRUE)
df$wind_mean_7d[rowSums(!is.na(df[, paste0("wind_speed_lag", 0:7)])) == 0] <- NA
df$wind_max_7d  <- safe_row_max(df[, paste0("wind_speed_lag", 0:7)])

# Wind gust
df$i10fg_mean_7d <- rowMeans(df[, paste0("i10fg_lag", 0:7)], na.rm = TRUE)
df$i10fg_mean_7d[rowSums(!is.na(df[, paste0("i10fg_lag", 0:7)])) == 0] <- NA
df$i10fg_max_7d  <- safe_row_max(df[, paste0("i10fg_lag", 0:7)])

# Mean wave period
df$mwp_mean_7d <- rowMeans(df[, paste0("mwp_lag", 0:7)], na.rm = TRUE)
df$mwp_mean_7d[rowSums(!is.na(df[, paste0("mwp_lag", 0:7)])) == 0] <- NA

# SST
df$sst_mean_7d <- rowMeans(df[, paste0("sst_lag", 0:7)], na.rm = TRUE)
df$sst_mean_7d[rowSums(!is.na(df[, paste0("sst_lag", 0:7)])) == 0] <- NA

# Precipitation
df$tp_sum_7d <- rowSums(df[, paste0("tp_lag", 0:7)], na.rm = TRUE)
df$tp_sum_7d[rowSums(!is.na(df[, paste0("tp_lag", 0:7)])) == 0] <- NA

# Hmax
df$hmax_mean_7d <- rowMeans(df[, paste0("hmax_lag", 0:7)], na.rm = TRUE)
df$hmax_mean_7d[rowSums(!is.na(df[, paste0("hmax_lag", 0:7)])) == 0] <- NA
df$hmax_max_7d  <- safe_row_max(df[, paste0("hmax_lag", 0:7)])

# Calendar variables for FE
df$month_fac <- factor(df$month)
df$year_fac <- factor(df$year)
df$grid_fac <- factor(df$grid_id)

cat("Derived variables computed.\n\n")

# ============================================================
# 7. Diagnostics
# ============================================================
cat("============================================================\n")
cat("5. DIAGNOSTICS\n")
cat("============================================================\n\n")

# Weather ranges across all windows
report_var <- function(label, day0, mean3, max3, mean7, max7) {
  cat(sprintf("  %-6s day0: [%5.2f, %5.2f] mean=%.2f NAs=%d\n",
      label, min(day0, na.rm=TRUE), max(day0, na.rm=TRUE),
      mean(day0, na.rm=TRUE), sum(is.na(day0))))
  cat(sprintf("  %-6s 3d mean: [%5.2f, %5.2f] NAs=%d | 3d max: [%5.2f, %5.2f] NAs=%d\n",
      "", min(mean3, na.rm=TRUE), max(mean3, na.rm=TRUE), sum(is.na(mean3)),
      min(max3, na.rm=TRUE), max(max3, na.rm=TRUE), sum(is.na(max3))))
  cat(sprintf("  %-6s 7d mean: [%5.2f, %5.2f] NAs=%d | 7d max: [%5.2f, %5.2f] NAs=%d\n",
      "", min(mean7, na.rm=TRUE), max(mean7, na.rm=TRUE), sum(is.na(mean7)),
      min(max7, na.rm=TRUE), max(max7, na.rm=TRUE), sum(is.na(max7))))
}

cat("--- Weather ranges: day0 / 3-day / 7-day ---\n")
report_var("SWH", df$swh_day0, df$swh_mean_3d, df$swh_max_3d,
           df$swh_mean_7d, df$swh_max_7d)
report_var("Wind", df$wind_day0, df$wind_mean_3d, df$wind_max_3d,
           df$wind_mean_7d, df$wind_max_7d)
report_var("Gust", df$i10fg_day0, df$i10fg_mean_3d, df$i10fg_max_3d,
           df$i10fg_mean_7d, df$i10fg_max_7d)
report_var("Hmax", df$hmax_day0, df$hmax_mean_3d, df$hmax_max_3d,
           df$hmax_mean_7d, df$hmax_max_7d)

cat(sprintf("  SST    day0: [%5.2f, %5.2f] NAs=%d | 3d: NAs=%d | 7d: NAs=%d\n",
    min(df$sst_day0, na.rm=TRUE), max(df$sst_day0, na.rm=TRUE),
    sum(is.na(df$sst_day0)), sum(is.na(df$sst_mean_3d)), sum(is.na(df$sst_mean_7d))))
cat(sprintf("  MWP    day0: [%5.2f, %5.2f] NAs=%d | 3d: NAs=%d | 7d: NAs=%d\n",
    min(df$mwp_day0, na.rm=TRUE), max(df$mwp_day0, na.rm=TRUE),
    sum(is.na(df$mwp_day0)), sum(is.na(df$mwp_mean_3d)), sum(is.na(df$mwp_mean_7d))))

cat(sprintf("\n  Correlations (day0 vs 3d mean / 7d mean):\n"))
cat(sprintf("    SWH:  %.3f / %.3f\n",
    cor(df$swh_day0, df$swh_mean_3d, use="complete"),
    cor(df$swh_day0, df$swh_mean_7d, use="complete")))
cat(sprintf("    Wind: %.3f / %.3f\n",
    cor(df$wind_day0, df$wind_mean_3d, use="complete"),
    cor(df$wind_day0, df$wind_mean_7d, use="complete")))
cat(sprintf("    Gust: %.3f / %.3f\n",
    cor(df$i10fg_day0, df$i10fg_mean_3d, use="complete"),
    cor(df$i10fg_day0, df$i10fg_mean_7d, use="complete")))

# Beaufort classification (daily wind speed should reach higher values)
bf_breaks <- c(0, 0.3, 1.6, 3.4, 5.5, 8.0, 10.8, 13.9, 17.2, Inf)
bf_labels <- 0:8
df$beaufort <- as.integer(as.character(
  cut(df$wind_day0, breaks = bf_breaks, labels = bf_labels,
      right = FALSE, include.lowest = TRUE)))

cat("\n--- Beaufort scale (incident day) ---\n")
print(table(df$beaufort, useNA = "ifany"))

# WMO sea state
wmo_breaks <- c(0, 0.001, 0.1, 0.5, 1.25, 2.5, 4.0, 6.0, Inf)
wmo_labels <- 0:7
df$wmo_state <- as.integer(as.character(
  cut(df$swh_day0, breaks = wmo_breaks, labels = wmo_labels,
      right = FALSE, include.lowest = TRUE)))

cat("\n--- WMO sea state (incident day) ---\n")
print(table(df$wmo_state, useNA = "ifany"))

# Grid cell distribution
cell_counts <- table(df$grid_id)
cat("\n--- Grid cell distribution ---\n")
cat(sprintf("  Unique cells: %d\n", length(cell_counts)))
cat(sprintf("  Singletons: %d (%.1f%%)\n",
    sum(cell_counts == 1), 100 * mean(cell_counts == 1)))
cat(sprintf("  Incidents per cell: median=%d, mean=%.1f, max=%d\n",
    median(cell_counts), mean(cell_counts), max(cell_counts)))

# Pre/post split
cat("\n--- Pre/post MoU ---\n")
cat(sprintf("  Pre:  %d incidents, mean dead=%.1f\n",
    sum(df$post_mou == 0), mean(df$dead[df$post_mou == 0])))
cat(sprintf("  Post: %d incidents, mean dead=%.1f\n",
    sum(df$post_mou == 1), mean(df$dead[df$post_mou == 1])))

# ============================================================
# 8. Save
# ============================================================
cat("\n============================================================\n")
cat("6. SAVING\n")
cat("============================================================\n\n")

# Select columns for output
out_cols <- c(
  # IOM identifiers
  "Main.ID", "date", "year", "month",
  # Location
  "lat", "lon", "grid_id", "grid_lon", "grid_lat",
  "flag_egypt", "Location.of.death",
  # Outcome
  "dead", "dead_missing", "No..survivors",
  # Covariates
  "Cause.of.death..category.", "Source.Quality",
  # Treatment
  "post_mou",
  # Weather: incident day
  "swh_day0", "wind_day0", "sst_day0", "i10fg_day0", "hmax_day0",
  "mwp_day0", "tp_day0",
  # Weather: 3-day summaries (lag 0-2, transit/distress window)
  "swh_mean_3d", "swh_max_3d", "wind_mean_3d", "wind_max_3d",
  "i10fg_mean_3d", "i10fg_max_3d", "mwp_mean_3d", "sst_mean_3d",
  "tp_sum_3d", "hmax_mean_3d", "hmax_max_3d",
  # Weather: 7-day summaries (lag 0-7, full transit window)
  "swh_mean_7d", "swh_max_7d", "wind_mean_7d", "wind_max_7d",
  "i10fg_mean_7d", "i10fg_max_7d", "mwp_mean_7d", "sst_mean_7d",
  "tp_sum_7d", "hmax_mean_7d", "hmax_max_7d",
  # Weather: all lags (for flexibility)
  paste0("swh_lag", 0:7),
  paste0("wind_speed_lag", 0:7),
  paste0("sst_lag", 0:7),
  paste0("mwp_lag", 0:7),
  paste0("tp_lag", 0:7),
  paste0("i10fg_lag", 0:7),
  paste0("hmax_lag", 0:7),
  # Beaufort/WMO
  "beaufort", "wmo_state",
  # FE variables
  "month_fac", "year_fac", "grid_fac"
)

# Keep only columns that exist
out_cols <- out_cols[out_cols %in% names(df)]

df_out <- df[, out_cols]
saveRDS(df_out, OUTPUT_PATH)
cat("Saved:", OUTPUT_PATH, "\n")
cat("Rows:", nrow(df_out), ", Cols:", ncol(df_out), "\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
