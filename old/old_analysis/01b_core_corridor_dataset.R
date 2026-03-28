# 01b_core_corridor_dataset.R
# ===========================
# Build the clean incident-level dataset for the core CMR corridor.
#
# This script:
#   1. Loads IOM MMP incident data (all routes, Incident type only)
#   2. Filters to Central Mediterranean route
#   3. Restricts to the core corridor [11,15]x[32,36] where weather
#      is measured — the Libya-to-Lampedusa crossing channel
#   4. Preserves ALL IOM MMP fields (no column dropping)
#   5. Merges ERA5 weather at each incident's nearest grid cell
#   6. Produces a map of the core corridor with incident locations
#   7. Saves a clean dataset
#
# Core corridor: lon [10.5, 15.5], lat [32.3, 36.2]
#   This box captures the Central Mediterranean Route crossing
#   channel from the Libyan/Tunisian coast to Lampedusa, where
#   the spatial-mean weather is computed.
#
# Input:
#   data/processed/iom_mmp_incidents_2014_2025_reg.csv
#   data/raw/era5/era5_daily_cmr_*.nc
#
# Output:
#   data/processed/core_corridor_incidents.RDS
#   data/processed/core_corridor_incidents.csv
#   output/figures/core_corridor_map.pdf
#   output/figures/core_corridor_map.png

library(dplyr)
library(lubridate)
library(ncdf4)
library(ggplot2)
library(sf)
library(rnaturalearth)

BASE_DIR <- here::here()

# ============================================================
# 0. Define the core corridor
# ============================================================
CORE <- list(lon_min = 10.5, lon_max = 15.5, lat_min = 32.3, lat_max = 36.2)

MOU_DATE <- as.Date("2017-07-01")

# ============================================================
# 1. Load IOM MMP data
# ============================================================
cat("============================================================\n")
cat("1. LOADING IOM MMP DATA\n")
cat("============================================================\n\n")

iom_path <- file.path(BASE_DIR, "data", "processed",
                       "iom_mmp_incidents_2014_2025_reg.csv")
raw <- read.csv(iom_path, stringsAsFactors = FALSE)
cat("Total IOM MMP records:", nrow(raw), "\n")

# ============================================================
# 2. Filter to CMR route
# ============================================================
df <- raw %>% filter(Route == "Central Mediterranean")
cat("Central Mediterranean:", nrow(df), "\n")

# Parse coordinates and dates
df <- df %>%
  mutate(
    lat = as.numeric(Latitude),
    lon = as.numeric(Longitude),
    date = as.Date(incident_date_clean)
  ) %>%
  filter(!is.na(lat) & !is.na(lon) & !is.na(date))
cat("With valid coordinates + dates:", nrow(df), "\n")

# ============================================================
# 3. Restrict to core corridor
# ============================================================
cat("\n============================================================\n")
cat("2. GEOGRAPHIC RESTRICTION\n")
cat("============================================================\n\n")

n_cmr <- nrow(df)
df <- df %>%
  filter(lon >= CORE$lon_min & lon <= CORE$lon_max &
         lat >= CORE$lat_min & lat <= CORE$lat_max)

cat(sprintf("Core corridor [%.1f,%.1f]x[%.1f,%.1f]:\n",
    CORE$lon_min, CORE$lon_max, CORE$lat_min, CORE$lat_max))
cat(sprintf("  CMR incidents total:    %d\n", n_cmr))
cat(sprintf("  In core corridor:       %d (%.1f%%)\n",
    nrow(df), 100 * nrow(df) / n_cmr))
cat(sprintf("  Dropped (outside core): %d\n\n", n_cmr - nrow(df)))

# ============================================================
# 4. Clean variables (preserve all IOM columns)
# ============================================================
cat("============================================================\n")
cat("3. VARIABLE PREPARATION\n")
cat("============================================================\n\n")

df <- df %>%
  mutate(
    dead = as.numeric(No..dead),
    dead = ifelse(is.na(dead), 0, dead),
    missing = as.numeric(No..missing),
    missing = ifelse(is.na(missing), 0, missing),
    dead_missing = as.numeric(No..dead.missing),
    dead_missing = ifelse(is.na(dead_missing), 0, dead_missing),
    survivors = as.numeric(No..survivors),
    n_female = as.numeric(No..Female),
    n_male = as.numeric(No..Male),
    n_minors = as.numeric(No..minors),
    source_quality = as.numeric(Source.Quality),
    year = year(date),
    month = month(date),
    post_mou = as.integer(date >= MOU_DATE),
    cause_category = Cause.of.death..category.,
    cause_reported = Cause.of.death..reported.,
    incident_type = Incident.Type,
    is_drowning = as.integer(cause_category == "Drowning"),
    # Fatality rate (where survivors known)
    total_on_board = dead_missing + ifelse(is.na(survivors), NA_real_, survivors),
    fatality_rate = ifelse(!is.na(total_on_board) & total_on_board > 0,
                           dead_missing / total_on_board, NA_real_)
  )

cat("--- Sample composition ---\n")
cat(sprintf("  Total incidents: %d\n", nrow(df)))
cat(sprintf("  Pre-MoU:  %d\n", sum(df$post_mou == 0)))
cat(sprintf("  Post-MoU: %d\n\n", sum(df$post_mou == 1)))

cat("  Cause of death:\n")
cause_tab <- sort(table(df$cause_category), decreasing = TRUE)
for (i in seq_along(cause_tab)) {
  cat(sprintf("    %-55s %4d (%4.1f%%)\n",
      names(cause_tab)[i], cause_tab[i],
      100 * cause_tab[i] / nrow(df)))
}

cat(sprintf("\n  Dead/missing: mean=%.1f, median=%.0f, max=%d\n",
    mean(df$dead_missing), median(df$dead_missing), max(df$dead_missing)))
cat(sprintf("  == 1: %d (%.1f%%)\n",
    sum(df$dead_missing == 1), 100 * mean(df$dead_missing == 1)))
cat(sprintf("  Survivors known: %d (%.1f%%)\n",
    sum(!is.na(df$survivors)), 100 * mean(!is.na(df$survivors))))
cat(sprintf("  Source quality: median=%.0f\n", median(df$source_quality, na.rm = TRUE)))

# ============================================================
# 5. Merge ERA5 weather
# ============================================================
cat("\n============================================================\n")
cat("4. MERGING ERA5 WEATHER\n")
cat("============================================================\n\n")

ERA5_DIR <- file.path(BASE_DIR, "data", "raw", "era5")

# --- Grid setup ---
FILE_TYPES <- c("atm_instant", "atm_accum", "atm_gust", "wave", "wave_hmax")
VAR_FILE_MAP <- list(
  u10   = "atm_instant", v10 = "atm_instant", sst = "atm_instant",
  tp    = "atm_accum",
  i10fg = "atm_gust",
  swh   = "wave",         mwp = "wave",
  hmax  = "wave_hmax"
)
WEATHER_VARS <- names(VAR_FILE_MAP)

# Read grids
nc_atm <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_atm_instant_2014.nc"))
atm_lon <- ncvar_get(nc_atm, "longitude")
atm_lat <- ncvar_get(nc_atm, "latitude")
nc_close(nc_atm)

nc_wav <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
wave_lon <- ncvar_get(nc_wav, "longitude")
wave_lat <- ncvar_get(nc_wav, "latitude")
nc_close(nc_wav)

cat(sprintf("Atm grid: %d x %d (0.25 deg)\n", length(atm_lon), length(atm_lat)))
cat(sprintf("Wave grid: %d x %d (0.50 deg)\n", length(wave_lon), length(wave_lat)))

# --- Helpers ---
find_nearest <- function(target_lon, target_lat, grid_lon, grid_lat) {
  ix <- which.min(abs(grid_lon - target_lon))
  iy <- which.min(abs(grid_lat - target_lat))
  list(ix = ix, iy = iy, grid_lon = grid_lon[ix], grid_lat = grid_lat[iy])
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

nc_cache <- list()
time_cache <- list()

open_nc <- function(year, file_type) {
  key <- paste0(year, "_", file_type)
  if (is.null(nc_cache[[key]])) {
    f <- file.path(ERA5_DIR, paste0("era5_daily_cmr_", file_type, "_", year, ".nc"))
    if (!file.exists(f)) return(NULL)
    nc_cache[[key]] <<- nc_open(f)
  }
  nc_cache[[key]]
}

get_dates <- function(year, file_type) {
  key <- paste0(year, "_", file_type)
  if (is.null(time_cache[[key]])) {
    nc <- open_nc(year, file_type)
    if (is.null(nc)) return(NULL)
    time_cache[[key]] <<- get_nc_time(nc)
  }
  time_cache[[key]]
}

# --- Nearest grid cells ---
df$atm_ix  <- NA_integer_; df$atm_iy  <- NA_integer_
df$wave_ix <- NA_integer_; df$wave_iy <- NA_integer_
df$grid_lon <- NA_real_;   df$grid_lat <- NA_real_

for (i in seq_len(nrow(df))) {
  nn_a <- find_nearest(df$lon[i], df$lat[i], atm_lon, atm_lat)
  nn_w <- find_nearest(df$lon[i], df$lat[i], wave_lon, wave_lat)
  df$atm_ix[i]  <- nn_a$ix;  df$atm_iy[i]  <- nn_a$iy
  df$wave_ix[i] <- nn_w$ix;  df$wave_iy[i] <- nn_w$iy
  df$grid_lon[i] <- nn_w$grid_lon
  df$grid_lat[i] <- nn_w$grid_lat
}

df$grid_1deg <- paste0(round(df$grid_lat), "_", round(df$grid_lon))

# --- Extract weather at lags 0, 1, 2, 3, 7 ---
LAGS <- c(0, 1, 2, 3, 7)

for (v in WEATHER_VARS) {
  for (lag in LAGS) df[[paste0(v, "_lag", lag)]] <- NA_real_
}

years_in_data <- sort(unique(df$year))
cat("Extracting weather (lags", paste(LAGS, collapse = ","),
    ") for years:", paste(years_in_data, collapse = ", "), "\n")

for (yr in years_in_data) {
  for (ft in FILE_TYPES) {
    open_nc(yr, ft)
    open_nc(yr - 1, ft)  # needed for lags crossing year boundary
  }
  idx_yr <- which(df$year == yr)
  cat("  ", yr, ":", length(idx_yr), "incidents... ")

  for (i in idx_yr) {
    incident_date <- df$date[i]
    for (lag in LAGS) {
      target_date <- incident_date - lag
      target_year <- year(target_date)
      if (target_year != yr && target_year != yr - 1) next

      for (v in WEATHER_VARS) {
        ft <- VAR_FILE_MAP[[v]]
        nc_use <- open_nc(target_year, ft)
        if (is.null(nc_use)) next
        dates_use <- get_dates(target_year, ft)
        t_idx <- which(dates_use == target_date)
        if (length(t_idx) == 0) next
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

# Close NetCDF files
for (nc in nc_cache) nc_close(nc)

# --- Derived variables ---
# Wind speed and SST at each lag
for (lag in LAGS) {
  u_col <- paste0("u10_lag", lag)
  v_col <- paste0("v10_lag", lag)
  df[[paste0("wind_lag", lag)]] <- sqrt(df[[u_col]]^2 + df[[v_col]]^2)
  sst_col <- paste0("sst_lag", lag)
  df[[sst_col]] <- df[[sst_col]] - 273.15  # K -> C
}

# Convenience aliases for day-0
df$swh_day0  <- df$swh_lag0
df$wind_day0 <- df$wind_lag0
df$sst_day0  <- df$sst_lag0
df$i10fg_day0 <- df$i10fg_lag0
df$hmax_day0 <- df$hmax_lag0
df$mwp_day0  <- df$mwp_lag0

cat(sprintf("\nWeather coverage: %d/%d incidents with SWH day-0 (%.1f%%)\n",
    sum(!is.na(df$swh_day0)), nrow(df),
    100 * mean(!is.na(df$swh_day0))))
cat(sprintf("  lag1: %d, lag2: %d, lag3: %d, lag7: %d\n",
    sum(!is.na(df$swh_lag1)), sum(!is.na(df$swh_lag2)),
    sum(!is.na(df$swh_lag3)), sum(!is.na(df$swh_lag7))))

# ============================================================
# 6. Core corridor map
# ============================================================
cat("\n============================================================\n")
cat("5. GENERATING MAP\n")
cat("============================================================\n\n")

world <- ne_countries(scale = "medium", returnclass = "sf")

# Bounding box for the map — wider view for context
map_bbox <- c(xmin = 5, xmax = 22, ymin = 29, ymax = 40)

corridor_rect <- data.frame(
  xmin = CORE$lon_min, xmax = CORE$lon_max,
  ymin = CORE$lat_min, ymax = CORE$lat_max
)

p_map <- ggplot() +
  geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.3) +
  # Core corridor box
  geom_rect(data = corridor_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "#2166AC", linewidth = 1.2, linetype = "solid") +
  # Incidents
  geom_point(data = df,
             aes(x = lon, y = lat, colour = factor(post_mou),
                 size = dead_missing),
             alpha = 0.5) +
  scale_colour_manual(
    values = c("0" = "grey50", "1" = "#B2182B"),
    labels = c("0" = "Pre-MoU", "1" = "Post-MoU"),
    name = "Period"
  ) +
  scale_size_continuous(
    range = c(0.5, 5), name = "Dead + missing",
    breaks = c(1, 10, 50, 200)
  ) +
  coord_sf(xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
           ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
           expand = FALSE) +
  # Labels
  annotate("text", x = 13, y = 36.5,
           label = sprintf("Core corridor\n[%.1f,%.1f] x [%.1f,%.1f]",
                           CORE$lon_min, CORE$lon_max,
                           CORE$lat_min, CORE$lat_max),
           colour = "#2166AC", size = 3, fontface = "bold") +
  annotate("text", x = 12.6, y = 35.6, label = "Lampedusa",
           size = 2.5, fontface = "italic") +
  annotate("text", x = 15, y = 32.5, label = "Libya",
           size = 3, fontface = "italic") +
  annotate("text", x = 10, y = 37, label = "Tunisia",
           size = 2.5, fontface = "italic") +
  annotate("text", x = 14.5, y = 38, label = "Sicily",
           size = 2.5, fontface = "italic") +
  labs(
    title = "Core CMR corridor: incident locations (2014-2025)",
    subtitle = sprintf("%d incidents within core corridor. Blue box = weather measurement area.",
                       nrow(df)),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave(file.path(BASE_DIR, "output", "figures", "core_corridor_map.pdf"),
       p_map, width = 10, height = 7)
ggsave(file.path(BASE_DIR, "output", "figures", "core_corridor_map.png"),
       p_map, width = 10, height = 7, dpi = 200)
cat("Saved: output/figures/core_corridor_map.pdf + .png\n")

# ============================================================
# 7. Save
# ============================================================
cat("\n============================================================\n")
cat("6. SAVING DATASET\n")
cat("============================================================\n\n")

# Drop internal grid-matching columns
df <- df %>% select(-atm_ix, -atm_iy, -wave_ix, -wave_iy)

saveRDS(df, file.path(BASE_DIR, "data", "processed",
                       "core_corridor_incidents.RDS"))
write.csv(df, file.path(BASE_DIR, "data", "processed",
                         "core_corridor_incidents.csv"),
          row.names = FALSE)

cat("Saved: data/processed/core_corridor_incidents.RDS\n")
cat("Saved: data/processed/core_corridor_incidents.csv\n")
cat(sprintf("  Rows: %d, Cols: %d\n", nrow(df), ncol(df)))

# ============================================================
# 8. Summary
# ============================================================
cat("\n============================================================\n")
cat("DATASET SUMMARY\n")
cat("============================================================\n\n")

cat(sprintf("Period: %s to %s\n", min(df$date), max(df$date)))
cat(sprintf("Incidents: %d (pre-MoU: %d, post-MoU: %d)\n",
    nrow(df), sum(df$post_mou == 0), sum(df$post_mou == 1)))
cat(sprintf("Dead+missing: %d total\n", sum(df$dead_missing)))
cat(sprintf("Core corridor: [%.1f,%.1f] x [%.1f,%.1f]\n",
    CORE$lon_min, CORE$lon_max, CORE$lat_min, CORE$lat_max))

cat("\nBy cause of death:\n")
for (i in seq_along(cause_tab)) {
  nm <- names(cause_tab)[i]
  n <- cause_tab[i]
  dm <- sum(df$dead_missing[df$cause_category == nm])
  cat(sprintf("  %-50s %4d incidents, %5d dead+missing\n", nm, n, dm))
}

cat(sprintf("\nDrowning incidents: %d\n", sum(df$is_drowning)))
cat(sprintf("  with dead/missing > 1: %d\n",
    sum(df$is_drowning & df$dead_missing > 1)))
cat(sprintf("  with survivors known:  %d (%.1f%%)\n",
    sum(df$is_drowning & !is.na(df$survivors)),
    100 * mean(!is.na(df$survivors[df$is_drowning == 1]))))

cat(sprintf("\nSWH day-0: mean=%.2f, range=[%.2f, %.2f]\n",
    mean(df$swh_day0, na.rm = TRUE),
    min(df$swh_day0, na.rm = TRUE),
    max(df$swh_day0, na.rm = TRUE)))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
