################################################################################
# Build Extended Mortality Dataset
# ================================
# Extension of Rodriguez Sanchez et al. (2023)
#
# This script:
#   1. Reads ERA5 NetCDF files (downloaded by download_era5.py)
#   2. Computes spatial means over the Central Med and N. Africa coast
#   3. Constructs derived variables (wind speed, fog index)
#   4. Creates lagged variables (lags 01-06, matching original convention)
#   5. Merges new sea condition variables with the existing df.RDS
#   6. Saves the extended dataset as df_extended.RDS
#
# PREREQUISITES:
#   - ERA5 NetCDF files in ../data/era5/ (run download_era5.py first)
#   - df.RDS from the replication (in Original code and data/)
#   - R packages: ncdf4, terra, tidyverse, lubridate
#
# USAGE:
#   Set working directory to the thesis root, then:
#   source("Extension-1-BSTS-mortality/code/build_mortality_dataset.R")
#   OR: Rscript Extension-1-BSTS-mortality/code/build_mortality_dataset.R
################################################################################

cat("=", rep("=", 70), "\n", sep = "")
cat("BUILD EXTENDED MORTALITY DATASET\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 71), "\n\n", sep = "")

# ==============================================================================
# 0. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(tidyverse)
  library(lubridate)
})

# --- Paths (relative to thesis root) ---
# Detect script location and set paths accordingly
THESIS_ROOT <- tryCatch({
  # When sourced
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, "..", ".."))
}, error = function(e) {
  # When run interactively or via Rscript — assume working dir is thesis root
  getwd()
})

ERA5_DIR <- file.path(THESIS_ROOT, "Extension-1-BSTS-mortality", "data", "era5")
REPLICATION_DIR <- file.path(THESIS_ROOT, "Original code and data")
OUTPUT_DIR <- file.path(THESIS_ROOT, "Extension-1-BSTS-mortality", "data")

cat("Thesis root:", THESIS_ROOT, "\n")
cat("ERA5 data dir:", ERA5_DIR, "\n")
cat("Replication dir:", REPLICATION_DIR, "\n\n")

# ERA5 file paths
FILE_ATMOS <- file.path(ERA5_DIR, "era5_central_med_atmos_monthly.nc")
FILE_WAVES <- file.path(ERA5_DIR, "era5_central_med_waves_monthly.nc")
# Coast data may be split into two files (instantaneous vs accumulated variables)
FILE_COAST_INST <- file.path(ERA5_DIR, "data_stream-moda_stepType-avgua.nc")  # t2m, u10, v10, tcc
FILE_COAST_ACCUM <- file.path(ERA5_DIR, "data_stream-moda_stepType-avgad.nc") # tp

# Check files exist
for (f in c(FILE_ATMOS, FILE_WAVES, FILE_COAST_INST, FILE_COAST_ACCUM)) {
  if (!file.exists(f)) {
    stop("ERA5 file not found: ", f,
         "\nRun download_era5.py first to download the data.")
  }
}
cat("All ERA5 files found.\n\n")

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

#' Extract spatial mean time series from an ERA5 NetCDF file
#'
#' Reads a single variable from a NetCDF file, computes the spatial mean
#' (unweighted, appropriate for small regions) for each time step, and returns
#' a data frame with columns: date, value.
#'
#' @param nc_file Path to NetCDF file
#' @param var_name Short name of the variable in the NetCDF (e.g., "sst", "u10")
#' @param col_name Name for the output column
#' @return data.frame with columns: date, [col_name]
extract_era5_variable <- function(nc_file, var_name, col_name) {
  nc <- nc_open(nc_file)
  on.exit(nc_close(nc))

  # Read the variable (lon x lat x time)
  vals <- ncvar_get(nc, var_name)

  # Detect time dimension name (ERA5 uses "valid_time" or "time")
  time_dim_name <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  if (is.na(time_dim_name)) stop("No time dimension found in ", nc_file)

  time_vals <- ncvar_get(nc, time_dim_name)
  time_units <- ncatt_get(nc, time_dim_name, "units")$value

  # Parse time based on units string
  if (grepl("seconds since 1970", time_units)) {
    dates <- as.POSIXct(time_vals, origin = "1970-01-01", tz = "UTC")
  } else if (grepl("hours since 1900", time_units)) {
    dates <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC") + time_vals * 3600
  } else if (grepl("seconds since", time_units)) {
    ref <- sub("seconds since ", "", time_units)
    dates <- as.POSIXct(ref, tz = "UTC") + time_vals
  } else {
    stop("Unexpected time units in ERA5 file: ", time_units)
  }

  # Convert to Date (first of each month)
  dates <- as.Date(dates)

  # Compute spatial mean for each time step
  n_time <- length(time_vals)
  spatial_means <- numeric(n_time)

  if (length(dim(vals)) == 3) {
    # 3D array: lon x lat x time
    for (t in seq_len(n_time)) {
      spatial_means[t] <- mean(vals[, , t], na.rm = TRUE)
    }
  } else if (length(dim(vals)) == 2) {
    # 2D array: space x time (if lon/lat are collapsed)
    for (t in seq_len(n_time)) {
      spatial_means[t] <- mean(vals[, t], na.rm = TRUE)
    }
  } else {
    stop("Unexpected array dimensions for variable: ", var_name)
  }

  # Build output data frame
  result <- data.frame(
    date = floor_date(dates, "month"),  # Ensure first-of-month
    value = spatial_means
  )
  names(result)[2] <- col_name

  # Deduplicate (in case of multiple time steps per month)
  result <- result %>%
    group_by(date) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    ungroup()

  return(result)
}


#' Create lagged variables (lags 01-06), matching original dataset convention
#'
#' @param df Data frame with a "date" column and variable columns
#' @param vars Character vector of column names to lag
#' @return Data frame with original + lagged columns
create_lags <- function(df, vars) {
  df <- df %>% arrange(date)
  for (v in vars) {
    for (lag_i in 1:6) {
      lag_name <- paste0(v, "_lag_", formatC(lag_i, width = 2, flag = "0"))
      df[[lag_name]] <- dplyr::lag(df[[v]], lag_i)
    }
  }
  return(df)
}


# ==============================================================================
# 2. EXTRACT ERA5 VARIABLES
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Extracting ERA5 variables...\n")

# --- Request 1: Atmospheric variables — Central Mediterranean ---
cat("  Reading Central Med atmospheric variables...\n")

# List available variables in the file
nc_atmos <- nc_open(FILE_ATMOS)
atmos_vars <- names(nc_atmos$var)
cat("    Variables in atmospheric file:", paste(atmos_vars, collapse = ", "), "\n")
nc_close(nc_atmos)

# Extract each variable with descriptive column names
# ERA5 short names: u10, v10, sst, t2m, tcc, lcc, d2m
df_u10 <- extract_era5_variable(FILE_ATMOS, "u10", "wind_u_central_med")
df_v10 <- extract_era5_variable(FILE_ATMOS, "v10", "wind_v_central_med")
df_sst <- extract_era5_variable(FILE_ATMOS, "sst", "sst_central_med")
df_t2m <- extract_era5_variable(FILE_ATMOS, "t2m", "temperature_central_med")
df_tcc <- extract_era5_variable(FILE_ATMOS, "tcc", "cloud_cover_central_med")
df_lcc <- extract_era5_variable(FILE_ATMOS, "lcc", "low_cloud_central_med")
df_d2m <- extract_era5_variable(FILE_ATMOS, "d2m", "dewpoint_central_med")

cat("  Done.\n")

# --- Request 2: Wave variables — Central Mediterranean ---
cat("  Reading Central Med wave variables...\n")

nc_waves <- nc_open(FILE_WAVES)
wave_vars <- names(nc_waves$var)
cat("    Variables in wave file:", paste(wave_vars, collapse = ", "), "\n")
nc_close(nc_waves)

df_swh <- extract_era5_variable(FILE_WAVES, "swh", "wave_height_central_med")
df_mwp <- extract_era5_variable(FILE_WAVES, "mwp", "wave_period_central_med")
df_mwd <- extract_era5_variable(FILE_WAVES, "mwd", "wave_direction_central_med")

cat("  Done.\n")

# --- Request 3: Atmospheric variables — North Africa departure coast ---
# Coast data is split into two files: instantaneous (t2m, u10, v10, tcc)
# and accumulated (tp) variables
cat("  Reading North Africa coast variables...\n")

nc_coast_inst <- nc_open(FILE_COAST_INST)
cat("    Variables in coast inst file:", paste(names(nc_coast_inst$var), collapse = ", "), "\n")
nc_close(nc_coast_inst)

nc_coast_accum <- nc_open(FILE_COAST_ACCUM)
cat("    Variables in coast accum file:", paste(names(nc_coast_accum$var), collapse = ", "), "\n")
nc_close(nc_coast_accum)

df_t2m_coast <- extract_era5_variable(FILE_COAST_INST, "t2m", "temperature_departure_coast")
df_u10_coast <- extract_era5_variable(FILE_COAST_INST, "u10", "wind_u_departure_coast")
df_v10_coast <- extract_era5_variable(FILE_COAST_INST, "v10", "wind_v_departure_coast")
df_tcc_coast <- extract_era5_variable(FILE_COAST_INST, "tcc", "cloud_cover_departure_coast")
df_tp_coast  <- extract_era5_variable(FILE_COAST_ACCUM, "tp", "precipitation_departure_coast")

cat("  Done.\n\n")

# ==============================================================================
# 3. MERGE AND DERIVE VARIABLES
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Merging and deriving variables...\n")

# Merge all ERA5 extractions into one data frame
df_era5 <- df_u10 %>%
  left_join(df_v10, by = "date") %>%
  left_join(df_sst, by = "date") %>%
  left_join(df_t2m, by = "date") %>%
  left_join(df_tcc, by = "date") %>%
  left_join(df_lcc, by = "date") %>%
  left_join(df_d2m, by = "date") %>%
  left_join(df_swh, by = "date") %>%
  left_join(df_mwp, by = "date") %>%
  left_join(df_mwd, by = "date") %>%
  left_join(df_t2m_coast, by = "date") %>%
  left_join(df_tp_coast, by = "date") %>%
  left_join(df_u10_coast, by = "date") %>%
  left_join(df_v10_coast, by = "date") %>%
  left_join(df_tcc_coast, by = "date")

# Unit conversions
df_era5 <- df_era5 %>%
  mutate(
    # Temperature: Kelvin to Celsius
    sst_central_med = sst_central_med - 273.15,
    temperature_central_med = temperature_central_med - 273.15,
    dewpoint_central_med = dewpoint_central_med - 273.15,
    temperature_departure_coast = temperature_departure_coast - 273.15,

    # Precipitation: meters to mm (ERA5 monthly accumulated)
    precipitation_departure_coast = precipitation_departure_coast * 1000,

    # Derived: wind speed magnitude (Central Med)
    wind_speed_central_med = sqrt(wind_u_central_med^2 + wind_v_central_med^2),

    # Derived: wind speed magnitude (departure coast)
    wind_speed_departure_coast = sqrt(wind_u_departure_coast^2 + wind_v_departure_coast^2),

    # Derived: fog index (dewpoint depression — small values indicate fog risk)
    # Dewpoint depression = T - Td; when near 0, fog is likely
    dewpoint_depression_central_med = temperature_central_med - dewpoint_central_med
  )

# Select final variable set (drop u/v components, keep derived speed)
era5_base_vars <- c(
  # Central Med sea conditions
  "wave_height_central_med",
  "wave_period_central_med",
  "wave_direction_central_med",
  "wind_speed_central_med",
  "sst_central_med",
  "temperature_central_med",
  "cloud_cover_central_med",
  "low_cloud_central_med",
  "dewpoint_depression_central_med",
  # Departure coast conditions
  "temperature_departure_coast",
  "precipitation_departure_coast",
  "wind_speed_departure_coast",
  "cloud_cover_departure_coast"
)

df_era5_final <- df_era5 %>%
  dplyr::select(date, all_of(era5_base_vars))

cat("  ERA5 variables (base, before lags):", length(era5_base_vars), "\n")
cat("  Date range:", as.character(min(df_era5_final$date)),
    "to", as.character(max(df_era5_final$date)), "\n")
cat("  Observations:", nrow(df_era5_final), "\n\n")

# ==============================================================================
# 4. CREATE LAGGED VARIABLES
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating lagged variables (lags 01-06)...\n")

df_era5_lags <- create_lags(df_era5_final, era5_base_vars)

total_new_cols <- length(era5_base_vars) * 7  # base + 6 lags each
cat("  Total new columns (base + lags):", total_new_cols, "\n\n")

# ==============================================================================
# 5. LOAD AND MERGE WITH EXISTING DATASET
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Loading existing dataset...\n")

df_original <- readRDS(file.path(REPLICATION_DIR, "df.RDS"))
cat("  Original df.RDS: ", nrow(df_original), " rows x ",
    ncol(df_original), " columns\n", sep = "")

# Merge ERA5 variables into the original dataset
df_extended <- left_join(df_original, df_era5_lags, by = "date")

cat("  Extended dataset: ", nrow(df_extended), " rows x ",
    ncol(df_extended), " columns\n", sep = "")
cat("  New columns added: ", ncol(df_extended) - ncol(df_original), "\n\n")

# ==============================================================================
# 6. CONSTRUCT MORTALITY OUTCOME VARIABLES
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Constructing mortality outcome variables...\n")

df_extended <- df_extended %>%
  mutate(
    # Impute NAs in key variables (matching original code)
    LCG_pushbacks_count = as.numeric(ifelse(is.na(LCG_pushbacks_count), 0, LCG_pushbacks_count)),
    TCG_pushbacks_count = as.numeric(ifelse(is.na(TCG_pushbacks_count), 0, TCG_pushbacks_count)),
    dead_and_missing_Central_Mediterranean = as.numeric(
      ifelse(is.na(dead_and_missing_Central_Mediterranean), 0,
             dead_and_missing_Central_Mediterranean)
    ),
    # Total crossing attempts (same as original)
    crossings_CMR = arrivals_CMR + LCG_pushbacks_count + TCG_pushbacks_count +
      dead_and_missing_Central_Mediterranean,
    # Mortality rate per 100 attempted crossings
    mortality_rate_100 = (dead_and_missing_Central_Mediterranean / crossings_CMR) * 100
  )

cat("  Outcome variables constructed:\n")
cat("    - dead_and_missing_Central_Mediterranean (raw death count)\n")
cat("    - crossings_CMR (total crossing attempts)\n")
cat("    - mortality_rate_100 (deaths per 100 crossings)\n\n")

# ==============================================================================
# 7. SAVE EXTENDED DATASET
# ==============================================================================

output_file <- file.path(OUTPUT_DIR, "df_extended.RDS")
saveRDS(df_extended, output_file)

cat("[", format(Sys.time(), "%H:%M:%S"), "] Extended dataset saved to:\n")
cat("  ", output_file, "\n\n")

# Print summary of new ERA5 variables
cat("--- SUMMARY OF NEW ERA5 VARIABLES ---\n\n")

era5_summary <- df_era5_final %>%
  filter(date >= "2011-02-01" & date < "2021-10-01") %>%
  summarise(across(-date, list(
    min = ~min(.x, na.rm = TRUE),
    mean = ~mean(.x, na.rm = TRUE),
    max = ~max(.x, na.rm = TRUE),
    na_count = ~sum(is.na(.x))
  )))

for (v in era5_base_vars) {
  cat(sprintf("  %-40s  min=%7.2f  mean=%7.2f  max=%7.2f  NAs=%d\n",
              v,
              era5_summary[[paste0(v, "_min")]],
              era5_summary[[paste0(v, "_mean")]],
              era5_summary[[paste0(v, "_max")]],
              era5_summary[[paste0(v, "_na_count")]]))
}

cat("\n")
cat(rep("=", 71), "\n", sep = "")
cat("DATASET BUILD COMPLETE\n")
cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 71), "\n", sep = "")
