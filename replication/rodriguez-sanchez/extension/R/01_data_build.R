# 01_data_build.R
# ===============
# Functions for building the analysis dataset from raw NetCDF + CSV inputs.
# Extracted from Extension-1/code/build_mortality_dataset.R and
# Extension-2/code/build_enhanced_dataset.R.


#' Extract spatial mean time series from a NetCDF file
#'
#' Reads a single variable, computes unweighted spatial mean per time step,
#' aggregates to monthly, and returns a two-column data frame.
#'
#' @param nc_file Path to NetCDF file
#' @param var_name Short name of the variable (e.g., "sst", "u10")
#' @param col_name Name for the output column
#' @return data.frame with columns: date, [col_name]
extract_nc_spatial_mean <- function(nc_file, var_name, col_name) {
  nc <- ncdf4::nc_open(nc_file)
  on.exit(ncdf4::nc_close(nc))

  vals <- ncdf4::ncvar_get(nc, var_name)

  time_dim_name <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  if (is.na(time_dim_name)) stop("No time dimension in ", nc_file)

  time_vals  <- ncdf4::ncvar_get(nc, time_dim_name)
  time_units <- ncdf4::ncatt_get(nc, time_dim_name, "units")$value

  if (grepl("seconds since 1970", time_units)) {
    dates <- as.POSIXct(time_vals, origin = "1970-01-01", tz = "UTC")
  } else if (grepl("hours since 1900", time_units)) {
    dates <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC") + time_vals * 3600
  } else if (grepl("minutes since", time_units)) {
    ref <- sub("minutes since ", "", time_units)
    dates <- as.POSIXct(ref, tz = "UTC") + time_vals * 60
  } else if (grepl("seconds since", time_units)) {
    ref <- sub("seconds since ", "", time_units)
    dates <- as.POSIXct(ref, tz = "UTC") + time_vals
  } else {
    stop("Unexpected time units: ", time_units)
  }
  dates <- as.Date(dates)

  n_time <- length(time_vals)
  spatial_means <- numeric(n_time)

  ndim <- length(dim(vals))
  if (ndim == 3) {
    for (t in seq_len(n_time)) spatial_means[t] <- mean(vals[, , t], na.rm = TRUE)
  } else if (ndim == 4) {
    for (t in seq_len(n_time)) spatial_means[t] <- mean(vals[, , 1, t], na.rm = TRUE)
  } else if (ndim == 2) {
    for (t in seq_len(n_time)) spatial_means[t] <- mean(vals[, t], na.rm = TRUE)
  }

  result <- data.frame(
    date  = lubridate::floor_date(dates, "month"),
    value = spatial_means
  )
  names(result)[2] <- col_name

  result %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(),
                                   \(x) mean(x, na.rm = TRUE)),
                     .groups = "drop")
}


#' Create lagged variables (lags 01-06)
#'
#' @param df Data frame with a "date" column
#' @param vars Character vector of column names to lag
#' @param n_lags Number of lags (default 6)
#' @return Data frame with original + lagged columns
create_lags <- function(df, vars, n_lags = 6L) {
  df <- dplyr::arrange(df, date)
  for (v in vars) {
    for (lag_i in seq_len(n_lags)) {
      lag_name <- paste0(v, "_lag_", formatC(lag_i, width = 2, flag = "0"))
      df[[lag_name]] <- dplyr::lag(df[[v]], lag_i)
    }
  }
  df
}


#' Build the base dataset from df.RDS + ERA5 monthly NetCDF files
#'
#' Replicates Extension-1/code/build_mortality_dataset.R: reads ERA5 atmospheric,
#' wave, and coast files; computes spatial means; derives wind speed, fog index;
#' creates lags 01-06; merges into the original authors' dataset; constructs
#' mortality outcome variables.
#'
#' @param file_df Path to the original df.RDS
#' @param file_atmos Path to ERA5 Central Med atmospheric NetCDF
#' @param file_waves Path to ERA5 Central Med waves NetCDF
#' @param file_coast_inst Path to departure coast instantaneous NetCDF
#' @param file_coast_accum Path to departure coast accumulated NetCDF
#' @return tibble (df_extended equivalent, ~14,933 columns)
build_base_dataset <- function(file_df, file_atmos, file_waves,
                               file_coast_inst, file_coast_accum) {
  message("[data_build] Extracting ERA5 variables...")

  # Central Med atmospheric
  df_u10 <- extract_nc_spatial_mean(file_atmos, "u10", "wind_u_central_med")
  df_v10 <- extract_nc_spatial_mean(file_atmos, "v10", "wind_v_central_med")
  df_sst <- extract_nc_spatial_mean(file_atmos, "sst", "sst_central_med")
  df_t2m <- extract_nc_spatial_mean(file_atmos, "t2m", "temperature_central_med")
  df_tcc <- extract_nc_spatial_mean(file_atmos, "tcc", "cloud_cover_central_med")
  df_lcc <- extract_nc_spatial_mean(file_atmos, "lcc", "low_cloud_central_med")
  df_d2m <- extract_nc_spatial_mean(file_atmos, "d2m", "dewpoint_central_med")

  # Central Med waves
  df_swh <- extract_nc_spatial_mean(file_waves, "swh", "wave_height_central_med")
  df_mwp <- extract_nc_spatial_mean(file_waves, "mwp", "wave_period_central_med")
  df_mwd <- extract_nc_spatial_mean(file_waves, "mwd", "wave_direction_central_med")

  # Departure coast
  df_t2m_coast <- extract_nc_spatial_mean(file_coast_inst, "t2m",
                                          "temperature_departure_coast")
  df_u10_coast <- extract_nc_spatial_mean(file_coast_inst, "u10",
                                          "wind_u_departure_coast")
  df_v10_coast <- extract_nc_spatial_mean(file_coast_inst, "v10",
                                          "wind_v_departure_coast")
  df_tcc_coast <- extract_nc_spatial_mean(file_coast_inst, "tcc",
                                          "cloud_cover_departure_coast")
  df_tp_coast  <- extract_nc_spatial_mean(file_coast_accum, "tp",
                                          "precipitation_departure_coast")

  # Merge all ERA5 extractions
  df_era5 <- df_u10 %>%
    dplyr::left_join(df_v10, by = "date") %>%
    dplyr::left_join(df_sst, by = "date") %>%
    dplyr::left_join(df_t2m, by = "date") %>%
    dplyr::left_join(df_tcc, by = "date") %>%
    dplyr::left_join(df_lcc, by = "date") %>%
    dplyr::left_join(df_d2m, by = "date") %>%
    dplyr::left_join(df_swh, by = "date") %>%
    dplyr::left_join(df_mwp, by = "date") %>%
    dplyr::left_join(df_mwd, by = "date") %>%
    dplyr::left_join(df_t2m_coast, by = "date") %>%
    dplyr::left_join(df_tp_coast,  by = "date") %>%
    dplyr::left_join(df_u10_coast, by = "date") %>%
    dplyr::left_join(df_v10_coast, by = "date") %>%
    dplyr::left_join(df_tcc_coast, by = "date")

  # Unit conversions and derived variables
  df_era5 <- df_era5 %>%
    dplyr::mutate(
      sst_central_med              = sst_central_med - 273.15,
      temperature_central_med      = temperature_central_med - 273.15,
      dewpoint_central_med         = dewpoint_central_med - 273.15,
      temperature_departure_coast  = temperature_departure_coast - 273.15,
      precipitation_departure_coast = precipitation_departure_coast * 1000,
      wind_speed_central_med       = sqrt(wind_u_central_med^2 +
                                          wind_v_central_med^2),
      wind_speed_departure_coast   = sqrt(wind_u_departure_coast^2 +
                                          wind_v_departure_coast^2),
      dewpoint_depression_central_med = temperature_central_med -
                                        dewpoint_central_med
    )

  era5_base_vars <- c(
    "wave_height_central_med", "wave_period_central_med",
    "wave_direction_central_med", "wind_speed_central_med",
    "sst_central_med", "temperature_central_med",
    "cloud_cover_central_med", "low_cloud_central_med",
    "dewpoint_depression_central_med",
    "temperature_departure_coast", "precipitation_departure_coast",
    "wind_speed_departure_coast", "cloud_cover_departure_coast"
  )

  df_era5_final <- df_era5 %>%
    dplyr::select(date, dplyr::all_of(era5_base_vars))

  # Create lags 01-06
  df_era5_lags <- create_lags(df_era5_final, era5_base_vars)

  # Load original dataset and merge
  message("[data_build] Merging ERA5 into df.RDS...")
  df_original <- readRDS(file_df)
  df_extended <- dplyr::left_join(df_original, df_era5_lags, by = "date")

  # Construct mortality outcome variables
  df_extended <- df_extended %>%
    dplyr::mutate(
      LCG_pushbacks_count = as.numeric(
        ifelse(is.na(LCG_pushbacks_count), 0, LCG_pushbacks_count)),
      TCG_pushbacks_count = as.numeric(
        ifelse(is.na(TCG_pushbacks_count), 0, TCG_pushbacks_count)),
      dead_and_missing_Central_Mediterranean = as.numeric(
        ifelse(is.na(dead_and_missing_Central_Mediterranean), 0,
               dead_and_missing_Central_Mediterranean)),
      crossings_CMR = arrivals_CMR + LCG_pushbacks_count +
        TCG_pushbacks_count + dead_and_missing_Central_Mediterranean,
      mortality_rate_100 = (dead_and_missing_Central_Mediterranean /
                              crossings_CMR) * 100
    )

  message("[data_build] Base dataset: ", nrow(df_extended), " rows x ",
          ncol(df_extended), " cols")
  df_extended
}


#' Build the enhanced dataset by adding Extension-2 variables
#'
#' Adds moon illumination, ocean currents, extreme wave statistics, and
#' SST anomaly to the base dataset (from build_enhanced_dataset.R).
#'
#' @param df_base Base dataset (output of build_base_dataset)
#' @param file_moon Path to moon_illumination.csv
#' @param file_daily_waves Path to ERA5 daily waves NetCDF
#' @param file_currents Path to ocean currents NetCDF
#' @return tibble (df_enhanced equivalent, ~14,982 columns)
build_enhanced_dataset <- function(df_base, file_moon, file_daily_waves,
                                   file_currents) {
  df <- df_base

  # --- 1. Moon illumination ---
  message("[data_build] Adding moon illumination...")
  df_moon <- readr::read_csv(file_moon, show_col_types = FALSE) %>%
    dplyr::mutate(date = as.Date(date))
  df_moon <- create_lags(df_moon, "moon_illumination_frac")
  df <- dplyr::left_join(df, df_moon, by = "date")

  # --- 2. Ocean surface currents ---
  message("[data_build] Adding ocean currents...")
  df_uo <- extract_nc_spatial_mean(file_currents, "uo",
                                   "current_eastward_central_med")
  df_vo <- extract_nc_spatial_mean(file_currents, "vo",
                                   "current_northward_central_med")

  df_currents <- df_uo %>%
    dplyr::left_join(df_vo, by = "date") %>%
    dplyr::mutate(
      current_speed_central_med = sqrt(current_eastward_central_med^2 +
                                       current_northward_central_med^2),
      current_against_route = -current_northward_central_med
    ) %>%
    dplyr::select(date, current_speed_central_med, current_against_route)

  df_currents <- create_lags(df_currents,
                             c("current_speed_central_med",
                               "current_against_route"))
  df <- dplyr::left_join(df, df_currents, by = "date")

  # --- 3. Extreme wave statistics from daily data ---
  message("[data_build] Adding extreme wave statistics...")
  nc <- ncdf4::nc_open(file_daily_waves)
  swh <- ncdf4::ncvar_get(nc, "swh")

  time_dim_name <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  time_vals  <- ncdf4::ncvar_get(nc, time_dim_name)
  time_units <- ncdf4::ncatt_get(nc, time_dim_name, "units")$value

  if (grepl("seconds since 1970", time_units)) {
    dates <- as.POSIXct(time_vals, origin = "1970-01-01", tz = "UTC")
  } else if (grepl("hours since 1900", time_units)) {
    dates <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC") + time_vals * 3600
  } else {
    ref <- sub("seconds since ", "", time_units)
    dates <- as.POSIXct(ref, tz = "UTC") + time_vals
  }
  dates <- as.Date(dates)
  ncdf4::nc_close(nc)

  n_time <- length(dates)
  spatial_means <- numeric(n_time)
  ndim <- length(dim(swh))
  if (ndim == 3) {
    for (t in seq_len(n_time)) spatial_means[t] <- mean(swh[, , t], na.rm = TRUE)
  } else if (ndim == 2) {
    for (t in seq_len(n_time)) spatial_means[t] <- mean(swh[, t], na.rm = TRUE)
  }

  df_daily <- data.frame(
    date       = dates,
    month_date = lubridate::floor_date(dates, "month"),
    swh_spatial_mean = spatial_means
  )

  df_extreme_waves <- df_daily %>%
    dplyr::group_by(month_date) %>%
    dplyr::summarise(
      wave_max_central_med   = max(swh_spatial_mean, na.rm = TRUE),
      wave_sd_central_med    = sd(swh_spatial_mean, na.rm = TRUE),
      wave_days_above_2m     = sum(swh_spatial_mean > 2.0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(date = month_date)

  df_extreme_waves <- create_lags(
    df_extreme_waves,
    c("wave_max_central_med", "wave_sd_central_med", "wave_days_above_2m")
  )
  df <- dplyr::left_join(df, df_extreme_waves, by = "date")

  # --- 4. SST anomaly (full-sample climatology; corrected to pre-period in prepare step) ---
  message("[data_build] Computing SST anomaly (full-sample, will be corrected later)...")
  if ("sst_central_med" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(month_num = lubridate::month(date)) %>%
      dplyr::group_by(month_num) %>%
      dplyr::mutate(sst_climatology = mean(sst_central_med, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sst_anomaly_central_med = sst_central_med - sst_climatology) %>%
      dplyr::select(-month_num, -sst_climatology)

    df <- dplyr::arrange(df, date)
    for (lag_i in 1:6) {
      lag_name <- paste0("sst_anomaly_central_med_lag_",
                         formatC(lag_i, width = 2, flag = "0"))
      df[[lag_name]] <- dplyr::lag(df$sst_anomaly_central_med, lag_i)
    }
  }

  message("[data_build] Enhanced dataset: ", nrow(df), " rows x ",
          ncol(df), " cols")
  df
}
