# 02_data_prepare.R
# =================
# Functions for preparing the model-ready dataset from df_enhanced.
# Extracted from run_new_model.R lines 84-250.
#
# Key corrections applied:
#   1. SST anomaly recomputed using pre-period climatology only
#   2. Seasonality via 11 month dummies (ref = January)
#   3. Lag discipline: all variables at lag 0 + lag 1


#' Define intervention specifications
#'
#' @return tibble with columns: name, label, pre_end, post_start,
#'         intervention_date, outcome
define_interventions <- function() {
  tibble::tibble(
    name = c("mare_nostrum", "ngo_sar", "mou",
             "mare_nostrum", "ngo_sar", "mou"),
    label = c("A_mortality", "B_mortality", "C_mortality",
              "A_deaths",    "B_deaths",    "C_deaths"),
    treatment = c("Mare Nostrum", "NGO SAR", "EU-Libya MoU",
                  "Mare Nostrum", "NGO SAR", "EU-Libya MoU"),
    outcome = c("mortality_rate", "mortality_rate", "mortality_rate",
                "death_count",    "death_count",    "death_count"),
    pre_end = as.Date(c("2013-09-01", "2014-10-01", "2017-01-01",
                        "2013-09-01", "2014-10-01", "2017-01-01")),
    post_start = as.Date(c("2013-10-01", "2014-11-01", "2017-02-01",
                           "2013-10-01", "2014-11-01", "2017-02-01"))
  )
}


#' Fix SST anomaly to use pre-period climatology only
#'
#' Recomputes sst_anomaly_central_med using climatology from dates
#' before the MoU (Feb 2017) to avoid post-period data leakage.
#'
#' @param df Data frame with sst_central_med and sst_anomaly_central_med
#' @param pre_period_end Date before which to compute climatology
#' @return df with corrected sst_anomaly_central_med
compute_sst_anomaly_preperiod <- function(df,
                                          pre_period_end = as.Date("2017-02-01")) {
  if (!"sst_central_med" %in% names(df)) {
    warning("sst_central_med not found, skipping SST anomaly correction.")
    return(df)
  }

  pre_sst <- df %>%
    dplyr::filter(date < pre_period_end) %>%
    dplyr::mutate(month_num = lubridate::month(date)) %>%
    dplyr::group_by(month_num) %>%
    dplyr::summarise(sst_clim_pre = mean(sst_central_med, na.rm = TRUE),
                     .groups = "drop")

  df %>%
    dplyr::mutate(
      sst_anomaly_ORIG = sst_anomaly_central_med,
      month_num = lubridate::month(date)
    ) %>%
    dplyr::left_join(
      pre_sst %>% dplyr::select(month_num, sst_clim_pre),
      by = "month_num"
    ) %>%
    dplyr::mutate(sst_anomaly_central_med = sst_central_med - sst_clim_pre) %>%
    dplyr::select(-month_num, -sst_clim_pre, -sst_anomaly_ORIG)
}


#' Select predictors with lag discipline
#'
#' All variables at lag 0 + lag 1.
#'
#' @param df Data frame (must contain all candidate columns)
#' @return Character vector of available predictor column names
select_predictors <- function(df) {
  # Base variables: only per-crossing danger predictors.
  # Volume-channel predictors (arrival-side weather, oil price) are excluded
  # because they predict crossing attempts, not per-crossing danger. Including
  # a few volume predictors would be inconsistent — either include all of them
  # (full predictor model) or none (this curated specification).
  base_vars <- c(
    # Central Mediterranean sea/atmosphere
    "wave_height_central_med", "wave_period_central_med",
    "wave_direction_central_med", "wind_speed_central_med",
    "sst_central_med", "sst_anomaly_central_med",
    "cloud_cover_central_med", "low_cloud_central_med",
    "dewpoint_depression_central_med", "temperature_central_med",
    # Departure coast conditions
    "wind_speed_departure_coast", "cloud_cover_departure_coast",
    "temperature_departure_coast", "precipitation_departure_coast",
    # Extreme wave statistics (from daily data)
    "wave_max_central_med", "wave_sd_central_med",
    "wave_days_above_2m",
    # Ocean currents
    "current_speed_central_med", "current_against_route"
  )

  # All variables at lag 0 + lag 1
  all_vars <- c(base_vars, paste0(base_vars, "_lag_01"))

  # Filter to columns actually available in the data
  all_vars[all_vars %in% names(df)]
}


#' Prepare model-ready datasets
#'
#' Applies SST anomaly correction, predictor selection, outcome construction,
#' month dummies, and time window filtering.
#'
#' @param df_enhanced Enhanced dataset (output of build_enhanced_dataset)
#' @return list with: df_mort (mortality rate), df_deaths (death count),
#'         cov_cols (predictor names), month_dummy_cols (month dummy names),
#'         interventions (tibble from define_interventions)
prepare_model_data <- function(df_enhanced) {
  EPSILON <- 0.01

  # Ensure outcome variables exist
  df <- df_enhanced %>%
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

  # Fix SST anomaly: pre-period climatology only
  message("[data_prepare] Correcting SST anomaly to pre-period climatology...")
  df <- compute_sst_anomaly_preperiod(df)

  # Select predictors with lag discipline
  message("[data_prepare] Selecting predictors with lag discipline...")
  cov_cols <- select_predictors(df)

  # Filter time window
  df_filtered <- df %>%
    dplyr::filter(date >= "2011-02-01" & date < "2021-10-01")

  # Re-check predictor availability after filtering
  cov_cols <- cov_cols[cov_cols %in% names(df_filtered)]

  # Month dummy column names
  month_dummy_cols <- paste0("month_", 2:12)

  # --- Mortality rate outcome (primary) ---
  df_mort <- df_filtered %>%
    dplyr::transmute(
      date = date,
      mortality_rate = log(mortality_rate_100 + EPSILON),
      dplyr::across(dplyr::all_of(cov_cols))
    ) %>%
    stats::na.omit()

  # Add month dummies (ref = January)
  for (m in 2:12) {
    df_mort[[paste0("month_", m)]] <- as.integer(
      lubridate::month(df_mort$date) == m)
  }

  # --- Death count outcome (robustness) ---
  df_deaths <- df_filtered %>%
    dplyr::transmute(
      date = date,
      deaths_cmr = log(dead_and_missing_Central_Mediterranean + 1),
      dplyr::across(dplyr::all_of(cov_cols))
    ) %>%
    stats::na.omit()

  for (m in 2:12) {
    df_deaths[[paste0("month_", m)]] <- as.integer(
      lubridate::month(df_deaths$date) == m)
  }

  # Store the full (unfiltered) df for descriptive figure
  df_full <- df

  message("[data_prepare] Mortality rate: ", nrow(df_mort), " rows x ",
          ncol(df_mort), " cols (", ncol(df_mort) - 2, " predictors)")
  message("[data_prepare] Death count: ", nrow(df_deaths), " rows x ",
          ncol(df_deaths), " cols")

  list(
    df_mort         = df_mort,
    df_deaths       = df_deaths,
    df_full         = df_full,
    cov_cols        = cov_cols,
    month_dummy_cols = month_dummy_cols,
    interventions   = define_interventions()
  )
}
