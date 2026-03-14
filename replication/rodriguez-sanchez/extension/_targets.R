# _targets.R
# ==========
# Pipeline definition for Extension-2: Minimal Exogenous BSTS/CausalImpact
# mortality model for the Central Mediterranean migration route.
#
# Usage:
#   targets::tar_make()           # Run pipeline
#   targets::tar_visnetwork()     # Visualize DAG
#   targets::tar_read(model_a)    # Read cached target
#
# Expected runtime: ~6-8 hours (MCMC-dominated)

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    "dplyr", "tidyr", "readr", "tibble", "lubridate", "stringr",
    "ggplot2", "scales", "gridExtra", "grid",
    "bsts", "CausalImpact", "ncdf4", "zoo"
  ),
  seed = 270488
)

# Auto-source all function files in R/
tar_source("R/")


list(
  # =========================================================================
  # Input files (format = "file" — targets tracks file modification times)
  # =========================================================================
  tar_target(file_df_original,
             "data/raw/df.RDS",
             format = "file"),
  tar_target(file_era5_atmos,
             "data/raw/era5_central_med_atmos_monthly.nc",
             format = "file"),
  tar_target(file_era5_waves,
             "data/raw/era5_central_med_waves_monthly.nc",
             format = "file"),
  tar_target(file_era5_coast_inst,
             "data/raw/data_stream-moda_stepType-avgua.nc",
             format = "file"),
  tar_target(file_era5_coast_accum,
             "data/raw/data_stream-moda_stepType-avgad.nc",
             format = "file"),
  tar_target(file_moon,
             "data/raw/moon_illumination.csv",
             format = "file"),
  tar_target(file_daily_waves,
             "data/raw/era5_daily_waves_central_med.nc",
             format = "file"),
  tar_target(file_currents,
             "data/raw/medsea_currents_central_med.nc",
             format = "file"),

  # =========================================================================
  # Data building
  # =========================================================================
  tar_target(
    df_base,
    build_base_dataset(
      file_df_original, file_era5_atmos, file_era5_waves,
      file_era5_coast_inst, file_era5_coast_accum
    )
  ),

  tar_target(
    df_enhanced,
    build_enhanced_dataset(
      df_base, file_moon, file_daily_waves, file_currents
    )
  ),

  tar_target(
    df_model,
    prepare_model_data(df_enhanced)
  ),

  # =========================================================================
  # CausalImpact models (each takes ~1-2 hours)
  # =========================================================================
  tar_target(
    model_a,
    fit_causal_impact(df_model, "A_mortality")
  ),

  tar_target(
    model_b,
    fit_causal_impact(df_model, "B_mortality")
  ),

  tar_target(
    model_c_rate,
    fit_causal_impact(df_model, "C_mortality")
  ),

  tar_target(
    model_a_deaths,
    fit_causal_impact(df_model, "A_deaths")
  ),

  tar_target(
    model_b_deaths,
    fit_causal_impact(df_model, "B_deaths")
  ),

  tar_target(
    model_c_deaths,
    fit_causal_impact(df_model, "C_deaths")
  ),

  # =========================================================================
  # Diagnostics
  # =========================================================================
  tar_target(
    placebos_a,
    run_placebo_tests(df_model, "A_mortality",
                      min_pre_months = 12, min_post_months = 6,
                      step_months = 2)
  ),

  tar_target(
    placebos_b,
    run_placebo_tests(df_model, "B_mortality",
                      min_pre_months = 12, min_post_months = 6,
                      step_months = 3)
  ),

  tar_target(
    placebos_c_rate,
    run_placebo_tests(df_model, "C_mortality",
                      min_pre_months = 12, min_post_months = 6,
                      step_months = 3)
  ),

  tar_target(
    placebos_a_deaths,
    run_placebo_tests(df_model, "A_deaths",
                      min_pre_months = 12, min_post_months = 6,
                      step_months = 2)
  ),

  tar_target(
    placebos_b_deaths,
    run_placebo_tests(df_model, "B_deaths",
                      min_pre_months = 12, min_post_months = 6,
                      step_months = 3)
  ),

  tar_target(
    placebos_c_deaths,
    run_placebo_tests(df_model, "C_deaths",
                      min_pre_months = 12, min_post_months = 6,
                      step_months = 3)
  ),

  tar_target(
    truncation_a,
    run_truncation_test(df_model)
  ),

  tar_target(
    forecast_diag,
    run_forecast_diagnostic(df_model)
  ),

  # =========================================================================
  # Figures (each returns file paths)
  # =========================================================================
  tar_target(
    fig_descriptive,
    plot_descriptive(df_model),
    format = "file"
  ),

  tar_target(
    fig_counterfactual,
    plot_counterfactual(model_a, model_b, model_c_rate,
                        model_a_deaths, model_b_deaths, model_c_deaths),
    format = "file"
  ),

  tar_target(
    fig_pointwise,
    plot_pointwise(model_a, model_b, model_c_rate,
                   model_a_deaths, model_b_deaths, model_c_deaths),
    format = "file"
  ),

  tar_target(
    fig_placebos,
    plot_placebo_distributions(
      placebos_a, placebos_b, placebos_c_rate,
      placebos_a_deaths, placebos_b_deaths, placebos_c_deaths,
      model_a, model_b, model_c_rate,
      model_a_deaths, model_b_deaths, model_c_deaths
    ),
    format = "file"
  ),

  tar_target(
    fig_truncation,
    plot_truncation(truncation_a),
    format = "file"
  ),

  tar_target(
    fig_forecast_diag,
    plot_forecast_diagnostic(forecast_diag),
    format = "file"
  ),

  # =========================================================================
  # Snapshot tables
  # =========================================================================
  tar_target(
    snapshot,
    write_snapshot(
      models    = list(model_a, model_b, model_c_rate,
                       model_a_deaths, model_b_deaths, model_c_deaths),
      placebos  = list(placebos_a, placebos_b, placebos_c_rate,
                       placebos_a_deaths, placebos_b_deaths,
                       placebos_c_deaths),
      truncation = truncation_a,
      forecast  = forecast_diag,
      output_dir = "output/tables"
    ),
    format = "file"
  )
)
