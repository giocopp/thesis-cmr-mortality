# 06_snapshot.R
# =============
# Functions for building consistency snapshot tables from model results.
# Extracted from build_extension2_snapshot.R.


#' Build model summary snapshot
#'
#' Extracts key statistics from all fitted models into a one-row-per-model
#' summary table.
#'
#' @param models Named list of model fits (output of fit_causal_impact)
#' @return tibble with one row per model
build_model_snapshot <- function(models) {
  dplyr::bind_rows(lapply(models, function(m) {
    extract_model_summary(m)
  }))
}


#' Build placebo summary snapshot
#'
#' Computes false positive rates and percentiles from placebo results.
#'
#' @param placebos Named list of placebo results (tibbles from run_placebo_tests)
#' @param models Named list of model fits (for actual p-values)
#' @return tibble with one row per model
build_placebo_snapshot <- function(placebos, models) {
  all_placebos <- dplyr::bind_rows(placebos)

  actual_results <- tibble::tibble(
    model    = sapply(models, function(m) m$label),
    actual_p = sapply(models, function(m) {
      as.numeric(m$impact$summary["Cumulative", "p"])
    })
  )

  dplyr::bind_rows(lapply(unique(all_placebos$model), function(mod) {
    psub <- all_placebos %>% dplyr::filter(model == mod)
    asub <- actual_results %>% dplyr::filter(model == mod)

    n_tests      <- nrow(psub)
    n_sig        <- sum(psub$significant, na.rm = TRUE)
    actual_p     <- as.numeric(asub$actual_p[1])
    more_extreme <- sum(psub$p_value <= actual_p, na.rm = TRUE)

    tibble::tibble(
      model                     = mod,
      n_placebo_tests           = n_tests,
      n_placebo_significant     = n_sig,
      false_positive_rate       = if (n_tests > 0) n_sig / n_tests else NA_real_,
      actual_p_value            = actual_p,
      placebos_p_le_actual      = more_extreme,
      share_placebos_p_le_actual = if (n_tests > 0) more_extreme / n_tests else NA_real_,
      placebo_p_min             = min(psub$p_value, na.rm = TRUE),
      placebo_p_median          = median(psub$p_value, na.rm = TRUE),
      placebo_p_max             = max(psub$p_value, na.rm = TRUE)
    )
  }))
}


#' Write snapshot tables to CSV + markdown
#'
#' @param models Named list of model fits
#' @param placebos Named list of placebo results
#' @param truncation Truncation test results (tibble)
#' @param forecast Forecast diagnostic results (list)
#' @param output_dir Output directory for tables
#' @return character vector of file paths written
write_snapshot <- function(models, placebos, truncation = NULL,
                           forecast = NULL, output_dir = "output/tables") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  models_df   <- build_model_snapshot(models)
  placebos_df <- build_placebo_snapshot(placebos, models)

  models_csv   <- file.path(output_dir, "consistency_snapshot_models.csv")
  placebos_csv <- file.path(output_dir, "consistency_snapshot_placebos.csv")
  snapshot_md  <- file.path(output_dir, "consistency_snapshot.md")

  readr::write_csv(models_df, models_csv)
  readr::write_csv(placebos_df, placebos_csv)

  # Build markdown
  fmt <- function(x, digits = 4) {
    if (is.na(x)) return("NA")
    format(round(x, digits), nsmall = digits, trim = TRUE)
  }

  md <- c(
    "# Extension-2 Consistency Snapshot",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "Pipeline: `targets::tar_make()`",
    "",
    "## Model Results",
    "",
    "Note: AbsEffect and CI are on the model's transformed outcome scale (log units).",
    "",
    "| Model | Treatment | Outcome | p-value | AbsEffect (log units) | 95% CI lower | 95% CI upper | RMSE/SD (pre) | Burn used (draws) |",
    "|---|---|---|---:|---:|---:|---:|---:|---:|"
  )

  for (i in seq_len(nrow(models_df))) {
    r <- models_df[i, ]
    md <- c(md, paste0(
      "| ", r$label,
      " | ", r$treatment,
      " | ", r$outcome_type,
      " | ", fmt(r$p_value),
      " | ", fmt(r$abs_effect),
      " | ", fmt(r$abs_effect_lower),
      " | ", fmt(r$abs_effect_upper),
      " | ", fmt(r$rmse_sd_pre),
      " | ", r$burn_used_draws,
      " |"
    ))
  }

  md <- c(
    md, "",
    "## Placebo Diagnostics", "",
    "| Model | # Placebos | # Significant | FPR | Actual p | Placebos p<=actual | Share p<=actual |",
    "|---|---:|---:|---:|---:|---:|---:|"
  )

  for (i in seq_len(nrow(placebos_df))) {
    r <- placebos_df[i, ]
    md <- c(md, paste0(
      "| ", r$model,
      " | ", r$n_placebo_tests,
      " | ", r$n_placebo_significant,
      " | ", fmt(r$false_positive_rate, 3),
      " | ", fmt(r$actual_p_value, 4),
      " | ", r$placebos_p_le_actual, "/", r$n_placebo_tests,
      " | ", fmt(r$share_placebos_p_le_actual, 3),
      " |"
    ))
  }

  if (!is.null(truncation) && nrow(truncation) > 0) {
    md <- c(md, "",
      "## Model A Truncation", "",
      "| Label | Post months | p-value | Significant |",
      "|---|---:|---:|---|"
    )
    for (i in seq_len(nrow(truncation))) {
      r <- truncation[i, ]
      md <- c(md, paste0(
        "| ", r$label,
        " | ", r$post_months,
        " | ", fmt(r$p_value, 4),
        " | ", ifelse(isTRUE(r$significant), "YES", "NO"),
        " |"
      ))
    }
  }

  if (!is.null(forecast)) {
    md <- c(md, "",
      "## Forecast Diagnostic (weighted averages)", "",
      "| Model | wt.RMSE | wt.MAE | RMSE/SD | vs M0a | vs M0c |",
      "|---|---:|---:|---:|---:|---:|"
    )
    for (i in seq_len(nrow(forecast$wavg))) {
      w <- forecast$wavg[i, ]
      md <- c(md, paste0(
        "| ", w$model,
        " | ", fmt(w$wavg_rmse, 3),
        " | ", fmt(w$wavg_mae, 3),
        " | ", fmt(w$wavg_rmse_sd, 3),
        " | ", sprintf("%+.1f%%", w$wavg_impr_m0a),
        " | ", sprintf("%+.1f%%", w$wavg_impr_m0c),
        " |"
      ))
    }
  }

  writeLines(md, snapshot_md)

  message("[snapshot] Written: ", models_csv)
  message("[snapshot] Written: ", placebos_csv)
  message("[snapshot] Written: ", snapshot_md)

  c(models_csv, placebos_csv, snapshot_md)
}
