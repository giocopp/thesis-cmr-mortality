# 03_model_fit.R
# ==============
# Functions for fitting CausalImpact models and extracting results.
# Extracted from run_new_model.R lines 357-492.


# --- Constants ---
NITER    <- 10000L
BURN     <- 2000L
EXP_SIZE <- 3L
EPSILON  <- 0.01
MOU_DATE <- as.Date("2017-02-01")


#' Fit a single CausalImpact model
#'
#' @param df_model Output of prepare_model_data() (list)
#' @param intervention_label One of: "A_mortality", "B_mortality",
#'        "C_mortality", "C_deaths"
#' @param seed Random seed (default 270488)
#' @param niter MCMC iterations (default 10000)
#' @return list with: impact (CausalImpact object), results (extracted df),
#'         label, outcome_type, elapsed_minutes
fit_causal_impact <- function(df_model, intervention_label,
                              seed = 270488, niter = NITER) {
  interventions <- df_model$interventions
  spec <- interventions %>% dplyr::filter(label == intervention_label)

  if (nrow(spec) != 1) {
    stop("Unknown intervention label: ", intervention_label)
  }

  # Select the right dataset
  if (spec$outcome == "death_count") {
    df <- df_model$df_deaths
  } else {
    df <- df_model$df_mort
  }

  pre_period  <- c(min(df$date), spec$pre_end)
  post_period <- c(spec$post_start, max(df$date))

  model_args <- list(
    dynamic.regression = FALSE,
    standardize.data   = TRUE,
    max.flips          = -1L,
    niter              = niter
  )

  message("[model_fit] Fitting ", intervention_label, "...")
  set.seed(seed)
  t1 <- Sys.time()

  impact <- CausalImpact::CausalImpact(
    df,
    pre.period  = pre_period,
    post.period = post_period,
    alpha       = 0.05,
    model.args  = model_args
  )

  t2 <- Sys.time()
  elapsed <- as.numeric(difftime(t2, t1, units = "mins"))
  message("[model_fit] ", intervention_label, " done in ",
          round(elapsed, 1), " minutes.")

  # Extract time series results
  results <- extract_results(impact, df$date)

  list(
    impact         = impact,
    results        = results,
    label          = intervention_label,
    outcome_type   = spec$outcome,
    treatment      = spec$treatment,
    pre_period     = pre_period,
    post_period    = post_period,
    elapsed_minutes = elapsed
  )
}


#' Extract time series results from a CausalImpact object
#'
#' @param impact_obj CausalImpact object
#' @param dates Date vector matching the data rows
#' @return data.frame with observed, predicted, effects, and cumulative effects
extract_results <- function(impact_obj, dates) {
  data.frame(
    date                     = dates,
    original                 = as.numeric(impact_obj$series$response),
    prediction               = as.numeric(impact_obj$series$point.pred),
    prediction_lower         = as.numeric(impact_obj$series$point.pred.lower),
    prediction_upper         = as.numeric(impact_obj$series$point.pred.upper),
    pointwise_effect         = as.numeric(impact_obj$series$point.effect),
    pointwise_effect_lower   = as.numeric(impact_obj$series$point.effect.lower),
    pointwise_effect_upper   = as.numeric(impact_obj$series$point.effect.upper),
    cumulative_effect        = as.numeric(impact_obj$series$cum.effect),
    cumulative_effect_lower  = as.numeric(impact_obj$series$cum.effect.lower),
    cumulative_effect_upper  = as.numeric(impact_obj$series$cum.effect.upper)
  )
}


#' Extract summary statistics from a fitted model
#'
#' @param model_fit Output of fit_causal_impact()
#' @return one-row tibble with key statistics
extract_model_summary <- function(model_fit) {
  impact <- model_fit$impact
  s <- impact$summary
  ser <- impact$series

  pre_idx <- which(ser$cum.effect == 0)
  pre_n   <- if (length(pre_idx) > 0) max(pre_idx) else NA_integer_
  post_n  <- nrow(ser) - pre_n

  rmse_pre <- NA_real_
  sd_pre   <- NA_real_
  rmse_sd_pre <- NA_real_

  if (!is.na(pre_n) && pre_n > 0) {
    resid    <- ser$response[seq_len(pre_n)] - ser$point.pred[seq_len(pre_n)]
    rmse_pre <- sqrt(mean(resid^2, na.rm = TRUE))
    sd_pre   <- sd(ser$response[seq_len(pre_n)], na.rm = TRUE)
    rmse_sd_pre <- if (!is.na(sd_pre) && sd_pre > 0) rmse_pre / sd_pre else NA_real_
  }

  # Inclusion probabilities
  inc_probs <- colMeans(impact$model$bsts.model$coefficients != 0)
  burn_used <- tryCatch(
    bsts::SuggestBurn(0.1, impact$model$bsts.model),
    error = function(e) NA_integer_
  )

  tibble::tibble(
    label           = model_fit$label,
    treatment       = model_fit$treatment,
    outcome_type    = model_fit$outcome_type,
    abs_effect      = as.numeric(s["Cumulative", "AbsEffect"]),
    abs_effect_lower = as.numeric(s["Cumulative", "AbsEffect.lower"]),
    abs_effect_upper = as.numeric(s["Cumulative", "AbsEffect.upper"]),
    rel_effect      = as.numeric(s["Cumulative", "RelEffect"]),
    p_value         = as.numeric(s["Cumulative", "p"]),
    n_pre           = pre_n,
    n_post          = post_n,
    rmse_pre        = rmse_pre,
    sd_pre          = sd_pre,
    rmse_sd_pre     = rmse_sd_pre,
    burn_used_draws = burn_used,
    n_inc_01        = sum(inc_probs > 0.01),
    n_inc_10        = sum(inc_probs > 0.10),
    elapsed_minutes = model_fit$elapsed_minutes
  )
}
