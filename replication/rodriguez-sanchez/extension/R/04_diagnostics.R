# 04_diagnostics.R
# ================
# Functions for placebo tests, truncation tests, and forecast diagnostics.
# Extracted from run_new_model.R lines 1018-1826.


# --- Constants (shared with 03_model_fit.R) ---
NITER    <- 10000L
BURN     <- 2000L
EXP_SIZE <- 3L


#' Run placebo / falsification tests for a single intervention
#'
#' Restricts data to the pre-period, then tests fake interventions at regular
#' intervals. Returns a tibble of placebo results.
#'
#' @param df_model Output of prepare_model_data() (list)
#' @param intervention_label Label from define_interventions()
#' @param min_pre_months Minimum pre-period for each placebo (default 12)
#' @param min_post_months Minimum post-period for each placebo (default 6)
#' @param step_months Step size between placebo dates (default 3)
#' @param seed Random seed
#' @param niter MCMC iterations
#' @return tibble with columns: model, placebo_date, pre_months, post_months,
#'         cum_effect, rel_effect, p_value, significant, direction
run_placebo_tests <- function(df_model, intervention_label,
                              min_pre_months = 12L, min_post_months = 6L,
                              step_months = 3L, seed = 270488,
                              niter = NITER) {
  interventions <- df_model$interventions
  spec <- interventions %>% dplyr::filter(label == intervention_label)

  if (spec$outcome == "death_count") {
    data <- df_model$df_deaths
  } else {
    data <- df_model$df_mort
  }

  actual_int_date <- spec$post_start
  dates      <- data$date
  data_start <- min(dates)
  pre_data   <- data %>% dplyr::filter(date < actual_int_date)
  pre_dates  <- pre_data$date
  n_pre      <- nrow(pre_data)

  if (n_pre < min_pre_months + min_post_months) {
    warning("Not enough pre-period data for placebo tests (",
            n_pre, " months, need ", min_pre_months + min_post_months, ")")
    return(tibble::tibble())
  }

  placebo_start_idx <- min_pre_months + 1
  placebo_end_idx   <- n_pre - min_post_months
  if (placebo_start_idx > placebo_end_idx) {
    warning("No valid placebo dates for ", intervention_label)
    return(tibble::tibble())
  }

  placebo_indices <- seq(placebo_start_idx, placebo_end_idx, by = step_months)
  placebo_dates   <- pre_dates[placebo_indices]

  model_args <- list(
    dynamic.regression = FALSE,
    standardize.data   = TRUE,
    max.flips          = -1L,
    niter              = niter
  )

  message("[diagnostics] Testing ", length(placebo_dates),
          " placebo dates for ", intervention_label)

  results <- list()
  set.seed(seed)

  for (i in seq_along(placebo_dates)) {
    pdate <- placebo_dates[i]
    message("  [", i, "/", length(placebo_dates), "] Placebo: ",
            as.character(pdate), "...")

    placebo_pre  <- c(data_start, pdate - months(1))
    placebo_post <- c(pdate, actual_int_date - months(1))

    pre_n  <- sum(pre_data$date >= placebo_pre[1] &
                    pre_data$date <= placebo_pre[2])
    post_n <- sum(pre_data$date >= placebo_post[1] &
                    pre_data$date <= placebo_post[2])

    tryCatch({
      impact <- CausalImpact::CausalImpact(
        pre_data,
        pre.period  = placebo_pre,
        post.period = placebo_post,
        alpha       = 0.05,
        model.args  = model_args
      )

      s       <- impact$summary
      p_val   <- s["Cumulative", "p"]
      cum_eff <- s["Cumulative", "AbsEffect"]
      rel_eff <- s["Cumulative", "RelEffect"]

      results[[i]] <- tibble::tibble(
        model       = intervention_label,
        placebo_date = pdate,
        pre_months  = pre_n,
        post_months = post_n,
        cum_effect  = cum_eff,
        rel_effect  = rel_eff,
        p_value     = p_val,
        significant = p_val < 0.05,
        direction   = ifelse(cum_eff > 0, "higher", "lower")
      )
    }, error = function(e) {
      message("    ERROR: ", conditionMessage(e))
    })
  }

  dplyr::bind_rows(results)
}


#' Run post-period truncation test for Model A (Mare Nostrum)
#'
#' Tests whether the Mare Nostrum effect holds when the post-period is
#' progressively restricted.
#'
#' @param df_model Output of prepare_model_data() (list)
#' @param seed Random seed
#' @param niter MCMC iterations
#' @return tibble with truncation results
run_truncation_test <- function(df_model, seed = 270488, niter = NITER) {
  df_mort <- df_model$df_mort

  truncation_specs <- tibble::tibble(
    post_end = as.Date(c("2014-10-01", "2015-10-01",
                         "2016-10-01", "2021-09-01")),
    label = c("Mare Nostrum window only (13 mo)", "+1yr after MN (25 mo)",
              "Pre-MoU window (37 mo)", "Full post-period (96 mo)")
  )

  model_args <- list(
    dynamic.regression = FALSE,
    standardize.data   = TRUE,
    max.flips          = -1L,
    niter              = niter
  )

  message("[diagnostics] Running truncation tests for Model A...")
  results <- list()
  set.seed(seed)

  for (i in seq_len(nrow(truncation_specs))) {
    end_d <- truncation_specs$post_end[i]
    lab   <- truncation_specs$label[i]

    message("  Truncation: ", lab, " (ends ", as.character(end_d), ")...")

    df_trunc <- df_mort %>% dplyr::filter(date <= end_d)

    tryCatch({
      impact <- CausalImpact::CausalImpact(
        df_trunc,
        pre.period  = c(min(df_trunc$date), as.Date("2013-09-01")),
        post.period = c(as.Date("2013-10-01"), end_d),
        alpha       = 0.05,
        model.args  = model_args
      )

      s      <- impact$summary
      post_n <- sum(df_trunc$date >= as.Date("2013-10-01"))

      results[[i]] <- tibble::tibble(
        post_end    = end_d,
        post_months = post_n,
        label       = lab,
        cum_effect  = s["Cumulative", "AbsEffect"],
        rel_effect  = s["Cumulative", "RelEffect"],
        p_value     = s["Cumulative", "p"],
        significant = s["Cumulative", "p"] < 0.05
      )
    }, error = function(e) {
      message("    ERROR: ", conditionMessage(e))
    })
  }

  dplyr::bind_rows(results)
}


#' Run pre-period forecasting diagnostic (3-fold rolling-origin CV)
#'
#' Cross-validates the BSTS specification on the pre-period (2011-02 to 2017-01)
#' to assess whether exogenous covariates add predictive signal beyond seasonal
#' baselines.
#'
#' Models compared:
#'   M0a: month-of-year mean
#'   M0c: BSTS AddLocalLevel + AddSeasonal(12), no covariates
#'   M1-A1: AddLocalLevel + month dummies + exogenous covariates
#'   M1-A2: AddLocalLinearTrend + month dummies + exogenous covariates
#'   M1-B:  AddLocalLevel + AddSeasonal(12) + exogenous covariates
#'
#' @param df_model Output of prepare_model_data() (list)
#' @param seed Random seed
#' @param niter MCMC iterations
#' @return list with: all_metrics, wavg, fold_details, sd_full_pre, settings
run_forecast_diagnostic <- function(df_model, seed = 270488, niter = NITER) {
  EPSILON <- 0.01
  PRE_START <- as.Date("2011-02-01")
  PRE_END   <- as.Date("2017-01-01")

  cov_cols        <- df_model$cov_cols
  month_dummy_cols <- df_model$month_dummy_cols

  # Build diagnostic dataset from the full (unflitered) data
  df_full <- df_model$df_full
  df_pre <- df_full %>%
    dplyr::filter(date >= PRE_START & date <= PRE_END)

  # Re-check available covariates
  cov_cols_diag <- cov_cols[cov_cols %in% names(df_pre)]

  df_diag <- df_pre %>%
    dplyr::transmute(
      date = date,
      y    = log(mortality_rate_100 + EPSILON),
      dplyr::across(dplyr::all_of(cov_cols_diag)),
      month_num = lubridate::month(date)
    ) %>%
    stats::na.omit()

  # Add month dummies
  for (m in 2:12) {
    df_diag[[paste0("month_", m)]] <- as.integer(df_diag$month_num == m)
  }

  sd_full_pre <- sd(df_diag$y)

  # Define folds
  diag_folds <- list(
    list(name = "Fold 1",
         train_end  = as.Date("2014-12-01"),
         test_start = as.Date("2015-01-01"),
         test_end   = as.Date("2015-12-01")),
    list(name = "Fold 2",
         train_end  = as.Date("2015-12-01"),
         test_start = as.Date("2016-01-01"),
         test_end   = as.Date("2016-12-01")),
    list(name = "Fold 3",
         train_end  = as.Date("2016-06-01"),
         test_start = as.Date("2016-07-01"),
         test_end   = as.Date("2017-01-01"))
  )

  # Helper: calculate metrics
  calc_metrics <- function(actual, predicted, sd_full) {
    resid    <- actual - predicted
    rmse     <- sqrt(mean(resid^2))
    mae      <- mean(abs(resid))
    sd_test  <- sd(actual)
    list(rmse = rmse, mae = mae, sd_test = sd_test,
         rmse_sd_test = rmse / sd_test,
         rmse_sd_full = rmse / sd_full)
  }

  # BSTS fitting helpers
  fit_bsts_with_cov <- function(y_train, X_train, X_test, n_test,
                                state_fn, label) {
    ss <- state_fn(list(), y_train)
    train_df <- data.frame(y = y_train, X_train)

    message("  Fitting ", label, " (K=", ncol(X_train),
            ", niter=", niter, ")...")

    model <- tryCatch({
      bsts::bsts(y ~ ., state.specification = ss, data = train_df,
                 niter = niter, expected.model.size = EXP_SIZE,
                 max.flips = -1, seed = seed, ping = 0)
    }, error = function(e) {
      message("    FAILED: ", conditionMessage(e))
      return(NULL)
    })

    if (is.null(model)) {
      return(list(pred = rep(mean(y_train), n_test),
                  inc_probs = NULL, ok = FALSE))
    }

    pred <- predict(model, newdata = as.data.frame(X_test),
                    horizon = n_test, burn = BURN)$mean
    coef_mat  <- model$coefficients[-(1:BURN), , drop = FALSE]
    inc_probs <- colMeans(coef_mat != 0)

    list(pred = as.numeric(pred), inc_probs = inc_probs, ok = TRUE)
  }

  fit_bsts_state_only <- function(y_train, n_test, state_fn, label) {
    ss <- state_fn(list(), y_train)
    message("  Fitting ", label, " (state-only, niter=", niter, ")...")

    model <- tryCatch({
      bsts::bsts(y_train, state.specification = ss, niter = niter,
                 seed = seed, ping = 0)
    }, error = function(e) {
      message("    FAILED: ", conditionMessage(e))
      return(NULL)
    })

    if (is.null(model)) {
      return(list(pred = rep(mean(y_train), n_test), ok = FALSE))
    }

    pred <- predict(model, horizon = n_test, burn = BURN)$mean
    list(pred = as.numeric(pred), ok = TRUE)
  }

  # State specification factories
  state_ll <- function(ss, y) bsts::AddLocalLevel(ss, y)
  state_llt <- function(ss, y) bsts::AddLocalLinearTrend(ss, y)
  state_ll_seas <- function(ss, y) {
    ss <- bsts::AddLocalLevel(ss, y)
    bsts::AddSeasonal(ss, y, nseasons = 12)
  }

  # Run diagnostic across folds
  message("[diagnostics] Running rolling-origin CV (5 models x 3 folds)...")
  set.seed(seed)

  diag_results_list <- list()

  for (fi in seq_along(diag_folds)) {
    f <- diag_folds[[fi]]
    message("--- ", f$name, " ---")

    train_data <- df_diag %>% dplyr::filter(date <= f$train_end)
    test_data  <- df_diag %>% dplyr::filter(date >= f$test_start &
                                              date <= f$test_end)

    y_train <- train_data$y
    y_test  <- test_data$y
    n_train <- length(y_train)
    n_test  <- length(y_test)

    # M0a: month-of-year mean
    month_means <- train_data %>%
      dplyr::group_by(month_num) %>%
      dplyr::summarise(pred = mean(y, na.rm = TRUE), .groups = "drop")
    pred_m0a <- test_data %>%
      dplyr::left_join(month_means, by = "month_num") %>%
      dplyr::pull(pred)
    pred_m0a[is.na(pred_m0a)] <- mean(y_train)

    # M0c: BSTS state-only
    res_m0c <- fit_bsts_state_only(y_train, n_test, state_ll_seas, "M0c")

    # M1-A1: LL + dummies + covariates
    X_train_a <- train_data %>%
      dplyr::select(dplyr::all_of(c(cov_cols_diag, month_dummy_cols)))
    X_test_a  <- test_data %>%
      dplyr::select(dplyr::all_of(c(cov_cols_diag, month_dummy_cols)))
    res_a1 <- fit_bsts_with_cov(y_train, X_train_a, X_test_a, n_test,
                                state_ll, "M1-A1")

    # M1-A2: LLT + dummies + covariates
    res_a2 <- fit_bsts_with_cov(y_train, X_train_a, X_test_a, n_test,
                                state_llt, "M1-A2")

    # M1-B: LL + Seasonal(12) + covariates (no dummies)
    X_train_b <- train_data %>%
      dplyr::select(dplyr::all_of(cov_cols_diag))
    X_test_b  <- test_data %>%
      dplyr::select(dplyr::all_of(cov_cols_diag))
    res_b <- fit_bsts_with_cov(y_train, X_train_b, X_test_b, n_test,
                               state_ll_seas, "M1-B")

    # Metrics
    models_vec <- c("M0a (month mean)", "M0c (state-only)",
                    "M1-A1 (LL+dummies)", "M1-A2 (LLT+dummies)",
                    "M1-B (LL+Seasonal)")
    metrics_list <- list(
      calc_metrics(y_test, pred_m0a, sd_full_pre),
      calc_metrics(y_test, res_m0c$pred, sd_full_pre),
      calc_metrics(y_test, res_a1$pred, sd_full_pre),
      calc_metrics(y_test, res_a2$pred, sd_full_pre),
      calc_metrics(y_test, res_b$pred, sd_full_pre)
    )
    preds_list <- list(pred_m0a, res_m0c$pred, res_a1$pred,
                       res_a2$pred, res_b$pred)

    m_m0a <- metrics_list[[1]]
    m_m0c <- metrics_list[[2]]

    fold_df <- data.frame(
      fold         = f$name,
      model        = models_vec,
      n_train      = n_train,
      n_test       = n_test,
      rmse         = sapply(metrics_list, `[[`, "rmse"),
      mae          = sapply(metrics_list, `[[`, "mae"),
      rmse_sd_test = sapply(metrics_list, `[[`, "rmse_sd_test"),
      rmse_sd_full = sapply(metrics_list, `[[`, "rmse_sd_full"),
      stringsAsFactors = FALSE
    )
    fold_df$impr_vs_m0a <- (m_m0a$rmse - fold_df$rmse) / m_m0a$rmse * 100
    fold_df$impr_vs_m0c <- (m_m0c$rmse - fold_df$rmse) / m_m0c$rmse * 100

    diag_results_list[[fi]] <- list(
      fold     = f$name,
      metrics  = fold_df,
      actuals  = y_test,
      preds    = stats::setNames(preds_list, models_vec),
      dates    = test_data$date,
      inc_a1   = res_a1$inc_probs,
      inc_a2   = res_a2$inc_probs,
      inc_b    = res_b$inc_probs
    )
  }

  # Aggregate results
  all_diag_metrics <- dplyr::bind_rows(
    lapply(diag_results_list, function(r) r$metrics)
  )

  diag_wavg <- all_diag_metrics %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      wavg_rmse     = weighted.mean(rmse, n_test),
      wavg_mae      = weighted.mean(mae, n_test),
      wavg_rmse_sd  = weighted.mean(rmse, n_test) / sd_full_pre,
      wavg_impr_m0a = NA_real_,
      wavg_impr_m0c = NA_real_,
      total_months  = sum(n_test),
      .groups = "drop"
    )

  diag_rmse_m0a <- diag_wavg$wavg_rmse[diag_wavg$model == "M0a (month mean)"]
  diag_rmse_m0c <- diag_wavg$wavg_rmse[diag_wavg$model == "M0c (state-only)"]

  diag_wavg$wavg_impr_m0a <- (diag_rmse_m0a - diag_wavg$wavg_rmse) /
    diag_rmse_m0a * 100
  diag_wavg$wavg_impr_m0c <- (diag_rmse_m0c - diag_wavg$wavg_rmse) /
    diag_rmse_m0c * 100

  list(
    all_metrics  = all_diag_metrics,
    wavg         = diag_wavg,
    fold_details = diag_results_list,
    sd_full_pre  = sd_full_pre,
    settings = list(
      version           = "v2 -- corrected seasonality",
      sst_climatology   = "pre-period only (dates < 2017-02-01)",
      lag_discipline    = "hazard vars lag 0 only; other vars lag 0-1",
      n_cov             = length(cov_cols_diag),
      n_month_dummies   = length(month_dummy_cols),
      predictor_list    = cov_cols_diag,
      bsts_niter        = niter,
      bsts_burn         = BURN,
      bsts_expected_size = EXP_SIZE,
      bsts_max_flips    = -1L,
      seed              = seed,
      outcome           = "log(mortality_rate_100 + 0.01)",
      pre_period        = c(PRE_START, PRE_END)
    )
  )
}
