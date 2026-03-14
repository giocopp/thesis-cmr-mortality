# run_forecast_diagnostic_deaths.R
# ================================
# Runs the same 5-specification rolling-origin CV as the mortality rate
# diagnostic, but for log(deaths + 1). Standalone script that loads
# cached data from the targets store.

library(dplyr)
library(lubridate)
library(bsts)

# --- Constants (match 04_diagnostics.R) ---
NITER    <- 10000L
BURN     <- 2000L
EXP_SIZE <- 3L
SEED     <- 270488L

# --- Load cached data from targets ---
message("[deaths_diag] Loading df_model from targets store...")
df_model <- targets::tar_read(df_model,
  store = "_targets")

cov_cols        <- df_model$cov_cols
month_dummy_cols <- df_model$month_dummy_cols
df_full         <- df_model$df_full

# --- Build diagnostic dataset (death count outcome) ---
PRE_START <- as.Date("2011-02-01")
PRE_END   <- as.Date("2017-01-01")

df_pre <- df_full %>%
  filter(date >= PRE_START & date <= PRE_END) %>%
  mutate(
    dead_and_missing_Central_Mediterranean = as.numeric(
      ifelse(is.na(dead_and_missing_Central_Mediterranean), 0,
             dead_and_missing_Central_Mediterranean))
  )

cov_cols_diag <- cov_cols[cov_cols %in% names(df_pre)]

df_diag <- df_pre %>%
  transmute(
    date = date,
    y    = log(dead_and_missing_Central_Mediterranean + 1),
    across(all_of(cov_cols_diag)),
    month_num = month(date)
  ) %>%
  na.omit()

for (m in 2:12) {
  df_diag[[paste0("month_", m)]] <- as.integer(df_diag$month_num == m)
}

sd_full_pre <- sd(df_diag$y)
message("[deaths_diag] SD of outcome (full pre-period): ", round(sd_full_pre, 4))
message("[deaths_diag] N obs: ", nrow(df_diag))

# --- Define folds (same as mortality rate diagnostic) ---
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

# --- Helpers ---
calc_metrics <- function(actual, predicted, sd_full) {
  resid <- actual - predicted
  rmse  <- sqrt(mean(resid^2))
  mae   <- mean(abs(resid))
  list(rmse = rmse, mae = mae,
       rmse_sd_full = rmse / sd_full)
}

fit_bsts_with_cov <- function(y_train, X_train, X_test, n_test,
                              state_fn, label) {
  ss <- state_fn(list(), y_train)
  train_df <- data.frame(y = y_train, X_train)

  message("  Fitting ", label, " (K=", ncol(X_train),
          ", niter=", NITER, ")...")

  model <- tryCatch({
    bsts::bsts(y ~ ., state.specification = ss, data = train_df,
               niter = NITER, expected.model.size = EXP_SIZE,
               max.flips = -1, seed = SEED, ping = 0)
  }, error = function(e) {
    message("    FAILED: ", conditionMessage(e))
    return(NULL)
  })

  if (is.null(model)) {
    return(list(pred = rep(mean(y_train), n_test), ok = FALSE))
  }

  pred <- predict(model, newdata = as.data.frame(X_test),
                  horizon = n_test, burn = BURN)$mean
  coef_mat  <- model$coefficients[-(1:BURN), , drop = FALSE]
  inc_probs <- colMeans(coef_mat != 0)

  list(pred = as.numeric(pred), inc_probs = inc_probs, ok = TRUE)
}

fit_bsts_state_only <- function(y_train, n_test, state_fn, label) {
  ss <- state_fn(list(), y_train)
  message("  Fitting ", label, " (state-only, niter=", NITER, ")...")

  model <- tryCatch({
    bsts::bsts(y_train, state.specification = ss, niter = NITER,
               seed = SEED, ping = 0)
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

# --- Run CV ---
message("[deaths_diag] Running rolling-origin CV (5 models x 3 folds)...")
set.seed(SEED)

all_fold_results <- list()

for (fi in seq_along(diag_folds)) {
  f <- diag_folds[[fi]]
  message("\n=== ", f$name, " ===")

  train_data <- df_diag %>% filter(date <= f$train_end)
  test_data  <- df_diag %>% filter(date >= f$test_start &
                                     date <= f$test_end)

  y_train <- train_data$y
  y_test  <- test_data$y
  n_train <- length(y_train)
  n_test  <- length(y_test)
  message("  Train: ", n_train, " months, Test: ", n_test, " months")

  # M0a: month-of-year mean
  month_means <- train_data %>%
    group_by(month_num) %>%
    summarise(pred = mean(y, na.rm = TRUE), .groups = "drop")
  pred_m0a <- test_data %>%
    left_join(month_means, by = "month_num") %>%
    pull(pred)
  pred_m0a[is.na(pred_m0a)] <- mean(y_train)

  # M0c: BSTS state-only (LL + Seasonal)
  res_m0c <- fit_bsts_state_only(y_train, n_test, state_ll_seas, "M0c")

  # A1-month-only: LL + month dummies only
  X_train_mo <- train_data %>% select(all_of(month_dummy_cols))
  X_test_mo  <- test_data  %>% select(all_of(month_dummy_cols))
  res_mo <- fit_bsts_with_cov(y_train, X_train_mo, X_test_mo, n_test,
                              state_ll, "A1-month-only")

  # A1-exog-only: LL + exogenous covariates only
  X_train_ex <- train_data %>% select(all_of(cov_cols_diag))
  X_test_ex  <- test_data  %>% select(all_of(cov_cols_diag))
  res_ex <- fit_bsts_with_cov(y_train, X_train_ex, X_test_ex, n_test,
                              state_ll, "A1-exog-only")

  # A1-full: LL + dummies + covariates
  X_train_full <- train_data %>%
    select(all_of(c(cov_cols_diag, month_dummy_cols)))
  X_test_full  <- test_data %>%
    select(all_of(c(cov_cols_diag, month_dummy_cols)))
  res_full <- fit_bsts_with_cov(y_train, X_train_full, X_test_full, n_test,
                                state_ll, "A1-full")

  # Collect metrics
  models_vec <- c("M0a (seasonal mean)", "M0c (state-only)",
                   "A1-month-only", "A1-exog-only", "A1-full")
  preds <- list(pred_m0a, res_m0c$pred, res_mo$pred,
                res_ex$pred, res_full$pred)
  metrics <- lapply(preds, function(p) calc_metrics(y_test, p, sd_full_pre))

  fold_df <- data.frame(
    fold    = f$name,
    model   = models_vec,
    n_train = n_train,
    n_test  = n_test,
    rmse    = sapply(metrics, `[[`, "rmse"),
    stringsAsFactors = FALSE
  )
  all_fold_results[[fi]] <- fold_df
}

# --- Aggregate (weighted by n_test) ---
all_metrics <- bind_rows(all_fold_results)
wavg <- all_metrics %>%
  group_by(model) %>%
  summarise(
    wt_rmse = weighted.mean(rmse, n_test),
    total_months = sum(n_test),
    .groups = "drop"
  )

rmse_m0a <- wavg$wt_rmse[wavg$model == "M0a (seasonal mean)"]
wavg$vs_m0a <- (rmse_m0a - wavg$wt_rmse) / rmse_m0a * 100

# --- Print results ---
message("\n======================================")
message("DEATH COUNT FORECAST DIAGNOSTIC RESULTS")
message("======================================")
message("Outcome: log(deaths + 1)")
message("SD (full pre-period): ", round(sd_full_pre, 4))
message("")

for (i in seq_len(nrow(wavg))) {
  message(sprintf("  %-20s  wt.RMSE = %.3f  RMSE/SD = %.3f  vs M0a = %+.1f%%",
                  wavg$model[i], wavg$wt_rmse[i],
                  wavg$wt_rmse[i] / sd_full_pre,
                  wavg$vs_m0a[i]))
}

# --- Save results ---
outfile <- "output/tables/forecast_diagnostic_deaths.csv"
readr::write_csv(wavg, outfile)
message("\nResults saved to: ", outfile)

# --- Also print per-fold details ---
message("\n--- Per-fold details ---")
for (fi in seq_along(all_fold_results)) {
  message("\n", diag_folds[[fi]]$name, ":")
  fd <- all_fold_results[[fi]]
  for (i in seq_len(nrow(fd))) {
    message(sprintf("  %-20s  RMSE = %.3f", fd$model[i], fd$rmse[i]))
  }
}
