# autocorrelation_diagnostics.R
# =============================
# Detailed exploration of residual autocorrelation in the monthly
# weather-mortality interaction regression.

library(dplyr)
library(lubridate)
library(MASS)

PROJECT_DIR <- here::here()

df <- readRDS("replication/rodriguez-sanchez/extension/data/df_extended.RDS") %>%
  filter(date >= as.Date("2011-02-01") & date <= as.Date("2021-09-01")) %>%
  mutate(log_rate = log(mortality_rate_100 + 0.01),
         deaths   = dead_and_missing_Central_Mediterranean,
         post_mou = as.integer(date >= as.Date("2017-07-01")),
         month_fac = factor(month(date)),
         year_fac  = factor(year(date)))

# Merge corridor/core weather from daily ERA5 panel
daily_panel <- readRDS(file.path(PROJECT_DIR, "data", "processed",
                                  "cmr_daily_weather_panel.RDS"))
monthly_wx <- daily_panel %>%
  mutate(ym_date = lubridate::floor_date(date, "month")) %>%
  group_by(ym_date) %>%
  summarise(
    swh_corridor  = mean(swh_mean, na.rm = TRUE),
    swh_core      = mean(swh_core, na.rm = TRUE),
    wind_corridor = mean(wind_mean, na.rm = TRUE),
    wind_core     = mean(wind_core, na.rm = TRUE),
    .groups = "drop"
  )
df <- df %>% left_join(monthly_wx, by = c("date" = "ym_date"))

cat("============================================================\n")
cat("AUTOCORRELATION DIAGNOSTICS (broad + core weather)\n")
cat("============================================================\n\n")

# --- Broad geography OLS (original) ---
m <- lm(log_rate ~ wave_height_central_med * post_mou +
         wind_speed_central_med * post_mou + month_fac, data = df)
e <- residuals(m)
n <- length(e)

# --- Core geography OLS (ERA5 subsample) ---
df_era5 <- df %>% filter(!is.na(swh_core))
m_core <- lm(log_rate ~ swh_core * post_mou + wind_core * post_mou + month_fac,
              data = df_era5)
e_core <- residuals(m_core)
n_core <- length(e_core)

# --- NegBin on death counts (ERA5 subsample) ---
m_nb <- glm.nb(deaths ~ swh_core * post_mou + wind_core + month_fac,
                data = df_era5)
e_nb <- residuals(m_nb, type = "pearson")
n_nb <- length(e_nb)

cat(sprintf("Models fitted:\n"))
cat(sprintf("  Broad OLS:  n=%d\n", n))
cat(sprintf("  Core OLS:   n=%d\n", n_core))
cat(sprintf("  Core NegBin: n=%d\n\n", n_nb))

# ============================================================
# 1. RESIDUAL ACF UP TO LAG 12 (all three models)
# ============================================================
cat("============================================================\n")
cat("1. RESIDUAL AUTOCORRELATION FUNCTION\n")
cat("============================================================\n\n")

acf_obj <- acf(e, lag.max = 12, plot = FALSE)
acf_vals <- acf_obj$acf[-1, 1, 1]  # drop lag 0
bound <- 1.96 / sqrt(n)

acf_core_obj <- acf(e_core, lag.max = 12, plot = FALSE)
acf_core_vals <- acf_core_obj$acf[-1, 1, 1]
bound_core <- 1.96 / sqrt(n_core)

acf_nb_obj <- acf(e_nb, lag.max = 12, plot = FALSE)
acf_nb_vals <- acf_nb_obj$acf[-1, 1, 1]
bound_nb <- 1.96 / sqrt(n_nb)

cat(sprintf("  %4s  %10s  %10s  %10s\n", "Lag", "Broad OLS", "Core OLS", "Core NegBin"))
cat(paste(rep("-", 45), collapse = ""), "\n")
for (i in 1:12) {
  flag_b <- ifelse(abs(acf_vals[i]) > bound, "*", " ")
  flag_c <- ifelse(i <= length(acf_core_vals) && abs(acf_core_vals[i]) > bound_core, "*", " ")
  flag_n <- ifelse(i <= length(acf_nb_vals) && abs(acf_nb_vals[i]) > bound_nb, "*", " ")
  cat(sprintf("  %4d  %+8.3f %s  %+8.3f %s  %+8.3f %s\n",
      i, acf_vals[i], flag_b,
      ifelse(i <= length(acf_core_vals), acf_core_vals[i], NA), flag_c,
      ifelse(i <= length(acf_nb_vals), acf_nb_vals[i], NA), flag_n))
}
cat(sprintf("\n  95%% bounds: broad=%.3f, core=%.3f, nb=%.3f\n",
    bound, bound_core, bound_nb))

# ============================================================
# 2. LJUNG-BOX TESTS (all three models)
# ============================================================
cat("\n============================================================\n")
cat("2. LJUNG-BOX TEST FOR RESIDUAL AUTOCORRELATION\n")
cat("============================================================\n\n")

for (model_label in list(list("Broad OLS", e), list("Core OLS", e_core),
                          list("Core NegBin", e_nb))) {
  cat(sprintf("  %s:\n", model_label[[1]]))
  for (lag in c(3, 6, 12)) {
    lb <- Box.test(model_label[[2]], lag = lag, type = "Ljung-Box")
    cat(sprintf("    Lags 1-%2d: Q = %7.3f, p = %.4f%s\n",
        lag, lb$statistic, lb$p.value,
        ifelse(lb$p.value < 0.05, " **", "")))
  }
  cat("\n")
}

# ============================================================
# 3. HAC SE SENSITIVITY TO BANDWIDTH
# ============================================================
cat("\n============================================================\n")
cat("3. HAC SE SENSITIVITY TO BANDWIDTH (Newey-West lag)\n")
cat("   Wind x PostMoU interaction\n")
cat("============================================================\n\n")

newey_west_vcov <- function(model, max_lag) {
  X <- model.matrix(model)
  e <- residuals(model)
  n <- length(e)
  k <- ncol(X)
  meat <- matrix(0, k, k)
  for (j in 0:max_lag) {
    w <- if (j == 0) 1 else 1 - j / (max_lag + 1)
    for (t in (j + 1):n) {
      if (j == 0) {
        meat <- meat + e[t]^2 * tcrossprod(X[t, ])
      } else {
        meat <- meat + w * e[t] * e[t - j] *
          (tcrossprod(X[t, ], X[t - j, ]) +
             tcrossprod(X[t - j, ], X[t, ]))
      }
    }
  }
  bread <- solve(crossprod(X))
  n / (n - k) * bread %*% meat %*% bread
}

wind_int <- grep("wind.*post|post.*wind", names(coef(m)), value = TRUE)
beta_wind <- coef(m)[wind_int]

cat(sprintf("  %-12s | %8s | %8s | %8s\n",
    "Bandwidth", "HAC SE", "HAC t", "HAC p"))
cat(paste(rep("-", 50), collapse = ""), "\n")

for (bw in c(1, 2, 3, 4, 6, 8, 10)) {
  vcov <- newey_west_vcov(m, bw)
  se <- sqrt(vcov[wind_int, wind_int])
  t_val <- beta_wind / se
  p_val <- 2 * pt(-abs(t_val), df = m$df.residual)
  cat(sprintf("  lag = %-5d | %8.4f | %8.3f | %8.4f%s\n",
      bw, se, t_val, p_val,
      ifelse(p_val < 0.05, " **", ifelse(p_val < 0.1, " *", ""))))
}

ols_se <- summary(m)$coefficients[wind_int, 2]
ols_p <- summary(m)$coefficients[wind_int, 4]
cat(sprintf("  OLS        | %8.4f | %8.3f | %8.4f%s\n",
    ols_se, beta_wind / ols_se, ols_p,
    ifelse(ols_p < 0.05, " **", "")))

opt_bw <- floor(n^(1/3))
cat(sprintf("\n  Plug-in optimal bandwidth (n^1/3): %d\n", opt_bw))

# Repeat for core weather model
cat("\n  --- Core Wind x PostMoU ---\n")
wind_int_core <- grep("wind_core.*post|post.*wind_core", names(coef(m_core)), value = TRUE)
if (length(wind_int_core) > 0) {
  beta_wind_core <- coef(m_core)[wind_int_core]
  cat(sprintf("  %-12s | %8s | %8s | %8s\n", "Bandwidth", "HAC SE", "HAC t", "HAC p"))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  for (bw in c(1, 2, 3, 4, 6, 8)) {
    vcov_c <- newey_west_vcov(m_core, bw)
    se_c <- sqrt(vcov_c[wind_int_core, wind_int_core])
    t_c <- beta_wind_core / se_c
    p_c <- 2 * pt(-abs(t_c), df = m_core$df.residual)
    cat(sprintf("  lag = %-5d | %8.4f | %8.3f | %8.4f%s\n",
        bw, se_c, t_c, p_c,
        ifelse(p_c < 0.05, " **", ifelse(p_c < 0.1, " *", ""))))
  }
}

# ============================================================
# 4. LAGGED DEPENDENT VARIABLE MODEL
# ============================================================
cat("\n============================================================\n")
cat("4. LAGGED DEPENDENT VARIABLE\n")
cat("============================================================\n\n")

df <- df %>% arrange(date) %>% mutate(log_rate_lag1 = lag(log_rate))

m_ldv <- lm(log_rate ~ log_rate_lag1 + wave_height_central_med * post_mou +
             wind_speed_central_med * post_mou + month_fac,
            data = df)

e_ldv <- residuals(m_ldv)
dw_orig <- sum(diff(e)^2) / sum(e^2)
dw_ldv <- sum(diff(e_ldv)^2) / sum(e_ldv^2)
acf1_ldv <- acf(e_ldv, lag.max = 1, plot = FALSE)$acf[2, 1, 1]

s_ldv <- summary(m_ldv)$coefficients

cat(sprintf("  LDV coef (lag1 of log_rate): %.3f (p = %.4f)\n",
    s_ldv["log_rate_lag1", 1], s_ldv["log_rate_lag1", 4]))
cat(sprintf("  DW: %.3f -> %.3f (with LDV)\n", dw_orig, dw_ldv))
cat(sprintf("  ACF(1): %.3f -> %.3f (with LDV)\n\n",
    acf_vals[1], acf1_ldv))

# Wind interaction in LDV model
wind_int_ldv <- grep("wind.*post|post.*wind", names(coef(m_ldv)), value = TRUE)
cat("  Wind x PostMoU in LDV model:\n")
cat(sprintf("    beta = %+.4f, OLS p = %.4f\n",
    s_ldv[wind_int_ldv, 1], s_ldv[wind_int_ldv, 4]))

vcov_ldv <- newey_west_vcov(m_ldv, 4)
se_ldv <- sqrt(vcov_ldv[wind_int_ldv, wind_int_ldv])
p_ldv <- 2 * pt(-abs(s_ldv[wind_int_ldv, 1] / se_ldv), df = m_ldv$df.residual)
cat(sprintf("    HAC(4) p = %.4f\n", p_ldv))

# Ljung-Box on LDV residuals
lb_ldv <- Box.test(e_ldv, lag = 6, type = "Ljung-Box")
cat(sprintf("    Ljung-Box(6) on LDV residuals: Q = %.3f, p = %.4f\n",
    lb_ldv$statistic, lb_ldv$p.value))

# Core weather LDV
cat("\n  --- Core weather LDV ---\n")
df_era5 <- df_era5 %>% arrange(date) %>% mutate(log_rate_lag1 = lag(log_rate))
m_ldv_core <- lm(log_rate ~ log_rate_lag1 + swh_core * post_mou +
                   wind_core * post_mou + month_fac, data = df_era5)
e_ldv_core <- residuals(m_ldv_core)
dw_ldv_core <- sum(diff(e_ldv_core)^2) / sum(e_ldv_core^2)
acf1_ldv_core <- acf(e_ldv_core, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
cat(sprintf("  DW: %.3f, ACF(1): %.3f\n", dw_ldv_core, acf1_ldv_core))

wind_int_ldv_core <- grep("wind_core.*post|post.*wind_core",
                            names(coef(m_ldv_core)), value = TRUE)
if (length(wind_int_ldv_core) > 0) {
  s_lc <- summary(m_ldv_core)$coefficients
  cat(sprintf("  Wind_core x Post in LDV: beta=%+.4f, p=%.4f\n",
      s_lc[wind_int_ldv_core, 1], s_lc[wind_int_ldv_core, 4]))
}

# ============================================================
# 5. RESIDUAL PATTERNS BY PERIOD
# ============================================================
cat("\n============================================================\n")
cat("5. RESIDUAL PATTERNS BY PERIOD\n")
cat("============================================================\n\n")

df$resid <- e
df_era5$resid_core <- e_core
df_era5$resid_nb <- e_nb

cat("  --- Broad OLS ---\n")
pre <- df %>% filter(post_mou == 0)
post <- df %>% filter(post_mou == 1)
acf1_pre <- acf(pre$resid, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
acf1_post <- acf(post$resid, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
cat(sprintf("  Pre-MoU:  mean = %+.4f, SD = %.4f, ACF(1) = %.3f (n=%d)\n",
    mean(pre$resid), sd(pre$resid), acf1_pre, nrow(pre)))
cat(sprintf("  Post-MoU: mean = %+.4f, SD = %.4f, ACF(1) = %.3f (n=%d)\n",
    mean(post$resid), sd(post$resid), acf1_post, nrow(post)))

cat("\n  --- Core OLS ---\n")
pre_c <- df_era5 %>% filter(post_mou == 0)
post_c <- df_era5 %>% filter(post_mou == 1)
acf1_pre_c <- acf(pre_c$resid_core, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
acf1_post_c <- acf(post_c$resid_core, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
cat(sprintf("  Pre-MoU:  mean = %+.4f, SD = %.4f, ACF(1) = %.3f (n=%d)\n",
    mean(pre_c$resid_core), sd(pre_c$resid_core), acf1_pre_c, nrow(pre_c)))
cat(sprintf("  Post-MoU: mean = %+.4f, SD = %.4f, ACF(1) = %.3f (n=%d)\n",
    mean(post_c$resid_core), sd(post_c$resid_core), acf1_post_c, nrow(post_c)))

cat("\n  --- Core NegBin (Pearson residuals) ---\n")
acf1_pre_nb <- acf(pre_c$resid_nb, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
acf1_post_nb <- acf(post_c$resid_nb, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
cat(sprintf("  Pre-MoU:  mean = %+.4f, SD = %.4f, ACF(1) = %.3f (n=%d)\n",
    mean(pre_c$resid_nb), sd(pre_c$resid_nb), acf1_pre_nb, nrow(pre_c)))
cat(sprintf("  Post-MoU: mean = %+.4f, SD = %.4f, ACF(1) = %.3f (n=%d)\n",
    mean(post_c$resid_nb), sd(post_c$resid_nb), acf1_post_nb, nrow(post_c)))

# Runs test
signs <- ifelse(df$resid > 0, "+", "-")
runs <- rle(signs)
n_runs <- length(runs$lengths)
n_pos <- sum(df$resid > 0)
n_neg <- sum(df$resid <= 0)
exp_runs <- 1 + 2 * n_pos * n_neg / n

cat(sprintf("\n  Runs: %d observed (%.1f expected). Mean length: %.1f. Max: %d months.\n",
    n_runs, exp_runs, mean(runs$lengths), max(runs$lengths)))

# ============================================================
# 6. EFFECTIVE SAMPLE SIZE
# ============================================================
cat("\n============================================================\n")
cat("6. EFFECTIVE SAMPLE SIZE\n")
cat("============================================================\n\n")

for (ml in list(
  list("Broad OLS", n, acf_vals[1]),
  list("Core OLS", n_core, acf_core_vals[1]),
  list("Core NegBin", n_nb, acf_nb_vals[1])
)) {
  rho1 <- ml[[3]]
  nn <- ml[[2]]
  n_eff <- nn * (1 - rho1) / (1 + rho1)
  cat(sprintf("  %s: n=%d, ACF(1)=%.3f, n_eff=%.0f, loss=%.0f%%\n",
      ml[[1]], nn, rho1, n_eff, 100 * (1 - n_eff / nn)))
}

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
