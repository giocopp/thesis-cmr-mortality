# autocorrelation_diagnostics.R
# =============================
# Detailed exploration of residual autocorrelation in the monthly
# weather-mortality interaction regression.

library(dplyr)
library(lubridate)

df <- readRDS("replication/rodriguez-sanchez/extension/data/df_extended.RDS") %>%
  filter(date >= as.Date("2011-02-01") & date <= as.Date("2021-09-01")) %>%
  mutate(log_rate = log(mortality_rate_100 + 0.01),
         post_mou = as.integer(date >= as.Date("2017-02-01")),
         month_fac = factor(month(date)))

m <- lm(log_rate ~ wave_height_central_med * post_mou +
         wind_speed_central_med * post_mou + month_fac, data = df)
e <- residuals(m)
n <- length(e)

# ============================================================
# 1. RESIDUAL ACF UP TO LAG 12
# ============================================================
cat("============================================================\n")
cat("1. RESIDUAL AUTOCORRELATION FUNCTION\n")
cat("============================================================\n\n")

acf_obj <- acf(e, lag.max = 12, plot = FALSE)
acf_vals <- acf_obj$acf[-1, 1, 1]  # drop lag 0
bound <- 1.96 / sqrt(n)

for (i in seq_along(acf_vals)) {
  cat(sprintf("  Lag %2d: %+.3f  %s\n", i, acf_vals[i],
      ifelse(abs(acf_vals[i]) > bound, " *", "")))
}
cat(sprintf("\n  95%% significance bound: +/- %.3f\n", bound))

# ============================================================
# 2. LJUNG-BOX TESTS
# ============================================================
cat("\n============================================================\n")
cat("2. LJUNG-BOX TEST FOR RESIDUAL AUTOCORRELATION\n")
cat("============================================================\n\n")

for (lag in c(3, 6, 12)) {
  lb <- Box.test(e, lag = lag, type = "Ljung-Box")
  cat(sprintf("  Lags 1-%d: Q = %.3f, p = %.4f%s\n",
      lag, lb$statistic, lb$p.value,
      ifelse(lb$p.value < 0.05, " **", "")))
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

# Optimal bandwidth: Newey-West (1994) plug-in rule ~ n^(1/3)
opt_bw <- floor(n^(1/3))
cat(sprintf("\n  Plug-in optimal bandwidth (n^1/3): %d\n", opt_bw))

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

# ============================================================
# 5. RESIDUAL PATTERNS BY PERIOD
# ============================================================
cat("\n============================================================\n")
cat("5. RESIDUAL PATTERNS BY PERIOD\n")
cat("============================================================\n\n")

df$resid <- e

pre <- df %>% filter(post_mou == 0)
post <- df %>% filter(post_mou == 1)

acf1_pre <- acf(pre$resid, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
acf1_post <- acf(post$resid, lag.max = 1, plot = FALSE)$acf[2, 1, 1]

cat(sprintf("  Pre-MoU:  mean = %+.4f, SD = %.4f, ACF(1) = %.3f (n=%d)\n",
    mean(pre$resid), sd(pre$resid), acf1_pre, nrow(pre)))
cat(sprintf("  Post-MoU: mean = %+.4f, SD = %.4f, ACF(1) = %.3f (n=%d)\n",
    mean(post$resid), sd(post$resid), acf1_post, nrow(post)))

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

rho1 <- acf_vals[1]
n_eff <- n * (1 - rho1) / (1 + rho1)
cat(sprintf("  Nominal n = %d\n", n))
cat(sprintf("  ACF(1) = %.3f\n", rho1))
cat(sprintf("  Effective n ~ n*(1-rho)/(1+rho) = %.0f\n", n_eff))
cat(sprintf("  Efficiency loss: %.0f%%\n", 100 * (1 - n_eff / n)))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
