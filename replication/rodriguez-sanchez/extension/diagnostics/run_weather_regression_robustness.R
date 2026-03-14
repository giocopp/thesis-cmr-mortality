# run_weather_regression_robustness.R
# ====================================
# Robustness checks for weather-mortality monthly regressions.
# Tests: HAC standard errors (manual), zero-inflation, functional form,
# mediator conditioning, degrees of freedom.

library(dplyr)
library(lubridate)

df_model <- targets::tar_read(df_model, store = "_targets")
df_full  <- df_model$df_full

PRE_START <- as.Date("2011-02-01")
END_DATE  <- as.Date("2021-09-01")
MOU_DATE  <- as.Date("2017-02-01")

df <- df_full %>%
  filter(date >= PRE_START & date <= END_DATE) %>%
  mutate(
    log_rate     = log(mortality_rate_100 + 0.01),
    log_rate_1   = log(mortality_rate_100 + 1),     # alternative constant
    log_deaths   = log(dead_and_missing_Central_Mediterranean + 1),
    log_cross    = log(crossings_CMR + 1),
    deaths_raw   = dead_and_missing_Central_Mediterranean,
    rate_raw     = mortality_rate_100,
    post_mou     = as.integer(date >= MOU_DATE),
    month_factor = factor(month(date))
  )

weather_vars <- c(
  "wave_height_central_med",
  "wind_speed_central_med",
  "current_speed_central_med",
  "wave_days_above_2m"
)

labels <- c(
  wave_height_central_med = "Wave height",
  wind_speed_central_med = "Wind speed",
  current_speed_central_med = "Current speed",
  wave_days_above_2m = "Days >2m"
)

# --- Manual Newey-West HAC SE (no external packages needed) ---
newey_west_se <- function(model, max_lag = 4) {
  X <- model.matrix(model)
  e <- residuals(model)
  n <- length(e)
  k <- ncol(X)

  # HC0 "meat"
  meat <- matrix(0, k, k)
  for (j in 0:max_lag) {
    w <- if (j == 0) 1 else 1 - j / (max_lag + 1)  # Bartlett kernel
    for (t in (j + 1):n) {
      if (j == 0) {
        meat <- meat + e[t]^2 * tcrossprod(X[t, ])
      } else {
        meat <- meat + w * e[t] * e[t - j] * (tcrossprod(X[t, ], X[t - j, ]) +
                                                 tcrossprod(X[t - j, ], X[t, ]))
      }
    }
  }

  bread <- solve(crossprod(X))
  vcov_hac <- n / (n - k) * bread %*% meat %*% bread
  sqrt(diag(vcov_hac))
}

# =========================================================================
# CHECK 1: Zero-death months and floor effects
# =========================================================================
message("========================================")
message("CHECK 1: Zero-death months")
message("========================================\n")

n_zero_deaths <- sum(df$deaths_raw == 0, na.rm = TRUE)
n_zero_rate   <- sum(df$rate_raw == 0, na.rm = TRUE)
message("Total months: ", nrow(df))
message("Months with zero deaths: ", n_zero_deaths,
  " (", round(100 * n_zero_deaths / nrow(df), 1), "%)")
message("Months with zero rate: ", n_zero_rate,
  " (", round(100 * n_zero_rate / nrow(df), 1), "%)")
message("Pre-MoU zero-death months: ",
  sum(df$deaths_raw[df$post_mou == 0] == 0, na.rm = TRUE))
message("Post-MoU zero-death months: ",
  sum(df$deaths_raw[df$post_mou == 1] == 0, na.rm = TRUE))

# =========================================================================
# CHECK 2: Autocorrelation in residuals
# =========================================================================
message("\n========================================")
message("CHECK 2: Autocorrelation (Durbin-Watson)")
message("========================================\n")

for (v in weather_vars) {
  fml <- as.formula(paste0("log_rate ~ ", v, " * post_mou + month_factor"))
  m <- lm(fml, data = df)
  dw <- sum(diff(residuals(m))^2) / sum(residuals(m)^2)
  message(sprintf("  %-20s  DW = %.3f %s",
    labels[v], dw,
    ifelse(dw < 1.5, " (positive autocorrelation)", "")))
}

# =========================================================================
# CHECK 3: OLS vs HAC (Newey-West) standard errors
# =========================================================================
message("\n========================================")
message("CHECK 3: OLS vs Newey-West HAC standard errors")
message("========================================")
message("Rate interaction: weather × post_mou\n")

message(sprintf("  %-15s | %8s %7s %2s | %8s %7s %2s",
  "Variable", "OLS SE", "OLS p", "", "HAC SE", "HAC p", ""))
message(paste0("  ", paste(rep("-", 60), collapse = "")))

for (v in weather_vars) {
  interact_name <- paste0(v, ":post_mou")
  fml <- as.formula(paste0("log_rate ~ ", v, " * post_mou + month_factor"))
  m <- lm(fml, data = df)
  coef_idx <- which(names(coef(m)) == interact_name)

  # OLS
  s_ols <- summary(m)$coefficients[interact_name, ]
  beta  <- s_ols[1]

  # HAC
  hac_ses <- newey_west_se(m, max_lag = 4)
  hac_se  <- hac_ses[coef_idx]
  hac_t   <- beta / hac_se
  hac_p   <- 2 * pt(-abs(hac_t), df = m$df.residual)

  message(sprintf("  %-15s | %8.4f %7.4f%s | %8.4f %7.4f%s",
    labels[v],
    s_ols[2], s_ols[4], ifelse(s_ols[4] < 0.05, "*", " "),
    hac_se, hac_p, ifelse(hac_p < 0.05, "*", " ")))
}

message("\nDeaths interaction (volume-controlled):\n")

message(sprintf("  %-15s | %8s %7s %2s | %8s %7s %2s",
  "Variable", "OLS SE", "OLS p", "", "HAC SE", "HAC p", ""))
message(paste0("  ", paste(rep("-", 60), collapse = "")))

for (v in weather_vars) {
  interact_name <- paste0(v, ":post_mou")
  fml <- as.formula(paste0("log_deaths ~ ", v, " * post_mou + log_cross + month_factor"))
  m <- lm(fml, data = df)
  coef_idx <- which(names(coef(m)) == interact_name)

  s_ols <- summary(m)$coefficients[interact_name, ]
  beta  <- s_ols[1]

  hac_ses <- newey_west_se(m, max_lag = 4)
  hac_se  <- hac_ses[coef_idx]
  hac_t   <- beta / hac_se
  hac_p   <- 2 * pt(-abs(hac_t), df = m$df.residual)

  message(sprintf("  %-15s | %8.4f %7.4f%s | %8.4f %7.4f%s",
    labels[v],
    s_ols[2], s_ols[4], ifelse(s_ols[4] < 0.05, "*", " "),
    hac_se, hac_p, ifelse(hac_p < 0.05, "*", " ")))
}

# =========================================================================
# CHECK 4: Sensitivity to functional form
# =========================================================================
message("\n========================================")
message("CHECK 4: Sensitivity to log(rate + c) constant")
message("========================================")
message("Rate interaction with c = 0.01 vs c = 1\n")

message(sprintf("  %-15s | %10s %7s %2s | %10s %7s %2s",
  "Variable", "c=0.01 β", "p", "", "c=1 β", "p", ""))
message(paste0("  ", paste(rep("-", 60), collapse = "")))

for (v in weather_vars) {
  interact_name <- paste0(v, ":post_mou")

  fml1 <- as.formula(paste0("log_rate ~ ", v, " * post_mou + month_factor"))
  fml2 <- as.formula(paste0("log_rate_1 ~ ", v, " * post_mou + month_factor"))

  s1 <- summary(lm(fml1, data = df))$coefficients[interact_name, ]
  s2 <- summary(lm(fml2, data = df))$coefficients[interact_name, ]

  message(sprintf("  %-15s | %+10.4f %7.4f%s | %+10.4f %7.4f%s",
    labels[v],
    s1[1], s1[4], ifelse(s1[4] < 0.05, "*", " "),
    s2[1], s2[4], ifelse(s2[4] < 0.05, "*", " ")))
}

# =========================================================================
# CHECK 5: Rate model (no mediator conditioning)
# =========================================================================
message("\n========================================")
message("CHECK 5: Mediator conditioning assessment")
message("========================================\n")

message("Rate model: log(rate) ~ weather × post_mou + month FE")
message("  → Rate = deaths/crossings. Volume enters via the outcome,")
message("    NOT as a covariate. No mediator conditioning problem.\n")
message("Deaths model: log(deaths) ~ weather × post_mou + log(cross) + month FE")
message("  → log(crossings) is a POST-TREATMENT covariate (MoU → fewer crossings).")
message("    Conditioning on it = conditioning on a mediator.")
message("    The interaction β may be biased.\n")
message("CONCLUSION: trust the RATE model for the interaction test.")
message("The deaths + volume model is informative but has the mediator issue.")

# =========================================================================
# CHECK 6: Degrees of freedom
# =========================================================================
message("\n========================================")
message("CHECK 6: Degrees of freedom")
message("========================================\n")

message("Pooled interaction model (n = ", nrow(df), "):")
for (v in weather_vars) {
  fml <- as.formula(paste0("log_rate ~ ", v, " * post_mou + month_factor"))
  m <- lm(fml, data = df)
  message(sprintf("  %-15s  df.residual = %d  (k = %d)",
    labels[v], m$df.residual, length(coef(m))))
}

message("\nPost-MoU only (n = ", sum(df$post_mou == 1), "):")
df_post_only <- df %>% filter(post_mou == 1)
for (v in weather_vars) {
  fml <- as.formula(paste0("log_rate ~ ", v, " + month_factor"))
  m <- lm(fml, data = df_post_only)
  message(sprintf("  %-15s  df.residual = %d  (k = %d)",
    labels[v], m$df.residual, length(coef(m))))
}

# =========================================================================
# CHECK 7: Quarter dummies instead of month dummies (save df)
# =========================================================================
message("\n========================================")
message("CHECK 7: Quarter dummies instead of month (conserve df)")
message("========================================\n")

df$quarter_factor <- factor(quarter(df$date))

message(sprintf("  %-15s | %12s %7s | %12s %7s",
  "Variable", "Month FE β", "p", "Quarter FE β", "p"))
message(paste0("  ", paste(rep("-", 64), collapse = "")))

for (v in weather_vars) {
  interact_name <- paste0(v, ":post_mou")

  fml_m <- as.formula(paste0("log_rate ~ ", v, " * post_mou + month_factor"))
  fml_q <- as.formula(paste0("log_rate ~ ", v, " * post_mou + quarter_factor"))

  s_m <- summary(lm(fml_m, data = df))$coefficients[interact_name, ]
  s_q <- summary(lm(fml_q, data = df))$coefficients[interact_name, ]

  message(sprintf("  %-15s | %+10.4f  %6.4f%s | %+10.4f  %6.4f%s",
    labels[v],
    s_m[1], s_m[4], ifelse(s_m[4] < 0.05, "*", " "),
    s_q[1], s_q[4], ifelse(s_q[4] < 0.05, "*", " ")))
}

# =========================================================================
# CHECK 8: Multiple testing (Bonferroni / Holm)
# =========================================================================
message("\n========================================")
message("CHECK 8: Multiple testing correction")
message("========================================\n")

# Collect all interaction p-values for rate model
rate_ps <- c()
for (v in weather_vars) {
  interact_name <- paste0(v, ":post_mou")
  fml <- as.formula(paste0("log_rate ~ ", v, " * post_mou + month_factor"))
  s <- summary(lm(fml, data = df))$coefficients[interact_name, ]
  rate_ps <- c(rate_ps, s[4])
}
names(rate_ps) <- labels[weather_vars]

bonf_ps <- p.adjust(rate_ps, method = "bonferroni")
holm_ps <- p.adjust(rate_ps, method = "holm")

message(sprintf("  %-15s | %8s | %8s | %8s",
  "Variable", "Raw p", "Bonf. p", "Holm p"))
message(paste0("  ", paste(rep("-", 50), collapse = "")))
for (i in seq_along(rate_ps)) {
  message(sprintf("  %-15s | %8.4f | %8.4f | %8.4f%s",
    names(rate_ps)[i], rate_ps[i], bonf_ps[i], holm_ps[i],
    ifelse(holm_ps[i] < 0.05, " *", "")))
}
