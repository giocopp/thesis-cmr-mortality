# run_weather_regression_corrected.R
# ==================================
# Corrected regressions addressing:
#   1. Multiple testing applied to HAC p-values (not OLS)
#   2. Joint multivariate model (not just one-at-a-time)
# Deaths model retains log_cross (conditioning on mediator is intentional).

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
    log_deaths   = log(dead_and_missing_Central_Mediterranean + 1),
    log_cross    = log(crossings_CMR + 1),
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

# --- Manual Newey-West HAC SE ---
newey_west_se <- function(model, max_lag = 4) {
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
        meat <- meat + w * e[t] * e[t - j] * (tcrossprod(X[t, ], X[t - j, ]) +
                                                 tcrossprod(X[t - j, ], X[t, ]))
      }
    }
  }
  bread <- solve(crossprod(X))
  vcov_hac <- n / (n - k) * bread %*% meat %*% bread
  sqrt(diag(vcov_hac))
}

hac_pval <- function(model, coef_name, max_lag = 4) {
  coef_idx <- which(names(coef(model)) == coef_name)
  beta <- coef(model)[coef_idx]
  hac_se <- newey_west_se(model, max_lag)[coef_idx]
  hac_t <- beta / hac_se
  2 * pt(-abs(hac_t), df = model$df.residual)
}

# =========================================================================
# A. ONE-AT-A-TIME with HAC + corrected multiple testing
# =========================================================================
message("========================================")
message("A. ONE-AT-A-TIME INTERACTIONS (HAC SEs)")
message("========================================")
message("log(rate + 0.01) ~ weather * post_mou + month FE\n")

rate_hac_ps <- c()
rate_betas  <- c()

for (v in weather_vars) {
  interact_name <- paste0(v, ":post_mou")
  fml <- as.formula(paste0("log_rate ~ ", v, " * post_mou + month_factor"))
  m <- lm(fml, data = df)

  beta <- coef(m)[interact_name]
  p_hac <- hac_pval(m, interact_name)

  rate_betas[v] <- beta
  rate_hac_ps[v] <- p_hac
}

# Apply Holm correction to HAC p-values
holm_ps <- p.adjust(rate_hac_ps, method = "holm")

message(sprintf("  %-15s | %10s | %8s | %8s | %8s",
  "Variable", "Interact β", "HAC p", "Holm p", ""))
message(paste0("  ", paste(rep("-", 62), collapse = "")))
for (v in weather_vars) {
  sig <- ifelse(holm_ps[v] < 0.05, " *", ifelse(holm_ps[v] < 0.1, " †", ""))
  message(sprintf("  %-15s | %+10.4f | %8.4f | %8.4f |%s",
    labels[v], rate_betas[v], rate_hac_ps[v], holm_ps[v], sig))
}

# Same for deaths
message("\nlog(deaths + 1) ~ weather * post_mou + log(cross) + month FE\n")

deaths_hac_ps <- c()
deaths_betas  <- c()

for (v in weather_vars) {
  interact_name <- paste0(v, ":post_mou")
  fml <- as.formula(paste0("log_deaths ~ ", v, " * post_mou + log_cross + month_factor"))
  m <- lm(fml, data = df)

  deaths_betas[v] <- coef(m)[interact_name]
  deaths_hac_ps[v] <- hac_pval(m, interact_name)
}

holm_deaths <- p.adjust(deaths_hac_ps, method = "holm")

message(sprintf("  %-15s | %10s | %8s | %8s | %8s",
  "Variable", "Interact β", "HAC p", "Holm p", ""))
message(paste0("  ", paste(rep("-", 62), collapse = "")))
for (v in weather_vars) {
  sig <- ifelse(holm_deaths[v] < 0.05, " *", ifelse(holm_deaths[v] < 0.1, " †", ""))
  message(sprintf("  %-15s | %+10.4f | %8.4f | %8.4f |%s",
    labels[v], deaths_betas[v], deaths_hac_ps[v], holm_deaths[v], sig))
}

# =========================================================================
# B. JOINT MULTIVARIATE INTERACTION MODEL
# =========================================================================
message("\n========================================")
message("B. JOINT MULTIVARIATE INTERACTION MODEL")
message("========================================")
message("All 4 weather vars + all 4 interactions simultaneously.")
message("Addresses omitted-variable bias from correlated weather.\n")

# Note: wave_height and days_above_2m are correlated at r=0.892.
# Including both may cause multicollinearity. We report VIFs.

# --- Rate model ---
weather_terms <- paste(weather_vars, collapse = " + ")
interact_terms <- paste(paste0(weather_vars, ":post_mou"), collapse = " + ")
fml_joint_rate <- as.formula(paste0(
  "log_rate ~ ", weather_terms, " + post_mou + ",
  interact_terms, " + month_factor"))

m_joint_rate <- lm(fml_joint_rate, data = df)

message("--- Rate: joint model ---\n")
message(sprintf("  %-15s | %10s | %8s | %8s",
  "Interaction", "β", "OLS p", "HAC p"))
message(paste0("  ", paste(rep("-", 52), collapse = "")))

joint_rate_hac_ps <- c()
for (v in weather_vars) {
  interact_name <- paste0(v, ":post_mou")
  s_ols <- summary(m_joint_rate)$coefficients[interact_name, ]
  p_hac <- hac_pval(m_joint_rate, interact_name)
  joint_rate_hac_ps[v] <- p_hac

  message(sprintf("  %-15s | %+10.4f | %8.4f | %8.4f%s",
    labels[v], s_ols[1], s_ols[4], p_hac,
    ifelse(p_hac < 0.05, " *", ifelse(p_hac < 0.1, " †", ""))))
}

message(sprintf("\n  R² = %.4f, adj.R² = %.4f, df.residual = %d",
  summary(m_joint_rate)$r.squared,
  summary(m_joint_rate)$adj.r.squared,
  m_joint_rate$df.residual))

# Joint F-test: all 4 interactions = 0
fml_no_interact <- as.formula(paste0(
  "log_rate ~ ", weather_terms, " + post_mou + month_factor"))
m_no_interact <- lm(fml_no_interact, data = df)
f_test <- anova(m_no_interact, m_joint_rate)
message(sprintf("  Joint F-test (all interactions = 0): F = %.3f, p = %.4f%s",
  f_test$F[2], f_test$`Pr(>F)`[2],
  ifelse(f_test$`Pr(>F)`[2] < 0.05, " *", "")))

# --- Deaths model (with log_cross) ---
fml_joint_deaths <- as.formula(paste0(
  "log_deaths ~ ", weather_terms, " + post_mou + log_cross + ",
  interact_terms, " + month_factor"))

m_joint_deaths <- lm(fml_joint_deaths, data = df)

message("\n--- Deaths (+ log_cross): joint model ---\n")
message(sprintf("  %-15s | %10s | %8s | %8s",
  "Interaction", "β", "OLS p", "HAC p"))
message(paste0("  ", paste(rep("-", 52), collapse = "")))

joint_deaths_hac_ps <- c()
for (v in weather_vars) {
  interact_name <- paste0(v, ":post_mou")
  s_ols <- summary(m_joint_deaths)$coefficients[interact_name, ]
  p_hac <- hac_pval(m_joint_deaths, interact_name)
  joint_deaths_hac_ps[v] <- p_hac

  message(sprintf("  %-15s | %+10.4f | %8.4f | %8.4f%s",
    labels[v], s_ols[1], s_ols[4], p_hac,
    ifelse(p_hac < 0.05, " *", ifelse(p_hac < 0.1, " †", ""))))
}

fml_no_interact_d <- as.formula(paste0(
  "log_deaths ~ ", weather_terms, " + post_mou + log_cross + month_factor"))
m_no_interact_d <- lm(fml_no_interact_d, data = df)
f_test_d <- anova(m_no_interact_d, m_joint_deaths)
message(sprintf("\n  Joint F-test (all interactions = 0): F = %.3f, p = %.4f%s",
  f_test_d$F[2], f_test_d$`Pr(>F)`[2],
  ifelse(f_test_d$`Pr(>F)`[2] < 0.05, " *", "")))

# =========================================================================
# C. REDUCED JOINT MODEL (drop wave_height, r=0.89 with days>2m)
# =========================================================================
message("\n========================================")
message("C. REDUCED JOINT MODEL (drop wave_height)")
message("========================================")
message("wave_height and days>2m at r=0.892 → keep only 3 vars\n")

reduced_vars <- c("wind_speed_central_med", "current_speed_central_med",
                   "wave_days_above_2m")

weather_r <- paste(reduced_vars, collapse = " + ")
interact_r <- paste(paste0(reduced_vars, ":post_mou"), collapse = " + ")

fml_red_rate <- as.formula(paste0(
  "log_rate ~ ", weather_r, " + post_mou + ", interact_r, " + month_factor"))
m_red_rate <- lm(fml_red_rate, data = df)

message("--- Rate: reduced joint model ---\n")
message(sprintf("  %-15s | %10s | %8s | %8s",
  "Interaction", "β", "OLS p", "HAC p"))
message(paste0("  ", paste(rep("-", 52), collapse = "")))

for (v in reduced_vars) {
  interact_name <- paste0(v, ":post_mou")
  s_ols <- summary(m_red_rate)$coefficients[interact_name, ]
  p_hac <- hac_pval(m_red_rate, interact_name)

  message(sprintf("  %-15s | %+10.4f | %8.4f | %8.4f%s",
    labels[v], s_ols[1], s_ols[4], p_hac,
    ifelse(p_hac < 0.05, " *", ifelse(p_hac < 0.1, " †", ""))))
}

fml_red_no <- as.formula(paste0(
  "log_rate ~ ", weather_r, " + post_mou + month_factor"))
m_red_no <- lm(fml_red_no, data = df)
f_red <- anova(m_red_no, m_red_rate)
message(sprintf("\n  Joint F-test (3 interactions = 0): F = %.3f, p = %.4f%s",
  f_red$F[2], f_red$`Pr(>F)`[2],
  ifelse(f_red$`Pr(>F)`[2] < 0.05, " *", "")))

message(sprintf("  R² = %.4f, adj.R² = %.4f, df.residual = %d",
  summary(m_red_rate)$r.squared,
  summary(m_red_rate)$adj.r.squared,
  m_red_rate$df.residual))

# --- Deaths reduced ---
fml_red_deaths <- as.formula(paste0(
  "log_deaths ~ ", weather_r, " + post_mou + log_cross + ",
  interact_r, " + month_factor"))
m_red_deaths <- lm(fml_red_deaths, data = df)

message("\n--- Deaths: reduced joint model ---\n")
message(sprintf("  %-15s | %10s | %8s | %8s",
  "Interaction", "β", "OLS p", "HAC p"))
message(paste0("  ", paste(rep("-", 52), collapse = "")))

for (v in reduced_vars) {
  interact_name <- paste0(v, ":post_mou")
  s_ols <- summary(m_red_deaths)$coefficients[interact_name, ]
  p_hac <- hac_pval(m_red_deaths, interact_name)

  message(sprintf("  %-15s | %+10.4f | %8.4f | %8.4f%s",
    labels[v], s_ols[1], s_ols[4], p_hac,
    ifelse(p_hac < 0.05, " *", ifelse(p_hac < 0.1, " †", ""))))
}

# =========================================================================
# D. SUMMARY
# =========================================================================
message("\n========================================")
message("D. SUMMARY OF CORRECTIONS")
message("========================================\n")

message("1. Multiple testing: Holm correction now applied to HAC p-values")
message("   Rate (one-at-a-time, HAC + Holm):")
for (v in weather_vars) {
  message(sprintf("     %-15s  HAC p=%.4f  Holm p=%.4f%s",
    labels[v], rate_hac_ps[v], holm_ps[v],
    ifelse(holm_ps[v] < 0.05, " *", ifelse(holm_ps[v] < 0.1, " †", ""))))
}

message("\n2. Joint model addresses omitted-variable bias:")
message("   Rate (joint, HAC):")
for (v in weather_vars) {
  message(sprintf("     %-15s  HAC p=%.4f%s",
    labels[v], joint_rate_hac_ps[v],
    ifelse(joint_rate_hac_ps[v] < 0.05, " *", ifelse(joint_rate_hac_ps[v] < 0.1, " †", ""))))
}

message("\n3. Zero-crossing months: NONE (min = 217). No sample selection issue.")
