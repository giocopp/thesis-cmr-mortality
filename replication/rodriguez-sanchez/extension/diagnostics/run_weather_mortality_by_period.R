# run_weather_mortality_by_period.R
# =================================
# Split-sample regressions: does weather predict mortality
# differently pre- vs post-MoU?

library(dplyr)
library(lubridate)

message("Loading data from targets store...")
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

df_pre  <- df %>% filter(post_mou == 0)
df_post <- df %>% filter(post_mou == 1)

message("N pre-MoU: ", nrow(df_pre), " months")
message("N post-MoU: ", nrow(df_post), " months")

weather_vars <- c(
  "wave_height_central_med",
  "wind_speed_central_med",
  "sst_central_med",
  "current_speed_central_med",
  "wave_max_central_med",
  "wave_days_above_2m"
)
available <- weather_vars[weather_vars %in% names(df)]

# =========================================================================
# A. SPLIT-SAMPLE: Rate ~ weather + month_dummies (per period)
# =========================================================================
message("\n================================================")
message("A. MORTALITY RATE: split-sample (season-controlled)")
message("================================================")
message("log(rate + 0.01) ~ weather + month_dummies\n")

message(sprintf("  %-28s | %12s %8s | %12s %8s",
  "Variable", "Pre-MoU β", "p", "Post-MoU β", "p"))
message(paste0("  ", paste(rep("-", 80), collapse = "")))

for (v in available) {
  fml <- as.formula(paste0("log_rate ~ ", v, " + month_factor"))

  m_pre  <- lm(fml, data = df_pre)
  m_post <- lm(fml, data = df_post)

  c_pre  <- summary(m_pre)$coefficients[v, ]
  c_post <- summary(m_post)$coefficients[v, ]

  message(sprintf("  %-28s | %+10.4f  %7.4f%s | %+10.4f  %7.4f%s",
    v,
    c_pre[1], c_pre[4], ifelse(c_pre[4] < 0.05, "*", " "),
    c_post[1], c_post[4], ifelse(c_post[4] < 0.05, "*", " ")))
}

# =========================================================================
# B. SPLIT-SAMPLE: Deaths ~ weather + log(cross) + month_dummies
# =========================================================================
message("\n================================================")
message("B. DEATHS: split-sample (season + volume controlled)")
message("================================================")
message("log(deaths + 1) ~ weather + log(crossings) + month_dummies\n")

message(sprintf("  %-28s | %12s %8s | %12s %8s",
  "Variable", "Pre-MoU β", "p", "Post-MoU β", "p"))
message(paste0("  ", paste(rep("-", 80), collapse = "")))

for (v in available) {
  fml <- as.formula(paste0("log_deaths ~ ", v, " + log_cross + month_factor"))

  m_pre  <- lm(fml, data = df_pre)
  m_post <- lm(fml, data = df_post)

  c_pre  <- summary(m_pre)$coefficients[v, ]
  c_post <- summary(m_post)$coefficients[v, ]

  message(sprintf("  %-28s | %+10.4f  %7.4f%s | %+10.4f  %7.4f%s",
    v,
    c_pre[1], c_pre[4], ifelse(c_pre[4] < 0.05, "*", " "),
    c_post[1], c_post[4], ifelse(c_post[4] < 0.05, "*", " ")))
}

# =========================================================================
# C. POOLED INTERACTION (formal test of difference)
# =========================================================================
message("\n================================================")
message("C. POOLED INTERACTION: formal test of period difference")
message("================================================")
message("Includes post_mou main effect + weather × post_mou interaction\n")

message("--- Rate model: log(rate) ~ weather * post_mou + month_dummies ---\n")

message(sprintf("  %-28s | %10s %7s | %12s %7s",
  "Variable", "Main β", "p", "Interact β", "p"))
message(paste0("  ", paste(rep("-", 76), collapse = "")))

for (v in available) {
  fml <- as.formula(paste0("log_rate ~ ", v, " * post_mou + month_factor"))
  m   <- lm(fml, data = df)
  s   <- summary(m)$coefficients
  interact_name <- paste0(v, ":post_mou")

  main <- s[v, ]
  inter <- s[interact_name, ]

  message(sprintf("  %-28s | %+8.4f  %6.4f%s | %+10.4f  %6.4f%s",
    v,
    main[1], main[4], ifelse(main[4] < 0.05, "*", " "),
    inter[1], inter[4], ifelse(inter[4] < 0.05, "*", " ")))
}

message("\n--- Deaths model: log(deaths) ~ weather * post_mou + log(cross) + month_dummies ---\n")

message(sprintf("  %-28s | %10s %7s | %12s %7s",
  "Variable", "Main β", "p", "Interact β", "p"))
message(paste0("  ", paste(rep("-", 76), collapse = "")))

for (v in available) {
  fml <- as.formula(paste0("log_deaths ~ ", v, " * post_mou + log_cross + month_factor"))
  m   <- lm(fml, data = df)
  s   <- summary(m)$coefficients
  interact_name <- paste0(v, ":post_mou")

  main <- s[v, ]
  inter <- s[interact_name, ]

  message(sprintf("  %-28s | %+8.4f  %6.4f%s | %+10.4f  %6.4f%s",
    v,
    main[1], main[4], ifelse(main[4] < 0.05, "*", " "),
    inter[1], inter[4], ifelse(inter[4] < 0.05, "*", " ")))
}

# =========================================================================
# D. R² comparison: how much does weather add per period?
# =========================================================================
message("\n================================================")
message("D. R² COMPARISON: weather's marginal contribution per period")
message("================================================\n")

for (period_label in c("Pre-MoU", "Post-MoU")) {
  d <- if (period_label == "Pre-MoU") df_pre else df_post

  # Rate
  null_r  <- lm(log_rate ~ month_factor, data = d)
  full_r  <- lm(as.formula(paste0("log_rate ~ ", paste(available, collapse = " + "),
    " + month_factor")), data = d)
  dr2_r <- summary(full_r)$r.squared - summary(null_r)$r.squared
  f_r   <- anova(null_r, full_r)

  # Deaths
  null_d  <- lm(log_deaths ~ log_cross + month_factor, data = d)
  full_d  <- lm(as.formula(paste0("log_deaths ~ ", paste(available, collapse = " + "),
    " + log_cross + month_factor")), data = d)
  dr2_d <- summary(full_d)$r.squared - summary(null_d)$r.squared
  f_d   <- anova(null_d, full_d)

  message(sprintf("  %s:", period_label))
  message(sprintf("    Rate:   ΔR²=%.4f  F=%.2f  p=%.4f %s",
    dr2_r, f_r$F[2], f_r$`Pr(>F)`[2],
    ifelse(f_r$`Pr(>F)`[2] < 0.05, " *", "")))
  message(sprintf("    Deaths: ΔR²=%.4f  F=%.2f  p=%.4f %s",
    dr2_d, f_d$F[2], f_d$`Pr(>F)`[2],
    ifelse(f_d$`Pr(>F)`[2] < 0.05, " *", "")))
  message("")
}
