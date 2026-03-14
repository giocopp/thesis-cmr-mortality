# 07_volume_channel_evidence.R
# ============================
# Empirical evidence for the volume channel claims:
#   1. Moon illumination predicts crossing volume
#   2. The pre-period death count structural break is driven by crossings
#
# This script uses the targets cache — run after tar_make().

library(targets)
library(dplyr)
library(ggplot2)

# =========================================================================
# Load data
# =========================================================================

df_model <- tar_read(df_model)
df_full  <- df_model$df_full

# Restrict to analysis window
df <- df_full %>%
  filter(date >= "2011-02-01" & date < "2021-10-01") %>%
  mutate(
    log_crossings = log(crossings_CMR),
    log_deaths    = log(dead_and_missing_Central_Mediterranean + 1),
    log_rate      = log(mortality_rate_100 + 0.01),
    moon          = moon_illumination_frac,
    moon_lag1     = moon_illumination_frac_lag_01,
    period = case_when(
      date < "2013-10-01" ~ "1. Pre-MN",
      date < "2014-11-01" ~ "2. Mare Nostrum",
      date < "2017-02-01" ~ "3. NGO SAR",
      TRUE                ~ "4. Post-MoU"
    ),
    pre_mou = date < "2017-02-01"
  )

cat("=======================================================\n")
cat("VOLUME CHANNEL EVIDENCE\n")
cat("=======================================================\n\n")


# =========================================================================
# 1. Does moon illumination predict crossing volume?
# =========================================================================

cat("=== 1. MOON ILLUMINATION → CROSSINGS ===\n\n")

# Simple correlations
cor_moon_cross <- cor(df$moon_lag1, df$log_crossings, use = "complete.obs")
cor_moon_death <- cor(df$moon_lag1, df$log_deaths, use = "complete.obs")
cor_moon_rate  <- cor(df$moon_lag1, df$log_rate, use = "complete.obs")
cor_cross_death <- cor(df$log_crossings, df$log_deaths, use = "complete.obs")

cat(sprintf("  cor(moon_lag1, log_crossings) = %+.3f\n", cor_moon_cross))
cat(sprintf("  cor(moon_lag1, log_deaths)    = %+.3f\n", cor_moon_death))
cat(sprintf("  cor(moon_lag1, log_rate)      = %+.3f\n", cor_moon_rate))
cat(sprintf("  cor(log_crossings, log_deaths)= %+.3f\n", cor_cross_death))

# Regression: moon → crossings
reg_moon_cross <- lm(log_crossings ~ moon_lag1, data = df)
cat("\n  Regression: log_crossings ~ moon_lag1\n")
cat(sprintf("    coef = %.3f, p = %.4f, R² = %.3f\n",
            coef(reg_moon_cross)[2],
            summary(reg_moon_cross)$coefficients[2, 4],
            summary(reg_moon_cross)$r.squared))

# Regression: moon → deaths (total)
reg_moon_death <- lm(log_deaths ~ moon_lag1, data = df)
cat("\n  Regression: log_deaths ~ moon_lag1\n")
cat(sprintf("    coef = %.3f, p = %.4f, R² = %.3f\n",
            coef(reg_moon_death)[2],
            summary(reg_moon_death)$coefficients[2, 4],
            summary(reg_moon_death)$r.squared))

# Regression: moon → deaths CONTROLLING for crossings
reg_moon_death_cond <- lm(log_deaths ~ moon_lag1 + log_crossings, data = df)
cat("\n  Regression: log_deaths ~ moon_lag1 + log_crossings\n")
s <- summary(reg_moon_death_cond)
cat(sprintf("    moon_lag1:      coef = %+.3f, p = %.4f\n",
            s$coefficients[2, 1], s$coefficients[2, 4]))
cat(sprintf("    log_crossings:  coef = %+.3f, p = %.4f\n",
            s$coefficients[3, 1], s$coefficients[3, 4]))
cat(sprintf("    R² = %.3f\n", s$r.squared))

cat("\n  Interpretation: if moon_lag1 becomes non-significant after\n")
cat("  controlling for crossings, its effect works THROUGH volume.\n")
cat("  If it remains significant, there is also a direct channel.\n\n")


# =========================================================================
# 2. Is the pre-period death structural break driven by crossings?
# =========================================================================

cat("=== 2. PRE-PERIOD STRUCTURAL BREAK ===\n\n")

df_pre <- df %>% filter(pre_mou)

# 2a. Unconditional means by sub-period
cat("  Unconditional sub-period means (pre-MoU):\n")
sub_means <- df_pre %>%
  group_by(period) %>%
  summarise(
    n = n(),
    mean_log_cross = mean(log_crossings, na.rm = TRUE),
    mean_log_death = mean(log_deaths, na.rm = TRUE),
    mean_log_rate  = mean(log_rate, na.rm = TRUE),
    .groups = "drop"
  )
for (i in seq_len(nrow(sub_means))) {
  r <- sub_means[i, ]
  cat(sprintf("    %-18s  n=%2d  log_cross=%.2f  log_death=%.2f  log_rate=%+.2f\n",
              r$period, r$n, r$mean_log_cross, r$mean_log_death, r$mean_log_rate))
}

# 2b. Regression of deaths on period dummies (unconditional break)
df_pre$period_f <- factor(df_pre$period)
reg_break_uncond <- lm(log_deaths ~ period_f, data = df_pre)
cat("\n  Unconditional break test: log_deaths ~ period\n")
s <- summary(reg_break_uncond)
cat(sprintf("    R² = %.3f, F-stat p = %.4e\n", s$r.squared,
            pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)))
for (i in 2:nrow(s$coefficients)) {
  cat(sprintf("    %s: coef = %+.3f, p = %.4f\n",
              rownames(s$coefficients)[i], s$coefficients[i, 1], s$coefficients[i, 4]))
}

# 2c. Regression of deaths on period dummies CONTROLLING for crossings
reg_break_cond <- lm(log_deaths ~ period_f + log_crossings, data = df_pre)
cat("\n  Conditional break test: log_deaths ~ period + log_crossings\n")
s <- summary(reg_break_cond)
cat(sprintf("    R² = %.3f\n", s$r.squared))
for (i in 2:nrow(s$coefficients)) {
  cat(sprintf("    %s: coef = %+.3f, p = %.4f\n",
              rownames(s$coefficients)[i], s$coefficients[i, 1], s$coefficients[i, 4]))
}

cat("\n  Interpretation: if period dummies become non-significant after\n")
cat("  controlling for crossings, the structural break is EXPLAINED\n")
cat("  by the volume change. If they remain significant, there are\n")
cat("  additional regime changes beyond volume.\n\n")

# 2d. Same for rate (should show no break even unconditionally)
reg_rate_uncond <- lm(log_rate ~ period_f, data = df_pre)
cat("  Rate break test (unconditional): log_rate ~ period\n")
s <- summary(reg_rate_uncond)
cat(sprintf("    R² = %.3f, F-stat p = %.4f\n", s$r.squared,
            pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)))
for (i in 2:nrow(s$coefficients)) {
  cat(sprintf("    %s: coef = %+.3f, p = %.4f\n",
              rownames(s$coefficients)[i], s$coefficients[i, 1], s$coefficients[i, 4]))
}


# =========================================================================
# 3. Partial correlations: moon → deaths with and without crossings
# =========================================================================

cat("\n=== 3. PARTIAL CORRELATION ANALYSIS ===\n\n")

# Pre-MoU only
df_pre_clean <- df_pre %>%
  select(moon_lag1, log_crossings, log_deaths, log_rate) %>%
  na.omit()

cat("  Pre-MoU period (n =", nrow(df_pre_clean), "):\n")

cor_pre <- cor(df_pre_clean)
cat(sprintf("    cor(moon_lag1, log_crossings) = %+.3f\n",
            cor_pre["moon_lag1", "log_crossings"]))
cat(sprintf("    cor(moon_lag1, log_deaths)    = %+.3f\n",
            cor_pre["moon_lag1", "log_deaths"]))
cat(sprintf("    cor(moon_lag1, log_rate)      = %+.3f\n",
            cor_pre["moon_lag1", "log_rate"]))

# Partial correlation: moon → deaths | crossings
# = cor(resid(deaths ~ crossings), resid(moon ~ crossings))
resid_deaths <- residuals(lm(log_deaths ~ log_crossings, data = df_pre_clean))
resid_moon   <- residuals(lm(moon_lag1 ~ log_crossings, data = df_pre_clean))
partial_cor  <- cor(resid_deaths, resid_moon)

cat(sprintf("\n    Partial cor(moon, deaths | crossings) = %+.3f\n", partial_cor))
cat("    (= moon's residual association with deaths after removing\n")
cat("     the volume channel)\n")


# =========================================================================
# 4. Summary
# =========================================================================

cat("\n=== 4. SUMMARY ===\n\n")

cat("  Claim 1: Moon illumination predicts crossing volume\n")
cat(sprintf("    Evidence: cor = %+.3f, regression p = %.4f\n",
            cor_moon_cross,
            summary(reg_moon_cross)$coefficients[2, 4]))

moon_uncond_p <- summary(reg_moon_death)$coefficients[2, 4]
moon_cond_p   <- summary(reg_moon_death_cond)$coefficients[2, 4]
cat(sprintf("    Moon → deaths unconditional: p = %.4f\n", moon_uncond_p))
cat(sprintf("    Moon → deaths | crossings:   p = %.4f\n", moon_cond_p))
if (moon_cond_p > 0.05 & moon_uncond_p < 0.05) {
  cat("    → Moon effect works ENTIRELY through volume channel\n")
} else if (moon_cond_p < 0.05) {
  cat("    → Moon has BOTH volume and direct channels\n")
} else {
  cat("    → Moon does not significantly predict deaths\n")
}

cat(sprintf("\n  Claim 2: Death structural break driven by crossings\n"))
period_uncond_p <- summary(reg_break_uncond)$coefficients[2, 4]
period_cond_p   <- summary(reg_break_cond)$coefficients[2, 4]
cat(sprintf("    MN dummy unconditional: p = %.4f\n", period_uncond_p))
cat(sprintf("    MN dummy | crossings:   p = %.4f\n", period_cond_p))
if (period_cond_p > 0.05 & period_uncond_p < 0.05) {
  cat("    → Structural break EXPLAINED by volume change\n")
} else if (period_cond_p < 0.05) {
  cat("    → Structural break only PARTIALLY explained by volume\n")
} else {
  cat("    → No significant structural break even unconditionally\n")
}

cat("\n=======================================================\n")
cat("DONE\n")
cat("=======================================================\n")
