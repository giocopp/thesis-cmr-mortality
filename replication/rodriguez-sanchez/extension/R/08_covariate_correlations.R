# 08_covariate_correlations.R
# ============================
# Correlations of exogenous covariates with:
#   1. log(deaths and missing + 1)
#   2. log(crossings)
#   3. log(mortality rate + 0.01)
#
# Key finding: covariates correlate with deaths in the WRONG direction
# (bad weather → fewer deaths) because they predict crossing volume,
# not per-crossing danger.
#
# This script uses the targets cache — run after tar_make().

library(targets)
library(dplyr)

# =========================================================================
# Load data
# =========================================================================

df_model <- tar_read(df_model)
cov_cols <- df_model$cov_cols

# Death count dataset
df_deaths <- df_model$df_deaths

# Mortality rate dataset
df_mort <- df_model$df_mort

# Full dataset for crossings
df_full <- df_model$df_full %>%
  filter(date >= "2011-02-01" & date < "2021-10-01") %>%
  mutate(log_crossings = log(crossings_CMR))


# =========================================================================
# 1. Correlations with log(deaths and missing + 1)
# =========================================================================

cat("=======================================================\n")
cat("COVARIATE CORRELATIONS\n")
cat("=======================================================\n\n")

cat("=== 1. CORRELATIONS WITH log(deaths+missing+1) ===\n\n")

cors_death <- sapply(cov_cols, function(v)
  cor(df_deaths[[v]], df_deaths$deaths_cmr, use = "complete.obs"))
pvals_death <- sapply(cov_cols, function(v)
  cor.test(df_deaths[[v]], df_deaths$deaths_cmr, use = "complete.obs")$p.value)

res_death <- data.frame(
  variable = names(cors_death),
  correlation = round(cors_death, 3),
  p_value = round(pvals_death, 4),
  abs_cor = round(abs(cors_death), 3),
  stringsAsFactors = FALSE
) %>% arrange(desc(abs_cor))

cat(sprintf("N = %d months\n\n", nrow(df_deaths)))

for (i in seq_len(nrow(res_death))) {
  r <- res_death[i, ]
  sig <- ifelse(r$p_value < 0.001, "***",
         ifelse(r$p_value < 0.01, "**",
         ifelse(r$p_value < 0.05, "*", "   ")))
  cat(sprintf("  %+.3f  p=%-8s %s  %s\n",
              r$correlation, format(r$p_value, digits = 4), sig, r$variable))
}

cat(sprintf("\n  # significant at 5%%: %d / %d\n\n",
            sum(res_death$p_value < 0.05), nrow(res_death)))


# =========================================================================
# 2. Correlations with log(crossings)
# =========================================================================

cat("=== 2. CORRELATIONS WITH log(crossings) ===\n\n")

cors_cross <- sapply(cov_cols, function(v) {
  if (v %in% names(df_full)) cor(df_full[[v]], df_full$log_crossings, use = "complete.obs") else NA
})
pvals_cross <- sapply(cov_cols, function(v) {
  if (v %in% names(df_full)) cor.test(df_full[[v]], df_full$log_crossings, use = "complete.obs")$p.value else NA
})

cat(sprintf("N = %d months\n\n", sum(!is.na(cors_cross))))

cat(sprintf("  %-42s  %8s  %8s  %8s\n",
            "Variable", "r(cross)", "p", "r(death)"))
cat(sprintf("  %-42s  %8s  %8s  %8s\n",
            paste(rep("-", 42), collapse = ""), "--------", "--------", "--------"))

res_combined <- data.frame(
  variable = names(cors_cross),
  cor_crossings = round(cors_cross, 3),
  p_crossings = round(pvals_cross, 4),
  cor_deaths = round(cors_death[names(cors_cross)], 3),
  stringsAsFactors = FALSE
) %>% arrange(desc(abs(cor_crossings)))

for (i in seq_len(nrow(res_combined))) {
  r <- res_combined[i, ]
  if (is.na(r$cor_crossings)) next
  sig <- ifelse(r$p_crossings < 0.001, "***",
         ifelse(r$p_crossings < 0.01, "**",
         ifelse(r$p_crossings < 0.05, "*", "   ")))
  cat(sprintf("  %-42s  %+.3f %s  %7s  %+.3f\n",
              r$variable, r$cor_crossings, sig,
              format(r$p_crossings, digits = 4), r$cor_deaths))
}

cat(sprintf("\n  # sig at 5%% with crossings: %d / %d\n",
            sum(res_combined$p_crossings < 0.05, na.rm = TRUE),
            sum(!is.na(res_combined$p_crossings))))


# =========================================================================
# 3. Correlations with log(mortality rate + 0.01)
# =========================================================================

cat("\n=== 3. CORRELATIONS WITH log(mortality_rate + 0.01) ===\n\n")

cors_rate <- sapply(cov_cols, function(v)
  cor(df_mort[[v]], df_mort$mortality_rate, use = "complete.obs"))
pvals_rate <- sapply(cov_cols, function(v)
  cor.test(df_mort[[v]], df_mort$mortality_rate, use = "complete.obs")$p.value)

res_rate <- data.frame(
  variable = names(cors_rate),
  correlation = round(cors_rate, 3),
  p_value = round(pvals_rate, 4),
  abs_cor = round(abs(cors_rate), 3),
  stringsAsFactors = FALSE
) %>% arrange(desc(abs_cor))

cat(sprintf("N = %d months\n\n", nrow(df_mort)))

for (i in seq_len(min(nrow(res_rate), 10))) {
  r <- res_rate[i, ]
  sig <- ifelse(r$p_value < 0.001, "***",
         ifelse(r$p_value < 0.01, "**",
         ifelse(r$p_value < 0.05, "*", "   ")))
  cat(sprintf("  %+.3f  p=%-8s %s  %s\n",
              r$correlation, format(r$p_value, digits = 4), sig, r$variable))
}

cat(sprintf("\n  # sig at 5%% with rate: %d / %d\n",
            sum(res_rate$p_value < 0.05), nrow(res_rate)))


# =========================================================================
# 4. Summary
# =========================================================================

cat("\n=== 4. SUMMARY ===\n\n")

cat(sprintf("  Covariates sig. with crossings:      %d / %d\n",
            sum(res_combined$p_crossings < 0.05, na.rm = TRUE),
            sum(!is.na(res_combined$p_crossings))))
cat(sprintf("  Covariates sig. with deaths:          %d / %d\n",
            sum(res_death$p_value < 0.05), nrow(res_death)))
cat(sprintf("  Covariates sig. with mortality rate:  %d / %d\n",
            sum(res_rate$p_value < 0.05), nrow(res_rate)))

cat(sprintf("\n  Max |r| with crossings:      %.3f\n",
            max(abs(res_combined$cor_crossings), na.rm = TRUE)))
cat(sprintf("  Max |r| with deaths:          %.3f\n",
            max(res_death$abs_cor)))
cat(sprintf("  Max |r| with mortality rate:  %.3f\n",
            max(res_rate$abs_cor)))

cat("\n  Conclusion: covariates predict crossing volume (30/38 sig.),\n")
cat("  correlate with deaths only through volume (wrong sign),\n")
cat("  and have near-zero signal for mortality rate (1/38 sig.).\n")

cat("\n=======================================================\n")
cat("DONE\n")
cat("=======================================================\n")
