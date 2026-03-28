# 04a_dml_robustness.R
# ====================
# Double Machine Learning robustness check for the SWH x Post interaction.
#
# Purpose: Test whether the NegBin log-link functional form drives the
#          primary result (core SWH x Post, week-year FE, IRR=2.55).
#          DML (Chernozhukov et al. 2018) uses ML to flexibly estimate
#          nuisance functions, then estimates the interaction coefficient
#          from orthogonalized residuals with cross-fitting.
#
# Approach: Partially Linear Model (Robinson 1988)
#   Y = theta * D + g(X) + epsilon
#   where D = swh_core * post_mou (the interaction)
#         X = {swh_core, post_mou, wind_core, temporal features}
#         g(X) is estimated by random forest (nonparametric)
#
#   If theta_DML ~ beta_3_NegBin in sign and rough magnitude,
#   the functional form is not driving the result.
#
# Input:  data/processed/cmr_daily_weather_panel.RDS
# Output: output/tables/dml_results.csv
#         printed comparison with NegBin baseline

library(DoubleML)
library(mlr3)
library(mlr3learners)
library(ranger)
library(data.table)
library(fixest)

BASE_DIR <- here::here()

# ============================================================
# 1. Load data and prepare variables
# ============================================================
cat("============================================================\n")
cat("DML ROBUSTNESS: Partially Linear Model\n")
cat("============================================================\n\n")

d <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                      "cmr_daily_weather_panel.RDS")))

# Create the interaction variable (this is "D" in the DML)
d[, swh_core_x_post := swh_core * post_mou]

cat(sprintf("Panel: %d days, outcome mean=%.3f, zeros=%.1f%%\n",
    nrow(d), mean(d$n_dead_missing), 100 * mean(d$n_dead_missing == 0)))

# ============================================================
# 2. NegBin baseline (for comparison)
# ============================================================
cat("\n--- NegBin baseline (from 03c) ---\n")
m_nb <- fenegbin(n_dead_missing ~ swh_core * post_mou + wind_core | week_year_fac,
                 data = d, vcov = "hetero")
ct_nb <- summary(m_nb, vcov = "hetero")$coeftable
int_row <- grep("swh_core:post_mou|post_mou:swh_core", rownames(ct_nb))
beta_nb <- ct_nb[int_row, 1]
se_nb   <- ct_nb[int_row, 2]
p_nb    <- ct_nb[int_row, 4]
cat(sprintf("  NegBin: beta=%+.4f, SE=%.4f, p=%.4f, IRR=%.4f\n\n",
    beta_nb, se_nb, p_nb, exp(beta_nb)))


# ============================================================
# 3. DML: Partially Linear Model with random forest
# ============================================================
cat("============================================================\n")
cat("3. DML PARTIALLY LINEAR MODEL\n")
cat("============================================================\n\n")

# Feature set for the ML learners.
# We include temporal features (week, year, month, dow) instead of
# 628 week-year dummies — the forest can learn these nonparametrically.
# Also include main effects of weather variables.
d[, week_num := as.integer(week)]
d[, year_num := as.integer(year)]
d[, month_num := as.integer(month)]
d[, dow_num := as.integer(dow)]

x_cols <- c("swh_core", "post_mou", "wind_core",
            "week_num", "year_num", "month_num", "dow_num")

# Complete cases only
d_dml <- d[complete.cases(d[, c("n_dead_missing", "swh_core_x_post", x_cols),
                             with = FALSE])]
cat(sprintf("DML sample: %d obs (dropped %d with NAs)\n",
    nrow(d_dml), nrow(d) - nrow(d_dml)))

# DoubleML data object
dml_data <- DoubleMLData$new(
  data = d_dml[, c("n_dead_missing", "swh_core_x_post", x_cols), with = FALSE],
  y_col = "n_dead_missing",
  d_cols = "swh_core_x_post",
  x_cols = x_cols
)

# --- 3a. Random forest learners ---
cat("\n--- DML with random forest (500 trees) ---\n")
set.seed(42)

ml_l_rf <- lrn("regr.ranger", num.trees = 500, min.node.size = 5,
                max.depth = NULL)
ml_m_rf <- lrn("regr.ranger", num.trees = 500, min.node.size = 5,
                max.depth = NULL)

dml_rf <- DoubleMLPLR$new(dml_data, ml_l = ml_l_rf, ml_m = ml_m_rf,
                           n_folds = 5, n_rep = 1)
dml_rf$fit()

cat("  Random Forest PLR:\n")
print(dml_rf$summary())
theta_rf <- dml_rf$coef
se_rf <- dml_rf$se
p_rf <- dml_rf$pval
ci_rf <- dml_rf$confint()

cat(sprintf("\n  theta=%+.4f, SE=%.4f, p=%.4f\n", theta_rf, se_rf, p_rf))
cat(sprintf("  95%% CI: [%+.4f, %+.4f]\n", ci_rf[1], ci_rf[2]))

# --- 3b. Multiple repetitions for stability ---
cat("\n--- DML with 5 repetitions (cross-fitting stability) ---\n")
set.seed(123)

ml_l_rf2 <- lrn("regr.ranger", num.trees = 500, min.node.size = 5)
ml_m_rf2 <- lrn("regr.ranger", num.trees = 500, min.node.size = 5)

dml_rf5 <- DoubleMLPLR$new(dml_data, ml_l = ml_l_rf2, ml_m = ml_m_rf2,
                            n_folds = 5, n_rep = 5)
dml_rf5$fit()

cat("  Random Forest PLR (5 reps, median aggregation):\n")
print(dml_rf5$summary())
theta_rf5 <- dml_rf5$coef
se_rf5 <- dml_rf5$se
p_rf5 <- dml_rf5$pval
ci_rf5 <- dml_rf5$confint()
cat(sprintf("\n  theta=%+.4f, SE=%.4f, p=%.4f\n", theta_rf5, se_rf5, p_rf5))
cat(sprintf("  95%% CI: [%+.4f, %+.4f]\n", ci_rf5[1], ci_rf5[2]))


# --- 3c. Boosted trees as alternative ML learner ---
cat("\n--- DML with larger forest (1000 trees, deeper) ---\n")
set.seed(42)

ml_l_rf3 <- lrn("regr.ranger", num.trees = 1000, min.node.size = 3,
                  max.depth = 20)
ml_m_rf3 <- lrn("regr.ranger", num.trees = 1000, min.node.size = 3,
                  max.depth = 20)

dml_deep <- DoubleMLPLR$new(dml_data, ml_l = ml_l_rf3, ml_m = ml_m_rf3,
                             n_folds = 5, n_rep = 1)
dml_deep$fit()

theta_deep <- dml_deep$coef
se_deep <- dml_deep$se
p_deep <- dml_deep$pval
cat(sprintf("  theta=%+.4f, SE=%.4f, p=%.4f\n", theta_deep, se_deep, p_deep))


# ============================================================
# 4. DML with richer feature set
# ============================================================
cat("\n============================================================\n")
cat("4. DML WITH RICHER FEATURES\n")
cat("============================================================\n\n")

# Add more weather variables as features
x_cols_rich <- c("swh_core", "post_mou", "wind_core",
                 "gust_core", "mwp_core",
                 "swh_mean", "wind_mean",
                 "week_num", "year_num", "month_num", "dow_num")

d_dml_rich <- d[complete.cases(d[, c("n_dead_missing", "swh_core_x_post",
                                      x_cols_rich), with = FALSE])]
cat(sprintf("Rich feature sample: %d obs\n", nrow(d_dml_rich)))

dml_data_rich <- DoubleMLData$new(
  data = d_dml_rich[, c("n_dead_missing", "swh_core_x_post", x_cols_rich),
                      with = FALSE],
  y_col = "n_dead_missing",
  d_cols = "swh_core_x_post",
  x_cols = x_cols_rich
)

set.seed(42)
ml_l_rich <- lrn("regr.ranger", num.trees = 500, min.node.size = 5)
ml_m_rich <- lrn("regr.ranger", num.trees = 500, min.node.size = 5)

dml_rich <- DoubleMLPLR$new(dml_data_rich, ml_l = ml_l_rich, ml_m = ml_m_rich,
                             n_folds = 5, n_rep = 5)
dml_rich$fit()

theta_rich <- dml_rich$coef
se_rich <- dml_rich$se
p_rich <- dml_rich$pval
ci_rich <- dml_rich$confint()
cat(sprintf("  Rich features: theta=%+.4f, SE=%.4f, p=%.4f\n",
    theta_rich, se_rich, p_rich))
cat(sprintf("  95%% CI: [%+.4f, %+.4f]\n", ci_rich[1], ci_rich[2]))


# ============================================================
# 5. OLS benchmark (linear model, same interaction)
# ============================================================
cat("\n============================================================\n")
cat("5. OLS BENCHMARK\n")
cat("============================================================\n\n")

# OLS with week-year dummies (linear probability-style for counts)
m_ols <- feols(n_dead_missing ~ swh_core * post_mou + wind_core | week_year_fac,
               data = d, vcov = "hetero")
ct_ols <- summary(m_ols)$coeftable
int_row_ols <- grep("swh_core:post_mou|post_mou:swh_core", rownames(ct_ols))
beta_ols <- ct_ols[int_row_ols, 1]
se_ols <- ct_ols[int_row_ols, 2]
p_ols <- ct_ols[int_row_ols, 4]
cat(sprintf("  OLS (feols, week-year FE): beta=%+.4f, SE=%.4f, p=%.4f\n",
    beta_ols, se_ols, p_ols))

# Poisson QMLE
m_poi <- fepois(n_dead_missing ~ swh_core * post_mou + wind_core | week_year_fac,
                data = d, vcov = "hetero")
ct_poi <- summary(m_poi, vcov = "hetero")$coeftable
int_row_poi <- grep("swh_core:post_mou|post_mou:swh_core", rownames(ct_poi))
beta_poi <- ct_poi[int_row_poi, 1]
se_poi <- ct_poi[int_row_poi, 2]
p_poi <- ct_poi[int_row_poi, 4]
cat(sprintf("  Poisson QMLE:              beta=%+.4f, SE=%.4f, p=%.4f\n",
    beta_poi, se_poi, p_poi))


# ============================================================
# 6. Comparison table
# ============================================================
cat("\n============================================================\n")
cat("6. COMPARISON: NegBin vs DML vs OLS vs Poisson\n")
cat("============================================================\n\n")

results <- data.table(
  model = c("NegBin (primary)", "DML RF (5-fold)", "DML RF (5x5 reps)",
            "DML RF (deep)", "DML RF (rich features)",
            "OLS (week-yr FE)", "Poisson QMLE"),
  beta = c(beta_nb, theta_rf, theta_rf5, theta_deep, theta_rich,
           beta_ols, beta_poi),
  se = c(se_nb, se_rf, se_rf5, se_deep, se_rich, se_ols, se_poi),
  p = c(p_nb, p_rf, p_rf5, p_deep, p_rich, p_ols, p_poi)
)

cat(sprintf("%-28s %10s %8s %8s\n", "Model", "Beta/Theta", "SE", "p"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i]
  stars <- ifelse(r$p < 0.01, "***", ifelse(r$p < 0.05, "**",
                  ifelse(r$p < 0.1, "*", "")))
  cat(sprintf("%-28s %+10.4f %8.4f %8.4f %s\n",
      r$model, r$beta, r$se, r$p, stars))
}

cat("\nNote: NegBin/Poisson betas are on log scale (IRR = exp(beta)).\n")
cat("DML/OLS thetas are on level scale (additional deaths per unit SWH*Post).\n")
cat("Sign agreement is the key comparison, not magnitude.\n")

# Save
fwrite(results, file.path(BASE_DIR, "output", "tables", "dml_results.csv"))
cat("\nSaved: output/tables/dml_results.csv\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
