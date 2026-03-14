# run_full_predictor_model.R
# ==========================
# "Full predictor model": Original paper's full predictor set + Extension-2
# improvements, applied to the MORTALITY outcome (not crossings).
#
# Combines:
#   - Original paper's ~1,000+ predictors (lags 01-06): economic indicators,
#     conflict data, airport flows, currency exchange, Google Trends, etc.
#   - Extension-2 curated environmental covariates: ERA5 sea conditions, moon
#     illumination, ocean currents, extreme wave statistics, SST anomaly
#   - 11 month-of-year dummies (replacing numeric month/quarter/semester)
#   - SST anomaly with pre-period climatology correction
#
# Improvements over the curated specification (Extension-2, 65 predictors):
#   1. Broader predictor set (~2,000 variables) giving the spike-and-slab
#      more candidates to build the counterfactual
#   2. Proper seasonality: month dummies instead of numeric month/quarter/semester
#   3. SST anomaly uses pre-period climatology only (no data leakage)
#   4. ERA5 sea variables added (not in original df.RDS)
#   5. max.flips = -1 (unlimited spike-and-slab proposals)
#
# This tests whether adding the original paper's rich predictor set
# improves the mortality counterfactual beyond Extension-2's curated specification.
#
# PREDICTOR COMPOSITION (2,006 predictors, 128 months)
# ────────────────────────────────────────────────────
# Source                             Columns  Base vars  Lags
# ─────────────────────────────────  ───────  ─────────  ────────
# Disaster counts (EM-DAT)              735       105    0-06
# Commodity prices                      602        86    0-06
# Currency exchange rates               266        38    0-06
# Google Trends (job search)            105        15    0-06
# Airport passenger flows                97        97    0 only
# ERA5 Central Med (sea/atmosphere)      91        13    0-06
# Weather Italy/Malta (original)         42         6    0-06
# ERA5 departure coast                   28         4    0-06
# Ocean currents                         14         2    0-06
# Month-of-year dummies                  11        11    ---
# Moon illumination                       7         1    0-06
# Syria Google Trends + year              8         2    0-06
# ─────────────────────────────────  ───────
# TOTAL                               2,006
#
# DROPPED DUE TO INCOMPLETE DATA (NAs in analysis window):
#   - ACLED conflict indicators:     7,175 cols (9 NAs / 128 months)
#   - Unemployment rates:              494 cols (3 NAs / 128 months)
#   - UCDP Syrian conflict deaths:       0 cols (not in df_enhanced)
#   Reason: the original df.RDS was compiled when ACLED data ended at
#   Dec 2020 (missing Jan-Sep 2021) and unemployment at Jun 2021 (missing
#   Jul-Sep 2021). Columns with any NAs are dropped to preserve the full
#   128-month window (Feb 2011 - Sep 2021), matching Extension-2 exactly.
#   Alternative: shorten the window to Dec 2020 for both this model and
#   Extension-2, which would retain ACLED and unemployment but lose 9
#   post-MoU months (47 -> 38 post-period observations).
#
# Usage:
#   cd Extension-2-new-data/
#   Rscript run_full_predictor_model.R
#
# Requires: targets pipeline has been run at least once (for df_enhanced cache).

library(dplyr)
library(lubridate)
library(CausalImpact)
library(bsts)
library(targets)

# --- Constants ---
SEED     <- 270488
NITER    <- 10000L
EPSILON  <- 0.01
MAX_FLIP <- -1L       # unlimited; reduce to 100 if too slow

cat("============================================================\n")
cat("  FULL PREDICTOR MODEL\n")
cat("============================================================\n\n")


# ---------------------------------------------------------------
# 0. Load cached df_enhanced from targets pipeline
# ---------------------------------------------------------------
cat("[", format(Sys.time(), "%H:%M:%S"), "] Loading df_enhanced from targets cache...\n")
df <- targets::tar_read(df_enhanced, store = "_targets")
cat("    Raw dimensions: ", nrow(df), " x ", ncol(df), "\n")


# ---------------------------------------------------------------
# 1. Apply SST anomaly correction (pre-period climatology only)
# ---------------------------------------------------------------
cat("[", format(Sys.time(), "%H:%M:%S"), "] Correcting SST anomaly...\n")
source("R/02_data_prepare.R")
df <- compute_sst_anomaly_preperiod(df)


# ---------------------------------------------------------------
# 2. Construct outcome variables
# ---------------------------------------------------------------
df <- df %>%
  mutate(
    LCG_pushbacks_count = as.numeric(
      ifelse(is.na(LCG_pushbacks_count), 0, LCG_pushbacks_count)),
    TCG_pushbacks_count = as.numeric(
      ifelse(is.na(TCG_pushbacks_count), 0, TCG_pushbacks_count)),
    dead_and_missing_Central_Mediterranean = as.numeric(
      ifelse(is.na(dead_and_missing_Central_Mediterranean), 0,
             dead_and_missing_Central_Mediterranean)),
    crossings_CMR = arrivals_CMR + LCG_pushbacks_count +
      TCG_pushbacks_count + dead_and_missing_Central_Mediterranean,
    mortality_rate_100 = (dead_and_missing_Central_Mediterranean /
                            crossings_CMR) * 100
  )


# ---------------------------------------------------------------
# 3. Filter time window
# ---------------------------------------------------------------
df <- df %>%
  filter(date >= "2011-02-01" & date < "2021-10-01")
cat("    After time filter: ", nrow(df), " x ", ncol(df), "\n")


# ---------------------------------------------------------------
# 4. Build predictor set
# ---------------------------------------------------------------
cat("[", format(Sys.time(), "%H:%M:%S"), "] Building predictor set...\n")

# 4a. Identify columns to REMOVE

# High-order lags (07-24) — keep only 01-06
high_lag_cols <- character(0)
for (lag_i in 7:24) {
  lag_str <- formatC(lag_i, width = 2, flag = "0")
  high_lag_cols <- c(high_lag_cols,
                     grep(paste0("_lag_", lag_str, "$"), names(df), value = TRUE))
}

# Outcome-related variables (would create endogeneity)
outcome_cols <- c(
  "arrivals_BSR", "arrivals_CMR", "arrivals_CRAG", "arrivals_EBR",
  "arrivals_EMR", "arrivals_OR", "arrivals_WAR", "arrivals_WBR", "arrivals_WMR",
  "dead_and_missing_Eastern_Mediterranean",
  "dead_and_missing_Western_Mediterranean",
  "dead_and_missing_Central_Mediterranean",
  "LCG_pushbacks_count", "TCG_pushbacks_count",
  "crossings_CMR", "mortality_rate_100"
)
# Also remove any lags of these that might exist
outcome_lag_cols <- character(0)
for (oc in c("arrivals_", "dead_and_missing_", "LCG_pushbacks", "TCG_pushbacks",
             "crossings_CMR", "mortality_rate")) {
  outcome_lag_cols <- c(outcome_lag_cols,
                        grep(paste0("^", oc), names(df), value = TRUE))
}

# Geographic dispersion and fraction indices (outcome-derived)
geo_cols <- c(grep("^sd_lat_", names(df), value = TRUE),
              grep("^sd_lon_", names(df), value = TRUE),
              grep("^frac_index_", names(df), value = TRUE))

# Problematic variables (same as original paper)
problem_cols <- c(grep("^airflow_Palestinian", names(df), value = TRUE),
                  grep("^asylum", names(df), value = TRUE))

# Numeric temporal indicators (REPLACED by month dummies)
temporal_cols <- intersect(c("month", "semester", "quarter"), names(df))

# Metadata
meta_cols <- c("date", "yearmonth")

# Intermediate variables from Extension-2 pipeline
intermediate_cols <- c(
  "wind_u_central_med", "wind_v_central_med",
  "wind_u_departure_coast", "wind_v_departure_coast",
  "dewpoint_central_med",
  "current_eastward_central_med", "current_northward_central_med"
)
# Also their lags
intermediate_lag_cols <- character(0)
for (ic in c("wind_u_central_med", "wind_v_central_med",
             "wind_u_departure_coast", "wind_v_departure_coast",
             "dewpoint_central_med",
             "current_eastward_central_med", "current_northward_central_med")) {
  intermediate_lag_cols <- c(intermediate_lag_cols,
                             grep(paste0("^", ic, "_lag_"), names(df), value = TRUE))
}

# Combine all exclusions
all_exclude <- unique(c(
  high_lag_cols, outcome_cols, outcome_lag_cols, geo_cols,
  problem_cols, temporal_cols, meta_cols,
  intermediate_cols, intermediate_lag_cols
))

# 4b. Get predictor columns
predictor_cols <- setdiff(names(df), all_exclude)

# Keep only numeric columns
predictor_cols <- predictor_cols[sapply(df[predictor_cols], is.numeric)]

# Remove any column that is constant (zero variance)
col_sds <- sapply(df[predictor_cols], function(x) sd(x, na.rm = TRUE))
zero_var_cols <- names(col_sds[is.na(col_sds) | col_sds == 0])
if (length(zero_var_cols) > 0) {
  cat("    Removing", length(zero_var_cols), "zero-variance columns\n")
  predictor_cols <- setdiff(predictor_cols, zero_var_cols)
}

cat("    Candidate predictor columns: ", length(predictor_cols), "\n")

# 4c. Drop columns with ANY NAs in the analysis window
#     (preserves all 128 months; matches Extension-2's time window exactly)
na_counts <- colSums(is.na(df[predictor_cols]))
incomplete_cols <- names(na_counts[na_counts > 0])
if (length(incomplete_cols) > 0) {
  cat("    Dropping", length(incomplete_cols),
      "predictors with missing values (preserving all months)\n")
  predictor_cols <- setdiff(predictor_cols, incomplete_cols)
}

cat("    Final predictor columns: ", length(predictor_cols), "\n")


# ---------------------------------------------------------------
# 5. Create model datasets
# ---------------------------------------------------------------
cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating model datasets...\n")

# Extension-2 reference window
EXT2_START <- as.Date("2011-02-01")
EXT2_END   <- as.Date("2021-09-01")   # last month in Extension-2

# Mortality rate dataset (no na.omit — columns are already complete)
df_mort <- df %>%
  transmute(
    date = date,
    mortality_rate = log(mortality_rate_100 + EPSILON),
    across(all_of(predictor_cols))
  )

# Add 11 month dummies (Feb-Dec, Jan = reference)
for (m in 2:12) {
  df_mort[[paste0("month_", m)]] <- as.integer(month(df_mort$date) == m)
}

# Death count dataset
df_deaths <- df %>%
  transmute(
    date = date,
    deaths_cmr = log(dead_and_missing_Central_Mediterranean + 1),
    across(all_of(predictor_cols))
  )

for (m in 2:12) {
  df_deaths[[paste0("month_", m)]] <- as.integer(month(df_deaths$date) == m)
}

n_pred_mort   <- ncol(df_mort) - 2    # exclude date and outcome
n_pred_deaths <- ncol(df_deaths) - 2

# Verify date range matches Extension-2
cat("    Mortality rate: ", nrow(df_mort), " obs x ", n_pred_mort, " predictors\n")
cat("    Death count:    ", nrow(df_deaths), " obs x ", n_pred_deaths, " predictors\n")
cat("    Date range:     ", as.character(min(df_mort$date)), " to ",
    as.character(max(df_mort$date)), "\n")
stopifnot(
  "Date range must start at 2011-02-01 (Extension-2 match)" =
    min(df_mort$date) == EXT2_START,
  "Date range must end at 2021-09-01 (Extension-2 match)" =
    max(df_mort$date) == EXT2_END
)
cat("    Date range matches Extension-2.\n\n")


# ---------------------------------------------------------------
# 6. Model configuration
# ---------------------------------------------------------------
model_args <- list(
  dynamic.regression = FALSE,
  standardize.data   = TRUE,
  max.flips          = MAX_FLIP,
  niter              = NITER
)

# Intervention periods (same as Extension-2)
interventions <- list(
  A = list(
    label     = "A_mortality (Mare Nostrum)",
    pre_end   = as.Date("2013-09-01"),
    post_start = as.Date("2013-10-01")
  ),
  B = list(
    label     = "B_mortality (NGO SAR)",
    pre_end   = as.Date("2014-10-01"),
    post_start = as.Date("2014-11-01")
  ),
  C = list(
    label     = "C (EU-Libya MoU)",
    pre_end   = as.Date("2017-01-01"),
    post_start = as.Date("2017-02-01")
  )
)


# ---------------------------------------------------------------
# 7. Helper: extract and report model results
# ---------------------------------------------------------------
extract_and_report <- function(impact, label) {
  s <- impact$summary
  ser <- impact$series

  # Pre-period fit (guarded against empty pre_idx)
  pre_idx <- which(ser$cum.effect == 0)
  if (length(pre_idx) == 0) {
    warning("Could not detect pre-period from cumulative effects for: ", label)
    pre_n  <- NA_integer_
    post_n <- nrow(ser)
    rmse   <- NA_real_
    sd_pre <- NA_real_
  } else {
    pre_n  <- max(pre_idx)
    post_n <- nrow(ser) - pre_n
    resid  <- ser$response[seq_len(pre_n)] - ser$point.pred[seq_len(pre_n)]
    rmse   <- sqrt(mean(resid^2, na.rm = TRUE))
    sd_pre <- sd(ser$response[seq_len(pre_n)], na.rm = TRUE)
  }

  # Inclusion probabilities
  inc_probs <- colMeans(impact$model$bsts.model$coefficients != 0)
  burn <- tryCatch(
    bsts::SuggestBurn(0.1, impact$model$bsts.model),
    error = function(e) NA_integer_
  )

  cat("\n========================================\n")
  cat(label, "\n")
  cat("========================================\n")
  cat(sprintf("p-value:         %.4f\n", s["Cumulative", "p"]))
  cat(sprintf("Abs effect:      %.2f\n", s["Cumulative", "AbsEffect"]))
  cat(sprintf("CI:              [%.2f, %.2f]\n",
              s["Cumulative", "AbsEffect.lower"],
              s["Cumulative", "AbsEffect.upper"]))
  rmse_sd <- if (!is.na(rmse) && !is.na(sd_pre) && sd_pre > 0) rmse / sd_pre else NA_real_
  cat(sprintf("RMSE/SD (pre):   %s\n",
              if (is.na(rmse_sd)) "NA" else sprintf("%.4f", rmse_sd)))
  cat(sprintf("Burn-in:         %d\n", burn))
  cat(sprintf("Pre months:      %d\n", pre_n))
  cat(sprintf("Post months:     %d\n", post_n))
  cat(sprintf("N predictors:    %d\n", length(inc_probs)))
  cat(sprintf("N inc > 1%%:      %d\n", sum(inc_probs > 0.01)))
  cat(sprintf("N inc > 5%%:      %d\n", sum(inc_probs > 0.05)))
  cat(sprintf("N inc > 10%%:     %d\n", sum(inc_probs > 0.10)))

  cat("\nTop 15 predictors by inclusion probability:\n")
  top <- sort(inc_probs, decreasing = TRUE)[1:min(15, length(inc_probs))]
  for (i in seq_along(top)) {
    cat(sprintf("  %2d. %-50s %.1f%%\n", i, names(top)[i], top[i] * 100))
  }

  invisible(list(
    p_value   = as.numeric(s["Cumulative", "p"]),
    abs_effect = as.numeric(s["Cumulative", "AbsEffect"]),
    abs_lower  = as.numeric(s["Cumulative", "AbsEffect.lower"]),
    abs_upper  = as.numeric(s["Cumulative", "AbsEffect.upper"]),
    rmse_sd    = rmse_sd,
    inc_probs  = inc_probs,
    burn       = burn,
    n_pre      = pre_n,
    n_post     = post_n
  ))
}


# ---------------------------------------------------------------
# 8. Fit all four models
# ---------------------------------------------------------------

# --- Model A: Mare Nostrum (mortality rate) ---
cat("\n=== Fitting Model A: Mare Nostrum (mortality rate) ===\n")
pre_A  <- c(EXT2_START, interventions$A$pre_end)
post_A <- c(interventions$A$post_start, EXT2_END)
set.seed(SEED)
t0 <- Sys.time()
impact_A <- CausalImpact(df_mort, pre.period = pre_A, post.period = post_A,
                          alpha = 0.05, model.args = model_args)
cat("  Elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
res_A <- extract_and_report(impact_A, "Model A (Mare Nostrum) - Full predictor")


# --- Model B: NGO SAR (mortality rate) ---
cat("\n=== Fitting Model B: NGO SAR (mortality rate) ===\n")
pre_B  <- c(EXT2_START, interventions$B$pre_end)
post_B <- c(interventions$B$post_start, EXT2_END)
set.seed(SEED)
t0 <- Sys.time()
impact_B <- CausalImpact(df_mort, pre.period = pre_B, post.period = post_B,
                          alpha = 0.05, model.args = model_args)
cat("  Elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
res_B <- extract_and_report(impact_B, "Model B (NGO SAR) - Full predictor")


# --- Model C: MoU (mortality rate) ---
cat("\n=== Fitting Model C: MoU (mortality rate) ===\n")
pre_C  <- c(EXT2_START, interventions$C$pre_end)
post_C <- c(interventions$C$post_start, EXT2_END)
set.seed(SEED)
t0 <- Sys.time()
impact_C_rate <- CausalImpact(df_mort, pre.period = pre_C, post.period = post_C,
                               alpha = 0.05, model.args = model_args)
cat("  Elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
res_C_rate <- extract_and_report(impact_C_rate, "Model C rate (MoU) - Full predictor")


# --- Model C: MoU (death count) ---
cat("\n=== Fitting Model C: MoU (death count) ===\n")
set.seed(SEED)
t0 <- Sys.time()
impact_C_deaths <- CausalImpact(df_deaths, pre.period = pre_C, post.period = post_C,
                                 alpha = 0.05, model.args = model_args)
cat("  Elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
res_C_deaths <- extract_and_report(impact_C_deaths, "Model C deaths (MoU) - Full predictor")


# ---------------------------------------------------------------
# 9. Placebo tests (Model C only)
# ---------------------------------------------------------------
# Placebo settings — matched to Extension-2 (_targets.R + R/04_diagnostics.R)
run_placebos <- function(df, actual_int_date, model_args,
                         min_pre_months = 12L, min_post_months = 6L,
                         step_months = 3L, seed = SEED) {
  # Restrict to pre-period only (same as Extension-2: date < actual_int_date)
  df_pre    <- df %>% filter(date < actual_int_date)
  n_pre     <- nrow(df_pre)
  pre_dates <- df_pre$date
  data_start <- min(pre_dates)

  if (n_pre < min_pre_months + min_post_months) {
    cat("  Too few pre-period observations for placebos.\n")
    return(data.frame(placebo_date = as.Date(character()),
                      p_value = numeric()))
  }

  # Placebo indices: every step_months, from min_pre+1 to n_pre - min_post
  placebo_start_idx <- min_pre_months + 1
  placebo_end_idx   <- n_pre - min_post_months
  if (placebo_start_idx > placebo_end_idx) {
    cat("  No valid placebo dates.\n")
    return(data.frame(placebo_date = as.Date(character()),
                      p_value = numeric()))
  }

  placebo_indices <- seq(placebo_start_idx, placebo_end_idx, by = step_months)
  placebo_dates   <- pre_dates[placebo_indices]

  cat("  Running ", length(placebo_dates), " placebos (step=", step_months,
      " months, min_pre=", min_pre_months, ")...\n")

  results <- vector("list", length(placebo_dates))

  # Set seed once before the loop (matches Extension-2 behavior)
  set.seed(seed)

  for (i in seq_along(placebo_dates)) {
    pdate <- placebo_dates[i]
    # Extension-2 period construction:
    #   pre:  [data_start, pdate - 1 month]
    #   post: [pdate, actual_int_date - 1 month]
    placebo_pre  <- c(data_start, pdate - months(1))
    placebo_post <- c(pdate, actual_int_date - months(1))

    tryCatch({
      imp <- CausalImpact(df_pre, pre.period = placebo_pre,
                           post.period = placebo_post,
                           alpha = 0.05, model.args = model_args)
      pval <- as.numeric(imp$summary["Cumulative", "p"])
      results[[i]] <- data.frame(placebo_date = pdate, p_value = pval)
      cat(sprintf("    Placebo %2d/%d (date=%s): p=%.4f\n",
                  i, length(placebo_dates), pdate, pval))
    }, error = function(e) {
      cat(sprintf("    Placebo %2d/%d FAILED: %s\n",
                  i, length(placebo_dates), e$message))
      results[[i]] <<- data.frame(placebo_date = pdate, p_value = NA_real_)
    })
  }

  do.call(rbind, results)
}


report_placebos <- function(placebos, actual_p, label) {
  placebos <- placebos[!is.na(placebos$p_value), ]
  n_total  <- nrow(placebos)

  cat("\n----------------------------------------\n")
  cat("Placebo results: ", label, "\n")
  cat("----------------------------------------\n")
  cat(sprintf("N placebos:      %d\n", n_total))

  if (n_total == 0) {
    cat("  No valid placebo results to report.\n")
    return(invisible(NULL))
  }

  n_sig     <- sum(placebos$p_value < 0.05)
  n_extreme <- sum(placebos$p_value <= actual_p)

  cat(sprintf("N significant:   %d\n", n_sig))
  cat(sprintf("FPR:             %.3f\n", n_sig / n_total))
  cat(sprintf("Actual p:        %.4f\n", actual_p))
  cat(sprintf("Placebos <= p:   %d / %d\n", n_extreme, n_total))
  cat(sprintf("Share <= p:      %.3f\n", n_extreme / n_total))
}


# MoU intervention date (same as Extension-2)
MOU_DATE <- as.Date("2017-02-01")

cat("\n\n=== Running placebo tests for Model C (mortality rate) ===\n")
placebos_C_rate <- run_placebos(df_mort, MOU_DATE, model_args)
report_placebos(placebos_C_rate, res_C_rate$p_value, "C rate - Full predictor")

cat("\n\n=== Running placebo tests for Model C (death count) ===\n")
placebos_C_deaths <- run_placebos(df_deaths, MOU_DATE, model_args)
report_placebos(placebos_C_deaths, res_C_deaths$p_value, "C deaths - Full predictor")


# ---------------------------------------------------------------
# 10. Comparison summary
# ---------------------------------------------------------------
cat("\n\n============================================================\n")
cat("  COMPARISON: Full predictor vs. Extension-2 (curated)\n")
cat("============================================================\n\n")

cat("Extension-2 results (from consistency_snapshot):\n")
cat("  C rate:   p=0.041, RMSE/SD=0.965, N_pred=65\n")
cat("  C deaths: p=0.014, RMSE/SD=0.757, N_pred=65\n")
cat("  C rate  placebos: FPR=0.500, Share<=p=0.500 (18 placebos, step=3, min_pre=12)\n")
cat("  C deaths placebos: FPR=0.444, Share<=p=0.167 (18 placebos, step=3, min_pre=12)\n\n")

fmt_rmse <- function(x) if (is.na(x)) "NA" else sprintf("%.4f", x)

cat("Full predictor results:\n")
cat(sprintf("  C rate:   p=%.4f, RMSE/SD=%s, N_pred=%d\n",
            res_C_rate$p_value, fmt_rmse(res_C_rate$rmse_sd),
            length(res_C_rate$inc_probs)))
cat(sprintf("  C deaths: p=%.4f, RMSE/SD=%s, N_pred=%d\n",
            res_C_deaths$p_value, fmt_rmse(res_C_deaths$rmse_sd),
            length(res_C_deaths$inc_probs)))

pl_rate <- placebos_C_rate[!is.na(placebos_C_rate$p_value), ]
if (nrow(pl_rate) > 0) {
  cat(sprintf("  C rate  placebos: FPR=%.3f, Share<=p=%.3f (%d placebos)\n",
              sum(pl_rate$p_value < 0.05) / nrow(pl_rate),
              sum(pl_rate$p_value <= res_C_rate$p_value) / nrow(pl_rate),
              nrow(pl_rate)))
}
pl_d <- placebos_C_deaths[!is.na(placebos_C_deaths$p_value), ]
if (nrow(pl_d) > 0) {
  cat(sprintf("  C deaths placebos: FPR=%.3f, Share<=p=%.3f (%d placebos)\n",
              sum(pl_d$p_value < 0.05) / nrow(pl_d),
              sum(pl_d$p_value <= res_C_deaths$p_value) / nrow(pl_d),
              nrow(pl_d)))
}

cat("\nAll four models:\n")
cat(sprintf("  A (Mare Nostrum): p=%.4f, RMSE/SD=%s\n",
            res_A$p_value, fmt_rmse(res_A$rmse_sd)))
cat(sprintf("  B (NGO SAR):     p=%.4f, RMSE/SD=%s\n",
            res_B$p_value, fmt_rmse(res_B$rmse_sd)))
cat(sprintf("  C (rate):        p=%.4f, RMSE/SD=%s\n",
            res_C_rate$p_value, fmt_rmse(res_C_rate$rmse_sd)))
cat(sprintf("  C (deaths):      p=%.4f, RMSE/SD=%s\n",
            res_C_deaths$p_value, fmt_rmse(res_C_deaths$rmse_sd)))


# ---------------------------------------------------------------
# 11. Save results
# ---------------------------------------------------------------
dir.create("output", showWarnings = FALSE)

save(
  impact_A, impact_B, impact_C_rate, impact_C_deaths,
  res_A, res_B, res_C_rate, res_C_deaths,
  placebos_C_rate, placebos_C_deaths,
  df_mort, df_deaths,
  file = "output/full_predictor_results.RData"
)

cat("\n[", format(Sys.time(), "%H:%M:%S"), "] Done! Results saved to output/full_predictor_results.RData\n")
