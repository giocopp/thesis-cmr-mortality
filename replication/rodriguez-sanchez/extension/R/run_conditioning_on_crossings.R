# run_conditioning_on_crossings.R
# ================================
# Sensitivity analysis: add log(crossings) as a covariate to the
# death count BSTS model. This absorbs the volume channel so the
# model tests whether deaths are higher than expected *at the
# observed crossing level*.
#
# Loads cached targets objects (df_model) to avoid re-building data.
# Runs CausalImpact for all three death-count interventions (A, B, C).

library(dplyr)
library(lubridate)
library(bsts)
library(CausalImpact)
library(targets)

# --- Load cached data from targets pipeline ---
message("Loading cached df_model from targets store...")
PROJECT_DIR <- "/Users/giocopp/Desktop/Uni/Hertie School/6th Semester/Thesis-MDS/Rodriguez-Sanchez-paper-replication/Extension-2-new-data"
withr::with_dir(PROJECT_DIR, {
  df_model <- tar_read(df_model)
})

# --- Build conditioned death-count dataset ---
# df_deaths has: date, deaths_cmr (= log(deaths+1)), covariates, month dummies
# df_full has: all original columns including crossings_CMR
# We need to add log(crossings_CMR) as a covariate to df_deaths

df_full <- df_model$df_full

# Create a lookup: date -> log(crossings)
crossings_lookup <- df_full %>%
  filter(!is.na(crossings_CMR), crossings_CMR > 0) %>%
  transmute(
    date = date,
    log_crossings = log(crossings_CMR)
  )

# Join to df_deaths
df_deaths_orig <- df_model$df_deaths
df_deaths_cond <- df_deaths_orig %>%
  left_join(crossings_lookup, by = "date") %>%
  filter(!is.na(log_crossings))

# Verify: log_crossings should be placed after the response (deaths_cmr)
# CausalImpact uses: col 1 = date (index), col 2 = response, cols 3+ = covariates
# Move log_crossings to be the third column (first covariate)
resp_col <- "deaths_cmr"
date_col <- "date"
other_cols <- setdiff(names(df_deaths_cond),
                      c(date_col, resp_col, "log_crossings"))
df_deaths_cond <- df_deaths_cond %>%
  select(all_of(c(date_col, resp_col, "log_crossings", other_cols)))

message("Conditioned dataset: ", nrow(df_deaths_cond), " rows x ",
        ncol(df_deaths_cond), " cols")
message("Original deaths dataset: ", nrow(df_deaths_orig), " rows x ",
        ncol(df_deaths_orig), " cols")
message("New covariate: log_crossings (range: ",
        round(min(df_deaths_cond$log_crossings), 2), " to ",
        round(max(df_deaths_cond$log_crossings), 2), ")")

# Verify correlation between log_crossings and deaths_cmr
r <- cor(df_deaths_cond$log_crossings, df_deaths_cond$deaths_cmr,
         use = "complete.obs")
message("Correlation log_crossings <-> deaths_cmr: ", round(r, 3))

# --- Model settings (same as main pipeline) ---
NITER <- 10000L
SEED  <- 270488

model_args <- list(
  dynamic.regression = FALSE,
  standardize.data   = TRUE,
  max.flips          = -1L,
  niter              = NITER
)

# --- Define interventions for death counts ---
interventions <- df_model$interventions %>%
  filter(outcome == "death_count")

# --- Fit models ---
results_list <- list()

for (i in seq_len(nrow(interventions))) {
  spec <- interventions[i, ]
  label <- spec$label

  pre_period  <- c(min(df_deaths_cond$date), spec$pre_end)
  post_period <- c(spec$post_start, max(df_deaths_cond$date))

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Fitting CONDITIONED model: ", label)
  message("  Treatment: ", spec$treatment)
  message("  Pre-period:  ", pre_period[1], " to ", pre_period[2])
  message("  Post-period: ", post_period[1], " to ", post_period[2])
  message(paste(rep("=", 60), collapse = ""))

  set.seed(SEED)
  t1 <- Sys.time()

  impact <- CausalImpact(
    df_deaths_cond,
    pre.period  = pre_period,
    post.period = post_period,
    alpha       = 0.05,
    model.args  = model_args
  )

  t2 <- Sys.time()
  elapsed <- as.numeric(difftime(t2, t1, units = "mins"))
  message("  Done in ", round(elapsed, 1), " minutes.")

  # --- Extract summary ---
  s <- impact$summary
  ser <- impact$series

  pre_idx <- which(ser$cum.effect == 0)
  pre_n   <- if (length(pre_idx) > 0) max(pre_idx) else NA_integer_
  post_n  <- nrow(ser) - pre_n

  # Pre-period fit
  resid    <- ser$response[seq_len(pre_n)] - ser$point.pred[seq_len(pre_n)]
  rmse_pre <- sqrt(mean(resid^2, na.rm = TRUE))
  sd_pre   <- sd(ser$response[seq_len(pre_n)], na.rm = TRUE)
  rmse_sd  <- rmse_pre / sd_pre

  # Inclusion probabilities
  inc_probs <- colMeans(impact$model$bsts.model$coefficients != 0)

  # Log crossings inclusion probability
  log_cross_inc <- if ("log_crossings" %in% names(inc_probs)) {
    inc_probs["log_crossings"]
  } else {
    NA_real_
  }

  # Top 10 included covariates
  top_inc <- sort(inc_probs[inc_probs > 0.01], decreasing = TRUE)
  top_10  <- head(top_inc, 10)

  result <- tibble::tibble(
    label           = label,
    treatment       = spec$treatment,
    p_value         = as.numeric(s["Cumulative", "p"]),
    rel_effect      = as.numeric(s["Cumulative", "RelEffect"]),
    abs_effect      = as.numeric(s["Cumulative", "AbsEffect"]),
    abs_lower       = as.numeric(s["Cumulative", "AbsEffect.lower"]),
    abs_upper       = as.numeric(s["Cumulative", "AbsEffect.upper"]),
    n_pre           = pre_n,
    n_post          = post_n,
    rmse_sd         = rmse_sd,
    n_inc_01        = sum(inc_probs > 0.01),
    n_inc_10        = sum(inc_probs > 0.10),
    log_cross_inc   = as.numeric(log_cross_inc),
    elapsed_min     = elapsed
  )

  results_list[[label]] <- list(
    summary = result,
    impact  = impact,
    top_inc = top_10
  )

  # Print summary
  message("\n  --- RESULTS ---")
  message("  p-value:     ", round(result$p_value, 4))
  message("  Rel. effect: ", round(result$rel_effect * 100, 1), "%")
  message("  Cum. effect: ", round(result$abs_effect, 1),
          " [", round(result$abs_lower, 1), ", ",
          round(result$abs_upper, 1), "]")
  message("  RMSE/SD:     ", round(result$rmse_sd, 3))
  message("  log_crossings inclusion: ",
          round(as.numeric(log_cross_inc) * 100, 1), "%")
  message("  Predictors > 1% incl.: ", result$n_inc_01)
  message("  Predictors > 10% incl.: ", result$n_inc_10)
  message("\n  Top included covariates:")
  for (nm in names(top_10)) {
    message("    ", sprintf("%-40s", nm), " ",
            round(top_10[nm] * 100, 1), "%")
  }
}

# --- Comparison table ---
message("\n\n", paste(rep("=", 70), collapse = ""))
message("COMPARISON: ORIGINAL vs CONDITIONED (death count models)")
message(paste(rep("=", 70), collapse = ""))

# Load original death count results from targets
withr::with_dir(PROJECT_DIR, {
  orig_a <- tar_read(model_a_deaths)
  orig_b <- tar_read(model_b_deaths)
  orig_c <- tar_read(model_c_deaths)
})

extract_orig <- function(m) {
  s <- m$impact$summary
  ser <- m$impact$series
  pre_idx <- which(ser$cum.effect == 0)
  pre_n   <- max(pre_idx)
  resid   <- ser$response[seq_len(pre_n)] - ser$point.pred[seq_len(pre_n)]
  rmse_pre <- sqrt(mean(resid^2, na.rm = TRUE))
  sd_pre   <- sd(ser$response[seq_len(pre_n)], na.rm = TRUE)
  tibble::tibble(
    p_value    = as.numeric(s["Cumulative", "p"]),
    rel_effect = as.numeric(s["Cumulative", "RelEffect"]),
    abs_effect = as.numeric(s["Cumulative", "AbsEffect"]),
    rmse_sd    = rmse_pre / sd_pre
  )
}

orig_results <- list(
  A_deaths = extract_orig(orig_a),
  B_deaths = extract_orig(orig_b),
  C_deaths = extract_orig(orig_c)
)

message("\n", sprintf("%-12s  %-8s %-10s %-10s %-8s  |  %-8s %-10s %-10s %-8s",
                      "Model", "p_orig", "rel_orig", "cum_orig", "RMSE/SD",
                      "p_cond", "rel_cond", "cum_cond", "RMSE/SD"))
message(paste(rep("-", 100), collapse = ""))

for (label in c("A_deaths", "B_deaths", "C_deaths")) {
  o <- orig_results[[label]]
  c_res <- results_list[[label]]$summary

  message(sprintf("%-12s  %-8s %-10s %-10s %-8s  |  %-8s %-10s %-10s %-8s",
                  label,
                  round(o$p_value, 4),
                  paste0(round(o$rel_effect * 100, 1), "%"),
                  round(o$abs_effect, 1),
                  round(o$rmse_sd, 3),
                  round(c_res$p_value, 4),
                  paste0(round(c_res$rel_effect * 100, 1), "%"),
                  round(c_res$abs_effect, 1),
                  round(c_res$rmse_sd, 3)))
}

# --- Save results ---
output_dir <- file.path(PROJECT_DIR, "output", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

comparison_df <- bind_rows(
  lapply(names(results_list), function(label) {
    results_list[[label]]$summary
  })
)

readr::write_csv(comparison_df,
                 file.path(output_dir, "conditioning_on_crossings.csv"))
message("\nResults saved to output/tables/conditioning_on_crossings.csv")

# Save the full R objects for later inspection
save(results_list,
     file = file.path(PROJECT_DIR, "output", "conditioning_results.RData"))
message("Full results saved to output/conditioning_results.RData")
