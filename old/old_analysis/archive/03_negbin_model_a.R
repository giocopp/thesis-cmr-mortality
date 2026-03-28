# 03_negbin_model_a.R
# ===================
# Event-level NegBin Model A: primary analysis.
#
# dead_missing_i ~ NegBin(mu_i)
# log(mu_i) = grid_FE + time_FE + beta(Weather_i x Post_i) + X'gamma
#
# Outcome: dead + missing per incident (>= 1 for all CMR incidents).
#
# Uses fixest::fenegbin for high-dimensional FE (proper singleton handling,
# no convergence issues from MASS::glm.nb with many grid cells).
#
# Weather windows:
#   day0:  incident-day value only
#   3-day: mean over lag 0-2 (transit/distress period)
#   7-day: mean over lag 0-7 (full transit window)
#
# Specifications per window:
#   A: SWH x Post interaction
#   B: Wind x Post interaction
#   C: Joint (both interactions)
#   D: Gust x Post interaction
#
# Sensitivity: grid FE resolution, outlier exclusion, FE alternatives,
#              clustered SEs, Poisson comparison
#
# Input:  data/processed/cmr_events_with_weather.RDS
# Output: printed results + saved model objects

library(fixest)
library(data.table)

BASE_DIR <- here::here()
DATA_PATH <- file.path(BASE_DIR, "data", "processed",
                        "cmr_events_with_weather.RDS")

df <- as.data.table(readRDS(DATA_PATH))
cat("Loaded:", nrow(df), "incidents\n\n")

# ============================================================
# 1. Sample description
# ============================================================
cat("============================================================\n")
cat("1. SAMPLE DESCRIPTION\n")
cat("============================================================\n\n")

cat(sprintf("Total incidents: %d\n", nrow(df)))
cat(sprintf("  Pre-MoU:  %d\n", sum(df$post_mou == 0)))
cat(sprintf("  Post-MoU: %d\n", sum(df$post_mou == 1)))
cat(sprintf("\nOutcome (dead + missing):\n"))
cat(sprintf("  Mean: %.2f, Median: %d, Max: %d, Min: %d\n",
    mean(df$dead_missing), median(df$dead_missing),
    max(df$dead_missing), min(df$dead_missing)))
cat(sprintf("  Zeros: %d (by construction, all incidents have dead+missing >= 1)\n",
    sum(df$dead_missing == 0)))
cat(sprintf("  Var/Mean ratio: %.1f (overdispersion if > 1)\n",
    var(df$dead_missing) / mean(df$dead_missing)))

# Estimation sample: rows with non-missing SWH and wind
est_mask <- !is.na(df$swh_day0) & !is.na(df$wind_day0)
cat(sprintf("\nEstimation sample (non-missing SWH + wind): %d / %d\n",
    sum(est_mask), nrow(df)))
cat(sprintf("  Dropped due to NA weather: %d\n", sum(!est_mask)))
cat(sprintf("  Pre-MoU:  %d\n", sum(est_mask & df$post_mou == 0)))
cat(sprintf("  Post-MoU: %d\n", sum(est_mask & df$post_mou == 1)))

# ============================================================
# 2. Grid variables
# ============================================================
cat("\n============================================================\n")
cat("2. GRID VARIABLES\n")
cat("============================================================\n\n")

# Create grid at multiple resolutions
df[, grid_1deg := paste0(sprintf("%.0f", round(grid_lat)), "_",
                          sprintf("%.0f", round(grid_lon)))]

cat(sprintf("0.25-degree grid: %d cells, %d singletons\n",
    length(unique(df$grid_id)),
    sum(table(df$grid_id) == 1)))
cat(sprintf("1-degree grid:    %d cells, %d singletons\n",
    length(unique(df$grid_1deg)),
    sum(table(df$grid_1deg) == 1)))
cat("Using 1-degree grid for primary FE (fixest drops singletons automatically).\n")

# ============================================================
# 3. Primary specifications across weather windows
# ============================================================
cat("\n============================================================\n")
cat("3. PRIMARY SPECIFICATIONS (day0, 3-day, 7-day)\n")
cat("============================================================\n\n")

# Helper: print interaction coefficients from fixest model
print_fixest_int <- function(label, model, vcov_type = "hetero") {
  ct <- summary(model, vcov = vcov_type)$coeftable
  int_rows <- grep(":", rownames(ct))
  for (r in int_rows) {
    b <- ct[r, 1]; se <- ct[r, 2]; p <- ct[r, 4]
    cat(sprintf("  %-40s  b=%+.4f  SE=%.4f  exp(b)=%.4f  p=%.4f %s\n",
        paste0(label, ": ", rownames(ct)[r]),
        b, se, exp(b), p,
        ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
               ifelse(p < 0.1, "*", "")))))
  }
}

# Helper: extract interaction row as data.table for summary
extract_int <- function(label, model, vcov_type = "hetero") {
  ct <- summary(model, vcov = vcov_type)$coeftable
  int_rows <- grep(":", rownames(ct))
  if (length(int_rows) == 0) return(NULL)
  rbindlist(lapply(int_rows, function(r) {
    data.table(
      spec = label, coef_name = rownames(ct)[r],
      beta = ct[r, 1], se = ct[r, 2], p = ct[r, 4],
      exp_beta = exp(ct[r, 1]), N = nobs(model), theta = model$theta
    )
  }))
}

# Define three weather windows
windows <- list(
  day0 = list(swh = "swh_day0",     wind = "wind_day0",     gust = "i10fg_day0"),
  d3   = list(swh = "swh_mean_3d",  wind = "wind_mean_3d",  gust = "i10fg_mean_3d"),
  d7   = list(swh = "swh_mean_7d",  wind = "wind_mean_7d",  gust = "i10fg_mean_7d")
)

all_models <- list()   # store all model objects
all_results <- list()  # store extracted coefficients

for (wname in names(windows)) {
  w <- windows[[wname]]
  wlabel <- switch(wname, day0 = "Day-0", d3 = "3-day mean", d7 = "7-day mean")

  cat(sprintf("========== Weather window: %s ==========\n\n", wlabel))

  # Build formulas with actual variable names
  # Spec A: SWH x Post (+ wind control)
  fml_a <- as.formula(sprintf(
    "dead_missing ~ %s * post_mou + %s | grid_1deg + year_fac + month_fac",
    w$swh, w$wind))

  # Spec B: Wind x Post (+ SWH control)
  fml_b <- as.formula(sprintf(
    "dead_missing ~ %s * post_mou + %s | grid_1deg + year_fac + month_fac",
    w$wind, w$swh))

  # Spec C: Joint (SWH + Wind both interacted)
  fml_c <- as.formula(sprintf(
    "dead_missing ~ %s * post_mou + %s * post_mou | grid_1deg + year_fac + month_fac",
    w$swh, w$wind))

  # Spec D: Gust x Post (+ SWH + wind controls)
  fml_d <- as.formula(sprintf(
    "dead_missing ~ %s * post_mou + %s + %s | grid_1deg + year_fac + month_fac",
    w$gust, w$swh, w$wind))

  cat(sprintf("--- Spec A [%s]: SWH x Post ---\n", wlabel))
  m_a <- fenegbin(fml_a, data = df, vcov = "hetero")
  cat(sprintf("  N: %d, theta: %.3f\n", nobs(m_a), m_a$theta))
  print_fixest_int(sprintf("A [%s]", wlabel), m_a)
  all_models[[paste0("a_", wname)]] <- m_a
  all_results[[length(all_results) + 1]] <- extract_int(
    sprintf("A: SWH×Post [%s]", wlabel), m_a)

  cat(sprintf("\n--- Spec B [%s]: Wind x Post ---\n", wlabel))
  m_b <- fenegbin(fml_b, data = df, vcov = "hetero")
  cat(sprintf("  N: %d, theta: %.3f\n", nobs(m_b), m_b$theta))
  print_fixest_int(sprintf("B [%s]", wlabel), m_b)
  all_models[[paste0("b_", wname)]] <- m_b
  all_results[[length(all_results) + 1]] <- extract_int(
    sprintf("B: Wind×Post [%s]", wlabel), m_b)

  cat(sprintf("\n--- Spec C [%s]: Joint ---\n", wlabel))
  m_c <- fenegbin(fml_c, data = df, vcov = "hetero")
  cat(sprintf("  N: %d, theta: %.3f\n", nobs(m_c), m_c$theta))
  print_fixest_int(sprintf("C [%s]", wlabel), m_c)
  all_models[[paste0("c_", wname)]] <- m_c
  all_results[[length(all_results) + 1]] <- extract_int(
    sprintf("C: Joint [%s]", wlabel), m_c)

  cat(sprintf("\n--- Spec D [%s]: Gust x Post ---\n", wlabel))
  m_d <- fenegbin(fml_d, data = df, vcov = "hetero")
  cat(sprintf("  N: %d, theta: %.3f\n", nobs(m_d), m_d$theta))
  print_fixest_int(sprintf("D [%s]", wlabel), m_d)
  all_models[[paste0("d_", wname)]] <- m_d
  all_results[[length(all_results) + 1]] <- extract_int(
    sprintf("D: Gust×Post [%s]", wlabel), m_d)

  cat("\n")
}

# ============================================================
# 4. Sensitivity checks (using 3-day mean as primary)
# ============================================================
cat("\n============================================================\n")
cat("4. SENSITIVITY CHECKS (3-day mean window)\n")
cat("============================================================\n\n")

# --- Grid FE sensitivity ---
cat("--- No grid FE (year + month only) ---\n")
m_a_nogrid <- fenegbin(dead_missing ~ swh_mean_3d * post_mou + wind_mean_3d |
                         year_fac + month_fac,
                       data = df, vcov = "hetero")
print_fixest_int("A no-grid [3d]", m_a_nogrid)
all_results[[length(all_results) + 1]] <- extract_int("A: no grid FE [3d]", m_a_nogrid)

cat("\n--- 0.25-degree grid FE ---\n")
m_a_fine <- fenegbin(dead_missing ~ swh_mean_3d * post_mou + wind_mean_3d |
                       grid_id + year_fac + month_fac,
                     data = df, vcov = "hetero")
cat(sprintf("  N: %d (after singleton drops)\n", nobs(m_a_fine)))
print_fixest_int("A 0.25deg [3d]", m_a_fine)
all_results[[length(all_results) + 1]] <- extract_int("A: 0.25° grid [3d]", m_a_fine)

# --- Temporal FE sensitivity ---
cat("\n--- Month-of-year only (no year FE) ---\n")
m_a_monthonly <- fenegbin(dead_missing ~ swh_mean_3d * post_mou + wind_mean_3d |
                            grid_1deg + month_fac,
                          data = df, vcov = "hetero")
print_fixest_int("A month-only [3d]", m_a_monthonly)
all_results[[length(all_results) + 1]] <- extract_int("A: month-only [3d]", m_a_monthonly)

# --- Outlier sensitivity ---
cat("\n--- Outlier sensitivity ---\n")
cat("Top 5 deadliest incidents (dead+missing):\n")
print(df[order(-dead_missing),
         .(date, dead_missing, dead, swh_mean_3d, wind_mean_3d, post_mou)][1:5])

df_trim <- df[dead_missing <= 100]
cat(sprintf("\nExcluding dead+missing > 100: dropped %d, kept %d\n",
    nrow(df) - nrow(df_trim), nrow(df_trim)))

m_a_trim <- fenegbin(dead_missing ~ swh_mean_3d * post_mou + wind_mean_3d |
                       grid_1deg + year_fac + month_fac,
                     data = df_trim, vcov = "hetero")
print_fixest_int("A trim [3d]", m_a_trim)
all_results[[length(all_results) + 1]] <- extract_int("A: trim ≤100 [3d]", m_a_trim)

m_d_trim <- fenegbin(dead_missing ~ i10fg_mean_3d * post_mou + swh_mean_3d + wind_mean_3d |
                       grid_1deg + year_fac + month_fac,
                     data = df_trim, vcov = "hetero")
print_fixest_int("D trim [3d]", m_d_trim)
all_results[[length(all_results) + 1]] <- extract_int("D: trim ≤100 [3d]", m_d_trim)

# --- Clustered SEs ---
cat("\n--- Clustered SEs (by grid cell) ---\n")
print_fixest_int("A cluster [3d]", all_models$a_d3, vcov_type = ~grid_1deg)
all_results[[length(all_results) + 1]] <- extract_int(
  "A: cluster SE [3d]", all_models$a_d3, vcov_type = ~grid_1deg)
print_fixest_int("D cluster [3d]", all_models$d_d3, vcov_type = ~grid_1deg)
all_results[[length(all_results) + 1]] <- extract_int(
  "D: cluster SE [3d]", all_models$d_d3, vcov_type = ~grid_1deg)

# --- Poisson comparison ---
cat("\n--- Poisson comparison ---\n")
m_a_pois <- feglm(dead_missing ~ swh_mean_3d * post_mou + wind_mean_3d |
                     grid_1deg + year_fac + month_fac,
                   data = df, family = poisson, vcov = "hetero")
m_d_pois <- feglm(dead_missing ~ i10fg_mean_3d * post_mou + swh_mean_3d + wind_mean_3d |
                     grid_1deg + year_fac + month_fac,
                   data = df, family = poisson, vcov = "hetero")
cat("NOTE: Poisson with Var/Mean >> 1; shown for comparison only.\n")
print_fixest_int("A Poisson [3d]", m_a_pois)
print_fixest_int("D Poisson [3d]", m_d_pois)
all_results[[length(all_results) + 1]] <- extract_int("A: Poisson [3d]", m_a_pois)
all_results[[length(all_results) + 1]] <- extract_int("D: Poisson [3d]", m_d_pois)

# ============================================================
# 5. Summary table
# ============================================================
cat("\n============================================================\n")
cat("5. SUMMARY TABLE\n")
cat("============================================================\n\n")

results_dt <- rbindlist(all_results, fill = TRUE)
results_dt[, sig := fifelse(p < 0.01, "***",
                    fifelse(p < 0.05, "**",
                    fifelse(p < 0.1, "*", "")))]

cat(sprintf("%-35s | %+8s %7s %7s | %8s %3s | %5s\n",
    "Specification", "Beta", "SE", "p", "exp(B)", "", "N"))
cat(paste(rep("-", 90), collapse = ""), "\n")

for (i in seq_len(nrow(results_dt))) {
  r <- results_dt[i]
  cat(sprintf("%-35s | %+8.4f %7.4f %7.4f | %8.4f %3s | %5d\n",
      r$spec, r$beta, r$se, r$p, r$exp_beta, r$sig, r$N))
}

# ============================================================
# 6. Save models and results table
# ============================================================
cat("\n============================================================\n")
cat("6. SAVING MODELS\n")
cat("============================================================\n\n")

model_dir <- file.path(BASE_DIR, "output", "models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

# Save all model objects
all_models[["a_nogrid"]]    <- m_a_nogrid
all_models[["a_fine"]]      <- m_a_fine
all_models[["a_monthonly"]]  <- m_a_monthonly
all_models[["a_trim"]]      <- m_a_trim
all_models[["d_trim"]]      <- m_d_trim
all_models[["a_pois"]]      <- m_a_pois
all_models[["d_pois"]]      <- m_d_pois

saveRDS(all_models, file.path(model_dir, "model_a_results.RDS"))
cat("Models saved to:", file.path(model_dir, "model_a_results.RDS"), "\n")

# Save results table
fwrite(results_dt, file.path(BASE_DIR, "output", "tables", "model_a_results.csv"))
cat("Results table saved to: output/tables/model_a_results.csv\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
