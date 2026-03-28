# 03d_event_model_revised.R
# ========================
# Complementary event-level analysis.
#
# The daily panel (03c) is the primary model. This script provides:
#   A. Event-level NegBin: does weather predict deaths per incident?
#   B. Fatality rate model: does weather predict the fraction who die?
#
# Key differences from daily panel:
#   - Unit = individual incident (not day)
#   - Weather = SWH at incident's nearest grid cell (not spatial mean)
#   - FE = grid + year + month (not week-year — too few per cell)
#   - Lag-1 is primary (same rationale: IOM date = reporting date)
#
# Uses the clean dataset from 01b_core_corridor_dataset.R:
#   - Core corridor [10.5, 15.5] x [32.3, 36.2]
#   - All IOM MMP fields preserved
#   - Weather at lags 0, 1, 2, 3, 7
#
# Input:  data/processed/core_corridor_incidents.RDS
# Output: output/tables/event_model_revised.csv

library(fixest)
library(data.table)

BASE_DIR <- here::here()

# ============================================================
# 0. Load clean dataset
# ============================================================
cat("============================================================\n")
cat("EVENT-LEVEL MODEL (complementary to daily panel)\n")
cat("============================================================\n\n")

df <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                       "core_corridor_incidents.RDS")))

# Filter to drowning + suspected drowning (same as daily panel)
sea_causes <- c("Drowning", "Mixed or unknown")
df_sea <- df[cause_category %in% sea_causes]

# FE variables
df_sea[, grid_fac := factor(grid_1deg)]
df_sea[, year_fac := factor(year)]
df_sea[, month_fac := factor(month)]

cat(sprintf("Core corridor incidents: %d\n", nrow(df)))
cat(sprintf("Drowning/suspected:      %d (pre=%d, post=%d)\n",
    nrow(df_sea), sum(df_sea$post_mou == 0), sum(df_sea$post_mou == 1)))
cat(sprintf("Grid cells (1deg):       %d\n", uniqueN(df_sea$grid_1deg)))
cat(sprintf("SWH lag-1 available:     %d\n", sum(!is.na(df_sea$swh_lag1))))
cat(sprintf("Outcome var/mean:        %.1f\n",
    var(df_sea$dead_missing) / mean(df_sea$dead_missing)))
cat(sprintf("With survivors known:    %d (%.1f%%)\n",
    sum(!is.na(df_sea$fatality_rate)),
    100 * mean(!is.na(df_sea$fatality_rate))))

# ============================================================
# Helpers
# ============================================================
report <- function(m, label, vcov_type = "hetero") {
  ct <- summary(m, vcov = vcov_type)$coeftable
  theta_str <- ""
  if (!is.null(m$theta)) theta_str <- sprintf(", theta=%.3f", m$theta)
  cat(sprintf("  %s (N=%d%s):\n", label, nobs(m), theta_str))
  for (i in seq_len(nrow(ct))) {
    stars <- ifelse(ct[i, 4] < 0.01, "***",
             ifelse(ct[i, 4] < 0.05, "**",
             ifelse(ct[i, 4] < 0.1, "*", "")))
    cat(sprintf("    %-35s %+8.4f (SE=%6.4f) p=%6.4f %s\n",
        rownames(ct)[i], ct[i, 1], ct[i, 2], ct[i, 4], stars))
  }
  cat("\n")
}

extract_int <- function(spec, model, vcov_type = "hetero") {
  ct <- summary(model, vcov = vcov_type)$coeftable
  int_rows <- grep(":post_mou|post_mou:", rownames(ct))
  if (length(int_rows) == 0) return(NULL)
  rbindlist(lapply(int_rows, function(r) {
    data.table(
      spec = spec, coef = rownames(ct)[r],
      beta = ct[r, 1], se = ct[r, 2], p = ct[r, 4],
      irr = exp(ct[r, 1]), n = nobs(model)
    )
  }))
}


# ============================================================
# A. EVENT-LEVEL NEGBIN: dead_missing ~ SWH x Post
# ============================================================
cat("\n============================================================\n")
cat("A. EVENT-LEVEL NEGBIN: dead_missing per incident\n")
cat("============================================================\n\n")

# A1. Primary: SWH lag-1 | grid + year + month FE
cat("--- A1. SWH lag-1 | grid + yr + mo FE ---\n")
m_lag1 <- fenegbin(dead_missing ~ swh_lag1 + swh_lag1:post_mou |
                     grid_fac + year_fac + month_fac,
                   data = df_sea[!is.na(swh_lag1)], vcov = "hetero")
report(m_lag1, "SWH lag-1 | grid + yr + mo")

# A2. SWH lag-1 | year + month FE (no grid)
cat("--- A2. SWH lag-1 | yr + mo FE (no grid) ---\n")
m_lag1_ng <- fenegbin(dead_missing ~ swh_lag1 + swh_lag1:post_mou |
                        year_fac + month_fac,
                      data = df_sea[!is.na(swh_lag1)], vcov = "hetero")
report(m_lag1_ng, "SWH lag-1 | yr + mo (no grid)")

# A3. SWH lag-2 for timing comparison
cat("--- A3. SWH lag-2 | grid + yr + mo FE ---\n")
m_lag2 <- fenegbin(dead_missing ~ swh_lag2 + swh_lag2:post_mou |
                     grid_fac + year_fac + month_fac,
                   data = df_sea[!is.na(swh_lag2)], vcov = "hetero")
report(m_lag2, "SWH lag-2 | grid + yr + mo")

# A4. SWH day-0 for comparison (old spec)
cat("--- A4. SWH day-0 | grid + yr + mo FE ---\n")
m_day0 <- fenegbin(dead_missing ~ swh_day0 + swh_day0:post_mou |
                     grid_fac + year_fac + month_fac,
                   data = df_sea[!is.na(swh_day0)], vcov = "hetero")
report(m_day0, "SWH day-0 | grid + yr + mo")


# ============================================================
# B. FATALITY RATE MODEL
# ============================================================
cat("\n============================================================\n")
cat("B. FATALITY RATE MODEL: fatality_rate ~ SWH x Post\n")
cat("============================================================\n\n")

# Subsample with known survivors
df_fr <- df_sea[!is.na(fatality_rate) & !is.na(swh_lag1)]

cat(sprintf("Fatality rate sample: %d incidents (pre=%d, post=%d)\n",
    nrow(df_fr), sum(df_fr$post_mou == 0), sum(df_fr$post_mou == 1)))
cat(sprintf("  Mean fatality rate: %.3f (pre=%.3f, post=%.3f)\n",
    mean(df_fr$fatality_rate),
    mean(df_fr$fatality_rate[df_fr$post_mou == 0]),
    mean(df_fr$fatality_rate[df_fr$post_mou == 1])))
cat(sprintf("  Median: %.3f\n", median(df_fr$fatality_rate)))
cat(sprintf("  Grid cells: %d\n\n", uniqueN(df_fr$grid_1deg)))

# Fatality rate is bounded [0,1]. Options:
#   - OLS (linear probability-style, simple, interpretable)
#   - Fractional logit (Papke & Wooldridge 1996)
#   - Beta regression
# We use OLS with FE as primary (transparent), fractional logit as robustness.

# B1. OLS: fatality_rate ~ SWH lag-1 x Post | grid + yr + mo
cat("--- B1. OLS: fatality_rate ~ SWH lag-1 x Post | grid + yr + mo ---\n")
m_fr_ols <- feols(fatality_rate ~ swh_lag1 + swh_lag1:post_mou |
                    grid_fac + year_fac + month_fac,
                  data = df_fr, vcov = "hetero")
report(m_fr_ols, "OLS fatality rate | grid + yr + mo")

# B2. OLS without grid FE
cat("--- B2. OLS: fatality_rate ~ SWH lag-1 x Post | yr + mo ---\n")
m_fr_ols_ng <- feols(fatality_rate ~ swh_lag1 + swh_lag1:post_mou |
                       year_fac + month_fac,
                     data = df_fr, vcov = "hetero")
report(m_fr_ols_ng, "OLS fatality rate | yr + mo (no grid)")

# B3. Fractional logit (Papke-Wooldridge): consistent for bounded outcomes
cat("--- B3. Fractional logit: fatality_rate ~ SWH lag-1 x Post ---\n")

# fepois with fractional response is the Papke-Wooldridge estimator
# when the outcome is in [0,1]. feglm with family=quasibinomial is equivalent.
# Use feglm for fractional logit.
m_fr_frac <- tryCatch(
  feglm(fatality_rate ~ swh_lag1 + swh_lag1:post_mou |
          year_fac + month_fac,
        data = df_fr, vcov = "hetero", family = quasibinomial(link = "logit")),
  error = function(e) { cat("  Fractional logit failed:", e$message, "\n"); NULL }
)
if (!is.null(m_fr_frac)) {
  report(m_fr_frac, "Fractional logit | yr + mo")
}

# B4. OLS with lag-2 for timing comparison
cat("--- B4. OLS: fatality_rate ~ SWH lag-2 x Post | grid + yr + mo ---\n")
df_fr2 <- df_sea[!is.na(fatality_rate) & !is.na(swh_lag2)]
m_fr_lag2 <- feols(fatality_rate ~ swh_lag2 + swh_lag2:post_mou |
                     grid_fac + year_fac + month_fac,
                   data = df_fr2, vcov = "hetero")
report(m_fr_lag2, "OLS fatality rate lag-2 | grid + yr + mo")


# ============================================================
# C. Summary table
# ============================================================
cat("\n============================================================\n")
cat("C. SUMMARY TABLE\n")
cat("============================================================\n\n")

results <- rbindlist(list(
  # Event-level NegBin
  extract_int("Dead/missing lag-1 | grid+yr+mo", m_lag1),
  extract_int("Dead/missing lag-1 | yr+mo", m_lag1_ng),
  extract_int("Dead/missing lag-2 | grid+yr+mo", m_lag2),
  extract_int("Dead/missing day-0 | grid+yr+mo", m_day0),
  # Fatality rate
  extract_int("Fatality rate OLS lag-1 | grid+yr+mo", m_fr_ols),
  extract_int("Fatality rate OLS lag-1 | yr+mo", m_fr_ols_ng),
  if (!is.null(m_fr_frac)) extract_int("Fatality rate frac.logit | yr+mo", m_fr_frac),
  extract_int("Fatality rate OLS lag-2 | grid+yr+mo", m_fr_lag2)
), fill = TRUE)

# Filter to SWH interaction only
results_swh <- results[grepl("swh|SWH", coef, ignore.case = TRUE)]

cat(sprintf("%-40s %+8s %7s %7s %8s %5s\n",
    "Specification", "Beta", "SE", "p", "IRR", "N"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in seq_len(nrow(results_swh))) {
  r <- results_swh[i]
  stars <- ifelse(r$p < 0.01, "***",
           ifelse(r$p < 0.05, "**",
           ifelse(r$p < 0.1, "*", "")))
  cat(sprintf("%-40s %+8.4f %7.4f %7.4f %8.4f %5d %s\n",
      r$spec, r$beta, r$se, r$p, r$irr, r$n, stars))
}

fwrite(results, file.path(BASE_DIR, "output", "tables", "event_model_revised.csv"))
cat("\nSaved: output/tables/event_model_revised.csv\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
