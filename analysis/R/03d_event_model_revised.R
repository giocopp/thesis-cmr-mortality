# 03d_event_model_revised.R
# ========================
# Event-level NegBin revised: informed by daily panel findings.
#
# Changes from 03_negbin_model_a.R:
#   1. Focus on day-0 weather (strongest timing per daily panel)
#   2. Compare incident-location vs corridor vs core geography
#   3. Test FE structures: grid+yr+mo, grid+quarter, grid+yr×mo, grid+mo
#   4. Gradient stability diagnostic (pre-period SWH × year)
#
# Input:  data/processed/cmr_events_with_weather.RDS
#         data/processed/cmr_daily_weather_panel.RDS (for corridor/core weather)
# Output: output/tables/event_model_revised.csv

library(fixest)
library(data.table)

BASE_DIR <- here::here()

# ============================================================
# 0. Load and merge
# ============================================================
df <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                       "cmr_events_with_weather.RDS")))
df[, grid_1deg := paste0(sprintf("%.0f", round(grid_lat)), "_",
                          sprintf("%.0f", round(grid_lon)))]

daily <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                          "cmr_daily_weather_panel.RDS")))

# Merge corridor/core spatial-mean weather onto events by date
merge_vars <- c("date", "swh_mean", "swh_core", "wind_mean", "wind_core")
df <- merge(df, daily[, ..merge_vars], by = "date", all.x = TRUE,
            suffixes = c("", "_panel"))

# Rename to avoid confusion with event-level variables
setnames(df,
         c("swh_mean", "swh_core", "wind_mean", "wind_core"),
         c("swh_corridor", "swh_core_geo", "wind_corridor", "wind_core_geo"))

# Create temporal FE variables
df[, quarter := paste0(year, "Q", ceiling(month(date) / 3))]
df[, quarter_fac := factor(quarter)]
df[, year_month := paste0(year, "_", sprintf("%02d", month(date)))]
df[, year_month_fac := factor(year_month)]
df[, half_year := paste0(year, ifelse(month(date) <= 6, "H1", "H2"))]
df[, half_year_fac := factor(half_year)]

# ============================================================
# Helpers
# ============================================================
report <- function(m, label) {
  ct <- summary(m, vcov = "hetero")$coeftable
  cat(sprintf("  %s:\n", label))
  for (i in seq_len(nrow(ct))) {
    stars <- ifelse(ct[i, 4] < 0.01, "***",
             ifelse(ct[i, 4] < 0.05, "**",
             ifelse(ct[i, 4] < 0.1, "*", "")))
    cat(sprintf("    %-35s %+8.4f (SE=%6.4f) p=%6.4f %s  IRR=%.4f\n",
        rownames(ct)[i], ct[i, 1], ct[i, 2], ct[i, 4], stars, exp(ct[i, 1])))
  }
  cat(sprintf("    N=%d, LogLik=%.1f\n\n", nobs(m), logLik(m)))
}

extract_int <- function(spec, model, vcov_type = "hetero") {
  ct <- summary(model, vcov = vcov_type)$coeftable
  int_rows <- grep(":", rownames(ct))
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
cat("============================================================\n")
cat("EVENT-LEVEL MODEL: Revised specs (informed by daily panel)\n")
cat("============================================================\n\n")

cat(sprintf("Events: %d (pre=%d, post=%d)\n",
    nrow(df), sum(df$post_mou == 0), sum(df$post_mou == 1)))
cat(sprintf("Corridor weather merged: %d matched\n", sum(!is.na(df$swh_corridor))))
cat(sprintf("Grid cells (1deg): %d\n", uniqueN(df$grid_1deg)))
cat(sprintf("Outcome var/mean: %.1f\n\n", var(df$dead_missing) / mean(df$dead_missing)))


# ============================================================
# 1. GEOGRAPHY: incident-location vs corridor vs core
# ============================================================
cat("============================================================\n")
cat("1. GEOGRAPHY COMPARISON (day-0 SWH × Post | grid+yr+mo FE)\n")
cat("============================================================\n\n")

# a) Incident-location SWH (0.5° nearest-neighbor)
m_geo_loc <- fenegbin(dead_missing ~ swh_day0 * post_mou + wind_day0 |
                        grid_1deg + year_fac + month_fac,
                      data = df, vcov = "hetero")
report(m_geo_loc, "Incident-location SWH")

# b) Corridor mean [10,18]x[31,38]
m_geo_corr <- fenegbin(dead_missing ~ swh_corridor * post_mou + wind_corridor |
                         grid_1deg + year_fac + month_fac,
                       data = df, vcov = "hetero")
report(m_geo_corr, "Corridor mean SWH [10-18,31-38]")

# c) Core mean [11,15]x[32,36]
m_geo_core <- fenegbin(dead_missing ~ swh_core_geo * post_mou + wind_core_geo |
                         grid_1deg + year_fac + month_fac,
                       data = df, vcov = "hetero")
report(m_geo_core, "Core mean SWH [11-15,32-36]")

# d) Core mean without grid FE (corridor weather is same for all locations)
m_geo_core_ng <- fenegbin(dead_missing ~ swh_core_geo * post_mou + wind_core_geo |
                            year_fac + month_fac,
                          data = df, vcov = "hetero")
report(m_geo_core_ng, "Core mean SWH | yr+mo (no grid FE)")


# ============================================================
# 2. FE STRUCTURE (incident-location day-0 SWH)
# ============================================================
cat("============================================================\n")
cat("2. FE STRUCTURE (incident-location SWH day-0 × Post)\n")
cat("============================================================\n\n")

# a) Current: grid + year + month (additive)
cat("Already estimated above (m_geo_loc)\n\n")

# b) Grid + quarter
m_fe_q <- fenegbin(dead_missing ~ swh_day0 * post_mou + wind_day0 |
                     grid_1deg + quarter_fac,
                   data = df, vcov = "hetero")
report(m_fe_q, "grid + quarter FE")

# c) Grid + year×month (interacted)
m_fe_ym <- tryCatch(
  fenegbin(dead_missing ~ swh_day0 * post_mou + wind_day0 |
             grid_1deg + year_month_fac,
           data = df, vcov = "hetero"),
  error = function(e) { cat("  grid+yr×mo FAILED:", e$message, "\n\n"); NULL }
)
if (!is.null(m_fe_ym)) {
  report(m_fe_ym, "grid + year×month FE")
}

# d) Grid + half-year
m_fe_hy <- fenegbin(dead_missing ~ swh_day0 * post_mou + wind_day0 |
                      grid_1deg + half_year_fac,
                    data = df, vcov = "hetero")
report(m_fe_hy, "grid + half-year FE")

# e) Grid + month only (no year)
m_fe_mo <- fenegbin(dead_missing ~ swh_day0 * post_mou + wind_day0 |
                      grid_1deg + month_fac,
                    data = df, vcov = "hetero")
report(m_fe_mo, "grid + month-only FE")

# f) No FE (baseline)
m_fe_none <- fenegbin(dead_missing ~ swh_day0 * post_mou + wind_day0,
                      data = df[!is.na(swh_day0)], vcov = "hetero")
report(m_fe_none, "No FE")


# ============================================================
# 3. BEST COMBO: core geography + FE variants
# ============================================================
cat("============================================================\n")
cat("3. CORE GEOGRAPHY + FE VARIANTS\n")
cat("============================================================\n\n")

# Core + grid + quarter
m_cq <- fenegbin(dead_missing ~ swh_core_geo * post_mou + wind_core_geo |
                   grid_1deg + quarter_fac,
                 data = df, vcov = "hetero")
report(m_cq, "Core SWH | grid + quarter")

# Core + grid + year×month
m_cym <- tryCatch(
  fenegbin(dead_missing ~ swh_core_geo * post_mou + wind_core_geo |
             grid_1deg + year_month_fac,
           data = df, vcov = "hetero"),
  error = function(e) { cat("  FAILED:", e$message, "\n\n"); NULL }
)
if (!is.null(m_cym)) {
  report(m_cym, "Core SWH | grid + year×month")
}

# Core + quarter only (no grid)
m_cq_ng <- fenegbin(dead_missing ~ swh_core_geo * post_mou + wind_core_geo |
                      quarter_fac,
                    data = df, vcov = "hetero")
report(m_cq_ng, "Core SWH | quarter only")


# ============================================================
# 4. GRADIENT STABILITY: Pre-period SWH × year
# ============================================================
cat("============================================================\n")
cat("4. GRADIENT STABILITY (pre-MoU only)\n")
cat("============================================================\n\n")

d_pre <- df[post_mou == 0 & !is.na(swh_day0)]
cat(sprintf("Pre-period: %d incidents (%s)\n\n",
    nrow(d_pre), paste(sort(unique(d_pre$year)), collapse = ", ")))

# a) Incident-location SWH × year
m_stab_loc <- tryCatch(
  fenegbin(dead_missing ~ swh_day0 * factor(year) + wind_day0 | grid_1deg + month_fac,
           data = d_pre, vcov = "hetero"),
  error = function(e) { cat("  Failed:", e$message, "\n"); NULL }
)

if (!is.null(m_stab_loc)) {
  ct <- summary(m_stab_loc, vcov = "hetero")$coeftable
  swh_rows <- grep("swh", rownames(ct))
  cat("Incident-location SWH gradient by pre-MoU year:\n")
  for (r in swh_rows) {
    stars <- ifelse(ct[r, 4] < 0.01, "***",
             ifelse(ct[r, 4] < 0.05, "**",
             ifelse(ct[r, 4] < 0.1, "*", "")))
    cat(sprintf("  %-40s %+8.4f (SE=%6.4f) p=%6.4f %s\n",
        rownames(ct)[r], ct[r, 1], ct[r, 2], ct[r, 4], stars))
  }

  # Wald test
  int_rows_idx <- grep("swh_day0:factor", rownames(ct))
  if (length(int_rows_idx) >= 2) {
    w_test <- tryCatch(
      wald(m_stab_loc, keep = "swh_day0:factor", vcov = "hetero"),
      error = function(e) NULL
    )
    if (!is.null(w_test)) {
      cat("\n  Wald test (H0: SWH gradient constant across pre-MoU years):\n")
      print(w_test)
    }
  }
}

# b) Year-by-year slopes (simple Poisson, no grid FE)
cat("\nYear-by-year SWH slopes (separate Poisson):\n")
for (yr in sort(unique(d_pre$year))) {
  dsub <- d_pre[year == yr]
  if (nrow(dsub) < 10) {
    cat(sprintf("  %d: n=%d (too few)\n", yr, nrow(dsub)))
    next
  }
  m_yr <- tryCatch(
    glm(dead_missing ~ swh_day0 + wind_day0, data = dsub, family = poisson),
    error = function(e) NULL)
  if (!is.null(m_yr)) {
    b <- coef(m_yr)["swh_day0"]
    se <- summary(m_yr)$coefficients["swh_day0", 2]
    cat(sprintf("  %d: n=%3d, beta=%+.4f (SE=%.4f), IRR=%.4f\n",
        yr, nrow(dsub), b, se, exp(b)))
  }
}

# c) Core-mean SWH gradient stability
cat("\n")
d_pre_core <- df[post_mou == 0 & !is.na(swh_core_geo)]
m_stab_core <- tryCatch(
  fenegbin(dead_missing ~ swh_core_geo * factor(year) + wind_core_geo | grid_1deg + month_fac,
           data = d_pre_core, vcov = "hetero"),
  error = function(e) { cat("  Core stability failed:", e$message, "\n"); NULL }
)

if (!is.null(m_stab_core)) {
  ct <- summary(m_stab_core, vcov = "hetero")$coeftable
  swh_rows <- grep("swh", rownames(ct))
  cat("Core-mean SWH gradient by pre-MoU year:\n")
  for (r in swh_rows) {
    stars <- ifelse(ct[r, 4] < 0.01, "***",
             ifelse(ct[r, 4] < 0.05, "**",
             ifelse(ct[r, 4] < 0.1, "*", "")))
    cat(sprintf("  %-40s %+8.4f (SE=%6.4f) p=%6.4f %s\n",
        rownames(ct)[r], ct[r, 1], ct[r, 2], ct[r, 4], stars))
  }

  int_rows_idx <- grep("swh_core_geo:factor", rownames(ct))
  if (length(int_rows_idx) >= 2) {
    w_test <- tryCatch(
      wald(m_stab_core, keep = "swh_core_geo:factor", vcov = "hetero"),
      error = function(e) NULL
    )
    if (!is.null(w_test)) {
      cat("\n  Wald test (H0: core SWH gradient constant pre-MoU):\n")
      print(w_test)
    }
  }
}


# ============================================================
# 5. SUMMARY TABLE
# ============================================================
cat("\n============================================================\n")
cat("5. SUMMARY TABLE\n")
cat("============================================================\n\n")

results <- rbindlist(list(
  # Geography comparison
  extract_int("Incident-loc | grid+yr+mo", m_geo_loc),
  extract_int("Corridor | grid+yr+mo", m_geo_corr),
  extract_int("Core | grid+yr+mo", m_geo_core),
  extract_int("Core | yr+mo (no grid)", m_geo_core_ng),
  # FE structure (incident-location)
  extract_int("Inc-loc | grid+quarter", m_fe_q),
  if (!is.null(m_fe_ym)) extract_int("Inc-loc | grid+yr×mo", m_fe_ym),
  extract_int("Inc-loc | grid+half-yr", m_fe_hy),
  extract_int("Inc-loc | grid+mo only", m_fe_mo),
  extract_int("Inc-loc | no FE", m_fe_none),
  # Core + FE variants
  extract_int("Core | grid+quarter", m_cq),
  if (!is.null(m_cym)) extract_int("Core | grid+yr×mo", m_cym),
  extract_int("Core | quarter (no grid)", m_cq_ng)
), fill = TRUE)

# Filter to SWH interaction only
results_swh <- results[grepl("swh", coef, ignore.case = TRUE)]

cat(sprintf("%-35s %+8s %7s %7s %8s %5s\n",
    "Specification", "Beta", "SE", "p", "IRR", "N"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (i in seq_len(nrow(results_swh))) {
  r <- results_swh[i]
  stars <- ifelse(r$p < 0.01, "***",
           ifelse(r$p < 0.05, "**",
           ifelse(r$p < 0.1, "*", "")))
  cat(sprintf("%-35s %+8.4f %7.4f %7.4f %8.4f %5d %s\n",
      r$spec, r$beta, r$se, r$p, r$irr, r$n, stars))
}

fwrite(results, file.path(BASE_DIR, "output", "tables", "event_model_revised.csv"))
cat("\nSaved: output/tables/event_model_revised.csv\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
