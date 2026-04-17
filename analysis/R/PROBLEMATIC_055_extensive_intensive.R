# 055_extensive_intensive.R
# =========================

library(glmmTMB)
library(sandwich)
library(lubridate)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
YEAR_START <- 2014
PERIODS <- c(2020, 2023)

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("055  EXTENSIVE vs INTENSIVE MARGIN\n")
cat("============================================================\n\n")

# ── 1. Load data ────────────────────────────────────────────
cat("--- 1. Loading panels ---\n")

# Drop the panel's broad n_dead_missing and replace with the analytical
# series via the shared helper. Default = incident-only, core corridor,
# all causes. Change the call to test sensitivity variants.
da <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_complete.RDS")) %>%
  select(-n_dead_missing) %>%
  left_join(build_iom_daily(), by = "date") %>%
  replace_na(list(n_dead_missing = 0))

zp <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_zone.RDS"))

# Collapse to 2 blocs
bloc <- zp %>%
  group_by(date, sar_bloc) %>%
  summarise(
    n_dead_missing = sum(n_dead_missing),
    swh_prevweek   = mean(swh_prevweek, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    post_mou   = as.integer(date >= MOU_DATE),
    year       = year(date),
    month_year = factor(format(date, "%Y-%m")),
    iso_week   = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))
  )
dim(bloc$date) <- NULL

cat(sprintf("  daily-agg:  %d rows\n", nrow(da)))
cat(sprintf("  2-bloc:     %d rows\n", nrow(bloc)))
cat(sprintf("  4-country:  %d rows\n", nrow(zp)))

# ── 2. Helpers ──────────────────────────────────────────────

fmt_row <- function(label, coef, se, n) {
  p <- 2 * pnorm(-abs(coef / se))
  sprintf("  %-40s b3 = %+.3f (SE = %.3f)  p = %.4f  (N = %d)",
          label, coef, se, p, n)
}

# glmmTMB reporting helper. Tries cluster-robust SE via sandwich::vcovCL;
# falls back to glmmTMB's model-based (Wald/Hessian) SE if that fails.
report_glmmtmb <- function(m, label, cluster_vec, n) {
  co <- fixef(m)$cond
  target <- "swh_prevweek:post_mou"
  if (!(target %in% names(co))) {
    cat("    target coef not found\n"); return(invisible())
  }
  coef_int <- co[[target]]

  # Try clustered SE
  vc <- tryCatch(
    sandwich::vcovCL(m, cluster = cluster_vec, type = "HC0"),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (!is.null(vc) && target %in% rownames(vc)) {
    se_int <- sqrt(vc[target, target])
    cat(fmt_row(paste(label, "(clust SE)"), coef_int, se_int, n), "\n")
    return(invisible())
  }

  # Fallback: model-based (Wald) SE from glmmTMB summary
  s <- tryCatch(summary(m), error = function(e) NULL)
  if (!is.null(s)) {
    se_int <- s$coefficients$cond[target, "Std. Error"]
    cat(fmt_row(paste(label, "(Wald SE)"), coef_int, se_int, n), "\n")
    return(invisible())
  }

  cat("    could not extract SE\n")
}

# ── 3. Run models ───────────────────────────────────────────
cat("\n--- 3. Estimation ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables", "055_extensive_intensive.txt")
sink(sink_file)

cat("055  EXTENSIVE vs INTENSIVE MARGIN\n")
cat("===================================\n")
cat("Target: swh_prevweek:post_mou coefficient\n\n")
cat("Extensive: logit via feglm on full panel, P(deaths > 0). FE: month_year.\n")
cat("           SE: NW(14).\n")
cat("Intensive: glmmTMB truncated_nbinom2 on days with deaths > 0.\n")
cat("           FE: factor(year) only (month_year caused non-pos-def Hessian\n")
cat("           in the small subset). SE: clustered by iso_week when possible,\n")
cat("           model-based Wald SE as fallback.\n\n")

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("\n=== %s ===\n", label))

  # ---- daily-agg ----
  d_da <- da %>%
    filter(year(date) >= YEAR_START, year(date) <= ye,
           !is.na(swh_prevweek)) %>%
    mutate(
      unit = 1L,
      any_death = as.integer(n_dead_missing > 0),
      year_fac = factor(year(date))
    )

  cat(sprintf("  [A] daily-agg N = %d | days with deaths = %d\n",
      nrow(d_da), sum(d_da$any_death)))

  # (1) Extensive: logit
  m_ext_da <- feglm(
    any_death ~ swh_prevweek + swh_prevweek:post_mou | month_year,
    data = d_da, family = binomial("logit"),
    vcov = NW(14), panel.id = ~unit + date
  )
  ct_ext_da <- coeftable(m_ext_da, vcov = NW(14))
  row_ext <- grep(":post_mou", rownames(ct_ext_da))

  cat(fmt_row("[A] extensive (logit, NW14)",
              ct_ext_da[row_ext, 1], ct_ext_da[row_ext, 2],
              nobs(m_ext_da)), "\n")

  # (2) Intensive: truncated NegBin via glmmTMB. Year FE only (month_year
  # overparameterizes the small subset and causes Hessian problems).
  d_da_int <- d_da %>% filter(any_death == 1)
  m_int_da <- tryCatch(
    glmmTMB(n_dead_missing ~ swh_prevweek + swh_prevweek:post_mou +
              year_fac,
            data = d_da_int,
            family = truncated_nbinom2),
    error = function(e) { cat("    glmmTMB err:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(m_int_da)) {
    report_glmmtmb(m_int_da, "[A] intensive (trunc NegBin)",
                    d_da_int$iso_week, nrow(d_da_int))
  }

  # ---- 2-bloc ----
  d_bl <- bloc %>%
    filter(year >= YEAR_START, year <= ye, !is.na(swh_prevweek)) %>%
    mutate(
      any_death = as.integer(n_dead_missing > 0),
      year_fac = factor(year)
    )

  cat(sprintf("\n  [B] 2-bloc N = %d | rows with deaths = %d\n",
      nrow(d_bl), sum(d_bl$any_death)))

  # (1) Extensive on 2-bloc
  m_ext_bl <- feglm(
    any_death ~ swh_prevweek + swh_prevweek:post_mou |
      month_year + sar_bloc,
    data = d_bl, family = binomial("logit"),
    vcov = NW(14), panel.id = ~sar_bloc + date
  )
  ct_ext_bl <- coeftable(m_ext_bl, vcov = NW(14))
  row_ext_bl <- grep(":post_mou", rownames(ct_ext_bl))

  cat(fmt_row("[B] extensive (logit, NW14)",
              ct_ext_bl[row_ext_bl, 1], ct_ext_bl[row_ext_bl, 2],
              nobs(m_ext_bl)), "\n")

  # (2) Intensive on 2-bloc subset (deaths > 0)
  d_bl_int <- d_bl %>% filter(any_death == 1)
  m_int_bl <- tryCatch(
    glmmTMB(n_dead_missing ~ swh_prevweek + swh_prevweek:post_mou +
              year_fac + sar_bloc,
            data = d_bl_int,
            family = truncated_nbinom2),
    error = function(e) { cat("    glmmTMB err:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(m_int_bl)) {
    report_glmmtmb(m_int_bl, "[B] intensive (trunc NegBin)",
                    d_bl_int$iso_week, nrow(d_bl_int))
  }

  # ---- 4-country ----
  d_zp <- zp %>%
    filter(year >= YEAR_START, year <= ye, !is.na(swh_prevweek)) %>%
    mutate(
      any_death = as.integer(n_dead_missing > 0),
      year_fac = factor(year)
    )

  cat(sprintf("\n  [C] 4-country N = %d | rows with deaths = %d\n",
      nrow(d_zp), sum(d_zp$any_death)))

  # (1) Extensive on 4-country
  m_ext_zp <- feglm(
    any_death ~ swh_prevweek + swh_prevweek:post_mou |
      month_year + country,
    data = d_zp, family = binomial("logit"),
    vcov = NW(14), panel.id = ~country + date
  )
  ct_ext_zp <- coeftable(m_ext_zp, vcov = NW(14))
  row_ext_zp <- grep(":post_mou", rownames(ct_ext_zp))

  cat(fmt_row("[C] extensive (logit, NW14)",
              ct_ext_zp[row_ext_zp, 1], ct_ext_zp[row_ext_zp, 2],
              nobs(m_ext_zp)), "\n")

  # (2) Intensive: truncated NegBin on 4-country subset (deaths > 0)
  d_zp_int <- d_zp %>% filter(any_death == 1)
  m_int_zp <- tryCatch(
    glmmTMB(n_dead_missing ~ swh_prevweek + swh_prevweek:post_mou +
              year_fac + country,
            data = d_zp_int,
            family = truncated_nbinom2),
    error = function(e) { cat("    glmmTMB err:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(m_int_zp)) {
    report_glmmtmb(m_int_zp, "[C] intensive (trunc NegBin)",
                    d_zp_int$iso_week, nrow(d_zp_int))
  }
}

sink()
cat(sprintf("\nSaved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
