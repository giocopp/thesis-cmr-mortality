# Primary count + volume-controlled model: SWH x post_MoU on daily CMR deaths.
# UNITED primary, IOM comparison. NegBin + Poisson, month_year FE, NW(14).

library(tidyverse)
library(lubridate)
library(fixest)
library(patchwork)
library(sf)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("20  SWH x POST-MOU GRADIENT ON DEATHS (UNITED primary + IOM comparison)\n")
cat("============================================================\n\n")

# ── 1. Load data ─────────────────────────────────────────────
cat("--- 1. Loading data ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

# IOM daily series via the shared helper (default primary filter).
iom_daily <- build_iom_daily()

# UNITED daily series via the shared builder. Defaults (corridor spatial
# join; country in CMR+Med; manner drowned/other_unknown) replicate the
# previous inline filter exactly — single source of truth, see _helpers.R.
united_daily <- build_united_daily()

panel <- panel |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  left_join(united_daily, by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) |>
  add_crossing_exposure() |>   # living_crossings used for lag controls
  mutate(
    swh_next1day  = dplyr::lead(swh, 1),
    swh_next3days = zoo::rollmean(dplyr::lead(swh, 1), k = 3,
                                  fill = NA, align = "left"),
    swh_next5days = zoo::rollmean(dplyr::lead(swh, 1), k = 5,
                                  fill = NA, align = "left"),
    swh_next7days = zoo::rollmean(dplyr::lead(swh, 1), k = 7,
                                  fill = NA, align = "left"),
    lc_lag7  = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 1),
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    log1p_lc_lag7  = log1p(lc_lag7),
    log1p_lc_lag14 = log1p(lc_lag14),
    unit           = 1L,
    year_fac       = factor(year),
    month_year_fac = factor(month_year)
  )

PANEL_END <- max(panel$date)

# Use the sample that has lag14 available (drops first 14 days)
d <- panel |> filter(!is.na(lc_lag14), !is.na(swh_prev5days))

cat(sprintf("  Panel: %s to %s (%d days)\n",
            min(d$date), max(d$date), nrow(d)))
cat(sprintf("  Deaths IOM:    %.0f over %d death-days (%.1f%% zeros)\n",
            sum(d$n_dead_iom), sum(d$n_dead_iom > 0),
            100 * mean(d$n_dead_iom == 0)))
cat(sprintf("  Deaths UNITED: %.0f over %d death-days (%.1f%% zeros)\n",
            sum(d$n_dead_united), sum(d$n_dead_united > 0),
            100 * mean(d$n_dead_united == 0)))

# ── 1b. Exposure sensitivity: past windows + future placebos ─────
# Past exposures: 1-day, 1-3 day, 1-5 day, 1-7 day lagged means.
# Future placebo exposures: next day, next 1-3, next 1-5, next 1-7 days.
# The primary analytical choice is swh_prev5days (1-5d lagged mean). The
# future windows should not predict deaths if the result is truly picking up
# prior sea-state exposure rather than generic seasonal/time patterns.
cat("\n--- 1b. Exposure sensitivity (past windows + future placebos) ---\n")

exposures <- c("swh_lag1", "swh_prev3days", "swh_prev5days", "swh_prevweek",
               "swh_next1day", "swh_next3days", "swh_next5days", "swh_next7days")
window_lookup <- c(
  swh_lag1      = "lag 1d",
  swh_prev3days = "lag 1-3d",
  swh_prev5days = "lag 1-5d",
  swh_prevweek  = "lag 1-7d",
  swh_next1day  = "lead 1d",
  swh_next3days = "lead 1-3d",
  swh_next5days = "lead 1-5d",
  swh_next7days = "lead 1-7d"
)
timing_lookup <- c(
  swh_lag1      = "past",
  swh_prev3days = "past",
  swh_prev5days = "past",
  swh_prevweek  = "past",
  swh_next1day  = "future_placebo",
  swh_next3days = "future_placebo",
  swh_next5days = "future_placebo",
  swh_next7days = "future_placebo"
)
source_lookup <- c(UNITED = "n_dead_united", IOM = "n_dead_iom")

missing_exp <- setdiff(exposures, names(d))
if (length(missing_exp) > 0) {
  stop("Missing exposure columns in panel: ",
       paste(missing_exp, collapse = ", "))
}

# Pin the sensitivity grid to a common row set across all past and future
# windows. This drops only the terminal lead-window days and keeps row-count
# changes from masquerading as window sensitivity.
d_sens <- d |> filter(if_all(all_of(exposures), ~ !is.na(.x)))

fit_sens_one <- function(x, outcome, source_label) {
  f <- as.formula(sprintf("%s ~ %s + %s:post_mou | month_year_fac",
                          outcome, x, x))
  models <- list(
    NegBin  = fenegbin(f, data = d_sens, vcov = NW(14), panel.id = ~unit + date),
    Poisson = fepois (f, data = d_sens, vcov = NW(14), panel.id = ~unit + date)
  )
  imap_dfr(models, function(m, fam) {
    ct <- coeftable(m, vcov = NW(14))
    r_main <- which(rownames(ct) == x)
    r_int  <- grep(":post_mou$", rownames(ct))
    r_int  <- r_int[grepl(x, rownames(ct)[r_int], fixed = TRUE)]
    tibble(
      source   = source_label,
      timing   = timing_lookup[[x]],
      window   = window_lookup[[x]],
      family   = fam,
      exposure = x,
      n_obs    = nobs(m),
      b1       = ct[r_main, 1],
      se1      = ct[r_main, 2],
      p1       = 2 * pnorm(-abs(ct[r_main, 1] / ct[r_main, 2])),
      b3       = ct[r_int,  1],
      se3      = ct[r_int,  2],
      p3       = 2 * pnorm(-abs(ct[r_int, 1] / ct[r_int, 2]))
    )
  })
}

sens_tbl <- map_dfr(names(source_lookup), function(src) {
  map_dfr(exposures, function(x) {
    fit_sens_one(x, source_lookup[[src]], src)
  })
}) |>
  mutate(
    source = factor(source, levels = c("UNITED", "IOM")),
    timing = factor(timing, levels = c("past", "future_placebo")),
    exposure_order = match(exposure, exposures)
  ) |>
  arrange(source, timing, exposure_order, family)

cat(sprintf("  Sensitivity sample: N = %d days (%s to %s)\n",
            nrow(d_sens), min(d_sens$date), max(d_sens$date)))
cat("\n  SWH × post_MoU interaction (b3) across windows:\n")
cat(sprintf("  %-7s %-14s %-9s %-8s  %+10s  %10s  %10s\n",
            "source", "timing", "window", "family", "b3", "SE", "p"))
for (i in seq_len(nrow(sens_tbl))) {
  r <- sens_tbl[i, ]
  star <- if (r$p3 < 0.05) " *" else ""
  cat(sprintf("  %-7s %-14s %-9s %-8s  %+10.3f  %10.3f  %10.4f%s\n",
              as.character(r$source), as.character(r$timing), r$window,
              r$family, r$b3, r$se3, r$p3, star))
}

# CSV for downstream tables/figures
sens_csv <- tbl_path("05_analysis", "01_exposure_sensitivity.csv")
write.csv(sens_tbl, sens_csv, row.names = FALSE)
cat(sprintf("  Saved: %s\n", sens_csv))

# ── 2. Count models: NegBin + Poisson, IOM + UNITED ──
cat("\n--- 2. Count models ---\n")

# IOM
m_nb   <- fenegbin(n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                   data = d, vcov = NW(14), panel.id = ~unit + date)
m_pois <- fepois  (n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                   data = d, vcov = NW(14), panel.id = ~unit + date)

# UNITED (same spec, same sample)
m_nb_u   <- fenegbin(n_dead_united ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                     data = d, vcov = NW(14), panel.id = ~unit + date)
m_pois_u <- fepois  (n_dead_united ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                     data = d, vcov = NW(14), panel.id = ~unit + date)

# ── 2b. Volume-controlled count model: log(crossing_attempts) as free covariate ──
d_rate <- d |>
  filter(crossing_attempts > 0) |>
  mutate(log_crossing_attempts = log(crossing_attempts))

m_rate_u <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  log_crossing_attempts | month_year_fac,
  data = d_rate, vcov = NW(14), panel.id = ~unit + date)

m_rate <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               log_crossing_attempts | month_year_fac,
  data = d_rate, vcov = NW(14), panel.id = ~unit + date)

# ── 2c. Elasticity check: Wald test on log(crossing_attempts) vs 1 ──

elast_test <- function(m, xname, label) {
  ct <- coeftable(m, vcov = NW(14))
  b  <- ct[xname, 1]
  se <- ct[xname, 2]
  z  <- (b - 1) / se
  tibble(
    source     = label,
    b_exposure = b,
    se         = se,
    z_vs_1     = z,
    p_vs_1     = 2 * pnorm(-abs(z))
  )
}

elast_tbl <- bind_rows(
  elast_test(m_rate_u, "log_crossing_attempts", "UNITED"),
  elast_test(m_rate,   "log_crossing_attempts", "IOM")
)

print_elast <- function(tbl) {
  cat(sprintf("  %-7s %12s %10s %10s %10s\n",
              "source", "b_exposure", "SE", "z_vs_1", "p_vs_1"))
  for (i in seq_len(nrow(tbl))) {
    r <- tbl[i, ]
    cat(sprintf("  %-7s %+12.3f %10.3f %+10.2f %10.4f\n",
                r$source, r$b_exposure, r$se, r$z_vs_1, r$p_vs_1))
  }
}

cat(sprintf("  Count sample: N = %d days (all)\n", nrow(d)))
cat(sprintf("  Rate sample (common): N = %d days (drops %d zero-crossing days)\n",
            nrow(d_rate), nrow(d) - nrow(d_rate)))

slope_summary <- function(m, label, family, estimand) {
  ct <- coeftable(m, vcov = NW(14))
  co <- coef(m)
  V  <- vcov(m, vcov = NW(14))
  b_pre <- unname(co["swh_prev5days"])
  b_shift <- unname(co["swh_prev5days:post_mou"])
  se_pre <- ct["swh_prev5days", 2]
  se_shift <- ct["swh_prev5days:post_mou", 2]
  b_post <- b_pre + b_shift
  v_post <- V["swh_prev5days", "swh_prev5days"] +
    V["swh_prev5days:post_mou", "swh_prev5days:post_mou"] +
    2 * V["swh_prev5days", "swh_prev5days:post_mou"]
  se_post <- sqrt(v_post)
  tibble(
    source = label,
    family = family,
    estimand = estimand,
    n_obs = nobs(m),
    b_pre = b_pre,
    se_pre = se_pre,
    p_pre = 2 * pnorm(-abs(b_pre / se_pre)),
    irr_pre = exp(b_pre),
    b_shift = b_shift,
    se_shift = se_shift,
    p_shift = 2 * pnorm(-abs(b_shift / se_shift)),
    irr_shift = exp(b_shift),
    b_post = b_post,
    se_post = se_post,
    p_post = 2 * pnorm(-abs(b_post / se_post)),
    irr_post = exp(b_post)
  )
}

count_slopes <- bind_rows(
  slope_summary(m_nb_u,   "UNITED", "NegBin",  "count"),
  slope_summary(m_pois_u, "UNITED", "Poisson", "count"),
  slope_summary(m_nb,     "IOM",    "NegBin",  "count"),
  slope_summary(m_pois,   "IOM",    "Poisson", "count")
)

rate_slopes <- bind_rows(
  slope_summary(m_rate_u, "UNITED", "Poisson", "rate-free-exposure"),
  slope_summary(m_rate,   "IOM",    "Poisson", "rate-free-exposure")
)

print_slope_summary <- function(tbl, show_irr = FALSE) {
  cat(sprintf("  %-7s %-8s %6s  %10s %10s %10s  %10s %10s %10s  %10s %10s %10s",
              "source", "family", "N", "b_pre", "SE", "p",
              "b_shift", "SE", "p", "b_post", "SE", "p"))
  if (show_irr) cat(sprintf("  %10s %10s %10s", "IRR_pre", "IRR_shift", "IRR_post"))
  cat("\n")
  for (i in seq_len(nrow(tbl))) {
    r <- tbl[i, ]
    cat(sprintf("  %-7s %-8s %6d  %+10.3f %10.3f %10.4f  %+10.3f %10.3f %10.4f  %+10.3f %10.3f %10.4f",
                r$source, r$family, r$n_obs,
                r$b_pre, r$se_pre, r$p_pre,
                r$b_shift, r$se_shift, r$p_shift,
                r$b_post, r$se_post, r$p_post))
    if (show_irr) {
      cat(sprintf("  %10.3f %10.3f %10.3f",
                  r$irr_pre, r$irr_shift, r$irr_post))
    }
    cat("\n")
  }
}

cat("\n  Count models, NW(14):\n")
print(etable(m_nb_u, m_pois_u, m_nb, m_pois,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED NB", "UNITED Poiss",
                         "IOM NB", "IOM Poiss")))

cat("\n  Count, slope decomposition, NW(14):\n")
print_slope_summary(count_slopes, show_irr = FALSE)

cat("\n  Volume-controlled model (free log exposure), NW(14):\n")
print(etable(m_rate_u, m_rate,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED rate", "IOM rate")))

cat("\n  Rate, slope decomposition, NW(14):\n")
print_slope_summary(rate_slopes, show_irr = FALSE)

cat("\n  Elasticity: log(crossing_attempts) vs 1, NW(14):\n")
print_elast(elast_tbl)

cat("\n  All models, cluster(month_year):\n")
print(etable(m_nb_u, m_pois_u, m_nb, m_pois, m_rate_u, m_rate,
             vcov = ~month_year_fac, se.below = TRUE,
             headers = c("UNITED NB", "UNITED Poiss",
                         "IOM NB", "IOM Poiss",
                         "UNITED rate", "IOM rate")))

# ── 3. Robustness: lagged crossing controls ──────────────────
cat("\n--- 3. Robustness: lagged crossing controls ---\n")

m_nb_u_lag7  <- fenegbin(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou + log1p_lc_lag7 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_nb_u_lag14 <- fenegbin(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou + log1p_lc_lag14 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_pois_u_lag7  <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou + log1p_lc_lag7 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_pois_u_lag14 <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou + log1p_lc_lag14 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_nb_lag7  <- fenegbin(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou + log1p_lc_lag7 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_nb_lag14 <- fenegbin(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou + log1p_lc_lag14 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_pois_lag7  <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou + log1p_lc_lag7 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_pois_lag14 <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou + log1p_lc_lag14 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  UNITED NegBin with crossing controls, NW(14):\n")
print(etable(m_nb_u, m_nb_u_lag7, m_nb_u_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))

cat("\n  UNITED Poisson with crossing controls, NW(14):\n")
print(etable(m_pois_u, m_pois_u_lag7, m_pois_u_lag14,
             vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))

cat("\n  IOM NegBin with crossing controls, NW(14):\n")
print(etable(m_nb, m_nb_lag7, m_nb_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))

cat("\n  IOM Poisson with crossing controls, NW(14):\n")
print(etable(m_pois, m_pois_lag7, m_pois_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))

# ── 3b. Robustness: ACLED Libya/Tunisia conflict + lag-14 ────
# ACLED columns come from daily_panel_complete.RDS; broadcast weekly, so
# month_year FE absorbs most monthly variation. Residual weekly variation
# proxies push-factor shocks orthogonal to weather.
cat("\n--- 3b. Robustness: ACLED Libya + Tunisia + lag-14 ---\n")

d <- d |>
  mutate(
    log1p_libya   = log1p(libya_conflict),
    log1p_tunisia = log1p(tunisia_conflict)
  )

m_nb_u_ctl <- fenegbin(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_pois_u_ctl <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_nb_ctl <- fenegbin(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_pois_ctl <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  UNITED NegBin with full controls, NW(14):\n")
print(etable(m_nb_u, m_nb_u_lag14, m_nb_u_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 14d", "Lag 14d + ACLED")))

cat("\n  IOM NegBin with full controls, NW(14):\n")
print(etable(m_nb, m_nb_lag14, m_nb_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 14d", "Lag 14d + ACLED")))

# ── 3c. Robustness: boundary leakage (drop +/-30d around MoU) ─
# Policy implementation is gradual; sharp 2017-02-02 cut may attribute
# leakage windows to the wrong regime. Dropping +/-30 days tests whether
# the post-MoU shift survives boundary noise.
cat("\n--- 3c. Robustness: boundary leakage (drop +/-30d around 2017-02-02) ---\n")

drop_window <- 30L
d_bl <- d |> filter(abs(as.integer(date - MOU_DATE)) > drop_window)

cat(sprintf("  Boundary-leakage sample: N = %d days (drops %d days around MoU)\n",
            nrow(d_bl), nrow(d) - nrow(d_bl)))

m_nb_u_bl   <- fenegbin(n_dead_united ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac, data = d_bl, vcov = NW(14),
                        panel.id = ~unit + date)
m_pois_u_bl <- fepois  (n_dead_united ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac, data = d_bl, vcov = NW(14),
                        panel.id = ~unit + date)
m_nb_bl     <- fenegbin(n_dead_iom    ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac, data = d_bl, vcov = NW(14),
                        panel.id = ~unit + date)
m_pois_bl   <- fepois  (n_dead_iom    ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac, data = d_bl, vcov = NW(14),
                        panel.id = ~unit + date)

m_nb_u_bl_ctl <- fenegbin(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d_bl, vcov = NW(14), panel.id = ~unit + date)
m_nb_bl_ctl <- fenegbin(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d_bl, vcov = NW(14), panel.id = ~unit + date)

cat("\n  Boundary leakage (drop +/-30d), NB, NW(14):\n")
print(etable(m_nb_u, m_nb_u_bl, m_nb_u_bl_ctl, m_nb, m_nb_bl, m_nb_bl_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED full", "UNITED -30d", "UNITED -30d + ctl",
                         "IOM full",    "IOM -30d",    "IOM -30d + ctl")))

# ── 3d. Robustness: Mare Nostrum exclusion ────────────────────
# Period 1 is internally heterogeneous: Mare Nostrum (state-led, very
# active) ran to 2014-10-31; Triton/Sophia (border-control + NGO SAR)
# from 2014-11-01. Restricting to date >= TRITON_START tests whether
# the pre-MoU baseline is an MN artefact.
cat("\n--- 3d. Robustness: Mare Nostrum exclusion (date >= 2014-11-01) ---\n")

d_post_mn <- d |> filter(date >= TRITON_START)
cat(sprintf("  Triton-onwards sample: N = %d days (drops %d MN days)\n",
            nrow(d_post_mn), nrow(d) - nrow(d_post_mn)))

m_nb_u_mn   <- fenegbin(n_dead_united ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac, data = d_post_mn, vcov = NW(14),
                        panel.id = ~unit + date)
m_pois_u_mn <- fepois  (n_dead_united ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac, data = d_post_mn, vcov = NW(14),
                        panel.id = ~unit + date)
m_nb_mn     <- fenegbin(n_dead_iom    ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac, data = d_post_mn, vcov = NW(14),
                        panel.id = ~unit + date)
m_pois_mn   <- fepois  (n_dead_iom    ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac, data = d_post_mn, vcov = NW(14),
                        panel.id = ~unit + date)

m_nb_u_mn_ctl <- fenegbin(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d_post_mn, vcov = NW(14), panel.id = ~unit + date)
m_nb_mn_ctl <- fenegbin(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d_post_mn, vcov = NW(14), panel.id = ~unit + date)

cat("\n  Mare Nostrum exclusion (post-MN sample), NB, NW(14):\n")
print(etable(m_nb_u, m_nb_u_mn, m_nb_u_mn_ctl, m_nb, m_nb_mn, m_nb_mn_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED full", "UNITED post-MN", "UNITED post-MN + ctl",
                         "IOM full",    "IOM post-MN",    "IOM post-MN + ctl")))

# ── 4. Save text output ──────────────────────────────────────
# No year-by-year gradients here; regime cuts are in 31_united_periods.R.
cat("\n--- 4. Saving results ---\n")

sink_file <- tbl_path("05_analysis", "01_primary_model.txt")
sink(sink_file)

cat("20  SWH x POST-MOU GRADIENT ON DEATHS (UNITED primary + IOM comparison)\n")
cat("=======================================================================\n")
cat(sprintf("Sample: %s to %s (N = %d days)\n",
            min(d$date), max(d$date), nrow(d)))
cat(sprintf("Deaths UNITED (corridor, drowned + other/unknown): %.0f\n",
            sum(d$n_dead_united)))
cat(sprintf("Deaths IOM    (corridor, drowning + mixed):       %.0f\n",
            sum(d$n_dead_iom)))
cat("\nModels (month_year FE; NW(14) SEs; UNITED primary, IOM comparison):\n")
cat("  count           : deaths ~ swh_prev5days + swh_prev5days:post_mou (NegBin, Poisson)\n")
cat("  volume-controlled: count spec + log(crossing_attempts) free       (Poisson)\n")
cat("  crossing_attempts = frx_persons + lcg_tcg_pushbacks + n_dead_united_for_ct\n")
cat("Slope decomposition: b_pre, b_shift (=SWH:post_mou), b_post = b_pre + b_shift (delta-method SE).\n")
cat("N = estimation N after fixest drops all-zero FE cells.\n\n")

cat("=== COUNT MODELS: UNITED PRIMARY, IOM COMPARISON ===\n")
cat(sprintf("Count sample: N = %d days.\n\n", nrow(d)))

cat("--- NW(14) SEs ---\n")
print(etable(m_nb_u, m_pois_u, m_nb, m_pois,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED NB", "UNITED Poiss",
                         "IOM NB", "IOM Poiss")))

cat("\n--- Slope decomposition, NW(14) ---\n")
print_slope_summary(count_slopes, show_irr = FALSE)

cat("\n--- Cluster(month_year) SEs ---\n")
print(etable(m_nb_u, m_pois_u, m_nb, m_pois,
             vcov = ~month_year_fac, se.below = TRUE,
             headers = c("UNITED NB", "UNITED Poiss",
                         "IOM NB", "IOM Poiss")))

cat("\n\n=== VOLUME-CONTROLLED MODEL (free log exposure) ===\n")
cat("crossing_attempts = frx_persons + lcg_tcg_pushbacks + n_dead_united_for_ct\n")
cat(sprintf(
  "Rate sample: N=%d days (drops %d zero-crossing days).\n\n",
  nrow(d_rate), nrow(d) - nrow(d_rate)))

cat("--- NW(14) SEs ---\n")
print(etable(m_rate_u, m_rate,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED rate", "IOM rate")))

cat("\n--- Slope decomposition, NW(14) ---\n")
print_slope_summary(rate_slopes, show_irr = FALSE)

cat("\n--- Elasticity: log(crossing_attempts) vs 1, NW(14) ---\n")
print_elast(elast_tbl)

cat("\n--- Cluster(month_year) SEs ---\n")
print(etable(m_rate_u, m_rate,
             vcov = ~month_year_fac, se.below = TRUE,
             headers = c("UNITED rate", "IOM rate")))

cat("\n\n=== EXPOSURE SENSITIVITY (PAST WINDOWS + FUTURE PLACEBOS) ===\n")
cat("Past exposures: lag 1d, 1-3d, 1-5d, 1-7d means.\n")
cat("Future placebo exposures: lead 1d, 1-3d, 1-5d, 1-7d means.\n")
cat(sprintf("Shared sensitivity sample: N = %d days (%s to %s).\n",
            nrow(d_sens), min(d_sens$date), max(d_sens$date)))
cat("b3 = SWH:post_mou; full coefficients in the CSV.\n\n")
cat(sprintf("  %-7s %-14s %-9s %-8s  %+10s  %10s  %10s\n",
            "source", "timing", "window", "family", "b3", "SE", "p"))
for (i in seq_len(nrow(sens_tbl))) {
  r <- sens_tbl[i, ]
  star <- if (r$p3 < 0.05) " *" else ""
  cat(sprintf("  %-7s %-14s %-9s %-8s  %+10.3f  %10.3f  %10.4f%s\n",
              as.character(r$source), as.character(r$timing), r$window,
              r$family, r$b3, r$se3, r$p3, star))
}
cat(sprintf("\n  Full table: %s\n", sens_csv))

cat("\n\n=== COUNT ROBUSTNESS: LAGGED CROSSING CONTROLS ===\n\n")
cat("UNITED NegBin:\n")
print(etable(m_nb_u, m_nb_u_lag7, m_nb_u_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))
cat("\nUNITED Poisson:\n")
print(etable(m_pois_u, m_pois_u_lag7, m_pois_u_lag14,
             vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))
cat("\nIOM NegBin:\n")
print(etable(m_nb, m_nb_lag7, m_nb_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))
cat("\nIOM Poisson:\n")
print(etable(m_pois, m_pois_lag7, m_pois_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))

# Summary table
cat("\n\n=== COUNT SUMMARY: b3 (SWH x post_mou) ===\n\n")
for (info in list(
  list(m_nb_u,         "UNITED NegBin, no control"),
  list(m_nb_u_lag7,    "UNITED NegBin, lag 7d"),
  list(m_nb_u_lag14,   "UNITED NegBin, lag 14d"),
  list(m_pois_u,       "UNITED Poisson, no control"),
  list(m_pois_u_lag7,  "UNITED Poisson, lag 7d"),
  list(m_pois_u_lag14, "UNITED Poisson, lag 14d"),
  list(m_nb,           "IOM NegBin, no control"),
  list(m_nb_lag7,      "IOM NegBin, lag 7d"),
  list(m_nb_lag14,     "IOM NegBin, lag 14d"),
  list(m_pois,         "IOM Poisson, no control"),
  list(m_pois_lag7,    "IOM Poisson, lag 7d"),
  list(m_pois_lag14,   "IOM Poisson, lag 14d")
)) {
  ct <- coeftable(info[[1]], vcov = NW(14))
  r <- grep(":post_mou", rownames(ct))
  p <- 2 * pnorm(-abs(ct[r, 1] / ct[r, 2]))
  ct_cl <- coeftable(info[[1]], vcov = ~month_year_fac)
  p_cl <- 2 * pnorm(-abs(ct_cl[r, 1] / ct_cl[r, 2]))
  cat(sprintf("  %-25s  b3=%+.3f  SE_NW=%.3f  p_NW=%.4f  SE_cl=%.3f  p_cl=%.4f\n",
              info[[2]], ct[r, 1], ct[r, 2], p, ct_cl[r, 2], p_cl))
}

cat("\n\n=== CREDIBILITY: FULL CONTROLS (lag-14 + ACLED Libya/Tunisia) ===\n\n")
cat("UNITED NegBin:\n")
print(etable(m_nb_u, m_nb_u_lag14, m_nb_u_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 14d", "Lag 14d + ACLED")))
cat("\nIOM NegBin:\n")
print(etable(m_nb, m_nb_lag14, m_nb_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 14d", "Lag 14d + ACLED")))
cat("\nUNITED Poisson:\n")
print(etable(m_pois_u, m_pois_u_lag14, m_pois_u_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 14d", "Lag 14d + ACLED")))
cat("\nIOM Poisson:\n")
print(etable(m_pois, m_pois_lag14, m_pois_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 14d", "Lag 14d + ACLED")))

cat("\n\n=== CREDIBILITY: BOUNDARY LEAKAGE (drop +/-30d around 2017-02-02) ===\n\n")
cat(sprintf("Sample: N=%d days (drops %d days around MoU).\n\n",
            nrow(d_bl), nrow(d) - nrow(d_bl)))
print(etable(m_nb_u, m_nb_u_bl, m_nb_u_bl_ctl, m_nb, m_nb_bl, m_nb_bl_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED full", "UNITED -30d", "UNITED -30d + ctl",
                         "IOM full",    "IOM -30d",    "IOM -30d + ctl")))

cat("\n\n=== CREDIBILITY: MARE NOSTRUM EXCLUSION (date >= 2014-11-01) ===\n\n")
cat(sprintf("Sample: N=%d days (drops %d MN days).\n\n",
            nrow(d_post_mn), nrow(d) - nrow(d_post_mn)))
print(etable(m_nb_u, m_nb_u_mn, m_nb_u_mn_ctl, m_nb, m_nb_mn, m_nb_mn_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED full", "UNITED post-MN", "UNITED post-MN + ctl",
                         "IOM full",    "IOM post-MN",    "IOM post-MN + ctl")))

cat("\n\n=== EXTENDED b3 SUMMARY: credibility checks ===\n\n")
for (info in list(
  list(m_nb_u_ctl,    "UNITED NegBin, lag-14 + ACLED"),
  list(m_nb_ctl,      "IOM NegBin, lag-14 + ACLED"),
  list(m_pois_u_ctl,  "UNITED Poisson, lag-14 + ACLED"),
  list(m_pois_ctl,    "IOM Poisson, lag-14 + ACLED"),
  list(m_nb_u_bl,     "UNITED NegBin, drop +/-30d (no ctl)"),
  list(m_nb_u_bl_ctl, "UNITED NegBin, drop +/-30d + full ctl"),
  list(m_nb_bl,       "IOM NegBin, drop +/-30d (no ctl)"),
  list(m_nb_bl_ctl,   "IOM NegBin, drop +/-30d + full ctl"),
  list(m_nb_u_mn,     "UNITED NegBin, post-MN (no ctl)"),
  list(m_nb_u_mn_ctl, "UNITED NegBin, post-MN + full ctl"),
  list(m_nb_mn,       "IOM NegBin, post-MN (no ctl)"),
  list(m_nb_mn_ctl,   "IOM NegBin, post-MN + full ctl")
)) {
  ct <- coeftable(info[[1]], vcov = NW(14))
  r <- grep(":post_mou", rownames(ct))
  p <- 2 * pnorm(-abs(ct[r, 1] / ct[r, 2]))
  ct_cl <- coeftable(info[[1]], vcov = ~month_year_fac)
  p_cl <- 2 * pnorm(-abs(ct_cl[r, 1] / ct_cl[r, 2]))
  cat(sprintf("  %-42s  b3=%+.3f  SE_NW=%.3f  p_NW=%.4f  SE_cl=%.3f  p_cl=%.4f\n",
              info[[2]], ct[r, 1], ct[r, 2], p, ct_cl[r, 2], p_cl))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

# ── 5. LaTeX tables (\input'd by paper/thesis.qmd) ───────────
cat("\n--- 5. Writing LaTeX tables ---\n")

sig_stars <- function(p) {
  ifelse(p < 0.001, "^{***}",
  ifelse(p < 0.01,  "^{**}",
  ifelse(p < 0.05,  "^{*}", "")))
}
fcoef <- function(b, p) sprintf("$%+.3f%s$", b, sig_stars(p))
fse   <- function(se)   sprintf("(%.3f)", se)
fint  <- function(x)    formatC(round(x), format = "d", big.mark = ",")
fr3   <- function(x)    formatC(x, format = "f", digits = 3)
fp    <- function(p)    ifelse(p < 0.001, "$< 0.001$", sprintf("$%.3f$", p))

pr2_nb_u   <- fixest::r2(m_nb_u,   "pr2")
pr2_pois_u <- fixest::r2(m_pois_u, "pr2")
pr2_nb     <- fixest::r2(m_nb,     "pr2")
pr2_pois   <- fixest::r2(m_pois,   "pr2")

# --- tab:primary -----------------------------------------------------
cs <- count_slopes
L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[h!]")
add("\\centering")
add("\\small")
add("\\caption{Primary count estimates of the SWH--mortality slope shift around the 2 February 2017 MoU.}")
add("\\label{tab:primary}")
add("\\begin{tabular}{lcccc}")
add("\\hline")
add("                              & \\multicolumn{2}{c}{UNITED} & \\multicolumn{2}{c}{IOM (comparison)} \\\\")
add("                              & NegBin & Poisson & NegBin & Poisson \\\\")
add("\\hline")
add("$\\beta_1$: SWH$_{t-1:t-5}$            & ",
    fcoef(cs$b_pre[1],   cs$p_pre[1]),   " & ", fcoef(cs$b_pre[2],   cs$p_pre[2]),   " & ",
    fcoef(cs$b_pre[3],   cs$p_pre[3]),   " & ", fcoef(cs$b_pre[4],   cs$p_pre[4]),   " \\\\")
add("                                      & ",
    fse(cs$se_pre[1]),   " & ", fse(cs$se_pre[2]),   " & ",
    fse(cs$se_pre[3]),   " & ", fse(cs$se_pre[4]),   " \\\\")
add("$\\beta_3$: SWH$_{t-1:t-5}$ $\\times$ Post-MoU & ",
    fcoef(cs$b_shift[1], cs$p_shift[1]), " & ", fcoef(cs$b_shift[2], cs$p_shift[2]), " & ",
    fcoef(cs$b_shift[3], cs$p_shift[3]), " & ", fcoef(cs$b_shift[4], cs$p_shift[4]), " \\\\")
add("                                      & ",
    fse(cs$se_shift[1]), " & ", fse(cs$se_shift[2]), " & ",
    fse(cs$se_shift[3]), " & ", fse(cs$se_shift[4]), " \\\\")
add("$\\beta_1+\\beta_3$: implied post-MoU slope    & ",
    fcoef(cs$b_post[1],  cs$p_post[1]),  " & ", fcoef(cs$b_post[2],  cs$p_post[2]),  " & ",
    fcoef(cs$b_post[3],  cs$p_post[3]),  " & ", fcoef(cs$b_post[4],  cs$p_post[4]),  " \\\\")
add("                                      & ",
    fse(cs$se_post[1]),  " & ", fse(cs$se_post[2]),  " & ",
    fse(cs$se_post[3]),  " & ", fse(cs$se_post[4]),  " \\\\")
add("\\hline")
add("Month--year fixed effects             & Yes            & Yes           & Yes            & Yes          \\\\")
add("Newey--West SEs (lag 14)              & Yes            & Yes           & Yes            & Yes          \\\\")
add("Observations (days)                   & ",
    fint(cs$n_obs[1]), " & ", fint(cs$n_obs[2]), " & ",
    fint(cs$n_obs[3]), " & ", fint(cs$n_obs[4]), " \\\\")
add("Pseudo $R^{2}$                        & ",
    fr3(pr2_nb_u),   " & ", fr3(pr2_pois_u), " & ",
    fr3(pr2_nb),     " & ", fr3(pr2_pois),   " \\\\")
add("\\hline")
add("\\multicolumn{5}{l}{\\footnotesize Stars denote two-sided $p$-values: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\multicolumn{5}{l}{\\footnotesize Newey--West standard errors in parentheses; delta-method SEs for $\\beta_1+\\beta_3$.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
out_primary <- tbl_path("05_analysis", "01_primary_model.tex")
writeLines(L, out_primary)
cat(sprintf("  Saved: %s\n", out_primary))

# --- tab:rate --------------------------------------------------------
rs <- rate_slopes  # 1: UNITED, 2: IOM (Poisson only)
ct_u <- coeftable(m_rate_u, vcov = NW(14))
ct_i <- coeftable(m_rate,   vcov = NW(14))
g_u  <- ct_u["log_crossing_attempts", 1]; g_se_u <- ct_u["log_crossing_attempts", 2]
g_i  <- ct_i["log_crossing_attempts", 1]; g_se_i <- ct_i["log_crossing_attempts", 2]
g_p_u <- 2 * pnorm(-abs(g_u  / g_se_u))
g_p_i <- 2 * pnorm(-abs(g_i  / g_se_i))

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[h!]")
add("\\centering")
add("\\small")
add("\\caption{Volume-controlled Poisson model with crossing-volume control.}")
add("\\label{tab:rate}")
add("\\begin{tabular}{lcc}")
add("\\hline")
add("                              & UNITED & IOM (comparison) \\\\")
add("\\hline")
add("$\\beta_1$: SWH$_{t-1:t-5}$            & ",
    fcoef(rs$b_pre[1],   rs$p_pre[1]),   " & ", fcoef(rs$b_pre[2],   rs$p_pre[2]),   " \\\\")
add("                                      & ",
    fse(rs$se_pre[1]),   " & ", fse(rs$se_pre[2]),   " \\\\")
add("$\\beta_3$: SWH$_{t-1:t-5}$ $\\times$ Post-MoU & ",
    fcoef(rs$b_shift[1], rs$p_shift[1]), " & ", fcoef(rs$b_shift[2], rs$p_shift[2]), " \\\\")
add("                                      & ",
    fse(rs$se_shift[1]), " & ", fse(rs$se_shift[2]), " \\\\")
add("$\\gamma$: $\\log C_t$                  & ",
    fcoef(g_u, g_p_u), " & ", fcoef(g_i, g_p_i), " \\\\")
add("                                      & ",
    fse(g_se_u), " & ", fse(g_se_i), " \\\\")
add("$\\beta_1+\\beta_3$: implied post-MoU slope    & ",
    fcoef(rs$b_post[1],  rs$p_post[1]),  " & ", fcoef(rs$b_post[2],  rs$p_post[2]),  " \\\\")
add("                                      & ",
    fse(rs$se_post[1]),  " & ", fse(rs$se_post[2]),  " \\\\")
add("\\hline")
add("Wald: $H_{0}: \\gamma = 1$, $z$        & ",
    sprintf("$%+.2f$", elast_tbl$z_vs_1[1]), " & ",
    sprintf("$%+.2f$", elast_tbl$z_vs_1[2]), " \\\\")
add("$p$-value                             & ", fp(elast_tbl$p_vs_1[1]),
    " & ", fp(elast_tbl$p_vs_1[2]), " \\\\")
add("Month--year fixed effects             & Yes            & Yes            \\\\")
add("Newey--West SEs (lag 14)              & Yes            & Yes            \\\\")
add("Observations (days, $C_t>0$)          & ",
    fint(rs$n_obs[1]), " & ", fint(rs$n_obs[2]), " \\\\")
add("\\hline")
add("\\multicolumn{3}{l}{\\footnotesize Stars denote two-sided $p$-values: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\multicolumn{3}{l}{\\footnotesize Newey--West standard errors in parentheses; delta-method SEs for $\\beta_1+\\beta_3$.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
out_rate <- tbl_path("05_analysis", "01_rate_model.tex")
writeLines(L, out_rate)
cat(sprintf("  Saved: %s\n", out_rate))

# --- tab:appx-exposure ----------------------------------------------
window_order <- c("lag 1d", "lag 1-3d", "lag 1-5d", "lag 1-7d",
                  "lead 1d", "lead 1-3d", "lead 1-5d", "lead 1-7d")
sens_wide <- sens_tbl |>
  dplyr::select(source, timing, window, family, b3, se3, p3) |>
  tidyr::pivot_wider(names_from = family,
                     values_from = c(b3, se3, p3)) |>
  dplyr::arrange(source, timing, match(window, window_order))

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[h!]")
add("\\centering")
add("\\small")
add("\\caption{SWH-window grid: SWH$\\times$Post-MoU shift coefficient ($\\beta_3$).}")
add("\\label{tab:appx-exposure}")
add("\\begin{tabular}{llcccc}")
add("\\hline")
add("Source & Window      & \\multicolumn{2}{c}{NegBin} & \\multicolumn{2}{c}{Poisson} \\\\")
add("       &             & $\\beta_3$ & SE    & $\\beta_3$ & SE \\\\")
add("\\hline")
write_block <- function(timing_label, header) {
  add(sprintf("\\multicolumn{6}{l}{\\textit{%s}} \\\\", header))
  rows <- sens_wide[sens_wide$timing == timing_label, ]
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    add(sprintf("%-7s & %-11s & %-15s & %.3f & %-15s & %.3f \\\\",
                as.character(r$source), r$window,
                fcoef(r$b3_NegBin,  r$p3_NegBin),  r$se3_NegBin,
                fcoef(r$b3_Poisson, r$p3_Poisson), r$se3_Poisson))
  }
}
write_block("past",           "Past (lagged) windows")
write_block("future_placebo", "Future (placebo) windows")
add("\\hline")
add("\\multicolumn{6}{l}{\\footnotesize SE columns report Newey--West standard errors. Stars: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add(sprintf("\\multicolumn{6}{l}{\\footnotesize Shared sample: %s days.}",
            fint(nrow(d_sens))))
add("\\end{tabular}")
add("\\end{table}")
out_exposure <- tbl_path("05_analysis", "01_exposure_sensitivity.tex")
writeLines(L, out_exposure)
cat(sprintf("  Saved: %s\n", out_exposure))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
