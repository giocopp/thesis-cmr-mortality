# 20_primary_model.R
# ===================
# Primary reduced-form model and robustness.
#
# Specification
# -------------
# Primary (no crossing control):
#   deaths ~ swh_prev5days + swh_prev5days:post_mou | month_year
#
# Interpretation
# --------------
# post_mou is a POLICY-EVENT INDICATOR (1 from 2017-07-01), not a measure
# of SAR capacity. The interaction swh_prev5days:post_mou is the REDUCED
# FORM for the guardrail argument: did the rate at which weather translates
# into deaths shift around the MoU? It does not by itself identify SAR as
# the channel -- many things changed at this date (Minniti Code, Salvini
# decrees, NGO targeting, LCG ramp-up, smuggler/boat composition).
#
# The MECHANISM check -- does the gradient track SAR capacity directly --
# lives in 23_mechanism_interactions.R, which uses three continuous SAR
# moderators (share, log absolute events, log absolute persons) and is
# the substantive complement to this script.
#
# So the division of labor is:
#   20 -> reduced form: gradient shift around the policy event
#   23 -> mechanism:    gradient covariation with SAR capacity
#
# Estimated on BOTH IOM MMP and UNITED death series (matched filter:
# country in CMR + Mediterranean; cause = drowning + other/unknown;
# spatial join to core corridor polygon). UNITED serves as an independent
# confirmation: same sample window, broader sourcing (press, NGOs,
# academic archives), typically ~13% more deaths in sample and tighter CIs.
#
# Both NegBin and Poisson QMLE, NW(14) SEs.
# Poisson QMLE is consistent under weaker assumptions (only needs correct
# conditional mean); NegBin adds a variance assumption that improves
# efficiency. Agreement between the two (and between IOM and UNITED)
# validates the result.
#
# Robustness:
#   - Lagged crossing controls (lag 7d, lag 14d) as covariates
#   - Cluster(month_year) SEs alongside NW(14)
#   - Boat composition controls for the rate model: 27_rate_with_boat_controls.R
#
# In:  analysis/data/daily_panel_complete.RDS
#      data/processed/iom_mmp_incidents.RDS
#      data/processed/united_incidents.RDS
#      data/processed/core_corridor.RDS
# Out: output/tables/20_primary_model.txt

library(tidyverse)
library(lubridate)
library(fixest)
library(patchwork)
library(sf)

BASE_DIR   <- here::here()
MOU_DATE   <- as.Date("2017-07-01")
START_DATE <- as.Date("2014-01-01")

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("05d  PRIMARY MODEL: SWH x POST-MOU -> DEATHS\n")
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

panel <- panel %>%
  left_join(iom_daily %>% rename(n_dead_iom = n_dead_missing), by = "date") %>%
  left_join(united_daily, by = "date") %>%
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) %>%
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
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
d <- panel %>% filter(!is.na(lc_lag14), !is.na(swh_prev5days))

cat(sprintf("  Panel: %s to %s (%d days)\n",
            min(d$date), max(d$date), nrow(d)))
cat(sprintf("  Deaths IOM:    %.0f over %d death-days (%.1f%% zeros)\n",
            sum(d$n_dead_iom), sum(d$n_dead_iom > 0),
            100 * mean(d$n_dead_iom == 0)))
cat(sprintf("  Deaths UNITED: %.0f over %d death-days (%.1f%% zeros)\n",
            sum(d$n_dead_united), sum(d$n_dead_united > 0),
            100 * mean(d$n_dead_united == 0)))

# ── 1b. Exposure sensitivity: window grid ────────────────────
# Four SWH exposures corresponding to lagged rolling-mean windows:
# 1-day, 1-3 day, 1-5 day, 1-7 day. Primary analytical choice is
# swh_prev5days (1-5d mean); these alternatives probe robustness to the
# window length.
# Sample is pinned to the primary sample (non-NA swh_prev5days) so that
# coefficient differences across rows reflect the exposure variable,
# not the sample.
cat("\n--- 1b. Exposure sensitivity (window length) ---\n")

exposures <- c("swh_lag1", "swh_prev3days", "swh_prev5days", "swh_prevweek")

# Verify all exposures exist in the panel
missing_exp <- setdiff(exposures, names(d))
if (length(missing_exp) > 0) {
  stop("Missing exposure columns in panel: ",
       paste(missing_exp, collapse = ", "))
}

fit_sens <- function(x) {
  f <- as.formula(sprintf(
    "n_dead_iom ~ %s + %s:post_mou | month_year_fac", x, x))
  list(
    nb   = fenegbin(f, data = d, vcov = NW(14), panel.id = ~unit + date),
    pois = fepois (f, data = d, vcov = NW(14), panel.id = ~unit + date)
  )
}

sens <- setNames(lapply(exposures, fit_sens), exposures)

extract_sens <- function(x) {
  get_int <- function(m, fam) {
    ct <- coeftable(m, vcov = NW(14))
    r_main <- which(rownames(ct) == x)
    r_int  <- grep(":post_mou$", rownames(ct))
    r_int  <- r_int[grepl(x, rownames(ct)[r_int], fixed = TRUE)]
    tibble(
      exposure = x,
      family   = fam,
      b1       = ct[r_main, 1],
      se1      = ct[r_main, 2],
      b3       = ct[r_int,  1],
      se3      = ct[r_int,  2],
      p3       = 2 * pnorm(-abs(ct[r_int, 1] / ct[r_int, 2]))
    )
  }
  bind_rows(
    get_int(sens[[x]]$nb,   "NegBin"),
    get_int(sens[[x]]$pois, "Poisson")
  )
}

sens_tbl <- bind_rows(lapply(exposures, extract_sens)) %>%
  mutate(
    window = case_when(
      grepl("lag1$",     exposure) ~ "1d",
      grepl("prev3days", exposure) ~ "1-3d",
      grepl("prev5days", exposure) ~ "1-5d",
      grepl("prevweek",  exposure) ~ "1-7d"
    )
  ) %>%
  select(window, family, exposure, b1, se1, b3, se3, p3)

cat("\n  SWH × post_MoU interaction (b3) across windows:\n")
cat(sprintf("  %-6s %-8s  %+10s  %10s  %10s\n",
            "window", "family", "b3", "SE", "p"))
for (i in seq_len(nrow(sens_tbl))) {
  r <- sens_tbl[i, ]
  star <- if (r$p3 < 0.05) " *" else ""
  cat(sprintf("  %-6s %-8s  %+10.3f  %10.3f  %10.4f%s\n",
              r$window, r$family, r$b3, r$se3, r$p3, star))
}

# Save a CSV for downstream figure-making / paper tables
write.csv(sens_tbl,
          file.path(BASE_DIR, "output", "tables",
                    "20_exposure_sensitivity.csv"),
          row.names = FALSE)
cat(sprintf("  Saved: output/tables/20_exposure_sensitivity.csv\n"))

# ── 2. Primary model: NegBin and Poisson (count), IOM + UNITED ────
cat("\n--- 2. Primary model (count) ---\n")

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

# ── 2b. Rate robustness: Poisson with offset(log(crossing_attempts)) ──
# Outcome reframed as deaths per crossing-attempt. Count model remains
# primary; this is reported alongside for estimand robustness, since the
# guardrail argument is fundamentally about fatality RATE.
# Days with crossing_attempts == 0 are dropped (log(0) undefined).
d_rate <- d %>%
  filter(crossing_attempts > 0) %>%
  mutate(log_attempts = log(crossing_attempts))

m_rate   <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               offset(log_attempts) | month_year_fac,
  data = d_rate, vcov = NW(14), panel.id = ~unit + date)

m_rate_u <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  offset(log_attempts) | month_year_fac,
  data = d_rate, vcov = NW(14), panel.id = ~unit + date)

cat(sprintf("  Count sample: N = %d days (all)\n", nrow(d)))
cat(sprintf("  Rate sample:  N = %d days (drops %d zero-attempt days)\n",
            nrow(d_rate), nrow(d) - nrow(d_rate)))

cat("\n  Primary (no crossing control), NW(14) -- IOM vs UNITED side by side:\n")
print(etable(m_nb, m_pois, m_rate, m_nb_u, m_pois_u, m_rate_u,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM NB", "IOM Poiss", "IOM rate",
                         "UNITED NB", "UNITED Poiss", "UNITED rate")))

cat("\n  Primary, cluster(month_year):\n")
print(etable(m_nb, m_pois, m_rate, m_nb_u, m_pois_u, m_rate_u,
             vcov = ~month_year_fac, se.below = TRUE,
             headers = c("IOM NB", "IOM Poiss", "IOM rate",
                         "UNITED NB", "UNITED Poiss", "UNITED rate")))

# ── 3. Robustness: lagged crossing controls ──────────────────
cat("\n--- 3. Robustness: lagged crossing controls ---\n")

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

cat("\n  NegBin with crossing controls, NW(14):\n")
print(etable(m_nb, m_nb_lag7, m_nb_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))

cat("\n  Poisson with crossing controls, NW(14):\n")
print(etable(m_pois, m_pois_lag7, m_pois_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))

# ── 4. Save text output ──────────────────────────────────────
# Year-by-year SWH gradients were estimated and reported in earlier versions
# of this script. They were dropped because calendar year is not a unit of
# policy or operational variation (MoU 2017-07-01 cuts through 2017, Salvini
# decrees through 2018, Meloni through 2022). The yearly design fragments
# the identifying variation into ~365-day windows that don't map to any
# regime change, producing noisy estimates with no substantive interpretation.
# The pre/post-MoU interaction here, and the 4-period gradient in
# 31_united_periods.R, are the appropriate units of inference.
cat("\n--- 4. Saving results ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables", "20_primary_model.txt")
sink(sink_file)

cat("20  PRIMARY MODEL: SWH x POST-MOU -> DEATHS (IOM + UNITED)\n")
cat("==========================================================\n")
cat(sprintf("Sample: %s to %s (N = %d days)\n",
            min(d$date), max(d$date), nrow(d)))
cat(sprintf("Deaths IOM    (corridor, drowning + mixed):       %.0f\n",
            sum(d$n_dead_iom)))
cat(sprintf("Deaths UNITED (corridor, drowned + other/unknown): %.0f\n",
            sum(d$n_dead_united)))
cat("\nPrimary specification:\n")
cat("  deaths ~ swh_prev5days + swh_prev5days:post_mou | month_year\n")
cat("  No crossing control. Month-year FE absorbs monthly confounders.\n")
cat("  NW(14) SEs for serial correlation.\n")
cat("  Estimated separately on IOM MMP and UNITED death series, same\n")
cat("  spatial+cause filter, same sample, same spec.\n\n")

cat("=== PRIMARY MODEL (count + rate): IOM vs UNITED side by side ===\n")
cat("Count models: deaths per day (NegBin, Poisson).\n")
cat(sprintf(
  "Rate model:   deaths per crossing-attempt (Poisson + offset(log_attempts); N=%d days, drops %d zero-attempt days).\n",
  nrow(d_rate), nrow(d) - nrow(d_rate)))
cat("The count specification is primary; the rate spec is estimand robustness.\n\n")

cat("--- NW(14) SEs ---\n")
print(etable(m_nb, m_pois, m_rate, m_nb_u, m_pois_u, m_rate_u,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM NB", "IOM Poiss", "IOM rate",
                         "UNITED NB", "UNITED Poiss", "UNITED rate")))

cat("\n--- Cluster(month_year) SEs ---\n")
print(etable(m_nb, m_pois, m_rate, m_nb_u, m_pois_u, m_rate_u,
             vcov = ~month_year_fac, se.below = TRUE,
             headers = c("IOM NB", "IOM Poiss", "IOM rate",
                         "UNITED NB", "UNITED Poiss", "UNITED rate")))

cat("\n\n=== EXPOSURE SENSITIVITY (WINDOW LENGTH) ===\n")
cat("Four SWH exposures: 1-day, 1-3d, 1-5d, 1-7d lagged rolling means.\n")
cat("All fits on the primary sample (same as the primary NegBin/Poisson above).\n")
cat("Reporting b3 (SWH:post_mou) only; see CSV for full coefficients.\n\n")
cat(sprintf("  %-6s %-8s  %+10s  %10s  %10s\n",
            "window", "family", "b3", "SE", "p"))
for (i in seq_len(nrow(sens_tbl))) {
  r <- sens_tbl[i, ]
  star <- if (r$p3 < 0.05) " *" else ""
  cat(sprintf("  %-6s %-8s  %+10.3f  %10.3f  %10.4f%s\n",
              r$window, r$family, r$b3, r$se3, r$p3, star))
}
cat(sprintf("\n  Full table: output/tables/20_exposure_sensitivity.csv\n"))

cat("\n\n=== ROBUSTNESS: LAGGED CROSSING CONTROLS ===\n\n")
cat("NegBin:\n")
print(etable(m_nb, m_nb_lag7, m_nb_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))
cat("\nPoisson:\n")
print(etable(m_pois, m_pois_lag7, m_pois_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))

# Summary table
cat("\n\n=== SUMMARY: b3 (SWH x post_mou) ===\n\n")
for (info in list(
  list(m_nb,         "NegBin, no control"),
  list(m_nb_lag7,    "NegBin, lag 7d"),
  list(m_nb_lag14,   "NegBin, lag 14d"),
  list(m_pois,       "Poisson, no control"),
  list(m_pois_lag7,  "Poisson, lag 7d"),
  list(m_pois_lag14, "Poisson, lag 14d")
)) {
  ct <- coeftable(info[[1]], vcov = NW(14))
  r <- grep(":post_mou", rownames(ct))
  p <- 2 * pnorm(-abs(ct[r, 1] / ct[r, 2]))
  ct_cl <- coeftable(info[[1]], vcov = ~month_year_fac)
  p_cl <- 2 * pnorm(-abs(ct_cl[r, 1] / ct_cl[r, 2]))
  cat(sprintf("  %-25s  b3=%+.3f  SE_NW=%.3f  p_NW=%.4f  SE_cl=%.3f  p_cl=%.4f\n",
              info[[2]], ct[r, 1], ct[r, 2], p, ct_cl[r, 2], p_cl))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
