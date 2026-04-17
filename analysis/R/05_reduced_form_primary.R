# 05_reduced_form_primary.R
# =========================
# Primary reduced-form: n_dead_missing ~ SWH + SWH × post_mou | month_year FE
# NegBin, Newey-West(14) SEs. Dual window (full panel and 2014–2020 symmetric).
#
# Death series are built via build_iom_daily() (analysis/R/_helpers.R) — not
# taken from the panel column, which is built from the same helper. The
# primary analytical filter is:
#   - Route == "Central Mediterranean"
#   - Incident Type in {"incident", "split incident"}
#   - Country of Incident in the 5 CMR countries
#   - Inside the core corridor polygon (spatial join)
#   - Cause of death in {Drowning, Mixed or unknown}
#
# Cause-of-death restriction — rationale:
#   Most CMR deaths in IOM MMP are sea-related already (the route is a sea
#   crossing). The point of the cause filter is NOT to exclude "land" deaths
#   but to keep the cases that map most directly to the act of crossing the
#   sea — where SWH is the relevant exposure. Drowning is unambiguous;
#   "Mixed or unknown" almost certainly contains many drownings that were
#   never confirmed. The other categories (violence, vehicle accident,
#   sickness, harsh exposure) include events from before/after the maritime
#   leg of the journey (LCG interception shootings, detention, land transit)
#   for which SWH is not the relevant treatment.
#
#   Either choice introduces some measurement error: the cause filter drops a
#   small number of misclassified drownings, the no-filter version adds a
#   small number of irrelevant non-crossing deaths. Drowning + Mixed-or-unknown
#   are the overwhelming majority of CMR deaths (~96% of categorised CMR
#   incident deaths), so we adopt them as the primary outcome and report the
#   no-filter version (`All causes`) as a robustness check.
#
# Robustness variants vary ONE dimension each from the primary:
#   (a) All causes              — drop the cause filter
#   (b) All CMR countries       — drop the spatial filter (still drown+mixed)
#   (c) Outliers (>100) removed — drop high-leverage days
# Falsification:
#   (d) Next-7d SWH             — future weather should be null
#
# In:  analysis/data/daily_panel_complete.RDS   (SWH, post_mou, FE ids)
#      data/processed/iom_mmp_incidents.RDS
#      data/processed/core_corridor.RDS
# Out: output/tables/05_reduced_form_results.txt
#      output/figures/05_reduced_form_coefplot.png
#      output/figures/05_reduced_form_yearly_gradient.png

library(tidyverse)
library(lubridate)
library(fixest)
library(patchwork)

BASE_DIR      <- here::here()
MOU_DATE      <- as.Date("2017-07-01")
START_DATE    <- as.Date("2014-01-01")
SYMMETRIC_END <- as.Date("2020-12-31")

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("05  REDUCED-FORM PRIMARY: SWH x POST-MOU -> DEATHS\n")
cat("============================================================\n\n")

# ── 1. Daily panel (SWH and exogenous controls only) ───────
cat("--- 1. Loading daily panel ---\n")

# Drop the panel's n_dead_missing and rebuild via build_iom_daily() below
# so the analytical series always reflects the helper's current default
# filters (incident + split incident, drowning + mixed, core corridor).
panel <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel_complete.RDS")) %>%
  select(-n_dead_missing) %>%
  arrange(date) %>%
  mutate(swh_next7avg = rowMeans(sapply(1:7, \(k) dplyr::lead(swh, k))),
         year         = year(date),
         unit         = 1L)

PANEL_END <- max(panel$date)
cat(sprintf("  panel: %d days, %s to %s\n", nrow(panel), min(panel$date), PANEL_END))

# ── 2. Build death series via shared helper ────────────────
# Each call applies a different filter combination — change one argument to
# swap in a sensitivity variant. The primary uses helper defaults
# (incident-only, core corridor, sea causes).
cat("\n--- 2. Building death series ---\n")

daily_primary  <- build_iom_daily()                                # Central + sea (default)
daily_allcause <- build_iom_daily(causes  = "all") %>%
  rename(n_dead_missing_allcause = n_dead_missing)                  # Central + all causes
daily_allcmr   <- build_iom_daily(spatial = "all_cmr") %>%
  rename(n_dead_missing_allcmr = n_dead_missing)                    # all CMR  + sea

panel <- panel %>%
  left_join(daily_primary,  by = "date") %>%
  left_join(daily_allcause, by = "date") %>%
  left_join(daily_allcmr,   by = "date") %>%
  replace_na(list(n_dead_missing          = 0,
                  n_dead_missing_allcause  = 0,
                  n_dead_missing_allcmr    = 0))

cat(sprintf("  primary  (Central, incident, sea causes): %.0f deaths\n",
            sum(panel$n_dead_missing)))
cat(sprintf("  all causes (Central, drop cause filter):  %.0f deaths\n",
            sum(panel$n_dead_missing_allcause)))
cat(sprintf("  all CMR  (drop spatial filter, sea):      %.0f deaths\n",
            sum(panel$n_dead_missing_allcmr)))

# ── 3. Estimation per period ────────────────────────────────
cat("\n--- 3. Estimation ---\n")

extract_b3 <- function(model, spec, period) {
  ct <- coeftable(model, vcov = NW(14))
  r  <- grep(":post_mou$", rownames(ct))
  tibble(spec   = spec,
         period = period,
         coef   = ct[r, 1],
         se     = ct[r, 2],
         p      = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

fit_period <- function(end_date, period_label) {
  d <- panel %>%
    filter(between(date, START_DATE, end_date)) %>%
    mutate(month_year = factor(format(date, "%Y-%m")))

  cat(sprintf("\n=== %s | N = %d | pre = %d | post = %d ===\n",
              period_label, nrow(d),
              sum(d$post_mou == 0), sum(d$post_mou == 1)))

  fe <- function(formula, data = d)
    fenegbin(formula, data = data, vcov = NW(14), panel.id = ~unit + date)

  m_primary   <- fe(n_dead_missing          ~ swh_prevweek + swh_prevweek:post_mou | month_year)
  m_nooutlier <- fe(n_dead_missing          ~ swh_prevweek + swh_prevweek:post_mou | month_year,
                    data = d %>% filter(n_dead_missing <= 100))
  m_allcause  <- fe(n_dead_missing_allcause ~ swh_prevweek + swh_prevweek:post_mou | month_year)
  m_allcmr    <- fe(n_dead_missing_allcmr   ~ swh_prevweek + swh_prevweek:post_mou | month_year)
  m_next7d    <- fe(n_dead_missing          ~ swh_next7avg + swh_next7avg:post_mou | month_year)

  models <- list("Primary (drowning + mixed)"  = m_primary,
                 "No outliers >100 (robust.)"  = m_nooutlier,
                 "All causes (robust.)"        = m_allcause,
                 "All CMR countries (robust.)" = m_allcmr,
                 "Next-7d (falsification)"     = m_next7d)

  print(etable(m_primary, m_nooutlier, m_allcause, m_allcmr, m_next7d,
               headers = names(models), se.below = TRUE, vcov = NW(14)))

  imap_dfr(models, \(m, nm) extract_b3(m, nm, period_label))
}

all_results <- bind_rows(
  fit_period(PANEL_END,     sprintf("2014..%s", PANEL_END)),
  fit_period(SYMMETRIC_END, sprintf("2014..%s", SYMMETRIC_END))
)

# ── 4. Coefficient plot ─────────────────────────────────────
cat("\n--- 4. Coefficient plot ---\n")

SPEC_LEVELS <- c("Primary (drowning + mixed)",
                 "No outliers >100 (robust.)",
                 "All causes (robust.)",
                 "All CMR countries (robust.)",
                 "Next-7d (falsification)")

plot_df <- all_results %>%
  mutate(spec   = factor(spec, levels = rev(SPEC_LEVELS)),
         ci_lo  = coef - 1.96 * se,
         ci_hi  = coef + 1.96 * se)

period_colours <- setNames(c("#2166AC", "#B2182B"),
                            sort(unique(plot_df$period)))

p_coef <- plot_df %>%
  ggplot(aes(coef, spec, colour = period)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = period_colours) +
  labs(title    = expression(paste("Reduced-form: SWH x post-MoU interaction (",
                                    beta[3], ")")),
       subtitle = "NegBin, Newey-West(14) SEs, 95% CI",
       x = "Coefficient (per 1 meter SWH)", y = NULL, colour = "Sample") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(BASE_DIR, "output", "figures", "05_reduced_form_coefplot.png"),
       p_coef, width = 10, height = 6, dpi = 200)
cat("Saved: output/figures/05_reduced_form_coefplot.png\n")

# ── 5. Year-by-year gradient ────────────────────────────────
cat("\n--- 5. Year-by-year gradient ---\n")

yearly_plot <- function(end_date, period_label) {
  d <- panel %>%
    filter(between(date, START_DATE, end_date), !is.na(swh_prevweek)) %>%
    mutate(year_fac   = factor(year),
           month_year = factor(format(date, "%Y-%m")))

  m_yr <- fenegbin(n_dead_missing ~ swh_prevweek:year_fac | month_year,
                   data = d, vcov = NW(14), panel.id = ~unit + date)

  co <- coef(m_yr)
  V  <- vcov(m_yr, vcov = NW(14))[seq_along(co), seq_along(co)]
  yr <- tibble(year = parse_number(names(co)),
               beta = co,
               se   = sqrt(diag(V)),
               ci_lo = beta - 1.96 * se,
               ci_hi = beta + 1.96 * se)

  cat(sprintf("  %s\n", period_label))
  walk(seq_len(nrow(yr)), \(i) {
    r <- yr[i, ]
    cat(sprintf("    %d: %+.3f (SE=%.3f)%s\n",
                r$year, r$beta, r$se, if (abs(r$beta/r$se) > 1.96) " *" else ""))
  })

  ggplot(yr, aes(year, beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 2017.5, linetype = "dotted",
               colour = "#D32F2F", linewidth = 0.5) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "grey40") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.5) +
    annotate("text", x = 2017.7, y = max(yr$ci_hi) * 0.95,
             label = "MoU", colour = "#D32F2F", size = 3.5, hjust = 0) +
    scale_x_continuous(breaks = seq(2014, year(end_date))) +
    labs(title    = sprintf("SWH-mortality gradient by year (%s)", period_label),
         subtitle = "NegBin on raw prev-week SWH (metres) | Month-year FE | NW(14) SEs | 95% CI",
         x = NULL, y = expression(beta[SWH])) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

p_full <- yearly_plot(PANEL_END,     sprintf("2014..%s", PANEL_END))
p_sym  <- yearly_plot(SYMMETRIC_END, sprintf("2014..%s", SYMMETRIC_END))

ggsave(file.path(BASE_DIR, "output", "figures", "05_reduced_form_yearly_gradient.png"),
       p_full / p_sym, width = 10, height = 8, dpi = 200)
cat("Saved: output/figures/05_reduced_form_yearly_gradient.png\n")

# ── 5b. Time-trend tests for β(SWH) ────────────────────────
# Two formal tests of whether the SWH coefficient trends over time.
#
# (1) Continuous linear trend:
#       n_dead_missing ~ swh_prevweek + swh_prevweek:year_index | month_year
#     year_index = year - 2014 so that `swh_prevweek` is the 2014 slope.
#     The interaction tests H0: slope is constant over time. A positive
#     coefficient means weather becomes a stronger predictor each year.
#
# (2) Yearly heterogeneity (joint Wald):
#       n_dead_missing ~ swh_prevweek + swh_prevweek:year_fac | month_year
#     with year_fac relevelled to 2014 as reference, so each interaction
#     coefficient is (slope_year - slope_2014). A joint Wald test that all
#     interactions = 0 is the non-parametric version of (1); it detects any
#     time-varying pattern, not just a linear trend.
cat("\n--- 5b. Trend tests for β(SWH) ---\n")

trend_test <- function(data, sample_label) {
  d <- data %>%
    filter(between(date, START_DATE, PANEL_END), !is.na(swh_prevweek)) %>%
    mutate(year_index = year - 2014,
           year_fac   = relevel(factor(year), ref = "2014"),
           month_year = factor(format(date, "%Y-%m")))

  cat(sprintf("\n  %s (N=%d, deaths=%.0f)\n",
              sample_label, nrow(d), sum(d$n_dead_missing)))

  # (1) Continuous linear trend in β(SWH)
  m_trend <- fenegbin(
    n_dead_missing ~ swh_prevweek + swh_prevweek:year_index | month_year,
    data = d, vcov = NW(14), panel.id = ~ unit + date
  )
  ct <- coeftable(m_trend, vcov = NW(14))
  r  <- grep(":year_index$", rownames(ct))
  b_trend  <- ct[r, 1]
  se_trend <- ct[r, 2]
  p_trend  <- 2 * pnorm(-abs(b_trend / se_trend))
  cat(sprintf("    (1) continuous trend swh_prevweek:year_index = %+.4f (SE=%.4f, p=%.4f)\n",
              b_trend, se_trend, p_trend))

  # (2) Joint Wald test on yearly heterogeneity (manual chi-square using NW(14)
  # variance-covariance to match the SE structure used throughout).
  m_heter <- fenegbin(
    n_dead_missing ~ swh_prevweek + swh_prevweek:year_fac | month_year,
    data = d, vcov = NW(14), panel.id = ~ unit + date
  )
  coef_nms <- names(coef(m_heter))
  idx <- grep("swh_prevweek:year_fac", coef_nms)
  b   <- coef(m_heter)[idx]
  V   <- vcov(m_heter, vcov = NW(14))[idx, idx, drop = FALSE]
  wald_stat <- as.numeric(t(b) %*% solve(V) %*% b)
  wald_df   <- length(idx)
  wald_p    <- pchisq(wald_stat, df = wald_df, lower.tail = FALSE)
  cat(sprintf("    (2) joint Wald H0: slope_year = slope_2014 for all years\n"))
  cat(sprintf("        chi2(%d) = %.3f  p = %.4f\n", wald_df, wald_stat, wald_p))

  tibble(
    sample     = sample_label,
    trend_coef = b_trend, trend_se = se_trend, trend_p = p_trend,
    wald_chi2  = wald_stat, wald_df = wald_df, wald_p = wald_p
  )
}

# AFR-only subset: pool Libya + Tunisia deaths from the zone panel, keep the
# same corridor-wide SWH. The zone panel uses the same filter as the primary
# spec (incident + split incident, drowning + mixed, core-corridor-intersected
# SAR polygons; see 03_build_zone_panel.R).
zone_panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                                 "daily_panel_zone.RDS"))
afr_panel <- zone_panel %>%
  filter(sar_bloc == "AFR") %>%
  group_by(date) %>%
  summarise(
    n_dead_missing    = sum(n_dead_missing, na.rm = TRUE),
    swh_prevweek      = first(swh_prevweek),
    crossing_attempts = first(crossing_attempts),
    .groups = "drop"
  ) %>%
  mutate(year = year(date), unit = 1L) %>%
  arrange(date)

trend_results <- bind_rows(
  trend_test(panel,     "Full CMR corridor"),
  trend_test(afr_panel, "African SAR zone (Libya+Tunisia)")
)

# ── 5c. Per-attempt trend test (offset model) ────────────
# Same tests as 5b but with offset(log(crossing_attempts)), so the
# coefficient is on the fatality RATE per attempt, not the raw count.
# Days with zero crossing attempts are dropped (log(0) undefined).
cat("\n--- 5c. Per-attempt trend tests (offset model) ---\n")

trend_test_rate <- function(data, sample_label) {
  d <- data %>%
    filter(between(date, START_DATE, PANEL_END),
           !is.na(swh_prevweek),
           crossing_attempts > 0) %>%
    mutate(year_index  = year - 2014,
           year_fac    = relevel(factor(year), ref = "2014"),
           month_year  = factor(format(date, "%Y-%m")),
           log_attempts = log(crossing_attempts))

  cat(sprintf("\n  %s (N=%d, deaths=%.0f, attempts=%.0f)\n",
              sample_label, nrow(d),
              sum(d$n_dead_missing), sum(d$crossing_attempts)))

  # (1) Continuous linear trend in per-attempt β(SWH)
  m_trend <- fenegbin(
    n_dead_missing ~ swh_prevweek + swh_prevweek:year_index +
                     offset(log_attempts) | month_year,
    data = d, vcov = NW(14), panel.id = ~ unit + date
  )
  ct <- coeftable(m_trend, vcov = NW(14))
  r  <- grep(":year_index$", rownames(ct))
  b_trend  <- ct[r, 1]
  se_trend <- ct[r, 2]
  p_trend  <- 2 * pnorm(-abs(b_trend / se_trend))
  cat(sprintf("    (1) continuous trend swh:year_index = %+.4f (SE=%.4f, p=%.4f)\n",
              b_trend, se_trend, p_trend))

  # (2) Joint Wald on yearly heterogeneity
  m_heter <- fenegbin(
    n_dead_missing ~ swh_prevweek + swh_prevweek:year_fac +
                     offset(log_attempts) | month_year,
    data = d, vcov = NW(14), panel.id = ~ unit + date
  )
  idx <- grep("swh_prevweek:year_fac", names(coef(m_heter)))
  b   <- coef(m_heter)[idx]
  V   <- vcov(m_heter, vcov = NW(14))[idx, idx, drop = FALSE]
  wald_stat <- as.numeric(t(b) %*% solve(V) %*% b)
  wald_df   <- length(idx)
  wald_p    <- pchisq(wald_stat, df = wald_df, lower.tail = FALSE)
  cat(sprintf("    (2) joint Wald chi2(%d) = %.3f  p = %.4f\n",
              wald_df, wald_stat, wald_p))

  tibble(
    sample     = sample_label,
    trend_coef = b_trend, trend_se = se_trend, trend_p = p_trend,
    wald_chi2  = wald_stat, wald_df = wald_df, wald_p = wald_p
  )
}

rate_results <- bind_rows(
  trend_test_rate(panel,     "Full CMR corridor (per-attempt)"),
  trend_test_rate(afr_panel, "African SAR (per-attempt)")
)

# ── 6. Results table ────────────────────────────────────────
cat("\n--- 6. Saving results table ---\n")

sink(file.path(BASE_DIR, "output", "tables", "05_reduced_form_results.txt"))
cat("REDUCED-FORM PRIMARY: n_dead_missing ~ SWH x post_mou | month-year FE\n")
cat("NegBin (fenegbin) | Newey-West(14) SEs | Raw prev-week SWH (per 1 metre)\n")
cat("Primary outcome: incident-only, core corridor polygon,\n")
cat("                 cause = Drowning or Mixed/unknown (cases most directly\n")
cat("                 tied to the act of crossing the sea).\n")
cat("Sample window: 2014-01-01 to ", as.character(PANEL_END), "\n\n", sep = "")

all_results %>%
  group_split(period) %>%
  walk(\(x) {
    cat(sprintf("=== %s ===\n", unique(x$period)))
    walk(seq_len(nrow(x)), \(i) {
      r <- x[i, ]
      cat(sprintf("  %-32s  b3 = %+.3f (SE=%.3f)  IRR = %.3f  p = %.4f\n",
                  r$spec, r$coef, r$se, exp(r$coef), r$p))
    })
    cat("\n")
  })

cat("TREND TESTS FOR β(SWH)\n")
cat("----------------------\n")
cat("(1) Continuous linear trend: slope of swh_prevweek:year_index\n")
cat("(2) Joint Wald on yearly heterogeneity (all (slope_year - slope_2014) = 0)\n")
cat("Sample window: 2014-01-01 to ", as.character(PANEL_END), "\n\n", sep = "")
walk(seq_len(nrow(trend_results)), \(i) {
  r <- trend_results[i, ]
  cat(sprintf("=== %s ===\n", r$sample))
  cat(sprintf("  (1) continuous  : %+.4f (SE=%.4f)  p = %.4f\n",
              r$trend_coef, r$trend_se, r$trend_p))
  cat(sprintf("  (2) joint Wald  : chi2(%d) = %.3f  p = %.4f\n\n",
              r$wald_df, r$wald_chi2, r$wald_p))
})

cat("TREND TESTS FOR β(SWH) — PER-ATTEMPT (offset model)\n")
cat("----------------------------------------------------\n")
cat("Same as above but with offset(log(crossing_attempts)).\n")
cat("Coefficients are on the fatality RATE per attempt, not raw counts.\n")
cat("Days with zero crossing attempts dropped.\n\n")
walk(seq_len(nrow(rate_results)), \(i) {
  r <- rate_results[i, ]
  cat(sprintf("=== %s ===\n", r$sample))
  cat(sprintf("  (1) continuous  : %+.4f (SE=%.4f)  p = %.4f\n",
              r$trend_coef, r$trend_se, r$trend_p))
  cat(sprintf("  (2) joint Wald  : chi2(%d) = %.3f  p = %.4f\n\n",
              r$wald_df, r$wald_chi2, r$wald_p))
})
sink()
cat("Saved: output/tables/05_reduced_form_results.txt\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
