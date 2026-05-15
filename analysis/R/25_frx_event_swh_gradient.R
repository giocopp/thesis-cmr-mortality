# 25_frx_event_swh_gradient.R
# ============================
# Frontex event-count gradient: how does the SWH-event slope shift around the
# MoU, separately for SAR-flagged and Not-SAR-flagged engagements?
#
# Specification
# -------------
#   frx_n_sar    ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac
#   frx_n_notsar ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac
# (plus stacked counterpart with SAR-flag interaction; year-by-year gradients)
#
# Why this exists
# ---------------
# The primary model (20_primary_model.R) shows the SWH-deaths gradient steepens
# at the MoU. The symmetric question on the Frontex side is whether the
# rescue side responds the way the guardrail argument predicts.
#
# Persons-per-event (intensive margin) does not move at the MoU: once Frontex
# engages a boat, crew sizes are roughly constant in SWH before and after.
# The retreat shows up in the COUNT of SAR engagements: how often Frontex
# engages in rough seas. That is the symmetric counterpart to deaths,
# because IOM deaths are observed unconditionally on engagement, while
# frx_persons is observed only conditional on engagement.
#
# Interpretation
# --------------
# A negative swh:post_mou on the SAR-event count = SAR engagements become
# less responsive to (or actively retreat from) rough seas after the MoU.
# Mirror coefficient on the Not-SAR count provides a placebo: if Not-SAR
# engagements show no comparable shift, the MoU break is concentrated in
# the SAR channel.
#
# Both NegBin (fenegbin) and Poisson QMLE (fepois), NW(14) SEs. Matches
# 20_primary_model.R modeling choices (daily aggregation, month-year FE,
# NW(14), dual family).
#
# In:  analysis/data/daily_panel_complete.RDS
# Out: output/tables/25_frx_event_swh_gradient.txt

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR   <- here::here()
MOU_DATE   <- as.Date("2017-07-01")

cat("============================================================\n")
cat("25  FRONTEX EVENT-COUNT SWH GRADIENT (SAR vs Not-SAR)\n")
cat("============================================================\n\n")

# ── 1. Load panel and build event counts ─────────────────────
cat("--- 1. Loading data ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

# SAR / Not-SAR daily event counts.
#   frx_n_sar      already exists (sum of sar_ops==TRUE per day).
#   frx_n_notsar   = frx_incidents - frx_n_sar
#                  (so SAR=NA events go to Not-SAR, consistent with
#                  01_build_daily_panel.R's sar_bucket definition).
panel <- panel %>%
  mutate(
    frx_n_notsar   = frx_incidents - frx_n_sar,
    unit           = 1L,
    year_fac       = factor(year),
    month_year_fac = factor(month_year)
  )

d <- panel %>% filter(!is.na(swh_prev5days))

cat(sprintf("  Panel: %s to %s (%d days)\n",
            min(d$date), max(d$date), nrow(d)))
cat(sprintf("  SAR events:    %d total over %d active days\n",
            sum(d$frx_n_sar), sum(d$frx_n_sar > 0)))
cat(sprintf("  Not-SAR events: %d total over %d active days\n",
            sum(d$frx_n_notsar), sum(d$frx_n_notsar > 0)))

# ── 2. Per-outcome gradient models (mirrors 20_primary_model.R) ──
cat("\n--- 2. Per-outcome SWH × post_MoU models ---\n")

# SAR-event count
m_nb_sar    <- fenegbin(frx_n_sar    ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac,
                        data = d, vcov = NW(14), panel.id = ~unit + date)
m_pois_sar  <- fepois  (frx_n_sar    ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac,
                        data = d, vcov = NW(14), panel.id = ~unit + date)

# Not-SAR-event count (placebo channel)
m_nb_not    <- fenegbin(frx_n_notsar ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac,
                        data = d, vcov = NW(14), panel.id = ~unit + date)
m_pois_not  <- fepois  (frx_n_notsar ~ swh_prev5days + swh_prev5days:post_mou |
                          month_year_fac,
                        data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  SAR & Not-SAR event counts, NW(14):\n")
print(etable(m_nb_sar, m_pois_sar, m_nb_not, m_pois_not,
             vcov = NW(14), se.below = TRUE,
             headers = c("SAR NB", "SAR Poiss", "NotSAR NB", "NotSAR Poiss")))

cat("\n  Same, cluster(month_year):\n")
print(etable(m_nb_sar, m_pois_sar, m_nb_not, m_pois_not,
             vcov = ~month_year_fac, se.below = TRUE,
             headers = c("SAR NB", "SAR Poiss", "NotSAR NB", "NotSAR Poiss")))

# ── 3. Stacked: SWH × SAR-flag interaction ───────────────────
# Long-format: each day contributes two rows (one for SAR count, one for
# Not-SAR count). Lets us test the differential gradient directly.
cat("\n--- 3. Stacked SAR-flag interaction ---\n")

d_long <- bind_rows(
  d %>% transmute(date, month_year_fac, post_mou, swh_prev5days, year_fac,
                  n_events = frx_n_sar,    sar = 1L),
  d %>% transmute(date, month_year_fac, post_mou, swh_prev5days, year_fac,
                  n_events = frx_n_notsar, sar = 0L)
)

m_nb_stack   <- fenegbin(n_events ~ swh_prev5days + swh_prev5days:sar + sar |
                           month_year_fac,
                         data = d_long, vcov = NW(14),
                         panel.id = ~sar + date)
m_pois_stack <- fepois  (n_events ~ swh_prev5days + swh_prev5days:sar + sar |
                           month_year_fac,
                         data = d_long, vcov = NW(14),
                         panel.id = ~sar + date)

cat("\n  Stacked, NW(14):\n")
print(etable(m_nb_stack, m_pois_stack, vcov = NW(14), se.below = TRUE,
             headers = c("NegBin", "Poisson")))

# ── 4. Save text output ──────────────────────────────────────
# Year-by-year SWH gradients were dropped: calendar year is not a unit of
# policy or operational variation, so per-year slopes confound regime
# changes (MoU, Salvini, Meloni) with noise. The pre/post-MoU interaction
# (Section 2) is the appropriate unit of inference for this script.
cat("\n--- 4. Saving results ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables",
                        "25_frx_event_swh_gradient.txt")
sink(sink_file)

cat("25  FRONTEX EVENT-COUNT SWH GRADIENT (SAR vs Not-SAR)\n")
cat("======================================================\n")
cat(sprintf("Sample: %s to %s (N = %d days)\n",
            min(d$date), max(d$date), nrow(d)))
cat(sprintf("SAR events:     %d (active on %d days)\n",
            sum(d$frx_n_sar), sum(d$frx_n_sar > 0)))
cat(sprintf("Not-SAR events: %d (active on %d days)\n",
            sum(d$frx_n_notsar), sum(d$frx_n_notsar > 0)))

cat("\nSpecifications:\n")
cat("  frx_n_sar    ~ swh_prev5days + swh_prev5days:post_mou | month_year\n")
cat("  frx_n_notsar ~ swh_prev5days + swh_prev5days:post_mou | month_year\n")
cat("  Daily aggregation, month-year FE, NW(14) SEs.\n")
cat("  Mirrors 20_primary_model.R applied to deaths.\n")
cat("  Not-SAR series serves as placebo channel.\n\n")

cat("=== SWH x post_MoU on event counts ===\n")
cat("--- NW(14) SEs ---\n")
print(etable(m_nb_sar, m_pois_sar, m_nb_not, m_pois_not,
             vcov = NW(14), se.below = TRUE,
             headers = c("SAR NB", "SAR Poiss", "NotSAR NB", "NotSAR Poiss")))

cat("\n--- Cluster(month_year) SEs ---\n")
print(etable(m_nb_sar, m_pois_sar, m_nb_not, m_pois_not,
             vcov = ~month_year_fac, se.below = TRUE,
             headers = c("SAR NB", "SAR Poiss", "NotSAR NB", "NotSAR Poiss")))

cat("\n=== Stacked: SWH x SAR-flag interaction ===\n")
print(etable(m_nb_stack, m_pois_stack, vcov = NW(14), se.below = TRUE,
             headers = c("NegBin", "Poisson")))

cat("\n=== SUMMARY: SWH x post_MoU coefficient (b3) ===\n\n")
for (info in list(
  list(m_nb_sar,   "SAR events, NegBin"),
  list(m_pois_sar, "SAR events, Poisson"),
  list(m_nb_not,   "Not-SAR events, NegBin"),
  list(m_pois_not, "Not-SAR events, Poisson")
)) {
  ct    <- coeftable(info[[1]], vcov = NW(14))
  r     <- grep(":post_mou", rownames(ct))
  p     <- 2 * pnorm(-abs(ct[r, 1] / ct[r, 2]))
  ct_cl <- coeftable(info[[1]], vcov = ~month_year_fac)
  p_cl  <- 2 * pnorm(-abs(ct_cl[r, 1] / ct_cl[r, 2]))
  cat(sprintf("  %-26s b3=%+.3f  SE_NW=%.3f  p_NW=%.4f  SE_cl=%.3f  p_cl=%.4f\n",
              info[[2]], ct[r, 1], ct[r, 2], p, ct_cl[r, 2], p_cl))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
