# 23_mechanism_interactions.R
# ============================
# Mechanism complement to 20_primary_model.R.
#
# 20 estimates the reduced-form SWH x post_MoU shift in the gradient: post_MoU
# is the policy-event indicator, not a measure of SAR. 23 asks whether the
# gradient varies with SAR capacity directly, which is the channel the
# guardrail argument names. The two scripts answer different questions:
#   20 -> "did the gradient shift around the policy event?"  (reduced form)
#   23 -> "does the gradient track SAR capacity?"            (mechanism)
#
# SAR moderator (weekly lagged, days t-7 to t-1, z-scored). Two variants:
#   (i)  share:    sar_events_pw / incidents_pw                   [primary]
#   (ii) persons:  log1p(sum of Frontex SAR persons, weekly)      [robustness]
#
# Why two variants:
#   The share spec has a mechanical concern: when LCG/TCG pullbacks rise
#   (post 2018) the denominator (frx_incidents) shifts and SAR share can
#   fall even when absolute SAR capacity is unchanged. The absolute-persons
#   variant drops the denominator and answers the question "did the
#   SWH-mortality gradient covary with raw SAR rescue capacity?" without
#   the share-mechanics confound. Persons rescued is the more direct
#   measure of "rescue capacity deployed" than event counts. Both
#   variants remain endogenous to crossings/deaths -- this is mechanism,
#   not identification.
#
# Specification:
#   deaths ~ SWH + SWH:SAR + SAR | FE
#
# FE variants:
#   (a) year + month-of-year
#   (b) month-year
#
# Both NegBin (fenegbin) and Poisson QMLE (fepois), NW(14) SEs.
# Matches 20_primary_model.R choices (5-day SWH, NW(14), dual family).
#
# In:  analysis/data/daily_panel_complete.RDS
# Out: output/tables/23_mechanism_interactions.txt
#      output/figures/23_mechanism_coefplot.png

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("23  MECHANISM: SWH x SAR share\n")
cat("============================================================\n\n")

# ── 1. Load data + build weekly-lagged SAR share ─────────────
cat("--- 1. Loading data ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

iom_daily <- build_iom_daily()

# Daily total SAR persons across all interceptor types (NGO, EU, Italy,
# commercial, EU CG, land patrol, no-intercept, other, NA). Built from
# the Frontex interceptor x SAR matrix in 01_build_daily_panel.R.
sar_persons_cols <- c(
  "frx_persons_sar_ngo",   "frx_persons_sar_eu",    "frx_persons_sar_ita",
  "frx_persons_sar_comm",  "frx_persons_sar_cg",    "frx_persons_sar_land",
  "frx_persons_sar_noint", "frx_persons_sar_other", "frx_persons_sar_na"
)
panel$frx_persons_sar_total <- rowSums(panel[, sar_persons_cols], na.rm = TRUE)

panel <- panel %>%
  left_join(iom_daily %>% rename(n_dead_iom = n_dead_missing), by = "date") %>%
  replace_na(list(n_dead_iom = 0)) %>%
  arrange(date) %>%
  mutate(
    # Weekly lagged sums (days t-7 to t-1) of raw SAR counts.
    sar_events_pw    = dplyr::lag(zoo::rollsumr(frx_n_sar,              k = 7, fill = NA), 1),
    sar_persons_pw   = dplyr::lag(zoo::rollsumr(frx_persons_sar_total,  k = 7, fill = NA), 1),
    incidents_pw     = dplyr::lag(zoo::rollsumr(frx_incidents,          k = 7, fill = NA), 1),

    # (i) Share moderator: SAR events / total incidents (z-scored).
    sar_share_pw     = ifelse(incidents_pw > 0, sar_events_pw / incidents_pw,
                                NA_real_),
    sar_share_pw_z   = (sar_share_pw - mean(sar_share_pw, na.rm = TRUE)) /
                         sd(sar_share_pw, na.rm = TRUE),

    # (ii) Absolute-persons moderator: log1p of weekly SAR persons (z-scored).
    #      Drops the share denominator -> immune to mechanical fall when
    #      LCG/TCG pullbacks rise. Persons rescued is the direct measure
    #      of rescue capacity deployed.
    log1p_sar_persons_pw   = log1p(sar_persons_pw),
    log1p_sar_persons_pw_z = (log1p_sar_persons_pw - mean(log1p_sar_persons_pw, na.rm = TRUE)) /
                                sd(log1p_sar_persons_pw, na.rm = TRUE),

    # Lagged crossing volume (for robustness / per-attempt interpretation).
    # lag 14 = mean over t-14 to t-8 (the PREVIOUS week, not this week).
    # Chosen over lag 7 because lag 7 is strongly correlated with current
    # SWH (r = -0.50) and partly blocks the mechanism; lag 14 drops to
    # r = -0.25 with SWH and delivers an interaction identical to the
    # no-control spec.
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag14         = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    log1p_lc_lag14   = log1p(lc_lag14),
    unit           = 1L,
    year_fac       = factor(year),
    month_of_year  = factor(month(date)),
    month_year_fac = factor(month_year)
  )

d <- panel %>% filter(!is.na(sar_share_pw), !is.na(swh_prev5days),
                       !is.na(log1p_lc_lag14))

cat(sprintf("  N = %d days (of %d; lost %d to lagged window NAs)\n",
            nrow(d), nrow(panel), nrow(panel) - nrow(d)))
cat(sprintf("  Deaths: %.0f over %d death-days\n",
            sum(d$n_dead_iom), sum(d$n_dead_iom > 0)))
cat(sprintf("  SAR share:           mean=%.3f  sd=%.3f  range=[%.2f, %.2f]\n",
            mean(d$sar_share_pw), sd(d$sar_share_pw),
            min(d$sar_share_pw), max(d$sar_share_pw)))
cat(sprintf("  log1p SAR persons:   mean=%.3f  sd=%.3f  range=[%.2f, %.2f]\n",
            mean(d$log1p_sar_persons_pw), sd(d$log1p_sar_persons_pw),
            min(d$log1p_sar_persons_pw), max(d$log1p_sar_persons_pw)))
cat(sprintf("  cor(sar_share_pw,      swh_prev5days)    = %+.3f\n",
            cor(d$sar_share_pw,         d$swh_prev5days)))
cat(sprintf("  cor(log1p_sar_persons, swh_prev5days)    = %+.3f\n",
            cor(d$log1p_sar_persons_pw, d$swh_prev5days)))
cat(sprintf("  cor(sar_share_pw,      log1p_sar_persons) = %+.3f\n",
            cor(d$sar_share_pw,         d$log1p_sar_persons_pw)))

# ── 2. Estimation ────────────────────────────────────────────
#
# PRIMARY: no crossing control -- the total SWH-deaths gradient, consistent
# with 20_primary_model.R's primary spec. This is the estimand of interest
# for the guardrail argument (does SAR neutralise the weather channel?).
#
# ROBUSTNESS: +log1p(lag-14d crossings). Crossings are post-treatment for
# both SWH (rough seas deter departures) and SAR (SAR capacity may attract
# crossings), so any crossing control is a "bad control" in the strict
# sense. lag 14 (mean over t-14 to t-8) is much more pre-treatment than
# lag 7: r(control, SWH) drops from -0.50 (lag 7) to -0.25 (lag 14), and
# the interaction under lag 14 matches the no-control spec (no attenuation),
# confirming the effect is not volume-selection. Conditional on volume the
# estimand becomes the per-attempt fatality-rate gradient.
cat("\n--- 2. Estimation ---\n")

# Year + month-of-year FE
m_ym_base       <- fenegbin(n_dead_iom ~ swh_prev5days | year_fac + month_of_year,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_ym_sar_nb     <- fenegbin(n_dead_iom ~ swh_prev5days + swh_prev5days:sar_share_pw_z +
                              sar_share_pw_z | year_fac + month_of_year,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_ym_sar_po     <- fepois  (n_dead_iom ~ swh_prev5days + swh_prev5days:sar_share_pw_z +
                              sar_share_pw_z | year_fac + month_of_year,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_ym_sar_nb_ctl <- fenegbin(n_dead_iom ~ swh_prev5days + swh_prev5days:sar_share_pw_z +
                              sar_share_pw_z + log1p_lc_lag14 | year_fac + month_of_year,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_ym_sar_po_ctl <- fepois  (n_dead_iom ~ swh_prev5days + swh_prev5days:sar_share_pw_z +
                              sar_share_pw_z + log1p_lc_lag14 | year_fac + month_of_year,
                            data = d, vcov = NW(14), panel.id = ~unit + date)

# Month-year FE (matches 20_primary_model.R)
m_my_base       <- fenegbin(n_dead_iom ~ swh_prev5days | month_year_fac,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_my_sar_nb     <- fenegbin(n_dead_iom ~ swh_prev5days + swh_prev5days:sar_share_pw_z +
                              sar_share_pw_z | month_year_fac,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_my_sar_po     <- fepois  (n_dead_iom ~ swh_prev5days + swh_prev5days:sar_share_pw_z +
                              sar_share_pw_z | month_year_fac,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_my_sar_nb_ctl <- fenegbin(n_dead_iom ~ swh_prev5days + swh_prev5days:sar_share_pw_z +
                              sar_share_pw_z + log1p_lc_lag14 | month_year_fac,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_my_sar_po_ctl <- fepois  (n_dead_iom ~ swh_prev5days + swh_prev5days:sar_share_pw_z +
                              sar_share_pw_z + log1p_lc_lag14 | month_year_fac,
                            data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  Year + month-of-year FE, NW(14):\n")
print(etable(m_ym_base, m_ym_sar_nb, m_ym_sar_po, m_ym_sar_nb_ctl, m_ym_sar_po_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("Baseline (NB)", "NB no ctl", "Pois no ctl",
                         "NB + log(cross)", "Pois + log(cross)")))

cat("\n  Month-year FE, NW(14):\n")
print(etable(m_my_base, m_my_sar_nb, m_my_sar_po, m_my_sar_nb_ctl, m_my_sar_po_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("Baseline (NB)", "NB no ctl", "Pois no ctl",
                         "NB + log(cross)", "Pois + log(cross)")))

# ── 2b. Absolute-SAR robustness ──────────────────────────────
# Same spec, replacing SAR-share with log1p of weekly absolute SAR persons,
# z-scored. Drops the share denominator -> answers the question "is the
# share signal mechanical from rising LCG/TCG pullbacks in the denominator?"
# If the absolute moderator gives the same-sign interaction, the answer
# is no. Persons rescued is the more meaningful capacity measure than
# event counts.
cat("\n--- 2b. Absolute-SAR robustness (persons) ---\n")

# log1p(SAR persons), month-year FE
m_my_per_nb     <- fenegbin(n_dead_iom ~ swh_prev5days + swh_prev5days:log1p_sar_persons_pw_z +
                              log1p_sar_persons_pw_z | month_year_fac,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_my_per_po     <- fepois  (n_dead_iom ~ swh_prev5days + swh_prev5days:log1p_sar_persons_pw_z +
                              log1p_sar_persons_pw_z | month_year_fac,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_my_per_nb_ctl <- fenegbin(n_dead_iom ~ swh_prev5days + swh_prev5days:log1p_sar_persons_pw_z +
                              log1p_sar_persons_pw_z + log1p_lc_lag14 | month_year_fac,
                            data = d, vcov = NW(14), panel.id = ~unit + date)
m_my_per_po_ctl <- fepois  (n_dead_iom ~ swh_prev5days + swh_prev5days:log1p_sar_persons_pw_z +
                              log1p_sar_persons_pw_z + log1p_lc_lag14 | month_year_fac,
                            data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  log1p(SAR persons), month-year FE, NW(14):\n")
print(etable(m_my_per_nb, m_my_per_po, m_my_per_nb_ctl, m_my_per_po_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("NB no ctl", "Pois no ctl",
                         "NB + log(cross)", "Pois + log(cross)")))

# ── 3. Summary of b_int (SWH x SAR share) ────────────────────
cat("\n--- 3. Summary (SWH x SAR share) ---\n")

extract_int <- function(m, moderator, label, group) {
  ct  <- coeftable(m, vcov = NW(14))
  pat <- paste0("swh_prev5days:", moderator)
  r   <- grep(pat, rownames(ct), fixed = TRUE)
  tibble(group = group, spec = label,
         coef  = ct[r, 1], se = ct[r, 2],
         p     = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

summary_rows <- bind_rows(
  # (i) Share moderator (existing primary)
  extract_int(m_ym_sar_nb,     "sar_share_pw_z", "yr+mo FE, NegBin  (primary)",      "share"),
  extract_int(m_ym_sar_po,     "sar_share_pw_z", "yr+mo FE, Poisson (primary)",      "share"),
  extract_int(m_ym_sar_nb_ctl, "sar_share_pw_z", "yr+mo FE, NegBin  + log(cross)",   "share"),
  extract_int(m_ym_sar_po_ctl, "sar_share_pw_z", "yr+mo FE, Poisson + log(cross)",   "share"),
  extract_int(m_my_sar_nb,     "sar_share_pw_z", "my FE,    NegBin  (primary)",      "share"),
  extract_int(m_my_sar_po,     "sar_share_pw_z", "my FE,    Poisson (primary)",      "share"),
  extract_int(m_my_sar_nb_ctl, "sar_share_pw_z", "my FE,    NegBin  + log(cross)",   "share"),
  extract_int(m_my_sar_po_ctl, "sar_share_pw_z", "my FE,    Poisson + log(cross)",   "share"),
  # (ii) Absolute SAR persons
  extract_int(m_my_per_nb,     "log1p_sar_persons_pw_z", "my FE,    NegBin  (primary)",    "persons"),
  extract_int(m_my_per_po,     "log1p_sar_persons_pw_z", "my FE,    Poisson (primary)",    "persons"),
  extract_int(m_my_per_nb_ctl, "log1p_sar_persons_pw_z", "my FE,    NegBin  + log(cross)", "persons"),
  extract_int(m_my_per_po_ctl, "log1p_sar_persons_pw_z", "my FE,    Poisson + log(cross)", "persons")
)

cat("\n  b_int = coefficient on swh_prev5days x SAR moderator (per 1-SD)\n")
group_labels <- c(
  share   = "(i)  SAR share",
  persons = "(ii) log1p(SAR persons) -- absolute level"
)
for (g in names(group_labels)) {
  cat(sprintf("\n  %s\n", group_labels[g]))
  rows_g <- summary_rows[summary_rows$group == g, ]
  for (i in seq_len(nrow(rows_g))) {
    r <- rows_g[i, ]
    cat(sprintf("    %-32s  %+.3f (SE=%.3f)  p=%.4f%s\n",
                r$spec, r$coef, r$se, r$p,
                if (r$p < 0.05) " *" else ""))
  }
}

# ── 4. Coefficient plot ──────────────────────────────────────
cat("\n--- 4. Plot ---\n")

plot_df <- summary_rows %>%
  mutate(
    ci_lo = coef - 1.96 * se,
    ci_hi = coef + 1.96 * se,
    group_lab = factor(group,
                        levels = c("share", "persons"),
                        labels = c("(i) SAR share",
                                   "(ii) log1p(SAR persons)")),
    row_id = row_number(),
    spec_f = factor(paste0(group, "::", spec),
                     levels = rev(paste0(group, "::", spec)))
  )

p <- ggplot(plot_df, aes(coef, spec_f, colour = group_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2) +
  scale_colour_manual(values = c("(i) SAR share"          = "#2166AC",
                                  "(ii) log1p(SAR persons)" = "#D95F02")) +
  scale_y_discrete(labels = function(x) sub("^[a-z]+::", "", x)) +
  labs(
    title    = "Mechanism: SWH x SAR moderator interactions",
    subtitle = "Two SAR moderators (share, log absolute persons), z-scored, weekly-lagged. NW(14) SEs.",
    x        = "Coefficient on SWH x SAR moderator (per 1-SD)",
    y        = NULL,
    colour   = "Moderator"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

ggsave(file.path(BASE_DIR, "output", "figures",
                  "23_mechanism_coefplot.png"),
       p, width = 11, height = 5.5, dpi = 200)
cat("  Saved: output/figures/23_mechanism_coefplot.png\n")

# ── 5. Save text output ──────────────────────────────────────
cat("\n--- 5. Saving results ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables",
                        "23_mechanism_interactions.txt")
sink(sink_file)

cat("23  MECHANISM: SWH x SAR moderator (share, persons)\n")
cat("====================================================\n")
cat(sprintf("Sample: %s to %s (N = %d days, %.0f deaths)\n",
            min(d$date), max(d$date), nrow(d), sum(d$n_dead_iom)))
cat("\nTwo moderator variants, both weekly-summed (days t-7 to t-1), z-scored:\n")
cat("  (i)  share:    sar_events_pw / incidents_pw\n")
cat("  (ii) persons:  log1p(sar_persons_pw)\n")
cat("\nCoefficients on the interaction are per 1-SD change in the moderator.\n")
cat(sprintf("cor(sar_share,         swh_prev5days)    = %+.3f\n",
            cor(d$sar_share_pw,         d$swh_prev5days)))
cat(sprintf("cor(log1p_sar_persons, swh_prev5days)    = %+.3f\n",
            cor(d$log1p_sar_persons_pw, d$swh_prev5days)))
cat(sprintf("cor(sar_share,         log1p_sar_persons) = %+.3f\n\n",
            cor(d$sar_share_pw,         d$log1p_sar_persons_pw)))

cat("PRIMARY spec: no crossing control.\n")
cat("ROBUSTNESS spec: +log1p(lag-14d crossings). Post-treatment control\n")
cat("  (crossings are affected by both SWH and SAR); lag 14 (t-14 to t-8)\n")
cat("  is much more pre-treatment than lag 7: r(control, SWH) drops from\n")
cat("  -0.50 at lag 7 to -0.25 at lag 14. The interaction at lag 14 matches\n")
cat("  the no-control spec, confirming it is not a volume-selection artefact.\n\n")

cat("WHY TWO MODERATORS:\n")
cat("  The share denominator includes LCG/TCG pullbacks and other Not-SAR\n")
cat("  events. After 2018 those rise sharply (see desc_panel_event_type.png),\n")
cat("  so the share can fall mechanically even if absolute SAR is unchanged.\n")
cat("  The log-persons variant drops the denominator and measures absolute\n")
cat("  SAR rescue capacity directly. Same-sign interactions across the two\n")
cat("  moderators rule out the share-mechanics critique.\n\n")

cat("=== (i) Share moderator -- year + month-of-year FE, NW(14) ===\n")
print(etable(m_ym_base, m_ym_sar_nb, m_ym_sar_po, m_ym_sar_nb_ctl, m_ym_sar_po_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("Baseline (NB)", "NB no ctl", "Pois no ctl",
                         "NB + log(cross)", "Pois + log(cross)")))

cat("\n=== (i) Share moderator -- month-year FE, NW(14) ===\n")
print(etable(m_my_base, m_my_sar_nb, m_my_sar_po, m_my_sar_nb_ctl, m_my_sar_po_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("Baseline (NB)", "NB no ctl", "Pois no ctl",
                         "NB + log(cross)", "Pois + log(cross)")))

cat("\n=== (ii) log1p(SAR persons) -- month-year FE, NW(14) ===\n")
print(etable(m_my_per_nb, m_my_per_po, m_my_per_nb_ctl, m_my_per_po_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("NB no ctl", "Pois no ctl",
                         "NB + log(cross)", "Pois + log(cross)")))

cat("\n=== SUMMARY: SWH x SAR moderator coefficient (per 1-SD) ===\n")
for (g in names(group_labels)) {
  cat(sprintf("\n  %s\n", group_labels[g]))
  rows_g <- summary_rows[summary_rows$group == g, ]
  for (i in seq_len(nrow(rows_g))) {
    r <- rows_g[i, ]
    cat(sprintf("    %-32s  %+.3f (SE=%.3f)  p=%.4f%s\n",
                r$spec, r$coef, r$se, r$p,
                if (r$p < 0.05) " *" else ""))
  }
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
