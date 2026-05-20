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
# UNITED is the primary source (consistent with 20_primary_model.R); IOM is
# the comparison. Both are estimated on the same sample with the same specs.
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
#      data/processed/{iom_mmp_incidents,united_incidents,core_corridor}.RDS
# Out: output/tables/23_mechanism_interactions.txt
#      output/figures/23_mechanism_coefplot.png

library(tidyverse)
library(lubridate)
library(fixest)
library(sf)

BASE_DIR <- here::here()

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("23  MECHANISM: SWH x SAR share (UNITED primary + IOM comparison)\n")
cat("============================================================\n\n")

# ── 1. Load data + build weekly-lagged SAR share ─────────────
cat("--- 1. Loading data ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

# IOM and UNITED daily series via shared corridor-joined builders.
iom_daily    <- build_iom_daily()
united_daily <- build_united_daily()

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
  left_join(iom_daily    %>% rename(n_dead_iom = n_dead_missing), by = "date") %>%
  left_join(united_daily,                                          by = "date") %>%
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) %>%
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
for (src in c("IOM", "UNITED")) {
  v <- if (src == "IOM") d$n_dead_iom else d$n_dead_united
  cat(sprintf("  %-7s deaths: %.0f over %d death-days\n",
              src, sum(v), sum(v > 0)))
}
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
cat("\n--- 2. Estimation ---\n")

# Helper: fit the full mechanism model set for one outcome variable.
# All models: year+month-of-year FE and month-year FE; share and persons
# moderators; with and without lag-14 crossing control.
fit_mech_set <- function(dep, d) {
  nb <- function(rhs) fenegbin(as.formula(sprintf("%s ~ %s", dep, rhs)),
                                data = d, vcov = NW(14), panel.id = ~unit + date)
  po <- function(rhs) fepois  (as.formula(sprintf("%s ~ %s", dep, rhs)),
                                data = d, vcov = NW(14), panel.id = ~unit + date)

  s0  <- "swh_prev5days | year_fac + month_of_year"
  s1  <- "swh_prev5days + swh_prev5days:sar_share_pw_z + sar_share_pw_z | year_fac + month_of_year"
  s1c <- "swh_prev5days + swh_prev5days:sar_share_pw_z + sar_share_pw_z + log1p_lc_lag14 | year_fac + month_of_year"
  m0  <- "swh_prev5days | month_year_fac"
  m1  <- "swh_prev5days + swh_prev5days:sar_share_pw_z + sar_share_pw_z | month_year_fac"
  m1c <- "swh_prev5days + swh_prev5days:sar_share_pw_z + sar_share_pw_z + log1p_lc_lag14 | month_year_fac"
  p1  <- "swh_prev5days + swh_prev5days:log1p_sar_persons_pw_z + log1p_sar_persons_pw_z | month_year_fac"
  p1c <- "swh_prev5days + swh_prev5days:log1p_sar_persons_pw_z + log1p_sar_persons_pw_z + log1p_lc_lag14 | month_year_fac"

  list(
    ym_base       = nb(s0),
    ym_sar_nb     = nb(s1),
    ym_sar_po     = po(s1),
    ym_sar_nb_ctl = nb(s1c),
    ym_sar_po_ctl = po(s1c),
    my_base       = nb(m0),
    my_sar_nb     = nb(m1),
    my_sar_po     = po(m1),
    my_sar_nb_ctl = nb(m1c),
    my_sar_po_ctl = po(m1c),
    my_per_nb     = nb(p1),
    my_per_po     = po(p1),
    my_per_nb_ctl = nb(p1c),
    my_per_po_ctl = po(p1c)
  )
}

cat("  Fitting IOM models...\n")
fits_iom    <- fit_mech_set("n_dead_iom",    d)
cat("  Fitting UNITED models...\n")
fits_united <- fit_mech_set("n_dead_united", d)

# Print etables for each source.
for (src in c("IOM", "UNITED")) {
  fits <- if (src == "IOM") fits_iom else fits_united
  cat(sprintf("\n  [%s] Year + month-of-year FE, NW(14):\n", src))
  print(etable(fits$ym_base, fits$ym_sar_nb, fits$ym_sar_po,
               fits$ym_sar_nb_ctl, fits$ym_sar_po_ctl,
               vcov = NW(14), se.below = TRUE,
               headers = c("Baseline (NB)", "NB no ctl", "Pois no ctl",
                           "NB + log(cross)", "Pois + log(cross)")))
  cat(sprintf("\n  [%s] Month-year FE, NW(14):\n", src))
  print(etable(fits$my_base, fits$my_sar_nb, fits$my_sar_po,
               fits$my_sar_nb_ctl, fits$my_sar_po_ctl,
               vcov = NW(14), se.below = TRUE,
               headers = c("Baseline (NB)", "NB no ctl", "Pois no ctl",
                           "NB + log(cross)", "Pois + log(cross)")))
  cat(sprintf("\n  [%s] log1p(SAR persons), month-year FE, NW(14):\n", src))
  print(etable(fits$my_per_nb, fits$my_per_po,
               fits$my_per_nb_ctl, fits$my_per_po_ctl,
               vcov = NW(14), se.below = TRUE,
               headers = c("NB no ctl", "Pois no ctl",
                           "NB + log(cross)", "Pois + log(cross)")))
}

# ── 3. Summary of b_int (SWH x SAR moderator) ────────────────
cat("\n--- 3. Summary (SWH x SAR moderator) ---\n")

extract_int <- function(m, moderator, label, group, source) {
  ct  <- coeftable(m, vcov = NW(14))
  pat <- paste0("swh_prev5days:", moderator)
  r   <- grep(pat, rownames(ct), fixed = TRUE)
  tibble(source = source, group = group, spec = label,
         coef  = ct[r, 1], se = ct[r, 2],
         p     = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

build_summary_rows <- function(fits, src) {
  bind_rows(
    extract_int(fits$ym_sar_nb,     "sar_share_pw_z", "yr+mo FE, NegBin  (primary)",    "share",   src),
    extract_int(fits$ym_sar_po,     "sar_share_pw_z", "yr+mo FE, Poisson (primary)",    "share",   src),
    extract_int(fits$ym_sar_nb_ctl, "sar_share_pw_z", "yr+mo FE, NegBin  + log(cross)", "share",   src),
    extract_int(fits$ym_sar_po_ctl, "sar_share_pw_z", "yr+mo FE, Poisson + log(cross)", "share",   src),
    extract_int(fits$my_sar_nb,     "sar_share_pw_z", "my FE,    NegBin  (primary)",    "share",   src),
    extract_int(fits$my_sar_po,     "sar_share_pw_z", "my FE,    Poisson (primary)",    "share",   src),
    extract_int(fits$my_sar_nb_ctl, "sar_share_pw_z", "my FE,    NegBin  + log(cross)", "share",   src),
    extract_int(fits$my_sar_po_ctl, "sar_share_pw_z", "my FE,    Poisson + log(cross)", "share",   src),
    extract_int(fits$my_per_nb,     "log1p_sar_persons_pw_z", "my FE,    NegBin  (primary)",    "persons", src),
    extract_int(fits$my_per_po,     "log1p_sar_persons_pw_z", "my FE,    Poisson (primary)",    "persons", src),
    extract_int(fits$my_per_nb_ctl, "log1p_sar_persons_pw_z", "my FE,    NegBin  + log(cross)", "persons", src),
    extract_int(fits$my_per_po_ctl, "log1p_sar_persons_pw_z", "my FE,    Poisson + log(cross)", "persons", src)
  )
}

summary_rows <- bind_rows(
  build_summary_rows(fits_united, "UNITED"),
  build_summary_rows(fits_iom,    "IOM")
)

group_labels <- c(
  share   = "(i)  SAR share",
  persons = "(ii) log1p(SAR persons) -- absolute level"
)

cat("\n  b_int = coefficient on swh_prev5days x SAR moderator (per 1-SD)\n")
for (src in c("UNITED", "IOM")) {
  for (g in names(group_labels)) {
    cat(sprintf("\n  [%s] %s\n", src, group_labels[g]))
    rows_g <- summary_rows[summary_rows$source == src & summary_rows$group == g, ]
    for (i in seq_len(nrow(rows_g))) {
      r <- rows_g[i, ]
      cat(sprintf("    %-32s  %+.3f (SE=%.3f)  p=%.4f%s\n",
                  r$spec, r$coef, r$se, r$p,
                  if (r$p < 0.05) " *" else ""))
    }
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
    source = factor(source, levels = c("UNITED", "IOM")),
    spec_f = factor(paste0(group, "::", spec),
                     levels = rev(unique(paste0(group, "::", spec))))
  )

p <- ggplot(plot_df, aes(coef, spec_f, colour = source, shape = source)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.8, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.25, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("UNITED" = "#B2182B", "IOM" = "#2166AC")) +
  scale_shape_manual(values  = c("UNITED" = 16,        "IOM" = 17)) +
  facet_wrap(~ group_lab, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[a-z]+::", "", x)) +
  labs(
    title    = "Mechanism: SWH x SAR moderator interactions",
    subtitle = "Two SAR moderators (share, log absolute persons), z-scored, weekly-lagged. NW(14) SEs.",
    x        = "Coefficient on SWH x SAR moderator (per 1-SD)",
    y        = NULL,
    colour   = "Source",
    shape    = "Source"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

ggsave(file.path(BASE_DIR, "output", "figures",
                  "23_mechanism_coefplot.png"),
       p, width = 12, height = 7, dpi = 200)
cat("  Saved: output/figures/23_mechanism_coefplot.png\n")

# ── 5. Save text output ──────────────────────────────────────
cat("\n--- 5. Saving results ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables",
                        "23_mechanism_interactions.txt")
sink(sink_file)

cat("23  MECHANISM: SWH x SAR moderator (UNITED primary + IOM comparison)\n")
cat("======================================================================\n")
cat(sprintf("Sample: %s to %s (N = %d days)\n",
            min(d$date), max(d$date), nrow(d)))
for (src in c("UNITED", "IOM")) {
  v <- if (src == "IOM") d$n_dead_iom else d$n_dead_united
  cat(sprintf("  %-7s deaths: %.0f over %d death-days\n",
              src, sum(v), sum(v > 0)))
}
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
cat("  -0.50 at lag 7 to -0.25 at lag 14.\n\n")

cat("WHY TWO MODERATORS:\n")
cat("  The share denominator includes LCG/TCG pullbacks and other Not-SAR\n")
cat("  events. After 2018 those rise sharply, so the share can fall\n")
cat("  mechanically even if absolute SAR is unchanged. The log-persons\n")
cat("  variant drops the denominator. Same-sign interactions across the two\n")
cat("  moderators rule out the share-mechanics critique.\n\n")

for (src in c("UNITED", "IOM")) {
  fits <- if (src == "IOM") fits_iom else fits_united

  cat(sprintf("=== [%s] (i) Share moderator -- year + month-of-year FE, NW(14) ===\n", src))
  print(etable(fits$ym_base, fits$ym_sar_nb, fits$ym_sar_po,
               fits$ym_sar_nb_ctl, fits$ym_sar_po_ctl,
               vcov = NW(14), se.below = TRUE,
               headers = c("Baseline (NB)", "NB no ctl", "Pois no ctl",
                           "NB + log(cross)", "Pois + log(cross)")))

  cat(sprintf("\n=== [%s] (i) Share moderator -- month-year FE, NW(14) ===\n", src))
  print(etable(fits$my_base, fits$my_sar_nb, fits$my_sar_po,
               fits$my_sar_nb_ctl, fits$my_sar_po_ctl,
               vcov = NW(14), se.below = TRUE,
               headers = c("Baseline (NB)", "NB no ctl", "Pois no ctl",
                           "NB + log(cross)", "Pois + log(cross)")))

  cat(sprintf("\n=== [%s] (ii) log1p(SAR persons) -- month-year FE, NW(14) ===\n", src))
  print(etable(fits$my_per_nb, fits$my_per_po,
               fits$my_per_nb_ctl, fits$my_per_po_ctl,
               vcov = NW(14), se.below = TRUE,
               headers = c("NB no ctl", "Pois no ctl",
                           "NB + log(cross)", "Pois + log(cross)")))
  cat("\n")
}

cat("=== SUMMARY: SWH x SAR moderator coefficient (per 1-SD) ===\n")
for (src in c("UNITED", "IOM")) {
  for (g in names(group_labels)) {
    cat(sprintf("\n  [%s] %s\n", src, group_labels[g]))
    rows_g <- summary_rows[summary_rows$source == src & summary_rows$group == g, ]
    for (i in seq_len(nrow(rows_g))) {
      r <- rows_g[i, ]
      cat(sprintf("    %-32s  %+.3f (SE=%.3f)  p=%.4f%s\n",
                  r$spec, r$coef, r$se, r$p,
                  if (r$p < 0.05) " *" else ""))
    }
  }
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
