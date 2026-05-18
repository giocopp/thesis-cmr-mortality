# 27_rate_with_boat_controls.R
# ============================
# Boat-composition robustness for the rate model (m_rate / m_rate_u) in
# 20_primary_model.R: does the SWH:post_mou shift hold with boat-composition
# controls (Deiana-style composition-mediation probe)?
#
# Sample: 20's source-specific rate sample (attempts_src > 0, lc_lag14 &
# swh_prev5days non-NA) further restricted to frx_incidents > 0 (boat shares
# defined only there); V1 here differs from 20 by sample restriction only.
#
# Specs (Poisson, source-specific offset, month_year FE, NW(14), IOM+UNITED):
#   V1: deaths ~ swh + swh:post_mou + offset
#   V2: V1 + frx_inflatable_share + frx_wooden_share
#   V3: V2 + swh:frx_inflatable_share
#
# In:  analysis/data/daily_panel_complete.RDS
#      data/processed/{iom_mmp_incidents,united_incidents,core_corridor}.RDS
# Out: output/tables/27_rate_with_boat_controls.txt
#      output/figures/27_rate_with_boat_controls_coefplot.png

library(tidyverse)
library(lubridate)
library(fixest)
library(sf)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("27  RATE MODEL WITH BOAT COMPOSITION CONTROLS\n")
cat("    Extension of m_rate / m_rate_u in 20_primary_model.R\n")
cat("============================================================\n\n")

# ── 1. Load data (mirror of 20's data prep) ──────────────────
cat("--- 1. Loading data + 20-style prep ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

# IOM daily (same helper / default filter as 20)
iom_daily <- build_iom_daily()

# UNITED daily via the shared builder. Defaults (corridor spatial join;
# country in CMR+Med; manner drowned/other_unknown) replicate the previous
# inline filter exactly — single source of truth, see _helpers.R.
united_daily <- build_united_daily()

panel <- panel %>%
  left_join(iom_daily %>% rename(n_dead_iom = n_dead_missing), by = "date") %>%
  left_join(united_daily, by = "date") %>%
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) %>%
  add_crossing_exposure() %>%   # living_crossings + source-specific attempts (shared, _helpers.R)
  mutate(
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    log1p_lc_lag14 = log1p(lc_lag14),
    unit           = 1L,
    year_fac       = factor(year),
    month_year_fac = factor(month_year)
  )

# 20's rate-model sample: attempts_source > 0 + lc_lag14/swh_prev5days non-NA
d_rate_full_iom <- panel %>%
  filter(!is.na(lc_lag14), !is.na(swh_prev5days),
         attempts_iom > 0) %>%
  mutate(log_attempts_iom = log(attempts_iom))

d_rate_full_united <- panel %>%
  filter(!is.na(lc_lag14), !is.na(swh_prev5days),
         attempts_united > 0) %>%
  mutate(log_attempts_united = log(attempts_united))

# 27 boat-observable sub-sample: additionally require frx_incidents > 0
# so that boat-composition shares are defined.
d_rate_boat_iom <- d_rate_full_iom %>%
  filter(frx_incidents > 0)

d_rate_boat_united <- d_rate_full_united %>%
  filter(frx_incidents > 0)

cat(sprintf("  Rate sample IOM (20 baseline):       N = %d days\n",
            nrow(d_rate_full_iom)))
cat(sprintf("  Rate sample UNITED (20 baseline):    N = %d days\n",
            nrow(d_rate_full_united)))
cat(sprintf("  Boat-observable IOM sub-sample:      N = %d days (%.1f%% of full)\n",
            nrow(d_rate_boat_iom),
            100 * nrow(d_rate_boat_iom) / nrow(d_rate_full_iom)))
cat(sprintf("  Boat-observable UNITED sub-sample:   N = %d days (%.1f%% of full)\n",
            nrow(d_rate_boat_united),
            100 * nrow(d_rate_boat_united) / nrow(d_rate_full_united)))
cat(sprintf("  Days dropped, IOM (no Frontex events):    %d\n",
            nrow(d_rate_full_iom) - nrow(d_rate_boat_iom)))
cat(sprintf("  Days dropped, UNITED (no Frontex events): %d\n",
            nrow(d_rate_full_united) - nrow(d_rate_boat_united)))

cat("\n  Boat composition summary (IOM boat-observable sample):\n")
cat(sprintf("    Inflatable share: mean = %.3f  sd = %.3f  range = [%.2f, %.2f]\n",
            mean(d_rate_boat_iom$frx_inflatable_share),
            sd(d_rate_boat_iom$frx_inflatable_share),
            min(d_rate_boat_iom$frx_inflatable_share),
            max(d_rate_boat_iom$frx_inflatable_share)))
cat(sprintf("    Wooden share:     mean = %.3f  sd = %.3f  range = [%.2f, %.2f]\n",
            mean(d_rate_boat_iom$frx_wooden_share),
            sd(d_rate_boat_iom$frx_wooden_share),
            min(d_rate_boat_iom$frx_wooden_share),
            max(d_rate_boat_iom$frx_wooden_share)))

cat("\n  Pre/post-MoU mean inflatable share (descriptive, matches Deiana Fig 9):\n")
print(d_rate_boat_iom %>%
        mutate(period = ifelse(date >= MOU_DATE, "Post-MoU", "Pre-MoU")) %>%
        group_by(period) %>%
        summarise(n_days = n(),
                  mean_inflatable = mean(frx_inflatable_share),
                  mean_wooden     = mean(frx_wooden_share),
                  .groups = "drop"))

# ── 2. Three specifications -- IOM ───────────────────────────
cat("\n--- 2. IOM rate model: V1, V2, V3 ---\n")

# V1: baseline rate model (matches 20 m_rate spec)
v1_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               offset(log_attempts_iom) | month_year_fac,
  data = d_rate_boat_iom, vcov = NW(14), panel.id = ~unit + date)

# V2: + boat composition (additive)
v2_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               frx_inflatable_share + frx_wooden_share +
               offset(log_attempts_iom) | month_year_fac,
  data = d_rate_boat_iom, vcov = NW(14), panel.id = ~unit + date)

# V3: + boat × SWH interaction (Deiana-style mediator probe)
v3_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               frx_inflatable_share + frx_wooden_share +
               swh_prev5days:frx_inflatable_share +
               offset(log_attempts_iom) | month_year_fac,
  data = d_rate_boat_iom, vcov = NW(14), panel.id = ~unit + date)

cat("\n  IOM, NW(14):\n")
print(etable(v1_iom, v2_iom, v3_iom,
             vcov = NW(14), se.below = TRUE,
             headers = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH")))

# ── 3. Three specifications -- UNITED ────────────────────────
cat("\n--- 3. UNITED rate model: V1, V2, V3 ---\n")

v1_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  offset(log_attempts_united) | month_year_fac,
  data = d_rate_boat_united, vcov = NW(14), panel.id = ~unit + date)

v2_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  frx_inflatable_share + frx_wooden_share +
                  offset(log_attempts_united) | month_year_fac,
  data = d_rate_boat_united, vcov = NW(14), panel.id = ~unit + date)

v3_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  frx_inflatable_share + frx_wooden_share +
                  swh_prev5days:frx_inflatable_share +
                  offset(log_attempts_united) | month_year_fac,
  data = d_rate_boat_united, vcov = NW(14), panel.id = ~unit + date)

cat("\n  UNITED, NW(14):\n")
print(etable(v1_united, v2_united, v3_united,
             vcov = NW(14), se.below = TRUE,
             headers = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH")))

# ── 4. Comparison summary ────────────────────────────────────
cat("\n--- 4. Summary: SWH x post_MoU across V1, V2, V3 ---\n")

extract_b3 <- function(m, label, source) {
  ct <- coeftable(m, vcov = NW(14))
  r <- which(rownames(ct) == "swh_prev5days:post_mou")
  if (length(r) == 0) {
    return(tibble(spec = label, source = source,
                  coef = NA_real_, se = NA_real_, p = NA_real_))
  }
  tibble(spec = label, source = source,
         coef = ct[r, 1], se = ct[r, 2],
         p = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

extract_swh_x_inflatable <- function(m, label, source) {
  ct <- coeftable(m, vcov = NW(14))
  r <- which(rownames(ct) == "swh_prev5days:frx_inflatable_share")
  if (length(r) == 0) {
    return(tibble(spec = label, source = source,
                  coef = NA_real_, se = NA_real_, p = NA_real_))
  }
  tibble(spec = label, source = source,
         coef = ct[r, 1], se = ct[r, 2],
         p = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

summary_b3 <- bind_rows(
  extract_b3(v1_iom,    "V1: Baseline",     "IOM"),
  extract_b3(v2_iom,    "V2: +boat ctrl",   "IOM"),
  extract_b3(v3_iom,    "V3: +boat x SWH",  "IOM"),
  extract_b3(v1_united, "V1: Baseline",     "UNITED"),
  extract_b3(v2_united, "V2: +boat ctrl",   "UNITED"),
  extract_b3(v3_united, "V3: +boat x SWH",  "UNITED")
)

cat("\n  SWH x post_MoU (recorded-death rate, NW(14)):\n\n")
for (i in seq_len(nrow(summary_b3))) {
  r <- summary_b3[i, ]
  star <- if (!is.na(r$p) && r$p < 0.05) " *" else ""
  cat(sprintf("    %-8s %-18s  b = %+.3f (SE=%.3f)  p = %.4f%s\n",
              r$source, r$spec, r$coef, r$se, r$p, star))
}

v3_inflatable <- bind_rows(
  extract_swh_x_inflatable(v3_iom,    "V3", "IOM"),
  extract_swh_x_inflatable(v3_united, "V3", "UNITED")
)

cat("\n  SWH x inflatable_share (V3 only, Deiana-style mediator probe):\n\n")
for (i in seq_len(nrow(v3_inflatable))) {
  r <- v3_inflatable[i, ]
  star <- if (!is.na(r$p) && r$p < 0.05) " *" else ""
  cat(sprintf("    %-8s %-18s  b = %+.3f (SE=%.3f)  p = %.4f%s\n",
              r$source, r$spec, r$coef, r$se, r$p, star))
}

# ── 5. Coefficient plot ──────────────────────────────────────
cat("\n--- 5. Plot ---\n")

plot_df <- summary_b3 %>%
  mutate(
    ci_lo = coef - 1.96 * se,
    ci_hi = coef + 1.96 * se,
    spec_f = factor(spec,
                    levels = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH"))
  )

p <- ggplot(plot_df, aes(coef, spec_f, colour = source)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                width = 0.2, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  labs(
    title    = "Recorded-death rate: SWH x post_MoU across boat-control specifications",
    subtitle = "Poisson with source-specific offsets. Month-year FE. NW(14) SEs.",
    x        = "SWH x post_MoU coefficient (per 1m SWH, log rate)",
    y        = NULL,
    colour   = "Source"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(BASE_DIR, "output", "figures",
                  "27_rate_with_boat_controls_coefplot.png"),
       p, width = 11, height = 5, dpi = 200)
cat("  Saved: output/figures/27_rate_with_boat_controls_coefplot.png\n")

# ── 6. Save text output ──────────────────────────────────────
cat("\n--- 6. Saving results ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables",
                        "27_rate_with_boat_controls.txt")
sink(sink_file)

cat("27  RATE MODEL WITH BOAT COMPOSITION CONTROLS\n")
cat("Extension of m_rate / m_rate_u in 20_primary_model.R\n")
cat("=====================================================\n\n")

cat(sprintf("Sample (boat-observable):  %s to %s\n",
            min(d_rate_boat_iom$date), max(d_rate_boat_iom$date)))
cat(sprintf("                           IOM N = %d days\n", nrow(d_rate_boat_iom)))
cat(sprintf("                           UNITED N = %d days\n", nrow(d_rate_boat_united)))
cat(sprintf("                           (vs. 20 rate sample: IOM N = %d, UNITED N = %d)\n",
            nrow(d_rate_full_iom), nrow(d_rate_full_united)))
cat("                           Sample is restricted to frx_incidents > 0 so boat shares are defined.\n")
cat("Outcome: recorded deaths per source-specific constructed crossing attempt\n")
cat("Standard errors: Newey-West (lag = 14)\n")
cat("IOM offset:    log(frx_persons + lcg_tcg_pushbacks + n_dead_iom)\n")
cat("UNITED offset: log(frx_persons + lcg_tcg_pushbacks + n_dead_united)\n")
cat("N = estimation N after fixest drops all-zero FE cells.\n\n")

cat("Three specifications:\n")
cat("  V1: Baseline (matches 20 m_rate spec, on the boat-observable sample)\n")
cat("       deaths ~ swh + swh:post_mou + offset | month_year_fac\n")
cat("  V2: + boat composition (additive)\n")
cat("       deaths ~ swh + swh:post_mou + inflatable_share + wooden_share\n")
cat("              + offset | month_year_fac\n")
cat("  V3: + boat x SWH interaction (Deiana-style mediator probe)\n")
cat("       deaths ~ swh + swh:post_mou + swh:inflatable_share +\n")
cat("              + inflatable_share + wooden_share + offset | month_year_fac\n\n")

cat("=== Boat composition: pre/post-MoU descriptive ===\n")
print(d_rate_boat_iom %>%
        mutate(period = ifelse(date >= MOU_DATE, "Post-MoU", "Pre-MoU")) %>%
        group_by(period) %>%
        summarise(n_days = n(),
                  mean_inflatable = mean(frx_inflatable_share),
                  mean_wooden     = mean(frx_wooden_share),
                  .groups = "drop"))

cat("\n=== IOM rate model -- V1, V2, V3 (NW(14)) ===\n")
print(etable(v1_iom, v2_iom, v3_iom,
             vcov = NW(14), se.below = TRUE,
             headers = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH")))

cat("\n=== UNITED rate model -- V1, V2, V3 (NW(14)) ===\n")
print(etable(v1_united, v2_united, v3_united,
             vcov = NW(14), se.below = TRUE,
             headers = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH")))

cat("\n=== SUMMARY: SWH x post_MoU coefficient across versions ===\n\n")
for (i in seq_len(nrow(summary_b3))) {
  r <- summary_b3[i, ]
  star <- if (!is.na(r$p) && r$p < 0.05) " *" else ""
  cat(sprintf("  %-8s %-18s  b = %+.3f (SE=%.3f)  p = %.4f%s\n",
              r$source, r$spec, r$coef, r$se, r$p, star))
}

cat("\n=== SWH x inflatable_share interaction (V3 only, Deiana-style) ===\n\n")
for (i in seq_len(nrow(v3_inflatable))) {
  r <- v3_inflatable[i, ]
  star <- if (!is.na(r$p) && r$p < 0.05) " *" else ""
  cat(sprintf("  %-8s %-18s  b = %+.3f (SE=%.3f)  p = %.4f%s\n",
              r$source, r$spec, r$coef, r$se, r$p, star))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
