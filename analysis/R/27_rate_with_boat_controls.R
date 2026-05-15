# 27_rate_with_boat_controls.R
# ============================
# Boat composition robustness for the per-person fatality rate model in
# 20_primary_model.R (m_rate / m_rate_u).
#
# Question: does the SWH x post_MoU shift in the per-person fatality rate
# survive controlling for boat composition?
#
# The Deiana et al. (2024) hypothesis: SAR availability -> smugglers
# substitute toward inflatables -> inflatables are more weather-sensitive
# -> the SWH-mortality slope steepens. Under this view, controlling for
# boat composition should ABSORB part of the SWH:post_mou interaction
# estimated by 20's rate model (since composition shifted post-MoU and
# composition mediates the weather effect).
#
# Three specifications, all Poisson with offset(log(crossing_attempts))
# matching 20's m_rate:
#   V1: Baseline (matches 20 m_rate, on the boat-observable sub-sample)
#         deaths ~ swh + swh:post_mou + offset | month_year_fac
#   V2: + boat composition as additive covariates
#         deaths ~ swh + swh:post_mou + inflatable_share + wooden_share + offset | month_year_fac
#   V3: + boat-share x SWH interaction (Deiana-style mediator probe)
#         deaths ~ swh + swh:post_mou + swh:inflatable_share + inflatable_share + wooden_share + offset | month_year_fac
#
# Interpretation guide (comparing SWH x post_MoU across versions):
#   - V1 ~ V2 ~ V3:           Composition isn't doing much; main result robust.
#   - V2 < V1 (in magnitude):  Composition absorbs part of SAR-removal effect
#                              (consistent with Deiana-style mediation).
#   - SWH x inflatable_share
#     significant in V3:       Independent evidence for the composition-
#                              as-mechanism story.
# We do NOT swap the controlled estimate in as the headline; we report all
# three and let the comparison itself be the evidence (see Theory).
#
# Sample note:
#   The boat-composition variables are defined only on days with
#   frx_incidents > 0. We therefore restrict the 27 sample to the
#   intersection of:
#     - 20's rate-sample filter (crossing_attempts > 0)
#     - frx_incidents > 0  (composition observable)
#   This is a slightly more restrictive sample than 20's m_rate. The V1
#   coefficients in 27 will therefore differ slightly from m_rate in 20
#   purely due to sample restriction; magnitude and sign should match.
#
# Estimated for both IOM and UNITED death series (matching 20's coverage).
#
# In:  analysis/data/daily_panel_complete.RDS
#      data/processed/iom_mmp_incidents.RDS
#      data/processed/united_incidents.RDS
#      data/processed/core_corridor.RDS
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
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    log1p_lc_lag14 = log1p(lc_lag14),
    unit           = 1L,
    year_fac       = factor(year),
    month_year_fac = factor(month_year)
  )

# 20's rate-model sample: crossing_attempts > 0 + lc_lag14/swh_prev5days non-NA
d_rate_full <- panel %>%
  filter(!is.na(lc_lag14), !is.na(swh_prev5days),
         crossing_attempts > 0) %>%
  mutate(log_attempts = log(crossing_attempts))

# 27 boat-observable sub-sample: additionally require frx_incidents > 0
# so that boat-composition shares are defined.
d_rate_boat <- d_rate_full %>%
  filter(frx_incidents > 0)

cat(sprintf("  Rate sample (20 baseline):     N = %d days\n", nrow(d_rate_full)))
cat(sprintf("  Boat-observable sub-sample:    N = %d days (%.1f%% of full)\n",
            nrow(d_rate_boat),
            100 * nrow(d_rate_boat) / nrow(d_rate_full)))
cat(sprintf("  Days dropped (no Frontex events): %d\n",
            nrow(d_rate_full) - nrow(d_rate_boat)))

cat("\n  Boat composition summary (boat-observable sample):\n")
cat(sprintf("    Inflatable share: mean = %.3f  sd = %.3f  range = [%.2f, %.2f]\n",
            mean(d_rate_boat$frx_inflatable_share), sd(d_rate_boat$frx_inflatable_share),
            min(d_rate_boat$frx_inflatable_share), max(d_rate_boat$frx_inflatable_share)))
cat(sprintf("    Wooden share:     mean = %.3f  sd = %.3f  range = [%.2f, %.2f]\n",
            mean(d_rate_boat$frx_wooden_share), sd(d_rate_boat$frx_wooden_share),
            min(d_rate_boat$frx_wooden_share), max(d_rate_boat$frx_wooden_share)))

cat("\n  Pre/post-MoU mean inflatable share (descriptive, matches Deiana Fig 9):\n")
print(d_rate_boat %>%
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
               offset(log_attempts) | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

# V2: + boat composition (additive)
v2_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               frx_inflatable_share + frx_wooden_share +
               offset(log_attempts) | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

# V3: + boat × SWH interaction (Deiana-style mediator probe)
v3_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               frx_inflatable_share + frx_wooden_share +
               swh_prev5days:frx_inflatable_share +
               offset(log_attempts) | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

cat("\n  IOM, NW(14):\n")
print(etable(v1_iom, v2_iom, v3_iom,
             vcov = NW(14), se.below = TRUE,
             headers = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH")))

# ── 3. Three specifications -- UNITED ────────────────────────
cat("\n--- 3. UNITED rate model: V1, V2, V3 ---\n")

v1_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  offset(log_attempts) | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

v2_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  frx_inflatable_share + frx_wooden_share +
                  offset(log_attempts) | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

v3_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  frx_inflatable_share + frx_wooden_share +
                  swh_prev5days:frx_inflatable_share +
                  offset(log_attempts) | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

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

cat("\n  SWH x post_MoU (per-person fatality rate, NW(14)):\n\n")
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
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2,
                 position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  labs(
    title    = "Per-person fatality rate: SWH x post_MoU across boat-control specifications",
    subtitle = "Poisson with offset(log(crossing_attempts)). Month-year FE. NW(14) SEs.",
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
            min(d_rate_boat$date), max(d_rate_boat$date)))
cat(sprintf("                           N = %d days\n", nrow(d_rate_boat)))
cat(sprintf("                           (vs. 20 rate sample N = %d; %d days dropped\n",
            nrow(d_rate_full), nrow(d_rate_full) - nrow(d_rate_boat)))
cat( "                            because frx_incidents = 0 -> boat shares undefined)\n")
cat("Outcome: per-person fatality rate via Poisson with offset(log(crossing_attempts))\n")
cat("Standard errors: Newey-West (lag = 14)\n")
cat("Crossing_attempts = frx_persons + lcg_tcg_pushbacks + n_dead_missing\n\n")

cat("Three specifications:\n")
cat("  V1: Baseline (matches 20 m_rate spec, on the boat-observable sample)\n")
cat("       deaths ~ swh + swh:post_mou + offset | month_year_fac\n")
cat("  V2: + boat composition (additive)\n")
cat("       deaths ~ swh + swh:post_mou + inflatable_share + wooden_share\n")
cat("              + offset | month_year_fac\n")
cat("  V3: + boat x SWH interaction (Deiana-style mediator probe)\n")
cat("       deaths ~ swh + swh:post_mou + swh:inflatable_share +\n")
cat("              + inflatable_share + wooden_share + offset | month_year_fac\n\n")

cat("Interpretation guide -- comparing SWH x post_MoU across V1, V2, V3:\n")
cat("  V1 ~ V2 ~ V3:           Composition isn't doing much; main result robust.\n")
cat("  V2 magnitude < V1:      Composition absorbs part of SAR-removal effect\n")
cat("                          (consistent with Deiana mediation).\n")
cat("  SWH x inflatable\n")
cat("    significant in V3:    Independent evidence for composition-as-mechanism.\n\n")

cat("=== Boat composition: pre/post-MoU descriptive ===\n")
print(d_rate_boat %>%
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
