# 051_zone_panel.R
# ================
# Extension to 05_reduced_form_primary.R: ZONE-LEVEL reduced-form regression.
#
# Splits the CMR corridor into SAR jurisdictions: African SAR (Libya + Tunisia,
# pooled) vs EU SAR (Italy + Malta, pooled), yielding a 2-bloc panel with
# within-day spatial variation.
#
# Primary spec (matches 05_reduced_form_primary.R but adds zone dimension):
#   n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou |
#                    month_year + sar_bloc
#
# Also estimates a 4-country robustness variant (| month_year + country).
#
# SE: NW(28), matching 05_reduced_form_primary.R. (Earlier a "dim" attribute
# on date from the zone build caused fixest's NW code to fail; that is now
# stripped in 03_build_zone_panel.R.)
#
# Input:  analysis/data/daily_panel_zone.RDS
# Output: output/tables/051_zone_panel.txt
#         output/figures/051_zone_panel_coef.png

library(tidyverse)
library(fixest)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
YEAR_START <- 2014
PERIODS <- c(2020, 2023)

cat("============================================================\n")
cat("051  ZONE-LEVEL REDUCED-FORM (AFR vs EU)\n")
cat("============================================================\n\n")

# ── 1. Load zone panel ──────────────────────────────────────
cat("--- 1. Loading zone panel ---\n")

p <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel_zone.RDS"))
cat(sprintf("  %d rows x %d columns\n", nrow(p), ncol(p)))

# ── 2. Collapse 4 zones -> 2 blocs (AFR, EU) ────────────────
cat("\n--- 2. Collapsing 4 zones to 2 blocs ---\n")

# Aggregate deaths by bloc (sum). SWH is identical across the 4 zones on
# any given day (under A2 — same corridor-wide SWH inherited from base),
# so a simple mean is exactly the same as the per-day SWH value.
bloc <- p %>%
  group_by(date, sar_bloc) %>%
  summarise(
    n_dead_missing = sum(n_dead_missing),
    swh            = mean(swh, na.rm = TRUE),
    swh_prev3days  = mean(swh_prev3days, na.rm = TRUE),
    swh_prevweek   = mean(swh_prevweek, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    post_mou   = as.integer(date >= MOU_DATE),
    year       = year(date),
    month_year = factor(format(date, "%Y-%m")),
    iso_week   = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))
  )

cat(sprintf("  Collapsed panel: %d rows (%d dates x 2 blocs)\n",
    nrow(bloc), length(unique(bloc$date))))

# ── 3. Estimate by period ───────────────────────────────────
cat("\n--- 3. Estimation ---\n")

all_results <- list()

sink_file <- file.path(BASE_DIR, "output", "tables", "051_zone_panel.txt")
sink(sink_file)

cat("051  ZONE-LEVEL REDUCED-FORM: AFR vs EU SAR\n")
cat("============================================\n")
cat("Outcome: n_dead_missing — incident-only, Cause = Drowning or Mixed/unknown\n")
cat("(matches 05_reduced_form_primary.R primary spec; built in 03_build_zone_panel.R)\n")
cat("Deaths assigned to SAR zones via spatial join (NOT country-of-incident)\n")
cat("Panel: 2 blocs x daily. FE: month_year + sar_bloc. NW(28) SEs.\n\n")

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("=== %s ===\n", label))
  cat(sprintf("--- sample: %d-%d ---\n", YEAR_START, ye))

  d_bloc <- bloc %>%
    filter(year >= YEAR_START, year <= ye, !is.na(swh_prevweek)) %>%
    mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)))

  d_zone <- p %>%
    filter(year >= YEAR_START, year <= ye, !is.na(swh_prevweek)) %>%
    mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)))

  cat(sprintf("  2-bloc panel rows: %d\n", nrow(d_bloc)))
  cat(sprintf("  4-zone panel rows: %d\n", nrow(d_zone)))

  # 2-bloc model: AFR vs EU, FE = month_year + sar_bloc, NW(28) SEs
  m_bloc <- fenegbin(
    n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou |
      month_year + sar_bloc,
    data = d_bloc, vcov = NW(28), panel.id = ~sar_bloc + date
  )

  # 4-country model: FE = month_year + country, NW(28) SEs
  m_zone <- fenegbin(
    n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou |
      month_year + country,
    data = d_zone, vcov = NW(28), panel.id = ~country + date
  )

  cat("\n  --- 2-bloc (AFR/EU) ---\n")
  print(etable(m_bloc, vcov = NW(28), se.below = TRUE))

  cat("\n  --- 4-country ---\n")
  print(etable(m_zone, vcov = NW(28), se.below = TRUE))

  # Extract the interaction term
  ct_bloc <- coeftable(m_bloc, vcov = NW(28))
  ct_zone <- coeftable(m_zone, vcov = NW(28))

  row_bloc <- which(grepl(":post_mou", rownames(ct_bloc)))
  row_zone <- which(grepl(":post_mou", rownames(ct_zone)))

  b3_bloc <- ct_bloc[row_bloc, 1]
  se_bloc <- ct_bloc[row_bloc, 2]
  p_bloc  <- 2 * pnorm(-abs(b3_bloc / se_bloc))

  b3_zone <- ct_zone[row_zone, 1]
  se_zone <- ct_zone[row_zone, 2]
  p_zone  <- 2 * pnorm(-abs(b3_zone / se_zone))

  cat(sprintf("\n  Summary:\n"))
  cat(sprintf("    2-bloc   (AFR/EU):   b3 = %+.3f (SE=%.3f)  IRR = %.3f  p = %.4f\n",
      b3_bloc, se_bloc, exp(b3_bloc), p_bloc))
  cat(sprintf("    4-country:           b3 = %+.3f (SE=%.3f)  IRR = %.3f  p = %.4f\n",
      b3_zone, se_zone, exp(b3_zone), p_zone))
  cat("\n")

  all_results[[label]] <- tibble(
    period = label,
    spec   = c("2-bloc (AFR/EU)", "4-country"),
    coef   = c(b3_bloc, b3_zone),
    se     = c(se_bloc, se_zone),
    p      = c(p_bloc, p_zone)
  )
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

# ── 4. Coefficient plot ─────────────────────────────────────
cat("\n--- 4. Coefficient plot ---\n")

combined <- bind_rows(all_results) %>%
  mutate(
    ci_lo = coef - 1.96 * se,
    ci_hi = coef + 1.96 * se,
    spec  = factor(spec, levels = c("4-country", "2-bloc (AFR/EU)"))
  )

p_coef <- ggplot(combined, aes(x = coef, y = spec, colour = period)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("2014-2020" = "#2166AC",
                                  "2014-2023" = "#B2182B")) +
  labs(
    title = expression(paste("Zone-level reduced-form: ",
                              beta[3], " (SWH x post-MoU)")),
    subtitle = "NegBin, NW(28) SEs, 95% CI. Deaths assigned via spatial join to SAR zones.",
    x = expression(paste(beta[3], " (per 1-SD SWH)")),
    y = NULL,
    colour = "Sample"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

fig_path <- file.path(BASE_DIR, "output", "figures", "051_zone_panel_coef.png")
ggsave(fig_path, p_coef, width = 9, height = 5, dpi = 200)
cat(sprintf("Saved: %s\n", fig_path))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
