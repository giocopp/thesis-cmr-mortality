# 21_zone_panel.R
# ================
# Zone-level reduced form, aligned with 20_primary_model.R.
#
# Extends the primary 05d spec (corridor-wide daily count) to zone-level
# variation by estimating the same model on a 2-bloc (AFR vs EU) panel.
# (4-country variant dropped — 2-bloc captures the spatial heterogeneity
# of interest for this project: African SAR vs EU SAR response.)
#
# Spec (identical to 05d, with an added FE for zone):
#   n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou |
#                month_year_fac + sar_bloc
#
# Both NegBin (fenegbin) and Poisson QMLE (fepois), NW(14) SEs — matching
# the dual-family approach in 05d.
#
# Sample: identical to 05d — requires !is.na(lc_lag14) & !is.na(swh_prev5days).
# The lc_lag14 filter drops the first 14 days of 2014 and is computed from
# living_crossings = frx_persons + lcg_tcg_pushbacks (see 05d header).
#
# Outcome naming: the zone panel's `n_dead_missing` column is renamed to
# `n_dead_iom` for consistency with 05d. The zone deaths use the SAME filter
# as build_iom_daily() defaults (incident only, split EXCLUDED; drowning + mixed,
# Central Mediterranean route), plus spatial assignment to the corridor-
# intersected SAR polygons (see 02_build_zone_panel.R).
#
# Input:  analysis/data/daily_panel_complete.RDS
#         analysis/data/daily_panel_zone.RDS
# Output: output/tables/21_zone_panel.txt

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("051  ZONE-LEVEL REDUCED-FORM (AFR vs EU)\n")
cat("     Aligned with 20_primary_model.R\n")
cat("============================================================\n\n")

# ── 1. 05d data prep on the daily panel ──────────────────────
cat("--- 1. Loading daily panel + 05d data prep ---\n")

panel_daily <- readRDS(file.path(BASE_DIR, "analysis", "data",
                                   "daily_panel_complete.RDS"))

iom_daily <- build_iom_daily()

panel_daily <- panel_daily %>%
  left_join(iom_daily %>% rename(n_dead_iom = n_dead_missing), by = "date") %>%
  replace_na(list(n_dead_iom = 0)) %>%
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag7  = dplyr::lag(zoo::rollmeanr(living_crossings, k = 7,
                                           fill = NA, align = "right"), 1),
    lc_lag14 = dplyr::lag(zoo::rollmeanr(living_crossings, k = 7,
                                           fill = NA, align = "right"), 8),
    log1p_lc_lag7  = log1p(lc_lag7),
    log1p_lc_lag14 = log1p(lc_lag14),
    unit           = 1L,
    year_fac       = factor(year),
    month_year_fac = factor(month_year)
  )

# Restrict to the 05d primary sample
sample_dates <- panel_daily %>%
  filter(!is.na(lc_lag14), !is.na(swh_prev5days)) %>%
  pull(date)

cat(sprintf("  05d sample: %s to %s (N = %d days)\n",
            min(sample_dates), max(sample_dates), length(sample_dates)))

# ── 2. Zone panel collapsed to 2 blocs, restricted to 05d sample ─
cat("\n--- 2. 2-bloc panel (AFR vs EU) ---\n")

zp <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_zone.RDS")) %>%
  filter(date %in% sample_dates) %>%
  rename(n_dead_iom = n_dead_missing) %>%
  mutate(month_year_fac = factor(format(date, "%Y-%m")))

# 2-bloc collapse. SWH, post_mou, crossing_attempts are identical across
# zones on any given day (A2: corridor-wide SWH; post_mou and attempts are
# corridor-wide by construction).
bloc <- zp %>%
  group_by(date, sar_bloc) %>%
  summarise(
    n_dead_iom        = sum(n_dead_iom),
    swh_prev5days      = first(swh_prev5days),
    post_mou          = first(post_mou),
    crossing_attempts = first(crossing_attempts),
    month_year_fac    = first(month_year_fac),
    .groups = "drop"
  )
dim(bloc$date) <- NULL

cat(sprintf("  2-bloc panel: %d rows (%d dates x 2 blocs)\n",
            nrow(bloc), length(unique(bloc$date))))
cat(sprintf("  Total deaths (IOM zone sum): %.0f\n", sum(bloc$n_dead_iom)))

# ── 3. 2-bloc (AFR vs EU): NegBin + Poisson ───────────────────
cat("\n--- 3. 2-bloc (AFR vs EU): NegBin + Poisson ---\n")

m_bloc_nb <- fenegbin(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou |
    month_year_fac + sar_bloc,
  data = bloc, vcov = NW(14), panel.id = ~sar_bloc + date
)

m_bloc_pois <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou |
    month_year_fac + sar_bloc,
  data = bloc, vcov = NW(14), panel.id = ~sar_bloc + date
)

print(etable(m_bloc_nb, m_bloc_pois, vcov = NW(14), se.below = TRUE,
             headers = c("NegBin", "Poisson")))

# ── 4. Extract b3 for summary and plot ────────────────────────
extract_b3 <- function(m, label) {
  ct <- coeftable(m, vcov = NW(14))
  r  <- grep(":post_mou$", rownames(ct))
  tibble(
    spec  = label,
    coef  = ct[r, 1],
    se    = ct[r, 2],
    p     = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])),
    ci_lo = ct[r, 1] - 1.96 * ct[r, 2],
    ci_hi = ct[r, 1] + 1.96 * ct[r, 2]
  )
}

summary_tbl <- bind_rows(
  extract_b3(m_bloc_nb,   "2-bloc (AFR/EU) - NegBin"),
  extract_b3(m_bloc_pois, "2-bloc (AFR/EU) - Poisson")
)

cat("\n  Summary (b3 = swh_prev5days:post_mou, NW(14) SEs):\n")
for (i in seq_len(nrow(summary_tbl))) {
  r <- summary_tbl[i, ]
  star <- if (r$p < 0.05) " *" else ""
  cat(sprintf("    %-30s  b3 = %+.3f (SE=%.3f)  IRR = %.3f  p = %.4f%s\n",
              r$spec, r$coef, r$se, exp(r$coef), r$p, star))
}

# ── 5. Save text output ───────────────────────────────────────
cat("\n--- 5. Saving results ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables", "21_zone_panel.txt")
sink(sink_file)

cat("051  ZONE-LEVEL REDUCED-FORM: AFR vs EU SAR\n")
cat("Aligned with 20_primary_model.R\n")
cat("===========================================\n")
cat("Spec: n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou |\n")
cat("      month_year_fac + sar_bloc\n")
cat("NegBin (fenegbin) AND Poisson QMLE (fepois), NW(14) SEs.\n\n")
cat(sprintf("Sample: %s to %s (N = %d days).\n",
            min(sample_dates), max(sample_dates), length(sample_dates)))
cat("Filter: !is.na(lc_lag14) & !is.na(swh_prev5days) (matches 05d exactly).\n")
cat("Outcome: n_dead_iom (zone-level deaths summed over 2 blocs,\n")
cat("         renamed from the zone panel's n_dead_missing).\n\n")

cat("=== 2-bloc (AFR vs EU) ===\n")
print(etable(m_bloc_nb, m_bloc_pois, vcov = NW(14), se.below = TRUE,
             headers = c("NegBin", "Poisson")))

cat("\n=== Summary (b3 = swh_prev5days:post_mou) ===\n")
for (i in seq_len(nrow(summary_tbl))) {
  r <- summary_tbl[i, ]
  star <- if (r$p < 0.05) " *" else ""
  cat(sprintf("  %-30s  b3 = %+.3f (SE=%.3f)  IRR = %.3f  p = %.4f%s\n",
              r$spec, r$coef, r$se, exp(r$coef), r$p, star))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
