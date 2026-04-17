# 05e_mechanism_interactions.R
# ============================
# Mechanism identification: does the SWH-mortality gradient vary with
# SAR capacity, boat composition, and boat size?
#
# Instead of a binary post_MoU moderator, we use CONTINUOUS daily measures
# of the operational environment as moderators. This tests the direct
# channels through which institutional change affects weather-related
# mortality.
#
# Moderators (all lagged weekly averages, days t-7 to t-1):
#   - SAR share:        frx_n_sar / frx_incidents
#   - Inflatable share: frx_n_inflatable / frx_incidents
#   - Persons per boat: frx_persons / frx_incidents
#
# Weekly averages are computed over RAW COUNTS then divided (not averages
# of daily shares) so that days with zero Frontex events don't produce
# undefined ratios.
#
# Specification:
#   deaths ~ SWH + SWH:sar_share_pw + SWH:inflatable_pw + SWH:ppb_pw | FE
#
# FE variants:
#   (a) year + month-of-year
#   (b) month-year
#
# In:  analysis/data/daily_panel_complete.RDS
# Out: output/tables/05e_mechanism_interactions.txt
#      output/figures/05e_mechanism_coefplot.png

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("05e  MECHANISM INTERACTIONS: SWH x SAR / BOATS / CROWDING\n")
cat("============================================================\n\n")

# ── 1. Load data and build weekly-lagged composition controls ─
cat("--- 1. Loading data and building controls ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

iom_daily <- build_iom_daily()

panel <- panel %>%
  left_join(iom_daily %>% rename(n_dead_iom = n_dead_missing), by = "date") %>%
  replace_na(list(n_dead_iom = 0)) %>%
  arrange(date) %>%
  mutate(
    # Weekly lagged rolling sums (days t-7 to t-1) of raw counts
    # Using sums over the window, then dividing, so zero-event days
    # don't produce NA ratios
    sar_events_pw     = dplyr::lag(zoo::rollsumr(frx_n_sar, k = 7, fill = NA), 1),
    inflatable_pw     = dplyr::lag(zoo::rollsumr(frx_n_inflatable, k = 7, fill = NA), 1),
    incidents_pw      = dplyr::lag(zoo::rollsumr(frx_incidents, k = 7, fill = NA), 1),
    persons_pw        = dplyr::lag(zoo::rollsumr(frx_persons, k = 7, fill = NA), 1),
    # Weekly shares (undefined if no incidents in the window)
    sar_share_pw      = ifelse(incidents_pw > 0, sar_events_pw / incidents_pw, NA_real_),
    inflatable_share_pw = ifelse(incidents_pw > 0, inflatable_pw / incidents_pw, NA_real_),
    ppb_pw            = ifelse(incidents_pw > 0, persons_pw / incidents_pw, NA_real_),
    # Standardize for interpretability (1-SD change)
    sar_share_pw_z    = (sar_share_pw - mean(sar_share_pw, na.rm = TRUE)) /
                         sd(sar_share_pw, na.rm = TRUE),
    inflatable_share_pw_z = (inflatable_share_pw - mean(inflatable_share_pw, na.rm = TRUE)) /
                             sd(inflatable_share_pw, na.rm = TRUE),
    ppb_pw_z          = (ppb_pw - mean(ppb_pw, na.rm = TRUE)) /
                         sd(ppb_pw, na.rm = TRUE),
    unit              = 1L,
    year_fac          = factor(year),
    month_of_year     = factor(month(date)),
    month_year_fac    = factor(month_year)
  )

d <- panel %>%
  filter(!is.na(sar_share_pw), !is.na(swh_prevweek))

cat(sprintf("  N = %d days (of %d; lost %d to lagged window NAs)\n",
            nrow(d), nrow(panel), nrow(panel) - nrow(d)))
cat(sprintf("  Deaths: %.0f over %d death-days\n",
            sum(d$n_dead_iom), sum(d$n_dead_iom > 0)))

cat("\n  Weekly-lagged composition controls (summary):\n")
cat(sprintf("    SAR share:        mean=%.3f  sd=%.3f  range=[%.2f, %.2f]\n",
            mean(d$sar_share_pw), sd(d$sar_share_pw),
            min(d$sar_share_pw), max(d$sar_share_pw)))
cat(sprintf("    Inflatable share: mean=%.3f  sd=%.3f  range=[%.2f, %.2f]\n",
            mean(d$inflatable_share_pw), sd(d$inflatable_share_pw),
            min(d$inflatable_share_pw), max(d$inflatable_share_pw)))
cat(sprintf("    Persons/boat:     mean=%.1f   sd=%.1f   range=[%.0f, %.0f]\n",
            mean(d$ppb_pw), sd(d$ppb_pw),
            min(d$ppb_pw), max(d$ppb_pw)))

# ── 2. Models: SWH interacted with each mechanism ────────────
cat("\n--- 2. Estimation ---\n")

# === (a) year + month FE ===
cat("\n  === year + month-of-year FE ===\n")

# Baseline (no moderator)
m_ym_base <- fenegbin(
  n_dead_iom ~ swh_prevweek | year_fac + month_of_year,
  data = d, vcov = NW(14), panel.id = ~unit + date)

# SAR share interaction only
m_ym_sar <- fenegbin(
  n_dead_iom ~ swh_prevweek + swh_prevweek:sar_share_pw_z +
    sar_share_pw_z | year_fac + month_of_year,
  data = d, vcov = NW(14), panel.id = ~unit + date)

# Inflatable share interaction only
m_ym_infl <- fenegbin(
  n_dead_iom ~ swh_prevweek + swh_prevweek:inflatable_share_pw_z +
    inflatable_share_pw_z | year_fac + month_of_year,
  data = d, vcov = NW(14), panel.id = ~unit + date)

# Persons per boat interaction only
m_ym_ppb <- fenegbin(
  n_dead_iom ~ swh_prevweek + swh_prevweek:ppb_pw_z +
    ppb_pw_z | year_fac + month_of_year,
  data = d, vcov = NW(14), panel.id = ~unit + date)

# All three together
m_ym_all <- fenegbin(
  n_dead_iom ~ swh_prevweek +
    swh_prevweek:sar_share_pw_z + sar_share_pw_z +
    swh_prevweek:inflatable_share_pw_z + inflatable_share_pw_z +
    swh_prevweek:ppb_pw_z + ppb_pw_z |
    year_fac + month_of_year,
  data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  NegBin, year + month FE, NW(14):\n")
print(etable(m_ym_base, m_ym_sar, m_ym_infl, m_ym_ppb, m_ym_all,
             vcov = NW(14), se.below = TRUE,
             headers = c("Baseline", "SAR", "Inflatable", "PPB", "All")))

# === (b) month-year FE ===
cat("\n  === month-year FE ===\n")

m_my_base <- fenegbin(
  n_dead_iom ~ swh_prevweek | month_year_fac,
  data = d, vcov = NW(14), panel.id = ~unit + date)

m_my_sar <- fenegbin(
  n_dead_iom ~ swh_prevweek + swh_prevweek:sar_share_pw_z +
    sar_share_pw_z | month_year_fac,
  data = d, vcov = NW(14), panel.id = ~unit + date)

m_my_infl <- fenegbin(
  n_dead_iom ~ swh_prevweek + swh_prevweek:inflatable_share_pw_z +
    inflatable_share_pw_z | month_year_fac,
  data = d, vcov = NW(14), panel.id = ~unit + date)

m_my_ppb <- fenegbin(
  n_dead_iom ~ swh_prevweek + swh_prevweek:ppb_pw_z +
    ppb_pw_z | month_year_fac,
  data = d, vcov = NW(14), panel.id = ~unit + date)

m_my_all <- fenegbin(
  n_dead_iom ~ swh_prevweek +
    swh_prevweek:sar_share_pw_z + sar_share_pw_z +
    swh_prevweek:inflatable_share_pw_z + inflatable_share_pw_z +
    swh_prevweek:ppb_pw_z + ppb_pw_z |
    month_year_fac,
  data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  NegBin, month-year FE, NW(14):\n")
print(etable(m_my_base, m_my_sar, m_my_infl, m_my_ppb, m_my_all,
             vcov = NW(14), se.below = TRUE,
             headers = c("Baseline", "SAR", "Inflatable", "PPB", "All")))

# === (c) Poisson robustness for the key specs ===
cat("\n  === Poisson robustness ===\n")

p_ym_sar <- fepois(
  n_dead_iom ~ swh_prevweek + swh_prevweek:sar_share_pw_z +
    sar_share_pw_z | year_fac + month_of_year,
  data = d, vcov = NW(14), panel.id = ~unit + date)

p_my_sar <- fepois(
  n_dead_iom ~ swh_prevweek + swh_prevweek:sar_share_pw_z +
    sar_share_pw_z | month_year_fac,
  data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  Poisson, SAR share interaction, NW(14):\n")
print(etable(p_ym_sar, p_my_sar, vcov = NW(14), se.below = TRUE,
             headers = c("year+month", "month-year")))

# ── 3. Summary table ─────────────────────────────────────────
cat("\n--- 3. Summary ---\n")

extract_int <- function(m, pattern, label) {
  ct <- coeftable(m, vcov = NW(14))
  r <- grep(pattern, rownames(ct))
  if (length(r) == 0) return(NULL)
  p <- 2 * pnorm(-abs(ct[r, 1] / ct[r, 2]))
  tibble(spec = label, coef = ct[r, 1], se = ct[r, 2], p = p)
}

summary_rows <- bind_rows(
  extract_int(m_ym_sar,  "swh_prevweek:sar_share",  "SAR share, yr+mo FE"),
  extract_int(m_ym_infl, "swh_prevweek:inflatable",  "Inflatable share, yr+mo FE"),
  extract_int(m_ym_ppb,  "swh_prevweek:ppb",         "Persons/boat, yr+mo FE"),
  extract_int(m_my_sar,  "swh_prevweek:sar_share",   "SAR share, my FE"),
  extract_int(m_my_infl, "swh_prevweek:inflatable",   "Inflatable share, my FE"),
  extract_int(m_my_ppb,  "swh_prevweek:ppb",          "Persons/boat, my FE"),
  extract_int(p_ym_sar,  "swh_prevweek:sar_share",   "SAR share, yr+mo FE (Poisson)"),
  extract_int(p_my_sar,  "swh_prevweek:sar_share",   "SAR share, my FE (Poisson)")
)

cat("\n  SWH x mechanism interactions (standardized, 1-SD change):\n")
for (i in seq_len(nrow(summary_rows))) {
  r <- summary_rows[i, ]
  cat(sprintf("    %-40s  %+.3f (SE=%.3f)  p=%.4f%s\n",
              r$spec, r$coef, r$se, r$p,
              if (r$p < 0.05) " *" else ""))
}

# ── 4. Coefficient plot ──────────────────────────────────────
cat("\n--- 4. Plot ---\n")

plot_df <- summary_rows %>%
  mutate(
    ci_lo = coef - 1.96 * se,
    ci_hi = coef + 1.96 * se,
    spec  = factor(spec, levels = rev(spec))
  )

p <- ggplot(plot_df, aes(coef, spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, colour = "#2166AC") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, colour = "#2166AC") +
  labs(
    title = "Mechanism interactions: SWH x operational environment",
    subtitle = "NegBin (+ Poisson for SAR). Standardized moderators (1-SD). NW(14) SEs.",
    x = "Coefficient on SWH x moderator (per 1-SD moderator change)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(BASE_DIR, "output", "figures",
                  "05e_mechanism_coefplot.png"),
       p, width = 10, height = 5, dpi = 200)
cat("  Saved: output/figures/05e_mechanism_coefplot.png\n")

# ── 5. Save text output ──────────────────────────────────────
cat("\n--- 5. Saving results ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables",
                        "05e_mechanism_interactions.txt")
sink(sink_file)

cat("05e  MECHANISM INTERACTIONS: SWH x SAR / BOATS / CROWDING\n")
cat("=========================================================\n")
cat(sprintf("Sample: %s to %s (N = %d days)\n",
            min(d$date), max(d$date), nrow(d)))
cat("\nModerators are WEEKLY LAGGED (days t-7 to t-1) averages of:\n")
cat("  SAR share:        frx_n_sar / frx_incidents\n")
cat("  Inflatable share: frx_n_inflatable / frx_incidents\n")
cat("  Persons per boat: frx_persons / frx_incidents\n")
cat("All standardized (mean=0, sd=1). Coefficients are per 1-SD change.\n")
cat("Weekly sums computed over raw counts before dividing.\n\n")

cat("=== NegBin, year + month FE, NW(14) ===\n")
print(etable(m_ym_base, m_ym_sar, m_ym_infl, m_ym_ppb, m_ym_all,
             vcov = NW(14), se.below = TRUE,
             headers = c("Baseline", "SAR", "Inflatable", "PPB", "All")))

cat("\n=== NegBin, month-year FE, NW(14) ===\n")
print(etable(m_my_base, m_my_sar, m_my_infl, m_my_ppb, m_my_all,
             vcov = NW(14), se.below = TRUE,
             headers = c("Baseline", "SAR", "Inflatable", "PPB", "All")))

cat("\n=== Poisson, SAR share, NW(14) ===\n")
print(etable(p_ym_sar, p_my_sar, vcov = NW(14), se.below = TRUE,
             headers = c("year+month", "month-year")))

cat("\n=== SUMMARY ===\n")
for (i in seq_len(nrow(summary_rows))) {
  r <- summary_rows[i, ]
  cat(sprintf("  %-40s  %+.3f (SE=%.3f)  p=%.4f%s\n",
              r$spec, r$coef, r$se, r$p,
              if (r$p < 0.05) " *" else ""))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
