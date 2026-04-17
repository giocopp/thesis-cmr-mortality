# 05d_primary_model.R
# ===================
# Primary reduced-form model and robustness.
#
# Specification
# -------------
# Primary (no crossing control):
#   deaths ~ swh_prevweek + swh_prevweek:post_mou | month_year
#
# Both NegBin and Poisson QMLE, NW(14) SEs.
# Poisson QMLE is consistent under weaker assumptions (only needs correct
# conditional mean); NegBin adds a variance assumption that improves
# efficiency. Agreement between the two validates the result.
#
# Also estimates year-by-year SWH gradients (swh_prevweek:year_fac) to
# show the trajectory is not a sharp break but a gradual shift.
#
# Robustness:
#   - Lagged crossing controls (lag 7d, lag 14d) as covariates
#   - Cluster(month_year) SEs alongside NW(14)
#
# In:  analysis/data/daily_panel_complete.RDS
#      data/processed/iom_mmp_incidents.RDS
#      data/processed/core_corridor.RDS
# Out: output/tables/05d_primary_model.txt
#      output/figures/05d_yearly_gradient.png

library(tidyverse)
library(lubridate)
library(fixest)
library(patchwork)

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

iom_daily <- build_iom_daily()

panel <- panel %>%
  left_join(iom_daily %>% rename(n_dead_iom = n_dead_missing), by = "date") %>%
  replace_na(list(n_dead_iom = 0)) %>%
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
d <- panel %>% filter(!is.na(lc_lag14), !is.na(swh_prevweek))

cat(sprintf("  Panel: %s to %s (%d days)\n",
            min(d$date), max(d$date), nrow(d)))
cat(sprintf("  Deaths (IOM primary): %.0f over %d death-days (%.1f%% zeros)\n",
            sum(d$n_dead_iom), sum(d$n_dead_iom > 0),
            100 * mean(d$n_dead_iom == 0)))

# ── 2. Primary model: NegBin and Poisson ─────────────────────
cat("\n--- 2. Primary model ---\n")

m_nb <- fenegbin(n_dead_iom ~ swh_prevweek + swh_prevweek:post_mou | month_year_fac,
                 data = d, vcov = NW(14), panel.id = ~unit + date)

m_pois <- fepois(n_dead_iom ~ swh_prevweek + swh_prevweek:post_mou | month_year_fac,
                 data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  Primary (no crossing control), NW(14):\n")
print(etable(m_nb, m_pois, vcov = NW(14), se.below = TRUE,
             headers = c("NegBin", "Poisson")))

cat("\n  Primary, cluster(month_year):\n")
print(etable(m_nb, m_pois, vcov = ~month_year_fac, se.below = TRUE,
             headers = c("NegBin", "Poisson")))

# ── 3. Robustness: lagged crossing controls ──────────────────
cat("\n--- 3. Robustness: lagged crossing controls ---\n")

m_nb_lag7  <- fenegbin(
  n_dead_iom ~ swh_prevweek + swh_prevweek:post_mou + log1p_lc_lag7 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_nb_lag14 <- fenegbin(
  n_dead_iom ~ swh_prevweek + swh_prevweek:post_mou + log1p_lc_lag14 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_pois_lag7  <- fepois(
  n_dead_iom ~ swh_prevweek + swh_prevweek:post_mou + log1p_lc_lag7 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

m_pois_lag14 <- fepois(
  n_dead_iom ~ swh_prevweek + swh_prevweek:post_mou + log1p_lc_lag14 |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  NegBin with crossing controls, NW(14):\n")
print(etable(m_nb, m_nb_lag7, m_nb_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))

cat("\n  Poisson with crossing controls, NW(14):\n")
print(etable(m_pois, m_pois_lag7, m_pois_lag14, vcov = NW(14), se.below = TRUE,
             headers = c("No control", "Lag 7d", "Lag 14d")))

# ── 4. Year-by-year SWH gradient ─────────────────────────────
cat("\n--- 4. Year-by-year SWH gradient ---\n")

extract_yearly <- function(m, label) {
  co <- coef(m)
  V  <- vcov(m, vcov = NW(14))
  idx <- grep("swh_prevweek:year_fac", names(co))
  tibble(
    year  = parse_number(names(co[idx])),
    beta  = co[idx],
    se    = sqrt(diag(V)[idx]),
    ci_lo = beta - 1.96 * se,
    ci_hi = beta + 1.96 * se,
    spec  = label
  )
}

# NegBin: no control and with lag 7d control
m_nb_yr <- fenegbin(n_dead_iom ~ swh_prevweek:year_fac | month_year_fac,
                    data = d, vcov = NW(14), panel.id = ~unit + date)

m_nb_yr_ctrl <- fenegbin(
  n_dead_iom ~ swh_prevweek:year_fac + log1p_lc_lag7 | month_year_fac,
  data = d, vcov = NW(14), panel.id = ~unit + date)

# Poisson: no control and with lag 7d control
m_pois_yr <- fepois(n_dead_iom ~ swh_prevweek:year_fac | month_year_fac,
                    data = d, vcov = NW(14), panel.id = ~unit + date)

m_pois_yr_ctrl <- fepois(
  n_dead_iom ~ swh_prevweek:year_fac + log1p_lc_lag7 | month_year_fac,
  data = d, vcov = NW(14), panel.id = ~unit + date)

yr_nb      <- extract_yearly(m_nb_yr,      "NegBin")
yr_nb_ctrl <- extract_yearly(m_nb_yr_ctrl, "NegBin + lag 7d control")
yr_pois    <- extract_yearly(m_pois_yr,    "Poisson")
yr_pois_ctrl <- extract_yearly(m_pois_yr_ctrl, "Poisson + lag 7d control")

cat("\n  NegBin, no control:\n")
walk(seq_len(nrow(yr_nb)), \(i) {
  r <- yr_nb[i, ]
  cat(sprintf("    %d: %+.3f (SE=%.3f)%s\n",
              r$year, r$beta, r$se, if (abs(r$beta/r$se) > 1.96) " *" else ""))
})

cat("\n  NegBin, with lag 7d crossing control:\n")
walk(seq_len(nrow(yr_nb_ctrl)), \(i) {
  r <- yr_nb_ctrl[i, ]
  cat(sprintf("    %d: %+.3f (SE=%.3f)%s\n",
              r$year, r$beta, r$se, if (abs(r$beta/r$se) > 1.96) " *" else ""))
})

cat("\n  Poisson, no control:\n")
walk(seq_len(nrow(yr_pois)), \(i) {
  r <- yr_pois[i, ]
  cat(sprintf("    %d: %+.3f (SE=%.3f)%s\n",
              r$year, r$beta, r$se, if (abs(r$beta/r$se) > 1.96) " *" else ""))
})

cat("\n  Poisson, with lag 7d crossing control:\n")
walk(seq_len(nrow(yr_pois_ctrl)), \(i) {
  r <- yr_pois_ctrl[i, ]
  cat(sprintf("    %d: %+.3f (SE=%.3f)%s\n",
              r$year, r$beta, r$se, if (abs(r$beta/r$se) > 1.96) " *" else ""))
})

# ── 5. Yearly gradient plots ─────────────────────────────────
cat("\n--- 5. Plots ---\n")

yr_all <- bind_rows(yr_nb, yr_nb_ctrl, yr_pois, yr_pois_ctrl) %>%
  mutate(
    model = ifelse(grepl("NegBin", spec), "NegBin", "Poisson"),
    control = ifelse(grepl("control", spec), "With lag 7d crossing control",
                     "No crossing control")
  )

p_yr <- ggplot(yr_all, aes(year, beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 2017.5, linetype = "dotted",
             colour = "#D32F2F", linewidth = 0.5) +
  annotate("text", x = 2017.7,
           y = max(yr_all$ci_hi, na.rm = TRUE) * 0.9,
           label = "MoU", colour = "#D32F2F", size = 3.5, hjust = 0) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "grey40") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2014:2023) +
  facet_grid(model ~ control) +
  labs(
    title = "Year-by-year SWH-mortality gradient: NegBin vs Poisson",
    subtitle = "Month-year FE, NW(14) SEs, 95% CI. IOM primary (corridor, sea causes).",
    x = NULL, y = expression(beta[SWH_prevweek])
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(BASE_DIR, "output", "figures", "05d_yearly_gradient.png"),
       p_yr, width = 12, height = 8, dpi = 200)
cat("  Saved: output/figures/05d_yearly_gradient.png\n")

# ── 6. Save text output ──────────────────────────────────────
cat("\n--- 6. Saving results ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables", "05d_primary_model.txt")
sink(sink_file)

cat("05d  PRIMARY MODEL: SWH x POST-MOU -> DEATHS\n")
cat("=============================================\n")
cat(sprintf("Sample: %s to %s (N = %d days)\n",
            min(d$date), max(d$date), nrow(d)))
cat(sprintf("Deaths (IOM primary — corridor, drowning+mixed): %.0f\n",
            sum(d$n_dead_iom)))
cat("\nPrimary specification:\n")
cat("  deaths ~ swh_prevweek + swh_prevweek:post_mou | month_year\n")
cat("  No crossing control. Month-year FE absorbs monthly confounders.\n")
cat("  NW(14) SEs for serial correlation.\n\n")

cat("=== PRIMARY MODEL ===\n\n")
cat("--- NW(14) SEs ---\n")
print(etable(m_nb, m_pois, vcov = NW(14), se.below = TRUE,
             headers = c("NegBin", "Poisson")))

cat("\n--- Cluster(month_year) SEs ---\n")
print(etable(m_nb, m_pois, vcov = ~month_year_fac, se.below = TRUE,
             headers = c("NegBin", "Poisson")))

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

cat("\n\n=== YEAR-BY-YEAR GRADIENT ===\n\n")

cat("NegBin, no control:\n")
walk(seq_len(nrow(yr_nb)), \(i) {
  r <- yr_nb[i, ]
  cat(sprintf("  %d: %+.3f (SE=%.3f)%s\n",
              r$year, r$beta, r$se, if (abs(r$beta/r$se) > 1.96) " *" else ""))
})

cat("\nNegBin, with lag 7d crossing control:\n")
walk(seq_len(nrow(yr_nb_ctrl)), \(i) {
  r <- yr_nb_ctrl[i, ]
  cat(sprintf("  %d: %+.3f (SE=%.3f)%s\n",
              r$year, r$beta, r$se, if (abs(r$beta/r$se) > 1.96) " *" else ""))
})

cat("\nPoisson, no control:\n")
walk(seq_len(nrow(yr_pois)), \(i) {
  r <- yr_pois[i, ]
  cat(sprintf("  %d: %+.3f (SE=%.3f)%s\n",
              r$year, r$beta, r$se, if (abs(r$beta/r$se) > 1.96) " *" else ""))
})

cat("\nPoisson, with lag 7d crossing control:\n")
walk(seq_len(nrow(yr_pois_ctrl)), \(i) {
  r <- yr_pois_ctrl[i, ]
  cat(sprintf("  %d: %+.3f (SE=%.3f)%s\n",
              r$year, r$beta, r$se, if (abs(r$beta/r$se) > 1.96) " *" else ""))
})

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
