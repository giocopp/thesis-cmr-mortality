# 05d_united_periods.R
# ====================
# UNITED data (2009-2025): SWH-mortality gradient by policy period.
#
# Periods:
#   1. Pre-Mare Nostrum:       2009-01-01 to 2013-10-17
#   2. Mare Nostrum:           2013-10-18 to 2014-10-31
#   3. NGO SAR:                2014-11-01 to 2017-01-31
#   4. EU-LCG Deal:           2017-02-01 to 2017-06-30
#   5. Libya SAR Zone:         2017-07-01 to 2025-12-31
#
# In:  data/processed/era5_swh_daily.RDS
#      data/processed/united_incidents.RDS
# Out: output/figures/05d_united_yearly_gradient.png
#      output/figures/05d_united_period_gradient.png
#      output/tables/05d_united_periods.txt

library(tidyverse)
library(fixest)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()

cat("============================================\n")
cat("UNITED: PERIOD-SPECIFIC SWH GRADIENTS (2009-2025)\n")
cat("============================================\n\n")

# ── 1. Build extended panel ──
cat("--- 1. Building extended panel ---\n")

era5 <- readRDS(file.path(BASE_DIR, "data", "processed",
                           "era5_swh_daily.RDS")) %>%
  select(date, swh, swh_prevweek) %>%
  filter(!is.na(swh_prevweek))

united_raw <- readRDS(file.path(BASE_DIR, "data", "processed",
                                 "united_incidents.RDS"))
CMR <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia", "Mediterranean")

united_daily <- united_raw %>%
  filter(country_of_death %in% CMR,
         (manner_of_death == "drowned" & !is.na(manner_of_death)) |
         (transport_means == "boat_ship_ferry" & !is.na(transport_means))) %>%
  group_by(date = incident_date_clean) %>%
  summarise(n_dead = sum(n_deaths, na.rm = TRUE), .groups = "drop")

panel <- era5 %>%
  filter(date >= as.Date("2009-01-01"), date <= as.Date("2025-12-31")) %>%
  left_join(united_daily, by = "date") %>%
  replace_na(list(n_dead = 0))

# Strip dim attribute from date (ERA5 artifact that breaks case_when)
dim(panel$date) <- NULL

panel <- panel %>%
  mutate(
    year       = year(date),
    year_fac   = factor(year),
    month_year = factor(format(date, "%Y-%m")),
    unit       = 1L,
    period = case_when(
      date < as.Date("2013-10-18")  ~ "1. Pre-Mare Nostrum",
      date <= as.Date("2014-10-31") ~ "2. Mare Nostrum",
      date <= as.Date("2017-01-31") ~ "3. NGO SAR",
      date <= as.Date("2017-06-30") ~ "4. EU-LCG Deal",
      TRUE                          ~ "5. Libya SAR Zone"
    ),
    period = factor(period)
  )

cat(sprintf("  Panel: %s to %s (%d days)\n",
            min(panel$date), max(panel$date), nrow(panel)))
cat(sprintf("  UNITED deaths: %.0f\n", sum(panel$n_dead)))

cat("\n  Period breakdown:\n")
panel %>%
  group_by(period) %>%
  summarise(from = min(date), to = max(date),
            n_days = n(), n_death_days = sum(n_dead > 0),
            total_dead = sum(n_dead), .groups = "drop") %>%
  print(n = Inf, width = Inf)

# ── 2. Year-by-year gradient ──
cat("\n--- 2. Year-by-year gradient ---\n")

m_yr <- fenegbin(n_dead ~ swh_prevweek:year_fac | month_year,
                 data = panel, vcov = NW(14), panel.id = ~unit + date)

co <- coef(m_yr)
V  <- vcov(m_yr, vcov = NW(14))
idx <- grep("swh_prevweek:year_fac", names(co))

yr_df <- tibble(
  year  = parse_number(names(co[idx])),
  beta  = co[idx],
  se    = sqrt(diag(V)[idx]),
  ci_lo = beta - 1.96 * se,
  ci_hi = beta + 1.96 * se
)

cat("\n  UNITED year-by-year SWH gradient (NegBin, month-year FE, NW(14)):\n")
for (i in seq_len(nrow(yr_df))) {
  r <- yr_df[i, ]
  cat(sprintf("    %d: %+.3f (SE=%.3f)%s\n",
              r$year, r$beta, r$se,
              if (abs(r$beta / r$se) > 1.96) " *" else ""))
}

# ── 3. Period-specific gradient ──
cat("\n--- 3. Period-specific gradient ---\n")

m_per <- fenegbin(n_dead ~ swh_prevweek:period | month_year,
                  data = panel, vcov = NW(14), panel.id = ~unit + date)

co_p <- coef(m_per)
V_p  <- vcov(m_per, vcov = NW(14))
idx_p <- grep("swh_prevweek:period", names(co_p))

per_df <- tibble(
  period = gsub("swh_prevweek:period", "", names(co_p[idx_p])),
  beta   = co_p[idx_p],
  se     = sqrt(diag(V_p)[idx_p]),
  ci_lo  = beta - 1.96 * se,
  ci_hi  = beta + 1.96 * se
)

cat("\n  NW(14) SEs:\n")
for (i in seq_len(nrow(per_df))) {
  r <- per_df[i, ]
  cat(sprintf("    %-30s: %+.3f (SE=%.3f)  p=%.4f%s\n",
              r$period, r$beta, r$se,
              2 * pnorm(-abs(r$beta / r$se)),
              if (abs(r$beta / r$se) > 1.96) " *" else ""))
}

ct_cl <- coeftable(m_per, vcov = ~month_year)
cat("\n  Cluster(month_year) SEs:\n")
for (nm in names(co_p[idx_p])) {
  r <- which(rownames(ct_cl) == nm)
  lab <- gsub("swh_prevweek:period", "", nm)
  cat(sprintf("    %-30s: %+.3f (SE=%.3f)  p=%.4f%s\n",
              lab, ct_cl[r, 1], ct_cl[r, 2], ct_cl[r, 4],
              if (abs(ct_cl[r, 1] / ct_cl[r, 2]) > 1.96) " *" else ""))
}

# ── 3b. Same periods on 2014-2023 subset WITH lagged crossing control ──
# Frontex data is only available 2014-2023, so the crossing control can
# only be computed for that window. We re-estimate the period model on
# this subset to show the effect of adding the control.
# Periods that start before 2014 are truncated; Pre-MN is dropped entirely.
cat("\n--- 3b. Period gradient WITH lagged crossing control (2014-2023) ---\n")

panel_frx <- readRDS(file.path(BASE_DIR, "analysis", "data",
                                "daily_panel_complete.RDS"))

# Merge UNITED deaths onto the Frontex panel
panel_frx <- panel_frx %>%
  left_join(united_daily, by = "date") %>%
  replace_na(list(n_dead = 0)) %>%
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag7 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 1),
    log1p_lc_lag7 = log1p(lc_lag7),
    unit       = 1L,
    month_year = factor(format(date, "%Y-%m")),
    period = case_when(
      date <= as.Date("2014-10-31") ~ "2. Mare Nostrum",
      date <= as.Date("2017-01-31") ~ "3. NGO SAR",
      date <= as.Date("2017-06-30") ~ "4. EU-LCG Deal",
      TRUE                          ~ "5. Libya SAR Zone"
    ),
    period = factor(period)
  ) %>%
  filter(!is.na(lc_lag7), !is.na(swh_prevweek))

cat(sprintf("  Frontex panel: %d days (%s to %s)\n",
            nrow(panel_frx), min(panel_frx$date), max(panel_frx$date)))

# Without crossing control (same subset for comparability)
m_per_no <- fenegbin(n_dead ~ swh_prevweek:period | month_year,
                     data = panel_frx, vcov = NW(14),
                     panel.id = ~unit + date)

# With lagged crossing control
m_per_lc <- fenegbin(n_dead ~ swh_prevweek:period + log1p_lc_lag7 | month_year,
                     data = panel_frx, vcov = NW(14),
                     panel.id = ~unit + date)

extract_per <- function(m, label) {
  co <- coef(m)
  V  <- vcov(m, vcov = NW(14))
  idx <- grep("swh_prevweek:period", names(co))
  df <- tibble(
    period = gsub("swh_prevweek:period", "", names(co[idx])),
    beta   = co[idx],
    se     = sqrt(diag(V)[idx]),
    ci_lo  = beta - 1.96 * se,
    ci_hi  = beta + 1.96 * se,
    spec   = label
  )
  df
}

per_no <- extract_per(m_per_no, "No crossing control")
per_lc <- extract_per(m_per_lc, "With lag 7d crossing control")

cat("\n  2014-2023, NO crossing control:\n")
for (i in seq_len(nrow(per_no))) {
  r <- per_no[i, ]
  cat(sprintf("    %-25s: %+.3f (SE=%.3f)  p=%.4f%s\n",
              r$period, r$beta, r$se,
              2 * pnorm(-abs(r$beta / r$se)),
              if (abs(r$beta / r$se) > 1.96) " *" else ""))
}

cat("\n  2014-2023, WITH lag 7d crossing control:\n")
for (i in seq_len(nrow(per_lc))) {
  r <- per_lc[i, ]
  cat(sprintf("    %-25s: %+.3f (SE=%.3f)  p=%.4f%s\n",
              r$period, r$beta, r$se,
              2 * pnorm(-abs(r$beta / r$se)),
              if (abs(r$beta / r$se) > 1.96) " *" else ""))
}

# Crossing control coefficient
co_lc <- coef(m_per_lc)
V_lc  <- vcov(m_per_lc, vcov = NW(14))
lc_i  <- which(names(co_lc) == "log1p_lc_lag7")
cat(sprintf("\n    log1p(lag7d_crossings): %+.3f (SE=%.3f)\n",
            co_lc[lc_i], sqrt(V_lc[lc_i, lc_i])))

# ── 4. Wald test ──
cat("\n--- 4. Wald test: all period gradients equal? ---\n")

wald_test <- function(m, label) {
  co <- coef(m)
  V  <- vcov(m, vcov = NW(14))
  idx <- grep("swh_prevweek:period", names(co))
  b <- co[idx]
  Vs <- V[idx, idx, drop = FALSE]
  k <- length(b)
  R <- matrix(0, nrow = k - 1, ncol = k)
  for (j in seq_len(k - 1)) {
    R[j, 1] <- -1
    R[j, j + 1] <- 1
  }
  Rb <- R %*% b
  RVR <- R %*% Vs %*% t(R)
  stat <- as.numeric(t(Rb) %*% solve(RVR) %*% Rb)
  p <- pchisq(stat, df = k - 1, lower.tail = FALSE)
  cat(sprintf("  %s: chi2(%d) = %.3f, p = %.4f\n", label, k - 1, stat, p))
  list(stat = stat, p = p, df = k - 1)
}

w_full <- wald_test(m_per, "Full 2009-2025 (no crossing control)")
w_no   <- wald_test(m_per_no, "2014-2023 (no crossing control)")
w_lc   <- wald_test(m_per_lc, "2014-2023 (with lag 7d crossing control)")

# Store for plot subtitle
wald_stat <- w_full$stat
wald_p    <- w_full$p

# ── 5. Plots ──
cat("\n--- 5. Plots ---\n")

period_breaks <- as.Date(c("2013-10-18", "2014-11-01",
                             "2017-02-01", "2017-07-01"))

p_yr <- ggplot(yr_df, aes(year, beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = decimal_date(period_breaks),
             linetype = "dotted", colour = "#D32F2F", linewidth = 0.4) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.15, fill = "grey40") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 2009:2025) +
  labs(
    title = "UNITED: year-by-year SWH-mortality gradient (2009-2025)",
    subtitle = "NegBin, month-year FE, NW(14) SEs. Dotted lines = policy regime changes.",
    x = NULL, y = expression(beta[SWH_prevweek])
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 9))

per_df <- per_df %>%
  mutate(
    period_short = c("Pre-MN\n(2009-Oct13)",
                     "Mare Nostrum\n(Oct13-Oct14)",
                     "NGO SAR\n(Nov14-Jan17)",
                     "EU-LCG Deal\n(Feb17-Jun17)",
                     "Libya SAR\n(Jul17-2025)"),
    period_short = factor(period_short, levels = period_short),
    spec = "Full 2009-2025 (no crossing control)"
  )

# Combined period plot: full sample + 2014-2023 with/without crossing control
period_labels <- c("Mare Nostrum\n(Oct13-Oct14)",
                   "NGO SAR\n(Nov14-Jan17)",
                   "EU-LCG Deal\n(Feb17-Jun17)",
                   "Libya SAR\n(Jul17-2023)")

per_no <- per_no %>%
  mutate(period_short = factor(period_labels, levels = period_labels))
per_lc <- per_lc %>%
  mutate(period_short = factor(period_labels, levels = period_labels))

per_compare <- bind_rows(per_no, per_lc)

p_per_full <- ggplot(per_df, aes(period_short, beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3.5, colour = "#2166AC") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.2, colour = "#2166AC", linewidth = 0.8) +
  labs(
    title = "Full 2009-2025, no crossing control",
    subtitle = sprintf("Wald test all equal: chi2(4)=%.1f, p=%.3f",
                       wald_stat, wald_p),
    x = NULL, y = expression(beta[SWH_prevweek])
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_per_compare <- ggplot(per_compare,
                         aes(period_short, beta, colour = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("No crossing control" = "#2166AC",
                                  "With lag 7d crossing control" = "#B2182B")) +
  labs(
    title = "2014-2023 subset: with vs without lagged crossing control",
    subtitle = sprintf("Wald (no ctrl): chi2(%d)=%.1f, p=%.3f | Wald (with ctrl): chi2(%d)=%.1f, p=%.3f",
                       w_no$df, w_no$stat, w_no$p,
                       w_lc$df, w_lc$stat, w_lc$p),
    x = NULL, y = expression(beta[SWH_prevweek]),
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_per <- p_per_full / p_per_compare

ggsave(file.path(BASE_DIR, "output", "figures",
                  "05d_united_yearly_gradient.png"),
       p_yr, width = 12, height = 5, dpi = 200)
cat("  Saved: output/figures/05d_united_yearly_gradient.png\n")

ggsave(file.path(BASE_DIR, "output", "figures",
                  "05d_united_period_gradient.png"),
       p_per, width = 10, height = 9, dpi = 200)
cat("  Saved: output/figures/05d_united_period_gradient.png\n")

ggsave(file.path(BASE_DIR, "output", "figures",
                  "05d_united_gradient_combined.png"),
       p_yr / p_per, width = 12, height = 14, dpi = 200)
cat("  Saved: output/figures/05d_united_gradient_combined.png\n")

# ── 6. Save text output ──
cat("\n--- 6. Saving text output ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables", "05d_united_periods.txt")
sink(sink_file)

cat("05d  UNITED: SWH-MORTALITY GRADIENT BY POLICY PERIOD (2009-2025)\n")
cat("================================================================\n")
cat(sprintf("Panel: %s to %s (%d days, %.0f UNITED deaths)\n",
            min(panel$date), max(panel$date), nrow(panel), sum(panel$n_dead)))
cat("Model: fenegbin(n_dead ~ swh_prevweek:period | month_year), NW(14)\n")
cat("No crossing control (Frontex unavailable before 2014).\n\n")

cat("=== Period breakdown ===\n")
panel %>%
  group_by(period) %>%
  summarise(from = min(date), to = max(date),
            n_days = n(), n_death_days = sum(n_dead > 0),
            total_dead = sum(n_dead), .groups = "drop") %>%
  print(n = Inf, width = Inf)

cat("\n=== Year-by-year gradient ===\n")
for (i in seq_len(nrow(yr_df))) {
  r <- yr_df[i, ]
  cat(sprintf("  %d: %+.3f (SE=%.3f)%s\n",
              r$year, r$beta, r$se,
              if (abs(r$beta / r$se) > 1.96) " *" else ""))
}

cat("\n=== Period-specific gradient, NW(14) ===\n")
for (i in seq_len(nrow(per_df))) {
  r <- per_df[i, ]
  cat(sprintf("  %-30s: %+.3f (SE=%.3f)  p=%.4f%s\n",
              r$period, r$beta, r$se,
              2 * pnorm(-abs(r$beta / r$se)),
              if (abs(r$beta / r$se) > 1.96) " *" else ""))
}

cat(sprintf("\nWald test H0 (all period gradients equal): chi2(4) = %.3f, p = %.4f\n",
            wald_stat, wald_p))

cat("\n=== Period gradient, 2014-2023 subset ===\n")
cat("(Frontex data available — can compute lagged crossing control)\n\n")

cat("--- No crossing control ---\n")
for (i in seq_len(nrow(per_no))) {
  r <- per_no[i, ]
  cat(sprintf("  %-25s: %+.3f (SE=%.3f)  p=%.4f%s\n",
              r$period, r$beta, r$se,
              2 * pnorm(-abs(r$beta / r$se)),
              if (abs(r$beta / r$se) > 1.96) " *" else ""))
}

cat("\n--- With lag 7d crossing control ---\n")
for (i in seq_len(nrow(per_lc))) {
  r <- per_lc[i, ]
  cat(sprintf("  %-25s: %+.3f (SE=%.3f)  p=%.4f%s\n",
              r$period, r$beta, r$se,
              2 * pnorm(-abs(r$beta / r$se)),
              if (abs(r$beta / r$se) > 1.96) " *" else ""))
}
cat(sprintf("\n  log1p(lag7d_crossings): %+.3f (SE=%.3f)\n",
            co_lc[lc_i], sqrt(V_lc[lc_i, lc_i])))

cat(sprintf("\nWald (no ctrl): chi2(%d) = %.3f, p = %.4f\n",
            w_no$df, w_no$stat, w_no$p))
cat(sprintf("Wald (with ctrl): chi2(%d) = %.3f, p = %.4f\n",
            w_lc$df, w_lc$stat, w_lc$p))

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================\n")
cat("DONE\n")
cat("============================================\n")
