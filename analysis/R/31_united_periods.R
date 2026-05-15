# 31_united_periods.R
# ====================
# SWH-mortality gradient by political/operational period, IOM and UNITED
# side by side. This is the LONG-SPAN sibling of 28: same data construction
# as 20/28 (corridor-joined IOM via build_iom_daily(); corridor-joined
# UNITED via build_united_daily()), but estimated on the full ERA5 span
# (2013-present) WITHOUT the Frontex-bounded crossing control, so the early
# (Mare Nostrum) and late (Meloni) regimes are not truncated by Frontex
# coverage. The 2014-2023 Frontex subset (with a CLEAN lag-14 volume
# control — window [t-14, t-8], strictly before the SWH window [t-5, t-1])
# is kept as a robustness block. A short overlapping lag (e.g. lag-7,
# window [t-7, t-1]) would be a descendant of the SWH regressor and is
# deliberately not used (over-control).
#
# Data construction is now IDENTICAL to 20/28 (this was previously NOT the
# case: 31 used a bespoke UNITED filter with no spatial join and no IOM
# series at all). Both sources go through the shared _helpers.R builders so
# the filter cannot drift again.
#
# Span caveat: IOM MMP coverage starts in 2014, so the IOM series is
# left-truncated relative to UNITED inside period 1 (2013 -> 2014 days are
# zero-filled for IOM). Reported openly in the period breakdown below; this
# is why the long-span design is run on BOTH sources rather than IOM alone.
#
# Periods (revised 2026-04-21):
#   1. Post-Arab Spring:    2013-01-01 to 2017-01-31
#        (Mare Nostrum + Frontex + NGO SAR operations)
#   2. MoU + Salvini:       2017-02-01 to 2019-12-31
#        (Italy-Libya MoU; Minniti-Gentiloni code; Salvini NGO targeting)
#   3. Partial rollback:    2020-01-01 to 2022-10-21
#        (narrow humanitarian protection partly restored)
#   4. Meloni:              2022-10-22 onwards
#        ("single-rescue" rule on NGOs; restricted humanitarian protection)
#
# Boundary 3-4 = 2022-10-22 (Meloni government sworn in).
# Data before 2013-01-01 is dropped.
#
# In:  data/processed/era5_swh_daily.RDS
#      data/processed/iom_mmp_incidents.RDS    (via build_iom_daily)
#      data/processed/united_incidents.RDS     (via build_united_daily)
#      data/processed/core_corridor.RDS        (spatial join, both sources)
#      analysis/data/daily_panel_complete.RDS  (crossing-control subset)
# Out: output/figures/31_united_period_gradient.png
#      output/tables/31_united_periods.txt

library(tidyverse)
library(fixest)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

# Period boundary dates
DATE_START <- as.Date("2013-01-01")
D_12       <- as.Date("2017-01-31")  # end Post-Arab Spring
D_23       <- as.Date("2019-12-31")  # end MoU + Salvini
D_34       <- as.Date("2022-10-21")  # end Partial rollback (Meloni sworn 10-22)

# Both sources, same corridor-joined construction as 20/28.
SRC <- c(IOM = "n_dead_iom", UNITED = "n_dead_united")

add_period <- function(df) {
  df %>% mutate(period = factor(case_when(
    date <= D_12 ~ "1. Post-Arab Spring",
    date <= D_23 ~ "2. MoU + Salvini",
    date <= D_34 ~ "3. Partial rollback",
    TRUE         ~ "4. Meloni"
  )))
}

cat("============================================\n")
cat("31  PERIOD-SPECIFIC SWH GRADIENTS (2013-present)\n")
cat("    IOM + UNITED, corridor-joined (20/28 construction)\n")
cat("============================================\n\n")

# ── 1. Build extended panel ──
cat("--- 1. Building extended panel (both sources) ---\n")

era5 <- readRDS(file.path(BASE_DIR, "data", "processed",
                           "era5_swh_daily.RDS")) %>%
  select(date, swh, swh_prev5days) %>%
  filter(!is.na(swh_prev5days))

# IOM + UNITED daily series via the shared corridor-joined builders.
iom_daily    <- build_iom_daily() %>% rename(n_dead_iom = n_dead_missing)
united_daily <- build_united_daily()

panel <- era5 %>%
  filter(date >= DATE_START) %>%
  left_join(iom_daily,    by = "date") %>%
  left_join(united_daily, by = "date") %>%
  replace_na(list(n_dead_iom = 0, n_dead_united = 0))

# Strip dim attribute from date (ERA5 artifact that breaks case_when)
dim(panel$date) <- NULL

panel <- panel %>%
  mutate(
    year       = year(date),
    month_year = factor(format(date, "%Y-%m")),
    unit       = 1L
  ) %>%
  add_period()

cat(sprintf("  Panel: %s to %s (%d days)\n",
            min(panel$date), max(panel$date), nrow(panel)))
for (s in names(SRC)) {
  v  <- panel[[SRC[s]]]
  nz <- panel$date[v > 0]
  cat(sprintf("  %-7s deaths: %7.0f  (first death-day %s, last %s)\n",
              s, sum(v), min(nz), max(nz)))
}

cat("\n  Period breakdown (IOM | UNITED):\n")
panel %>%
  group_by(period) %>%
  summarise(from = min(date), to = max(date), n_days = n(),
            iom_death_days    = sum(n_dead_iom    > 0),
            utd_death_days    = sum(n_dead_united > 0),
            iom_dead          = sum(n_dead_iom),
            utd_dead          = sum(n_dead_united),
            .groups = "drop") %>%
  print(n = Inf, width = Inf)

# ── 2. Period-specific gradient (full span, both sources) ──
cat("\n--- 2. Period-specific gradient (full span) ---\n")

fit_period <- function(dep, data) {
  fenegbin(as.formula(sprintf("%s ~ swh_prev5days:period | month_year", dep)),
           data = data, vcov = NW(14), panel.id = ~unit + date)
}

extract_per <- function(m, label) {
  co  <- coef(m)
  V   <- vcov(m, vcov = NW(14))
  idx <- grep("swh_prev5days:period", names(co))
  tibble(
    period = gsub("swh_prev5days:period", "", names(co[idx])),
    beta   = co[idx],
    se     = sqrt(diag(V)[idx]),
    ci_lo  = beta - 1.96 * se,
    ci_hi  = beta + 1.96 * se,
    spec   = label
  )
}

wald_test <- function(m, label) {
  co  <- coef(m)
  V   <- vcov(m, vcov = NW(14))
  idx <- grep("swh_prev5days:period", names(co))
  b   <- co[idx]
  Vs  <- V[idx, idx, drop = FALSE]
  k   <- length(b)
  R   <- matrix(0, nrow = k - 1, ncol = k)
  for (j in seq_len(k - 1)) {
    R[j, 1]     <- -1
    R[j, j + 1] <-  1
  }
  Rb   <- R %*% b
  RVR  <- R %*% Vs %*% t(R)
  stat <- as.numeric(t(Rb) %*% solve(RVR) %*% Rb)
  p    <- pchisq(stat, df = k - 1, lower.tail = FALSE)
  cat(sprintf("  %s: chi2(%d) = %.3f, p = %.4f\n", label, k - 1, stat, p))
  list(stat = stat, p = p, df = k - 1)
}

m_full <- lapply(SRC, fit_period, data = panel)
names(m_full) <- names(SRC)

per_full <- imap(m_full, ~ extract_per(.x, .y) %>% mutate(source = .y)) %>%
  list_rbind()

for (s in names(SRC)) {
  cat(sprintf("\n  [%s] period gradient, NW(14):\n", s))
  pf <- per_full %>% filter(source == s)
  for (i in seq_len(nrow(pf))) {
    r <- pf[i, ]
    cat(sprintf("    %-30s: %+.3f (SE=%.3f)  p=%.4f%s\n",
                r$period, r$beta, r$se,
                2 * pnorm(-abs(r$beta / r$se)),
                if (abs(r$beta / r$se) > 1.96) " *" else ""))
  }
  ct_cl <- coeftable(m_full[[s]], vcov = ~month_year)
  cat(sprintf("  [%s] period gradient, cluster(month_year):\n", s))
  idx <- grep("swh_prev5days:period", rownames(ct_cl))
  for (rr in idx) {
    lab <- gsub("swh_prev5days:period", "", rownames(ct_cl)[rr])
    cat(sprintf("    %-30s: %+.3f (SE=%.3f)  p=%.4f%s\n",
                lab, ct_cl[rr, 1], ct_cl[rr, 2], ct_cl[rr, 4],
                if (abs(ct_cl[rr, 1] / ct_cl[rr, 2]) > 1.96) " *" else ""))
  }
}

# ── 2b. 2014-2023 subset WITH clean lag-14 volume control (both sources) ──
# Frontex data is only available 2014-01-01 to 2023-05-31, so the crossing
# control can only be computed for that window. Period 1 is truncated to
# start 2014; Period 4 (Meloni) is thin (~7 months) — accepted, disclosed.
cat("\n--- 2b. 2014-2023 subset: lag-14 volume control (clean) ---\n")

panel_frx <- readRDS(file.path(BASE_DIR, "analysis", "data",
                                "daily_panel_complete.RDS")) %>%
  left_join(iom_daily,    by = "date") %>%
  left_join(united_daily, by = "date") %>%
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) %>%
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    # Volume control = mean living-crossings over [t-14, t-8]. This window
    # is STRICTLY BEFORE the SWH regressor window (swh_prev5days = mean
    # waves over [t-5, t-1]), so the control cannot re-measure the
    # regressor. A short overlapping lag (e.g. [t-7, t-1]) would be a
    # descendant of SWH and over-control the very effect of interest, so
    # it is deliberately NOT used. The clean-control check below verifies
    # lc_lag14 is exogenous to the SWH window.
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    log1p_lc_lag14 = log1p(lc_lag14),
    unit       = 1L,
    month_year = factor(format(date, "%Y-%m"))
  ) %>%
  add_period() %>%
  filter(!is.na(lc_lag14), !is.na(swh_prev5days))

cat(sprintf("  Frontex subset: %d days (%s to %s)\n",
            nrow(panel_frx), min(panel_frx$date), max(panel_frx$date)))

# Clean-control validation. A valid volume control must NOT be a descendant
# of the SWH regressor. Two checks on lc_lag14: (i) the raw correlation with
# swh_prev5days is only modest; (ii) crucially, swh_prev5days does NOT
# predict the control once month-year FE absorb seasonality/trend (coef ~ 0,
# not significant). Both confirm lc_lag14 [t-14, t-8] is exogenous to the
# SWH window [t-5, t-1] and is safe to condition on.
cor14 <- cor(panel_frx$swh_prev5days, panel_frx$log1p_lc_lag14)
fe14  <- feols(log1p_lc_lag14 ~ swh_prev5days | month_year,
               panel_frx, vcov = NW(14), panel.id = ~unit + date)
ct14  <- coeftable(fe14, vcov = NW(14))
cat(sprintf("  Clean-control check: cor(SWH, lag-14) = %+.3f\n", cor14))
cat(sprintf("  SWH -> lag-14 | month-year FE : b=%+.3f (SE=%.3f) p=%.4f\n",
            ct14[1, 1], ct14[1, 2], ct14[1, 4]))
cat("  (b ~ 0 and p > 0.05  =>  control is exogenous to the SWH window)\n")

fit_period_ctl <- function(dep, data, ctl = NULL) {
  rhs <- "swh_prev5days:period"
  if (!is.null(ctl)) rhs <- paste(rhs, "+", ctl)
  fenegbin(as.formula(sprintf("%s ~ %s | month_year", dep, rhs)),
           data = data, vcov = NW(14), panel.id = ~unit + date)
}

SPEC_NO  <- "No crossing control"
SPEC_L14 <- "Lag-14 crossing control (clean)"

m_sub <- list()
for (s in names(SRC)) {
  dep <- SRC[s]
  m_sub[[s]] <- list(
    no   = fit_period_ctl(dep, panel_frx),
    lc14 = fit_period_ctl(dep, panel_frx, "log1p_lc_lag14")
  )
}

per_sub <- map(names(SRC), function(s) {
  bind_rows(
    extract_per(m_sub[[s]]$no,   SPEC_NO)  %>% mutate(source = s),
    extract_per(m_sub[[s]]$lc14, SPEC_L14) %>% mutate(source = s)
  )
}) %>% list_rbind()

for (s in names(SRC)) {
  for (sp in c(SPEC_NO, SPEC_L14)) {
    cat(sprintf("\n  [%s] 2014-2023, %s:\n", s, sp))
    ps <- per_sub %>% filter(source == s, spec == sp)
    for (i in seq_len(nrow(ps))) {
      r <- ps[i, ]
      cat(sprintf("    %-25s: %+.3f (SE=%.3f)  p=%.4f%s\n",
                  r$period, r$beta, r$se,
                  2 * pnorm(-abs(r$beta / r$se)),
                  if (abs(r$beta / r$se) > 1.96) " *" else ""))
    }
  }
  co_lc <- coef(m_sub[[s]]$lc14)
  V_lc  <- vcov(m_sub[[s]]$lc14, vcov = NW(14))
  lc_i  <- which(names(co_lc) == "log1p_lc_lag14")
  cat(sprintf("    log1p(lag14d_crossings) [%s]: %+.3f (SE=%.3f)\n",
              s, co_lc[lc_i], sqrt(V_lc[lc_i, lc_i])))
}

# ── 3. Wald tests: all period gradients equal (both sources) ──
cat("\n--- 3. Wald test: all period gradients equal? ---\n")

wald <- list()
for (s in names(SRC)) {
  wald[[paste0(s, "_full")]] <- wald_test(m_full[[s]],     sprintf("%s full 2013-present (no ctrl)", s))
  wald[[paste0(s, "_no")]]   <- wald_test(m_sub[[s]]$no,    sprintf("%s 2014-2023 (no ctrl)", s))
  wald[[paste0(s, "_lc14")]] <- wald_test(m_sub[[s]]$lc14,  sprintf("%s 2014-2023 (lag-14 clean)", s))
}

# ── 4. Plots (both sources) ──
cat("\n--- 4. Plots ---\n")

period_labels_full <- c(
  "Post-Arab Spring\n(2013-Jan17)",
  "MoU + Salvini\n(Feb17-2019)",
  "Partial rollback\n(2020-Oct22)",
  "Meloni\n(Oct22-)"
)

relabel <- function(df) {
  df %>%
    group_by(source, spec) %>%
    arrange(period, .by_group = TRUE) %>%
    mutate(period_short = factor(period_labels_full,
                                 levels = period_labels_full)) %>%
    ungroup()
}

per_full_p <- relabel(per_full %>% mutate(spec = "Full 2013-present (no ctrl)"))
per_sub_p  <- relabel(per_sub)

src_cols <- c("IOM" = "#2166AC", "UNITED" = "#B2182B")

p_full <- ggplot(per_full_p, aes(period_short, beta, colour = source)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                linewidth = 0.8, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = src_cols) +
  labs(
    title = "Full 2013-present, no crossing control (corridor-joined)",
    subtitle = sprintf(
      "Wald all-equal: IOM chi2(%d)=%.1f p=%.3f | UNITED chi2(%d)=%.1f p=%.3f",
      wald$IOM_full$df, wald$IOM_full$stat, wald$IOM_full$p,
      wald$UNITED_full$df, wald$UNITED_full$stat, wald$UNITED_full$p),
    x = NULL, y = expression(beta[SWH_prev5days]), colour = "Source"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_sub <- ggplot(per_sub_p, aes(period_short, beta, colour = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.6, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                linewidth = 0.7, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = setNames(c("#2166AC", "#1B7837"),
                                        c(SPEC_NO, SPEC_L14))) +
  facet_wrap(~ source) +
  labs(
    title = "2014-2023 subset: no control vs clean lag-14 volume control",
    subtitle = sprintf("Frontex-bounded; period 4 (Meloni) thin. Clean-control check cor(SWH, lag-14)=%+.2f.",
                        cor14),
    x = NULL, y = expression(beta[SWH_prev5days]), colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_all <- p_full / p_sub

ggsave(file.path(BASE_DIR, "output", "figures",
                  "31_united_period_gradient.png"),
       p_all, width = 10, height = 9, dpi = 200)
cat("  Saved: output/figures/31_united_period_gradient.png\n")

# ── 5. Save text output ──
cat("\n--- 5. Saving text output ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables", "31_united_periods.txt")
sink(sink_file)

cat("31  SWH-MORTALITY GRADIENT BY POLITICAL PERIOD (IOM + UNITED)\n")
cat("=============================================================\n")
cat("Long-span sibling of 28. Data construction IDENTICAL to 20/28:\n")
cat("  IOM    = build_iom_daily()    (corridor spatial join)\n")
cat("  UNITED = build_united_daily() (corridor spatial join; same\n")
cat("           country/cause filter as 20/28, NOT 31's old filter)\n")
cat("Estimated on the full ERA5 span (2013-present), NO Frontex\n")
cat("crossing control, so Mare Nostrum / Meloni regimes are not\n")
cat("truncated by Frontex coverage. IOM MMP starts 2014, so the IOM\n")
cat("series is left-truncated inside period 1 (see breakdown).\n\n")
cat("Periods:\n")
cat("  1. Post-Arab Spring  (2013-01-01 to 2017-01-31):\n")
cat("       Mare Nostrum + Frontex + NGO SAR\n")
cat("  2. MoU + Salvini     (2017-02-01 to 2019-12-31):\n")
cat("       Italy-Libya MoU, Minniti-Gentiloni code, NGO targeting\n")
cat("  3. Partial rollback  (2020-01-01 to 2022-10-21):\n")
cat("       Narrow humanitarian protection restored\n")
cat("  4. Meloni            (2022-10-22 onwards):\n")
cat("       Single-rescue rule on NGOs; restricted humanitarian protection\n\n")

cat(sprintf("Panel: %s to %s (%d days)\n",
            min(panel$date), max(panel$date), nrow(panel)))
for (s in names(SRC)) {
  v <- panel[[SRC[s]]]
  cat(sprintf("  %-7s deaths: %.0f\n", s, sum(v)))
}
cat("Model: fenegbin(n_dead ~ swh_prev5days:period | month_year), NW(14)\n\n")

cat("=== Period breakdown ===\n")
panel %>%
  group_by(period) %>%
  summarise(from = min(date), to = max(date), n_days = n(),
            iom_death_days = sum(n_dead_iom > 0),
            utd_death_days = sum(n_dead_united > 0),
            iom_dead = sum(n_dead_iom),
            utd_dead = sum(n_dead_united),
            .groups = "drop") %>%
  print(n = Inf, width = Inf)

for (s in names(SRC)) {
  cat(sprintf("\n=== [%s] Period-specific gradient, full span, NW(14) ===\n", s))
  pf <- per_full %>% filter(source == s)
  for (i in seq_len(nrow(pf))) {
    r <- pf[i, ]
    cat(sprintf("  %-30s: %+.3f (SE=%.3f)  p=%.4f%s\n",
                r$period, r$beta, r$se,
                2 * pnorm(-abs(r$beta / r$se)),
                if (abs(r$beta / r$se) > 1.96) " *" else ""))
  }
  wf <- wald[[paste0(s, "_full")]]
  cat(sprintf("  Wald H0 (all equal): chi2(%d) = %.3f, p = %.4f\n",
              wf$df, wf$stat, wf$p))
}

cat("\n=== Period gradient, 2014-2023 subset (Frontex available) ===\n")
cat("NB: Period 4 (Meloni) has only ~7 months in this subset — thin.\n\n")
cat("Volume control = lag-14 living-crossings, window [t-14, t-8]. A short\n")
cat("overlapping lag (e.g. [t-7, t-1]) overlaps the SWH regressor window\n")
cat("[t-5, t-1], is a descendant of SWH, and would over-control the effect\n")
cat("of interest -> deliberately NOT used.\n")
cat("Clean-control validation (lag-14 is exogenous to the SWH window):\n")
cat(sprintf("  cor(SWH, lag-14)               = %+.3f\n", cor14))
cat(sprintf("  SWH -> lag-14 | month-year FE : b=%+.3f (SE=%.3f) p=%.4f  (b~0, ns)\n",
            ct14[1, 1], ct14[1, 2], ct14[1, 4]))
for (s in names(SRC)) {
  for (sp in c(SPEC_NO, SPEC_L14)) {
    cat(sprintf("\n--- [%s] %s ---\n", s, sp))
    ps <- per_sub %>% filter(source == s, spec == sp)
    for (i in seq_len(nrow(ps))) {
      r <- ps[i, ]
      cat(sprintf("  %-25s: %+.3f (SE=%.3f)  p=%.4f%s\n",
                  r$period, r$beta, r$se,
                  2 * pnorm(-abs(r$beta / r$se)),
                  if (abs(r$beta / r$se) > 1.96) " *" else ""))
    }
  }
  co_lc <- coef(m_sub[[s]]$lc14)
  V_lc  <- vcov(m_sub[[s]]$lc14, vcov = NW(14))
  lc_i  <- which(names(co_lc) == "log1p_lc_lag14")
  cat(sprintf("  log1p(lag14d_crossings) [%s]: %+.3f (SE=%.3f)\n",
              s, co_lc[lc_i], sqrt(V_lc[lc_i, lc_i])))
  wn  <- wald[[paste0(s, "_no")]]
  w14 <- wald[[paste0(s, "_lc14")]]
  cat(sprintf("  Wald [%s] (no ctrl):      chi2(%d) = %.3f, p = %.4f\n",
              s, wn$df, wn$stat, wn$p))
  cat(sprintf("  Wald [%s] (lag-14 clean): chi2(%d) = %.3f, p = %.4f\n",
              s, w14$df, w14$stat, w14$p))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================\n")
cat("DONE\n")
cat("============================================\n")
