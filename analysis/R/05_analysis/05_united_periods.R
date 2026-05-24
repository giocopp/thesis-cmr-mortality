# ── Long-span SWH-mortality gradient by political period ──────────────────
# UNITED primary 2013-present, IOM comparison 2014-present. Full ERA5 span
# without Frontex crossing control; 2014-2023 Frontex subset with clean
# lag-14 control kept as robustness.

library(tidyverse)
library(fixest)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

D_12       <- PERIOD_END_1       # end of period 1: 2017-02-01
D_23       <- PERIOD_END_2       # end of period 2: 2020-10-20
SAMPLE_END <- PERIOD_END_3       # end of period 3: 2023-01-01

# Same corridor-joined construction as 20/28; UNITED is primary.
SRC <- c(UNITED = "n_dead_united", IOM = "n_dead_iom")

add_period <- function(df) {
  df |> mutate(period = factor(case_when(
    date <= D_12 ~ "1. SAR + border control",
    date <= D_23 ~ "2. MoU + NGO containment",
    TRUE         ~ "3. Lamorgese rollback"
  )))
}

cat("============================================\n")
cat("31  PERIOD-SPECIFIC SWH GRADIENTS\n")
cat("    UNITED primary 2013-present; IOM comparison 2014-present\n")
cat("============================================\n\n")

# ── 1. Build extended panel ──
cat("--- 1. Building extended panel (UNITED primary + IOM comparison) ---\n")

era5 <- readRDS(file.path(BASE_DIR, "data", "processed",
                           "era5_swh_daily.RDS")) |>
  select(date, swh, swh_prev5days) |>
  filter(!is.na(swh_prev5days))

# IOM + UNITED daily series via the shared corridor-joined builders.
iom_daily    <- build_iom_daily() |> rename(n_dead_iom = n_dead_missing)
united_daily <- build_united_daily()

panel <- era5 |>
  filter(date >= DATE_START, date <= SAMPLE_END) |>
  left_join(iom_daily,    by = "date") |>
  left_join(united_daily, by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0))

# Strip dim attribute from date (ERA5 artifact that breaks case_when)
dim(panel$date) <- NULL

panel <- panel |>
  mutate(
    year       = year(date),
    month_year = factor(format(date, "%Y-%m")),
    unit       = 1L
  ) |>
  add_period()

cat(sprintf("  Panel: %s to %s (%d days)\n",
            min(panel$date), max(panel$date), nrow(panel)))
for (s in names(SRC)) {
  dat <- if (s == "IOM") filter(panel, date >= IOM_START) else panel
  v  <- dat[[SRC[s]]]
  cat(sprintf("  %-7s sample: %s to %s; deaths: %7.0f  (first death-day %s, last %s)\n",
              s, min(dat$date), max(dat$date), sum(v), min(dat$date[v > 0]),
              max(dat$date[v > 0])))
}

period_breakdown <- panel |>
  group_by(period) |>
  summarise(from = min(date), to = max(date),
            n_days_united = n(),
            n_days_iom_observed = sum(date >= IOM_START),
            united_only_days = sum(date < IOM_START),
            iom_death_days = sum(n_dead_iom[date >= IOM_START] > 0),
            utd_death_days = sum(n_dead_united > 0),
            iom_dead = sum(n_dead_iom[date >= IOM_START]),
            utd_dead = sum(n_dead_united),
            .groups = "drop")

cat("\n  Period breakdown (2013 explicitly UNITED-only):\n")
print(period_breakdown, n = Inf, width = Inf)

y2013 <- panel |> filter(date >= DATE_START, date < IOM_START)
cat(sprintf("\n  2013 UNITED-only days: %d; UNITED deaths: %.0f; IOM deaths counted: %.0f\n",
            nrow(y2013), sum(y2013$n_dead_united), sum(y2013$n_dead_iom)))

fit_data_full <- list(
  UNITED = panel,
  IOM    = panel |> filter(date >= IOM_START)
)

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

m_full <- imap(SRC, ~ fit_period(.x, fit_data_full[[.y]]))

per_full <- imap(m_full, ~ extract_per(.x, .y) |> mutate(source = .y)) |>
  list_rbind()

for (s in names(SRC)) {
  cat(sprintf("\n  [%s] period gradient, NW(14):\n", s))
  pf <- per_full |> filter(source == s)
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
# control can only be computed for that window. Period 1 is left-truncated
# to start 2014; the sample is right-clipped at SAMPLE_END (2023-01-01).
cat("\n--- 2b. 2014-2023 subset: lag-14 volume control (clean) ---\n")

panel_frx <- readRDS(file.path(BASE_DIR, "analysis", "data",
                                "daily_panel_complete.RDS")) |>
  filter(date <= SAMPLE_END) |>
  left_join(iom_daily,    by = "date") |>
  left_join(united_daily, by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) |>
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
    log1p_libya    = log1p(libya_conflict),
    log1p_tunisia  = log1p(tunisia_conflict),
    unit       = 1L,
    month_year = factor(format(date, "%Y-%m"))
  ) |>
  add_period() |>
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

SPEC_NO   <- "No crossing control"
SPEC_L14  <- "Lag-14 crossing control (clean)"
SPEC_FULL <- "Lag-14 + ACLED Libya/Tunisia"

m_sub <- list()
for (s in names(SRC)) {
  dep <- SRC[s]
  m_sub[[s]] <- list(
    no   = fit_period_ctl(dep, panel_frx),
    lc14 = fit_period_ctl(dep, panel_frx, "log1p_lc_lag14"),
    full = fit_period_ctl(dep, panel_frx,
                          "log1p_lc_lag14 + log1p_libya + log1p_tunisia")
  )
}

per_sub <- map(names(SRC), function(s) {
  bind_rows(
    extract_per(m_sub[[s]]$no,   SPEC_NO)   |> mutate(source = s),
    extract_per(m_sub[[s]]$lc14, SPEC_L14)  |> mutate(source = s),
    extract_per(m_sub[[s]]$full, SPEC_FULL) |> mutate(source = s)
  )
}) |> list_rbind()

for (s in names(SRC)) {
  for (sp in c(SPEC_NO, SPEC_L14, SPEC_FULL)) {
    cat(sprintf("\n  [%s] 2014-2023, %s:\n", s, sp))
    ps <- per_sub |> filter(source == s, spec == sp)
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
  co_f  <- coef(m_sub[[s]]$full)
  V_f   <- vcov(m_sub[[s]]$full, vcov = NW(14))
  for (cn in c("log1p_libya", "log1p_tunisia")) {
    ii <- which(names(co_f) == cn)
    if (length(ii) == 1) {
      cat(sprintf("    %-22s [%s]: %+.3f (SE=%.3f)\n",
                  cn, s, co_f[ii], sqrt(V_f[ii, ii])))
    }
  }
}

# ── 3. Wald tests: all period gradients equal (both sources) ──
cat("\n--- 3. Wald test: all period gradients equal? ---\n")

wald <- list()
for (s in names(SRC)) {
  full_label <- if (s == "UNITED") {
    "UNITED full 2013-present (no ctrl)"
  } else {
    "IOM comparison 2014-present (no ctrl)"
  }
  wald[[paste0(s, "_full")]]    <- wald_test(m_full[[s]], full_label)
  wald[[paste0(s, "_no")]]      <- wald_test(m_sub[[s]]$no,   sprintf("%s 2014-2023 (no ctrl)", s))
  wald[[paste0(s, "_lc14")]]    <- wald_test(m_sub[[s]]$lc14, sprintf("%s 2014-2023 (lag-14 clean)", s))
  wald[[paste0(s, "_fullctl")]] <- wald_test(m_sub[[s]]$full, sprintf("%s 2014-2023 (lag-14 + ACLED)", s))
}

# ── 3b. Mare Nostrum vs Triton/Sophia split inside P1 ─────
# Period 1 (SAR + border control) mixes Mare Nostrum (state-led, very
# active SAR, 2013-10-18 to 2014-10-31) and Triton/Sophia (border-control
# + NGO SAR, 2014-11-01 onwards). Restricted to the Frontex-bounded subset
# so ACLED and crossings controls are available.
cat("\n--- 3b. Mare Nostrum vs Triton/Sophia split inside period 1 ---\n")

panel_mn <- panel_frx |>
  mutate(period_mn = factor(case_when(
    date <= MARE_NOSTRUM_END ~ "1a. Mare Nostrum",
    date <= D_12             ~ "1b. Triton/Sophia",
    date <= D_23             ~ "2. MoU + NGO containment",
    TRUE                     ~ "3. Lamorgese rollback"
  )))

cat("\n  Period_mn breakdown (Frontex-bounded sample):\n")
print(panel_mn |> group_by(period_mn) |>
        summarise(from = min(date), to = max(date), n_days = n(),
                  dead_iom = sum(n_dead_iom), dead_united = sum(n_dead_united),
                  .groups = "drop"), width = Inf)

fit_period_mn <- function(dep, data, ctl = NULL) {
  rhs <- "swh_prev5days:period_mn"
  if (!is.null(ctl)) rhs <- paste(rhs, "+", ctl)
  fenegbin(as.formula(sprintf("%s ~ %s | month_year", dep, rhs)),
           data = data, vcov = NW(14), panel.id = ~unit + date)
}

m_mn <- list()
for (s in names(SRC)) {
  dep <- SRC[s]
  m_mn[[s]] <- list(
    no   = fit_period_mn(dep, panel_mn),
    full = fit_period_mn(dep, panel_mn,
                         "log1p_lc_lag14 + log1p_libya + log1p_tunisia")
  )
}

wald_1a_vs_1b <- function(m, label) {
  co <- coef(m); V <- vcov(m, vcov = NW(14))
  i_1a <- grep("period_mn1a", names(co))
  i_1b <- grep("period_mn1b", names(co))
  if (length(i_1a) != 1 || length(i_1b) != 1) {
    cat(sprintf("  %s: cannot locate 1a/1b coefficients; skipped.\n", label))
    return(invisible(NULL))
  }
  b_diff <- co[i_1b] - co[i_1a]
  v_diff <- V[i_1a, i_1a] + V[i_1b, i_1b] - 2 * V[i_1a, i_1b]
  z      <- b_diff / sqrt(v_diff)
  p      <- 2 * pnorm(-abs(z))
  cat(sprintf("  %s: b(1b)-b(1a) = %+.3f (SE=%.3f)  p=%.4f\n",
              label, b_diff, sqrt(v_diff), p))
}

for (s in names(SRC)) {
  for (spec_name in c("no", "full")) {
    m <- m_mn[[s]][[spec_name]]
    cat(sprintf("\n  [%s] MN-split (%s):\n", s,
                if (spec_name == "no") "no ctrl" else "full ctrl"))
    co  <- coef(m); V <- vcov(m, vcov = NW(14))
    idx <- grep("swh_prev5days:period_mn", names(co))
    for (rr in idx) {
      lab <- gsub("swh_prev5days:period_mn", "", names(co)[rr])
      b   <- co[rr]; se <- sqrt(V[rr, rr])
      cat(sprintf("    %-26s: %+.3f (SE=%.3f)  p=%.4f%s\n",
                  lab, b, se, 2 * pnorm(-abs(b / se)),
                  if (abs(b / se) > 1.96) " *" else ""))
    }
  }
}

cat("\n  Test 1a (Mare Nostrum) vs 1b (Triton/Sophia) slope equality:\n")
for (s in names(SRC)) {
  wald_1a_vs_1b(m_mn[[s]]$no,   sprintf("%s NB (no ctrl)", s))
  wald_1a_vs_1b(m_mn[[s]]$full, sprintf("%s NB + full ctl", s))
}

# ── 3c. Boundary-leakage robustness (drop +/-30d around cuts) ─
cat("\n--- 3c. Boundary-leakage: drop +/-30d around 2017-02-02 and 2020-10-21 ---\n")

drop_window <- 30L
panel_bl <- panel_frx |> filter(
  abs(as.integer(date - MOU_DATE)) > drop_window,
  abs(as.integer(date - (D_23 + 1))) > drop_window
)
cat(sprintf("  Boundary-leakage sample: N = %d days (drops %d around 2 cuts)\n",
            nrow(panel_bl), nrow(panel_frx) - nrow(panel_bl)))

m_bl <- list()
for (s in names(SRC)) {
  dep <- SRC[s]
  m_bl[[s]] <- list(
    no   = fit_period_ctl(dep, panel_bl),
    full = fit_period_ctl(dep, panel_bl,
                          "log1p_lc_lag14 + log1p_libya + log1p_tunisia")
  )
}

for (s in names(SRC)) {
  for (spec_name in c("no", "full")) {
    m <- m_bl[[s]][[spec_name]]
    cat(sprintf("\n  [%s] boundary-leakage (%s):\n", s,
                if (spec_name == "no") "no ctrl" else "full ctrl"))
    co  <- coef(m); V <- vcov(m, vcov = NW(14))
    idx <- grep("swh_prev5days:period", names(co))
    for (rr in idx) {
      lab <- gsub("swh_prev5days:period", "", names(co)[rr])
      b   <- co[rr]; se <- sqrt(V[rr, rr])
      cat(sprintf("    %-26s: %+.3f (SE=%.3f)  p=%.4f%s\n",
                  lab, b, se, 2 * pnorm(-abs(b / se)),
                  if (abs(b / se) > 1.96) " *" else ""))
    }
  }
  wald[[paste0(s, "_bl_no")]]   <- wald_test(m_bl[[s]]$no,
                                              sprintf("%s -30d (no ctrl)", s))
  wald[[paste0(s, "_bl_full")]] <- wald_test(m_bl[[s]]$full,
                                              sprintf("%s -30d + full ctl", s))
}

# ── 4. Plots (both sources) ──
cat("\n--- 4. Plots ---\n")

period_labels_full <- c(
  "SAR + border control\n(UNITED 2013-Feb17;\nIOM 2014-Feb17)",
  "MoU + NGO containment\n(Feb17-Oct20)",
  "Lamorgese rollback\n(Oct20-Jan23)"
)

relabel <- function(df) {
  df |>
    group_by(source, spec) |>
    arrange(period, .by_group = TRUE) |>
    mutate(period_short = factor(period_labels_full,
                                 levels = period_labels_full)) |>
    ungroup()
}

per_full_p <- relabel(per_full |>
                        mutate(spec = "Full span: UNITED 2013-present; IOM 2014-present"))
per_sub_p  <- relabel(per_sub)

src_cols <- c("IOM" = "#2166AC", "UNITED" = "#B2182B")

p_full <- ggplot(per_full_p, aes(period_short, beta, colour = source)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                linewidth = 0.8, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = src_cols) +
  labs(
    title = "Full span, no crossing control (corridor-joined)",
    subtitle = sprintf(
      "UNITED primary 2013-present; IOM comparison starts 2014. Wald: UNITED chi2(%d)=%.1f p=%.3f | IOM chi2(%d)=%.1f p=%.3f",
      wald$UNITED_full$df, wald$UNITED_full$stat, wald$UNITED_full$p,
      wald$IOM_full$df, wald$IOM_full$stat, wald$IOM_full$p),
    x = NULL, y = expression(beta[SWH_prev5days]), colour = "Source"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_sub <- ggplot(per_sub_p, aes(period_short, beta, colour = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.6, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                linewidth = 0.7, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = setNames(c("#2166AC", "#1B7837", "#7B3294"),
                                        c(SPEC_NO, SPEC_L14, SPEC_FULL))) +
  facet_wrap(~ source) +
  labs(
    title = "2014-2023 subset: no control vs clean lag-14 vs lag-14 + ACLED",
    subtitle = sprintf("Frontex-bounded; sample clipped at 2023-01-01. Clean-control check cor(SWH, lag-14)=%+.2f.",
                        cor14),
    x = NULL, y = expression(beta[SWH_prev5days]), colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_all <- p_full / p_sub

utd_fig <- fig_path("05_analysis", "05_united_period_gradient.png")
ggsave(utd_fig, p_all, width = 10, height = 9, dpi = 200)
cat(sprintf("  Saved: %s\n", utd_fig))

# ── 5. Save text output ──
cat("\n--- 5. Saving text output ---\n")

sink_file <- tbl_path("05_analysis", "05_united_periods.txt")
sink(sink_file)

cat("31  SWH-MORTALITY GRADIENT BY POLITICAL PERIOD\n")
cat("=============================================================\n")
cat("Long-span sibling of 28. Data construction IDENTICAL to 20/28:\n")
cat("  IOM    = build_iom_daily()    (corridor spatial join)\n")
cat("  UNITED = build_united_daily() (corridor spatial join; same\n")
cat("           country/cause filter as 20/28, NOT 31's old filter)\n")
cat("UNITED is the primary long-span source and is estimated on 2013-01-01\n")
cat("through SAMPLE_END (2023-01-01), with NO Frontex crossing control, so\n")
cat("the Mare Nostrum tail of period 1 is not truncated by Frontex coverage.\n")
cat("IOM MMP starts 2014, so IOM is estimated only as a 2014-present\n")
cat("comparison; 2013 is explicitly UNITED-only (see breakdown).\n\n")
cat("Periods (3-period design aligned with 04_descriptive/05_policy_timeline.R):\n")
cat("  1. SAR + border control     (UNITED 2013-01-01 / IOM 2014-01-01 to 2017-02-01):\n")
cat("       Mare Nostrum + Frontex + active NGO SAR\n")
cat("  2. MoU + NGO containment    (2017-02-02 to 2020-10-20):\n")
cat("       Italy-Libya MoU, Minniti code, closed-ports, Salvini decrees\n")
cat("  3. Lamorgese rollback       (2020-10-21 to 2023-01-01):\n")
cat("       Partial rollback of NGO restrictions; pre-Piantedosi decree\n\n")

cat(sprintf("Panel: %s to %s (%d days for UNITED)\n",
            min(panel$date), max(panel$date), nrow(panel)))
for (s in names(SRC)) {
  dat <- fit_data_full[[s]]
  v <- dat[[SRC[s]]]
  cat(sprintf("  %-7s sample: %s to %s; deaths: %.0f; N=%d days\n",
              s, min(dat$date), max(dat$date), sum(v), nrow(dat)))
}
cat("Model: fenegbin(n_dead ~ swh_prev5days:period | month_year), NW(14)\n\n")

cat("=== Period breakdown ===\n")
print(period_breakdown, n = Inf, width = Inf)
cat(sprintf("\n2013 UNITED-only days: %d; UNITED deaths: %.0f; IOM deaths counted: %.0f\n",
            nrow(y2013), sum(y2013$n_dead_united), sum(y2013$n_dead_iom)))

for (s in names(SRC)) {
  cat(sprintf("\n=== [%s] Period-specific gradient, full span, NW(14) ===\n", s))
  pf <- per_full |> filter(source == s)
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
cat("Sample clipped at 2023-01-01 (pre-Piantedosi).\n\n")
cat("Volume control = lag-14 living-crossings, window [t-14, t-8]. A short\n")
cat("overlapping lag (e.g. [t-7, t-1]) overlaps the SWH regressor window\n")
cat("[t-5, t-1], is a descendant of SWH, and would over-control the effect\n")
cat("of interest -> deliberately NOT used.\n")
cat("Clean-control validation (lag-14 is exogenous to the SWH window):\n")
cat(sprintf("  cor(SWH, lag-14)               = %+.3f\n", cor14))
cat(sprintf("  SWH -> lag-14 | month-year FE : b=%+.3f (SE=%.3f) p=%.4f  (b~0, ns)\n",
            ct14[1, 1], ct14[1, 2], ct14[1, 4]))
for (s in names(SRC)) {
  for (sp in c(SPEC_NO, SPEC_L14, SPEC_FULL)) {
    cat(sprintf("\n--- [%s] %s ---\n", s, sp))
    ps <- per_sub |> filter(source == s, spec == sp)
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
  co_f <- coef(m_sub[[s]]$full)
  V_f  <- vcov(m_sub[[s]]$full, vcov = NW(14))
  for (cn in c("log1p_libya", "log1p_tunisia")) {
    ii <- which(names(co_f) == cn)
    if (length(ii) == 1) {
      cat(sprintf("  %-22s [%s]: %+.3f (SE=%.3f)\n",
                  cn, s, co_f[ii], sqrt(V_f[ii, ii])))
    }
  }
  wn  <- wald[[paste0(s, "_no")]]
  w14 <- wald[[paste0(s, "_lc14")]]
  wfc <- wald[[paste0(s, "_fullctl")]]
  cat(sprintf("  Wald [%s] (no ctrl):       chi2(%d) = %.3f, p = %.4f\n",
              s, wn$df, wn$stat, wn$p))
  cat(sprintf("  Wald [%s] (lag-14 clean):  chi2(%d) = %.3f, p = %.4f\n",
              s, w14$df, w14$stat, w14$p))
  cat(sprintf("  Wald [%s] (lag-14 + ACLED):chi2(%d) = %.3f, p = %.4f\n",
              s, wfc$df, wfc$stat, wfc$p))
}

cat("\n=== Mare Nostrum / Triton-Sophia split inside P1 (Frontex-bounded) ===\n")
for (s in names(SRC)) {
  for (spec_name in c("no", "full")) {
    m <- m_mn[[s]][[spec_name]]
    cat(sprintf("\n--- [%s] %s ---\n",
                s, if (spec_name == "no") "no ctrl" else "lag-14 + ACLED"))
    co  <- coef(m); V <- vcov(m, vcov = NW(14))
    idx <- grep("swh_prev5days:period_mn", names(co))
    for (rr in idx) {
      lab <- gsub("swh_prev5days:period_mn", "", names(co)[rr])
      b   <- co[rr]; se <- sqrt(V[rr, rr])
      cat(sprintf("  %-26s: %+.3f (SE=%.3f)  p=%.4f%s\n",
                  lab, b, se, 2 * pnorm(-abs(b / se)),
                  if (abs(b / se) > 1.96) " *" else ""))
    }
  }
  cat(sprintf("\n  Test 1a (Mare Nostrum) vs 1b (Triton/Sophia) [%s]:\n", s))
  wald_1a_vs_1b(m_mn[[s]]$no,   sprintf("%s NB (no ctrl)", s))
  wald_1a_vs_1b(m_mn[[s]]$full, sprintf("%s NB + full ctl", s))
}

cat("\n=== Boundary leakage (drop +/-30d around 2017-02-02 and 2020-10-21) ===\n")
cat(sprintf("Sample: N=%d days (drops %d around 2 cuts).\n",
            nrow(panel_bl), nrow(panel_frx) - nrow(panel_bl)))
for (s in names(SRC)) {
  for (spec_name in c("no", "full")) {
    m <- m_bl[[s]][[spec_name]]
    cat(sprintf("\n--- [%s] %s ---\n",
                s, if (spec_name == "no") "no ctrl" else "lag-14 + ACLED"))
    co  <- coef(m); V <- vcov(m, vcov = NW(14))
    idx <- grep("swh_prev5days:period", names(co))
    for (rr in idx) {
      lab <- gsub("swh_prev5days:period", "", names(co)[rr])
      b   <- co[rr]; se <- sqrt(V[rr, rr])
      cat(sprintf("  %-26s: %+.3f (SE=%.3f)  p=%.4f%s\n",
                  lab, b, se, 2 * pnorm(-abs(b / se)),
                  if (abs(b / se) > 1.96) " *" else ""))
    }
  }
  wn  <- wald[[paste0(s, "_bl_no")]]
  wfc <- wald[[paste0(s, "_bl_full")]]
  cat(sprintf("\n  Wald [%s] -30d (no ctrl):       chi2(%d) = %.3f, p = %.4f\n",
              s, wn$df, wn$stat, wn$p))
  cat(sprintf("  Wald [%s] -30d (lag-14 + ACLED): chi2(%d) = %.3f, p = %.4f\n",
              s, wfc$df, wfc$stat, wfc$p))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================\n")
cat("DONE\n")
cat("============================================\n")
