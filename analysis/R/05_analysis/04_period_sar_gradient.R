# 28_period_sar_gradient.R
# =========================
# Regime-resolved SWH-mortality models, IOM and UNITED side by side.
#
# Consolidates three model families that the existing pipeline only covers
# partially:
#   F1  count, period interaction      (2-period + 4-period; fills IOM gap)
#   F2  count, continuous-SAR moderator (raw scale + UNITED; extends 23)
#   F3  deadly-day probability          (logit + LPM; new)
#
# Every model is estimated on BOTH n_dead_iom and n_dead_united, on the SAME
# sample, so the two sources are directly comparable.
#
# Date conventions (deliberately two different operationalisations of "the
# MoU" -- kept distinct, not conflated):
#   - 2-period cut  : post_mou = 1 from 2017-07-01 (MoU enforcement; the 20
#                     primary-model definition).
#   - 4-period cut  : boundary 1->2 at 2017-01-31 (MoU signed Feb 2017; the
#                     31_united_periods definition). Periods:
#       1 Post-Arab Spring  (.. 2017-01-31)  Mare Nostrum + Frontex + NGO SAR
#       2 MoU + Salvini      (2017-02-01 .. 2019-12-31)
#       3 Partial rollback   (2020-01-01 .. 2022-10-21)
#       4 Meloni             (2022-10-22 ..)
#
# Sample = the 20_primary_model sample (!is.na(lc_lag14) & !is.na(swh_prev5days)),
# so the F1 IOM 2-period gradient reconciles exactly with 20's b3. Because
# the panel is Frontex-bounded (~2014-01-15 to 2023-05-31), period 1 is
# left-truncated at 2014 and period 4 (Meloni) is thin (~7 months) -- as in
# 31's 2014-2023 subset. Stated openly in the output.
#
# UNITED filter = build_united_daily() defaults (country in CMR+Med, cause
# drowned/other_unknown, spatial join to core corridor). 20/27/28/31 all
# use the same shared builder, so UNITED numbers reconcile EXACTLY with 31
# on the shared (Frontex-bounded) sample.
#
# Scope: count/probability regime decomposition. The offset rate model is
# in 20_primary_model.R (boat-composition robustness in 27). Per-incident
# fatality rates are not computed (incident counts are not crossing exposure).
# No control group: estimates are descriptive regime-resolved, not causal.
#
# In:  analysis/data/daily_panel_complete.RDS
#      data/processed/iom_mmp_incidents.RDS (via build_iom_daily)
#      data/processed/united_incidents.RDS
#      data/processed/core_corridor.RDS
# Out: output/tables/28_period_sar_gradient.txt
#      output/figures/28_period_sar_gradient.png

library(tidyverse)
library(lubridate)
library(fixest)
library(patchwork)
library(sf)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")   # 2-period cut (20 definition)
D_12     <- as.Date("2017-01-31")   # 4-period: end Post-Arab Spring
D_23     <- as.Date("2019-12-31")   # 4-period: end MoU + Salvini
D_34     <- as.Date("2022-10-21")   # 4-period: end Partial rollback

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("28  REGIME-RESOLVED SWH-MORTALITY MODELS (IOM + UNITED)\n")
cat("    F1 period count | F2 continuous-SAR | F3 deadly-day prob\n")
cat("============================================================\n\n")

# ── 1. Load data + 20-style prep ─────────────────────────────
cat("--- 1. Loading data ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

# IOM daily via shared helper (default primary filter), as in 20.
iom_daily <- build_iom_daily()

# UNITED daily via the shared builder. Defaults (corridor spatial join;
# country in CMR+Med; manner drowned/other_unknown) replicate the previous
# inline filter exactly — single source of truth, see _helpers.R.
united_daily <- build_united_daily()

# Join + 20 sample prep + period factor + SAR moderators.
panel <- panel |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  left_join(united_daily, by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) |>
  arrange(date) |>
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    unit           = 1L,
    month_year_fac = factor(month_year),
    post_mou       = as.integer(date >= MOU_DATE),
    # 4-period factor (31 construction, verbatim).
    period = factor(case_when(
      date <= D_12 ~ "1. Post-Arab Spring",
      date <= D_23 ~ "2. MoU + Salvini",
      date <= D_34 ~ "3. Partial rollback",
      TRUE         ~ "4. Meloni"
    )),
    # Continuous SAR moderators (raw scale; 23 builds z-scored weekly).
    #   daily : frx_n_sar (count), frx_sar_share (= frx_n_sar/frx_incidents)
    #   weekly: 7-day rolling sum lagged 1 (zero-heavy daily -> weekly shown)
    sar_n_pw     = dplyr::lag(zoo::rollsumr(frx_n_sar,     k = 7, fill = NA), 1),
    sar_inc_pw   = dplyr::lag(zoo::rollsumr(frx_incidents, k = 7, fill = NA), 1),
    sar_share_pw = ifelse(sar_inc_pw > 0, sar_n_pw / sar_inc_pw, NA_real_)
  )

# Strip any dim attribute on date (ERA5 artefact) that breaks case_when /
# downstream — 31 does this defensively.
dim(panel$date) <- NULL

# 20 sample (makes F1 IOM 2-period reconcile exactly with 20's b3).
d <- panel |> filter(!is.na(lc_lag14), !is.na(swh_prev5days))

# I(deadly day) outcomes for F3.
d <- d |>
  mutate(dd_iom    = as.integer(n_dead_iom    > 0),
         dd_united = as.integer(n_dead_united > 0))

cat(sprintf("  Sample: %s to %s (N = %d days)\n",
            min(d$date), max(d$date), nrow(d)))
cat(sprintf("  Deaths IOM:    %.0f over %d death-days (%.1f%% deadly days)\n",
            sum(d$n_dead_iom), sum(d$dd_iom), 100 * mean(d$dd_iom)))
cat(sprintf("  Deaths UNITED: %.0f over %d death-days (%.1f%% deadly days)\n",
            sum(d$n_dead_united), sum(d$dd_united), 100 * mean(d$dd_united)))

cat("\n  Period breakdown (Frontex-bounded sample):\n")
print(d |> group_by(period) |>
        summarise(from = min(date), to = max(date), n_days = n(),
                  dead_iom = sum(n_dead_iom), dead_united = sum(n_dead_united),
                  .groups = "drop"), width = Inf)

cat(sprintf("\n  SAR moderators (sample):\n"))
cat(sprintf("    frx_n_sar (daily):  mean=%.2f  median=%.0f  max=%.0f  %%zero=%.1f\n",
            mean(d$frx_n_sar), median(d$frx_n_sar), max(d$frx_n_sar),
            100 * mean(d$frx_n_sar == 0)))
cat(sprintf("    frx_sar_share:      mean=%.3f  (NA on %d days w/o Frontex events)\n",
            mean(d$frx_sar_share, na.rm = TRUE), sum(is.na(d$frx_sar_share))))
cat(sprintf("    sar_n_pw (weekly):  mean=%.2f  median=%.0f  max=%.0f\n",
            mean(d$sar_n_pw, na.rm = TRUE), median(d$sar_n_pw, na.rm = TRUE),
            max(d$sar_n_pw, na.rm = TRUE)))

# ── 2. FAMILY 1: count, period interaction ───────────────────
cat("\n--- 2. FAMILY 1: count, period interaction ---\n")

# 2-period (= 20 primary spec): swh + swh:post_mou
f1_2p_iom_nb <- fenegbin(n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                         data = d, vcov = NW(14), panel.id = ~unit + date)
f1_2p_iom_po <- fepois  (n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                         data = d, vcov = NW(14), panel.id = ~unit + date)
f1_2p_utd_nb <- fenegbin(n_dead_united ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                         data = d, vcov = NW(14), panel.id = ~unit + date)
f1_2p_utd_po <- fepois  (n_dead_united ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                         data = d, vcov = NW(14), panel.id = ~unit + date)

# 4-period (31 construction): swh:period
f1_4p_iom_nb <- fenegbin(n_dead_iom ~ swh_prev5days:period | month_year_fac,
                         data = d, vcov = NW(14), panel.id = ~unit + date)
f1_4p_iom_po <- fepois  (n_dead_iom ~ swh_prev5days:period | month_year_fac,
                         data = d, vcov = NW(14), panel.id = ~unit + date)
f1_4p_utd_nb <- fenegbin(n_dead_united ~ swh_prev5days:period | month_year_fac,
                         data = d, vcov = NW(14), panel.id = ~unit + date)
f1_4p_utd_po <- fepois  (n_dead_united ~ swh_prev5days:period | month_year_fac,
                         data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  2-period (swh:post_mou), NW(14):\n")
print(etable(f1_2p_iom_nb, f1_2p_iom_po, f1_2p_utd_nb, f1_2p_utd_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM NB", "IOM Pois", "UNITED NB", "UNITED Pois")))

cat("\n  4-period (swh:period), NW(14):\n")
print(etable(f1_4p_iom_nb, f1_4p_iom_po, f1_4p_utd_nb, f1_4p_utd_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM NB", "IOM Pois", "UNITED NB", "UNITED Pois")))

# Wald test: all 4 period gradients equal (31's function, verbatim logic).
wald_periods <- function(m, label) {
  co <- coef(m); V <- vcov(m, vcov = NW(14))
  idx <- grep("swh_prev5days:period", names(co))
  b <- co[idx]; Vs <- V[idx, idx, drop = FALSE]; k <- length(b)
  R <- matrix(0, nrow = k - 1, ncol = k)
  for (j in seq_len(k - 1)) { R[j, 1] <- -1; R[j, j + 1] <- 1 }
  Rb <- R %*% b; RVR <- R %*% Vs %*% t(R)
  stat <- as.numeric(t(Rb) %*% solve(RVR) %*% Rb)
  p <- pchisq(stat, df = k - 1, lower.tail = FALSE)
  cat(sprintf("  %s: chi2(%d) = %.3f, p = %.4f\n", label, k - 1, stat, p))
  tibble(spec = label, chi2 = stat, df = k - 1, p = p)
}

cat("\n  Wald test (4 period gradients equal):\n")
wald_tbl <- bind_rows(
  wald_periods(f1_4p_iom_nb, "IOM NB"),
  wald_periods(f1_4p_iom_po, "IOM Pois"),
  wald_periods(f1_4p_utd_nb, "UNITED NB"),
  wald_periods(f1_4p_utd_po, "UNITED Pois")
)

# ── 3. FAMILY 2: count, continuous-SAR moderator ─────────────
cat("\n--- 3. FAMILY 2: count, continuous-SAR moderator ---\n")
# Spec: n_dead ~ swh + swh:SARmod + SARmod | month_year  (23 design, raw scale)
# Moderators: frx_n_sar (daily count), frx_sar_share (daily share),
#             sar_n_pw (weekly count), sar_share_pw (weekly share).

fit_f2 <- function(dep, mod) {
  f <- as.formula(sprintf(
    "%s ~ swh_prev5days + swh_prev5days:%s + %s | month_year_fac", dep, mod, mod))
  list(
    nb = fenegbin(f, data = d, vcov = NW(14), panel.id = ~unit + date),
    po = fepois  (f, data = d, vcov = NW(14), panel.id = ~unit + date)
  )
}

f2_mods <- c("frx_n_sar", "frx_sar_share", "sar_n_pw", "sar_share_pw")
f2_iom  <- setNames(lapply(f2_mods, function(m) fit_f2("n_dead_iom",    m)), f2_mods)
f2_utd  <- setNames(lapply(f2_mods, function(m) fit_f2("n_dead_united", m)), f2_mods)

extract_f2 <- function(fits, mod, source) {
  get1 <- function(m, fam) {
    ct <- coeftable(m, vcov = NW(14))
    r  <- grep(paste0("swh_prev5days:", mod), rownames(ct), fixed = TRUE)
    tibble(source = source, moderator = mod, family = fam,
           b_int = ct[r, 1], se_int = ct[r, 2],
           p_int = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
  }
  bind_rows(get1(fits$nb, "NegBin"), get1(fits$po, "Poisson"))
}

f2_tbl <- bind_rows(
  lapply(f2_mods, function(m) extract_f2(f2_iom[[m]], m, "IOM")),
  lapply(f2_mods, function(m) extract_f2(f2_utd[[m]], m, "UNITED"))
)

cat("\n  SWH x SAR-moderator interaction (b_int), NW(14):\n")
cat(sprintf("  %-7s %-14s %-8s %+10s %10s %10s\n",
            "source", "moderator", "family", "b_int", "SE", "p"))
for (i in seq_len(nrow(f2_tbl))) {
  r <- f2_tbl[i, ]
  star <- if (r$p_int < 0.05) " *" else ""
  cat(sprintf("  %-7s %-14s %-8s %+10.4f %10.4f %10.4f%s\n",
              r$source, r$moderator, r$family, r$b_int, r$se_int, r$p_int, star))
}

# ── 4. FAMILY 3: deadly-day probability ──────────────────────
cat("\n--- 4. FAMILY 3: deadly-day probability ---\n")
# Outcome I(n_dead>0). Logit (feglm) AND linear probability (feols).
# Period (2- and 4-) and continuous-SAR moderators.

# 4-period, LPM (clean predicted probabilities) + logit.
f3_4p_iom_lpm <- feols(dd_iom    ~ swh_prev5days:period | month_year_fac,
                       data = d, vcov = NW(14), panel.id = ~unit + date)
f3_4p_utd_lpm <- feols(dd_united ~ swh_prev5days:period | month_year_fac,
                       data = d, vcov = NW(14), panel.id = ~unit + date)
f3_4p_iom_lgt <- feglm(dd_iom    ~ swh_prev5days:period | month_year_fac,
                       data = d, family = binomial("logit"),
                       vcov = ~month_year_fac, panel.id = ~unit + date)
f3_4p_utd_lgt <- feglm(dd_united ~ swh_prev5days:period | month_year_fac,
                       data = d, family = binomial("logit"),
                       vcov = ~month_year_fac, panel.id = ~unit + date)

# 2-period, LPM.
f3_2p_iom_lpm <- feols(dd_iom    ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                       data = d, vcov = NW(14), panel.id = ~unit + date)
f3_2p_utd_lpm <- feols(dd_united ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                       data = d, vcov = NW(14), panel.id = ~unit + date)

# Continuous SAR (weekly share, the most stable moderator), LPM.
f3_sar_iom_lpm <- feols(dd_iom    ~ swh_prev5days + swh_prev5days:sar_share_pw + sar_share_pw | month_year_fac,
                        data = d, vcov = NW(14), panel.id = ~unit + date)
f3_sar_utd_lpm <- feols(dd_united ~ swh_prev5days + swh_prev5days:sar_share_pw + sar_share_pw | month_year_fac,
                        data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  Deadly-day, 4-period LPM, NW(14):\n")
print(etable(f3_4p_iom_lpm, f3_4p_utd_lpm, vcov = NW(14), se.below = TRUE,
             headers = c("IOM LPM", "UNITED LPM")))
cat("\n  Deadly-day, 4-period logit, cluster(month_year):\n")
print(etable(f3_4p_iom_lgt, f3_4p_utd_lgt, vcov = ~month_year_fac,
             se.below = TRUE, headers = c("IOM logit", "UNITED logit")))

# Descriptive predicted P(deadly day) at high SWH by period (LPM-based
# projection): baseline period mean P + period slope * (swh90 - period mean
# swh). Labelled as a descriptive projection, not a structural prediction.
swh90 <- quantile(d$swh_prev5days, 0.90, na.rm = TRUE)
pred_dd <- function(m, dep_dd, source) {
  ct <- coeftable(m, vcov = NW(14))
  base <- d |> group_by(period) |>
    summarise(p_bar = mean(.data[[dep_dd]]),
              swh_bar = mean(swh_prev5days), .groups = "drop")
  base$slope <- sapply(base$period, function(pp) {
    rn <- grep(paste0("swh_prev5days:period", pp), rownames(ct), fixed = TRUE)
    if (length(rn) == 0) NA_real_ else ct[rn, 1]
  })
  base |> mutate(source = source, swh90 = as.numeric(swh90),
                  p_at_swh90 = pmin(pmax(p_bar + slope * (swh90 - swh_bar), 0), 1))
}
pred_iom <- pred_dd(f3_4p_iom_lpm, "dd_iom",    "IOM")
pred_utd <- pred_dd(f3_4p_utd_lpm, "dd_united", "UNITED")

cat(sprintf("\n  Descriptive projected P(deadly day) at SWH=%.2f (90th pct) by period:\n",
            as.numeric(swh90)))
cat("  (LPM-based projection: period mean P + period slope x (swh90 - period mean swh))\n")
for (df in list(pred_iom, pred_utd)) {
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    cat(sprintf("    %-7s %-22s base P=%.3f  slope=%+.4f  -> P@swh90=%.3f\n",
                r$source, r$period, r$p_bar, r$slope, r$p_at_swh90))
  }
}

# ── 5. Figure ────────────────────────────────────────────────
cat("\n--- 5. Plot ---\n")

extract_period_betas <- function(m, source, family) {
  co <- coef(m); V <- vcov(m, vcov = NW(14))
  idx <- grep("swh_prev5days:period", names(co))
  tibble(period = gsub("swh_prev5days:period", "", names(co[idx])),
         beta = co[idx], se = sqrt(diag(V)[idx]),
         source = source, family = family) |>
    mutate(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)
}
p4_df <- bind_rows(
  extract_period_betas(f1_4p_iom_nb, "IOM", "NegBin"),
  extract_period_betas(f1_4p_iom_po, "IOM", "Poisson"),
  extract_period_betas(f1_4p_utd_nb, "UNITED", "NegBin"),
  extract_period_betas(f1_4p_utd_po, "UNITED", "Poisson")
)

p_panel1 <- ggplot(p4_df, aes(period, beta, colour = source)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.6, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  facet_wrap(~ family) +
  labs(title = "F1: SWH-mortality gradient by 4-period SAR regime",
       subtitle = "Count model, month-year FE, NW(14) 95% CI. Frontex-bounded sample.",
       x = NULL, y = expression(beta[SWH]), colour = "Source") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1))

p2_df <- f2_tbl |>
  mutate(ci_lo = b_int - 1.96 * se_int, ci_hi = b_int + 1.96 * se_int,
         lab = paste(moderator, family))
p_panel2 <- ggplot(p2_df, aes(b_int, lab, colour = source)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.6, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                width = 0.2, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  labs(title = "F2: SWH x SAR-moderator interaction (raw scale)",
       subtitle = "Count model, month-year FE, NW(14). Negative = more SAR flattens the SWH-death slope.",
       x = "b(SWH x SAR moderator)", y = NULL, colour = "Source") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

pred_all <- bind_rows(pred_iom, pred_utd)
p_panel3 <- ggplot(pred_all, aes(period, p_at_swh90, colour = source)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.4)) +
  geom_point(aes(y = p_bar), shape = 1, size = 2.6,
             position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  labs(title = "F3: deadly-day probability by period (LPM projection)",
       subtitle = sprintf("Filled = projected P at SWH=%.2f (90th pct); hollow = period mean P.",
                          as.numeric(swh90)),
       x = NULL, y = "P(deadly day)", colour = "Source") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1))

p_all <- p_panel1 / p_panel2 / p_panel3
per_fig <- fig_path("05_analysis", "04_period_sar_gradient.png")
ggsave(per_fig, p_all, width = 11, height = 14, dpi = 200)
cat(sprintf("  Saved: %s\n", per_fig))

# ── 6. Save text output ──────────────────────────────────────
cat("\n--- 6. Saving results ---\n")

sink_file <- tbl_path("05_analysis", "04_period_sar_gradient.txt")
sink(sink_file)

cat("28  REGIME-RESOLVED SWH-MORTALITY MODELS (IOM + UNITED)\n")
cat("=======================================================\n")
cat(sprintf("Sample: %s to %s (N = %d days; 20 primary sample)\n",
            min(d$date), max(d$date), nrow(d)))
cat("UNITED filter = build_united_daily() (corridor join); shared by 20/27/28/31.\n")
cat("2-period cut: post_mou from 2017-07-01. 4-period: 31 boundaries\n")
cat("(1->2 at 2017-01-31). Frontex-bounded -> period 1 left-truncated at\n")
cat("2014, period 4 (Meloni) thin (~7 months). Count/probability regime\n")
cat("decomposition; offset rate model is in 20 (boat controls in 27).\n")
cat("No control group -> descriptive.\n\n")

cat("=== FAMILY 1: count, 2-period (swh:post_mou) ===\n")
print(etable(f1_2p_iom_nb, f1_2p_iom_po, f1_2p_utd_nb, f1_2p_utd_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM NB", "IOM Pois", "UNITED NB", "UNITED Pois")))
cat("\n=== FAMILY 1: count, 4-period (swh:period) ===\n")
print(etable(f1_4p_iom_nb, f1_4p_iom_po, f1_4p_utd_nb, f1_4p_utd_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM NB", "IOM Pois", "UNITED NB", "UNITED Pois")))
cat("\nWald test (4 period gradients equal):\n")
for (i in seq_len(nrow(wald_tbl))) {
  r <- wald_tbl[i, ]
  cat(sprintf("  %-12s chi2(%d) = %.3f, p = %.4f\n",
              r$spec, r$df, r$chi2, r$p))
}

cat("\n=== FAMILY 2: count, continuous-SAR moderator (raw scale) ===\n")
cat("Spec: n_dead ~ swh + swh:SARmod + SARmod | month_year_fac (23 design)\n")
cat("Negative b_int = more SAR flattens the SWH-death slope (guardrail).\n\n")
cat(sprintf("  %-7s %-14s %-8s %+10s %10s %10s\n",
            "source", "moderator", "family", "b_int", "SE", "p"))
for (i in seq_len(nrow(f2_tbl))) {
  r <- f2_tbl[i, ]
  star <- if (r$p_int < 0.05) " *" else ""
  cat(sprintf("  %-7s %-14s %-8s %+10.4f %10.4f %10.4f%s\n",
              r$source, r$moderator, r$family, r$b_int, r$se_int, r$p_int, star))
}

cat("\n=== FAMILY 3: deadly-day probability ===\n")
cat("\n--- 4-period LPM, NW(14) ---\n")
print(etable(f3_4p_iom_lpm, f3_4p_utd_lpm, vcov = NW(14), se.below = TRUE,
             headers = c("IOM LPM", "UNITED LPM")))
cat("\n--- 4-period logit, cluster(month_year) ---\n")
print(etable(f3_4p_iom_lgt, f3_4p_utd_lgt, vcov = ~month_year_fac,
             se.below = TRUE, headers = c("IOM logit", "UNITED logit")))
cat("\n--- 2-period LPM, NW(14) ---\n")
print(etable(f3_2p_iom_lpm, f3_2p_utd_lpm, vcov = NW(14), se.below = TRUE,
             headers = c("IOM LPM", "UNITED LPM")))
cat("\n--- Continuous SAR (weekly share) LPM, NW(14) ---\n")
print(etable(f3_sar_iom_lpm, f3_sar_utd_lpm, vcov = NW(14), se.below = TRUE,
             headers = c("IOM LPM", "UNITED LPM")))

cat(sprintf("\n--- Descriptive projected P(deadly day) at SWH=%.2f (90th pct) ---\n",
            as.numeric(swh90)))
cat("LPM projection: period mean P + period slope x (swh90 - period mean swh).\n")
cat("Descriptive only (FE absorbed at sample mean), not a structural prediction.\n\n")
for (df in list(pred_iom, pred_utd)) {
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    cat(sprintf("  %-7s %-22s base P=%.3f  slope=%+.4f  P@swh90=%.3f\n",
                r$source, r$period, r$p_bar, r$slope, r$p_at_swh90))
  }
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
