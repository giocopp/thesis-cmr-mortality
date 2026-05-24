# ── Regime-resolved SWH-mortality gradient (IOM + UNITED) ─────────────────
# Three families: 2/3-period interaction, continuous SAR moderator,
# deadly-day probability. Sample matches the primary model.

library(tidyverse)
library(lubridate)
library(fixest)
library(patchwork)
library(sf)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

D_12        <- PERIOD_END_1     # end of period 1: 2017-02-01
D_23        <- PERIOD_END_2     # end of period 2: 2020-10-20
SAMPLE_END  <- PERIOD_END_3     # end of period 3: 2023-01-01

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
    # 3-period factor aligned with the policy timeline.
    period = factor(case_when(
      date <= D_12 ~ "1. SAR + border control",
      date <= D_23 ~ "2. MoU + NGO containment",
      TRUE         ~ "3. Lamorgese rollback"
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

# 20 sample (makes F1 IOM 2-period reconcile exactly with 20's b3), clipped
# at SAMPLE_END so the Piantedosi-era tail (post 2023-01-01) does not bleed
# into period 3. Build log1p transforms of the credibility-audit controls:
# lag-14 living-crossings (already used in 20 robustness) and ACLED conflict
# events for Libya + Tunisia (in panel from daily_panel_complete.RDS).
d <- panel |>
  filter(!is.na(lc_lag14), !is.na(swh_prev5days), date <= SAMPLE_END) |>
  mutate(
    log1p_lc_lag14 = log1p(lc_lag14),
    log1p_libya    = log1p(libya_conflict),
    log1p_tunisia  = log1p(tunisia_conflict)
  )

# I(deadly day) outcomes for F3.
d <- d |>
  mutate(dd_iom    = as.integer(n_dead_iom    > 0),
         dd_united = as.integer(n_dead_united > 0))

# Wider samples for the 3-period count decomposition (tab:period). The
# decomposition model has no crossing-volume control, so it does not
# require the Frontex-bounded daily_panel_complete sample (which is
# left-truncated at 2014-01-01). UNITED extends back to 2013-01-01 to
# take advantage of its longer recording span; IOM is left-floored at
# 2014-01-01 (when MMP coverage begins). Built directly from the
# ERA5 SWH series and the source-specific daily death builders so we
# are not bounded by the Frontex sample. Both end at SAMPLE_END.
swh_daily_full <- readRDS(file.path(BASE_DIR, "data", "processed",
                                     "era5_swh_daily.RDS"))
dim(swh_daily_full$date) <- NULL

build_period_sample <- function(deaths_daily, count_col, date_floor) {
  swh_daily_full |>
    dplyr::select(date, swh_prev5days) |>
    dplyr::left_join(deaths_daily |> dplyr::select(date, dplyr::all_of(count_col)),
                     by = "date") |>
    tidyr::replace_na(setNames(list(0), count_col)) |>
    dplyr::filter(date >= date_floor, date <= SAMPLE_END,
                  !is.na(swh_prev5days)) |>
    dplyr::mutate(
      unit           = 1L,
      month_year_fac = factor(format(date, "%Y-%m")),
      period = factor(dplyr::case_when(
        date <= D_12 ~ "1. SAR + border control",
        date <= D_23 ~ "2. MoU + NGO containment",
        TRUE         ~ "3. Lamorgese rollback"
      ))
    )
}

d_utd_3p <- build_period_sample(united_daily, "n_dead_united",
                                 as.Date("2013-01-01"))
d_iom_3p <- build_period_sample(iom_daily |> rename(n_dead_iom = n_dead_missing),
                                 "n_dead_iom",
                                 as.Date("2014-01-01"))

cat(sprintf("  3-period samples: UNITED N=%d (%s..%s), IOM N=%d (%s..%s)\n",
            nrow(d_utd_3p), min(d_utd_3p$date), max(d_utd_3p$date),
            nrow(d_iom_3p), min(d_iom_3p$date), max(d_iom_3p$date)))

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

# 3-period (policy-timeline cuts): swh:period. Fit on the source-specific
# wider samples so UNITED uses its 2013-onwards span and IOM uses
# 2014-onwards (see note above d_utd_3p / d_iom_3p).
f1_3p_iom_nb <- fenegbin(n_dead_iom ~ swh_prev5days:period | month_year_fac,
                         data = d_iom_3p, vcov = NW(14), panel.id = ~unit + date)
f1_3p_iom_po <- fepois  (n_dead_iom ~ swh_prev5days:period | month_year_fac,
                         data = d_iom_3p, vcov = NW(14), panel.id = ~unit + date)
f1_3p_utd_nb <- fenegbin(n_dead_united ~ swh_prev5days:period | month_year_fac,
                         data = d_utd_3p, vcov = NW(14), panel.id = ~unit + date)
f1_3p_utd_po <- fepois  (n_dead_united ~ swh_prev5days:period | month_year_fac,
                         data = d_utd_3p, vcov = NW(14), panel.id = ~unit + date)

cat("\n  2-period (swh:post_mou), NW(14):\n")
print(etable(f1_2p_iom_nb, f1_2p_iom_po, f1_2p_utd_nb, f1_2p_utd_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM NB", "IOM Pois", "UNITED NB", "UNITED Pois")))

cat("\n  3-period (swh:period), NW(14):\n")
print(etable(f1_3p_iom_nb, f1_3p_iom_po, f1_3p_utd_nb, f1_3p_utd_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM NB", "IOM Pois", "UNITED NB", "UNITED Pois")))

# Wald test: all 3 period gradients equal.
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

cat("\n  Wald test (3 period gradients equal):\n")
wald_tbl <- bind_rows(
  wald_periods(f1_3p_iom_nb, "IOM NB"),
  wald_periods(f1_3p_iom_po, "IOM Pois"),
  wald_periods(f1_3p_utd_nb, "UNITED NB"),
  wald_periods(f1_3p_utd_po, "UNITED Pois")
)

# ── 2b. FAMILY 1 with full controls (lag-14 + ACLED) ──────────
# Auxiliary to the 2-period primary model: tests whether the 3-period
# decomposition pattern survives ACLED Libya/Tunisia + lag-14 controls.
cat("\n--- 2b. FAMILY 1 with full controls (lag-14 + ACLED) ---\n")

f1_3p_iom_nb_ctl <- fenegbin(
  n_dead_iom ~ swh_prev5days:period + log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)
f1_3p_iom_po_ctl <- fepois(
  n_dead_iom ~ swh_prev5days:period + log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)
f1_3p_utd_nb_ctl <- fenegbin(
  n_dead_united ~ swh_prev5days:period + log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)
f1_3p_utd_po_ctl <- fepois(
  n_dead_united ~ swh_prev5days:period + log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d, vcov = NW(14), panel.id = ~unit + date)

cat("\n  3-period with full controls, NB, NW(14):\n")
print(etable(f1_3p_iom_nb, f1_3p_iom_nb_ctl, f1_3p_utd_nb, f1_3p_utd_nb_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM (no ctl)", "IOM + ctl",
                         "UNITED (no ctl)", "UNITED + ctl")))

cat("\n  Wald test (3 period gradients equal) with full controls:\n")
wald_tbl_ctl <- bind_rows(
  wald_periods(f1_3p_iom_nb_ctl, "IOM NB + ctl"),
  wald_periods(f1_3p_iom_po_ctl, "IOM Pois + ctl"),
  wald_periods(f1_3p_utd_nb_ctl, "UNITED NB + ctl"),
  wald_periods(f1_3p_utd_po_ctl, "UNITED Pois + ctl")
)

# ── 2c. Mare Nostrum vs Triton/Sophia split inside P1 ─────────
# Period 1 mixes Mare Nostrum (state-led, very active) and Triton/Sophia
# (border-control + NGO SAR). Splitting tests whether the negative beta_1
# is driven solely by MN, or also appears under Triton/Sophia.
cat("\n--- 2c. Mare Nostrum vs Triton/Sophia split inside period 1 ---\n")

d_mn <- d |>
  mutate(period_mn = factor(case_when(
    date <= MARE_NOSTRUM_END ~ "1a. Mare Nostrum",
    date <= D_12             ~ "1b. Triton/Sophia",
    date <= D_23             ~ "2. MoU + NGO containment",
    TRUE                     ~ "3. Lamorgese rollback"
  )))

cat("\n  Period_mn breakdown:\n")
print(d_mn |> group_by(period_mn) |>
        summarise(from = min(date), to = max(date), n_days = n(),
                  dead_iom = sum(n_dead_iom), dead_united = sum(n_dead_united),
                  .groups = "drop"), width = Inf)

f1_mn_iom_nb <- fenegbin(n_dead_iom    ~ swh_prev5days:period_mn | month_year_fac,
                         data = d_mn, vcov = NW(14), panel.id = ~unit + date)
f1_mn_utd_nb <- fenegbin(n_dead_united ~ swh_prev5days:period_mn | month_year_fac,
                         data = d_mn, vcov = NW(14), panel.id = ~unit + date)
f1_mn_iom_nb_ctl <- fenegbin(
  n_dead_iom ~ swh_prev5days:period_mn + log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d_mn, vcov = NW(14), panel.id = ~unit + date)
f1_mn_utd_nb_ctl <- fenegbin(
  n_dead_united ~ swh_prev5days:period_mn + log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d_mn, vcov = NW(14), panel.id = ~unit + date)

cat("\n  MN/Triton split, NB, NW(14):\n")
print(etable(f1_mn_iom_nb, f1_mn_iom_nb_ctl, f1_mn_utd_nb, f1_mn_utd_nb_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM (no ctl)", "IOM + ctl",
                         "UNITED (no ctl)", "UNITED + ctl")))

# Wald test: 1a vs 1b slopes equal?
wald_1a_vs_1b <- function(m, label) {
  co <- coef(m); V <- vcov(m, vcov = NW(14))
  i_1a <- grep("period_mn1a", names(co))
  i_1b <- grep("period_mn1b", names(co))
  if (length(i_1a) != 1 || length(i_1b) != 1) {
    cat(sprintf("  %s: cannot locate 1a/1b coefficients; skipped.\n", label))
    return(NULL)
  }
  b_diff  <- co[i_1b] - co[i_1a]
  v_diff  <- V[i_1a, i_1a] + V[i_1b, i_1b] - 2 * V[i_1a, i_1b]
  z       <- b_diff / sqrt(v_diff)
  p       <- 2 * pnorm(-abs(z))
  cat(sprintf("  %s: b(1b)-b(1a) = %+.3f (SE=%.3f)  p=%.4f\n",
              label, b_diff, sqrt(v_diff), p))
}

cat("\n  Test 1a (Mare Nostrum) vs 1b (Triton/Sophia) slope equality:\n")
wald_1a_vs_1b(f1_mn_iom_nb,     "IOM NB (no ctl)")
wald_1a_vs_1b(f1_mn_iom_nb_ctl, "IOM NB + ctl")
wald_1a_vs_1b(f1_mn_utd_nb,     "UNITED NB (no ctl)")
wald_1a_vs_1b(f1_mn_utd_nb_ctl, "UNITED NB + ctl")

# ── 2d. Boundary-leakage robustness (drop +/-30d around cuts) ─
# Policy implementation is gradual; sharp 2017-02-02 and 2020-10-21 cuts
# may attribute leakage windows to the wrong regime. Dropping +/-30 days
# around each boundary tests whether the 3-period decomposition survives.
cat("\n--- 2d. Boundary-leakage: drop +/-30d around 2017-02-02 and 2020-10-21 ---\n")

drop_window <- 30L
d_bl <- d |> filter(
  abs(as.integer(date - MOU_DATE)) > drop_window,
  abs(as.integer(date - (D_23 + 1))) > drop_window
)
cat(sprintf("  Boundary-leakage sample: N = %d days (drops %d days around 2 cuts)\n",
            nrow(d_bl), nrow(d) - nrow(d_bl)))

f1_3p_iom_nb_bl <- fenegbin(n_dead_iom    ~ swh_prev5days:period | month_year_fac,
                            data = d_bl, vcov = NW(14), panel.id = ~unit + date)
f1_3p_utd_nb_bl <- fenegbin(n_dead_united ~ swh_prev5days:period | month_year_fac,
                            data = d_bl, vcov = NW(14), panel.id = ~unit + date)
f1_3p_iom_nb_bl_ctl <- fenegbin(
  n_dead_iom ~ swh_prev5days:period + log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d_bl, vcov = NW(14), panel.id = ~unit + date)
f1_3p_utd_nb_bl_ctl <- fenegbin(
  n_dead_united ~ swh_prev5days:period + log1p_lc_lag14 + log1p_libya + log1p_tunisia |
    month_year_fac, data = d_bl, vcov = NW(14), panel.id = ~unit + date)

cat("\n  Boundary-leakage 3-period, NB, NW(14):\n")
print(etable(f1_3p_iom_nb_bl, f1_3p_iom_nb_bl_ctl, f1_3p_utd_nb_bl, f1_3p_utd_nb_bl_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM (no ctl)", "IOM + ctl",
                         "UNITED (no ctl)", "UNITED + ctl")))

cat("\n  Wald test (3 period gradients equal) under boundary leakage:\n")
wald_tbl_bl <- bind_rows(
  wald_periods(f1_3p_iom_nb_bl,     "IOM NB -30d (no ctl)"),
  wald_periods(f1_3p_iom_nb_bl_ctl, "IOM NB -30d + ctl"),
  wald_periods(f1_3p_utd_nb_bl,     "UNITED NB -30d (no ctl)"),
  wald_periods(f1_3p_utd_nb_bl_ctl, "UNITED NB -30d + ctl")
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

# 3-period, LPM (clean predicted probabilities) + logit.
f3_3p_iom_lpm <- feols(dd_iom    ~ swh_prev5days:period | month_year_fac,
                       data = d, vcov = NW(14), panel.id = ~unit + date)
f3_3p_utd_lpm <- feols(dd_united ~ swh_prev5days:period | month_year_fac,
                       data = d, vcov = NW(14), panel.id = ~unit + date)
f3_3p_iom_lgt <- feglm(dd_iom    ~ swh_prev5days:period | month_year_fac,
                       data = d, family = binomial("logit"),
                       vcov = ~month_year_fac, panel.id = ~unit + date)
f3_3p_utd_lgt <- feglm(dd_united ~ swh_prev5days:period | month_year_fac,
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

cat("\n  Deadly-day, 3-period LPM, NW(14):\n")
print(etable(f3_3p_iom_lpm, f3_3p_utd_lpm, vcov = NW(14), se.below = TRUE,
             headers = c("IOM LPM", "UNITED LPM")))
cat("\n  Deadly-day, 3-period logit, cluster(month_year):\n")
print(etable(f3_3p_iom_lgt, f3_3p_utd_lgt, vcov = ~month_year_fac,
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
pred_iom <- pred_dd(f3_3p_iom_lpm, "dd_iom",    "IOM")
pred_utd <- pred_dd(f3_3p_utd_lpm, "dd_united", "UNITED")

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
  extract_period_betas(f1_3p_iom_nb, "IOM", "NegBin"),
  extract_period_betas(f1_3p_iom_po, "IOM", "Poisson"),
  extract_period_betas(f1_3p_utd_nb, "UNITED", "NegBin"),
  extract_period_betas(f1_3p_utd_po, "UNITED", "Poisson")
)

p_panel1 <- ggplot(p4_df, aes(period, beta, colour = source)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.6, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  facet_wrap(~ family) +
  labs(title = "F1: SWH-mortality gradient by 3-period SAR regime",
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
cat("2-period cut: post_mou from 2017-02-02. 3-period cuts: 2017-02-01\n")
cat("(MoU) and 2020-10-20 (Lamorgese). Sample clipped at 2023-01-01 so\n")
cat("the Piantedosi-era tail does not contaminate period 3. Frontex-\n")
cat("bounded sample is left-truncated at 2014. Count/probability regime\n")
cat("decomposition; offset rate model is in 20 (boat controls in 27).\n")
cat("No control group -> descriptive.\n\n")

cat("=== FAMILY 1: count, 2-period (swh:post_mou) ===\n")
print(etable(f1_2p_iom_nb, f1_2p_iom_po, f1_2p_utd_nb, f1_2p_utd_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM NB", "IOM Pois", "UNITED NB", "UNITED Pois")))
cat("\n=== FAMILY 1: count, 3-period (swh:period) ===\n")
print(etable(f1_3p_iom_nb, f1_3p_iom_po, f1_3p_utd_nb, f1_3p_utd_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM NB", "IOM Pois", "UNITED NB", "UNITED Pois")))
cat("\nWald test (3 period gradients equal):\n")
for (i in seq_len(nrow(wald_tbl))) {
  r <- wald_tbl[i, ]
  cat(sprintf("  %-12s chi2(%d) = %.3f, p = %.4f\n",
              r$spec, r$df, r$chi2, r$p))
}

cat("\n=== FAMILY 1 with full controls (lag-14 + ACLED Libya/Tunisia) ===\n")
print(etable(f1_3p_iom_nb, f1_3p_iom_nb_ctl, f1_3p_utd_nb, f1_3p_utd_nb_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM (no ctl)", "IOM + ctl",
                         "UNITED (no ctl)", "UNITED + ctl")))
cat("\nWald test with full controls:\n")
for (i in seq_len(nrow(wald_tbl_ctl))) {
  r <- wald_tbl_ctl[i, ]
  cat(sprintf("  %-20s chi2(%d) = %.3f, p = %.4f\n",
              r$spec, r$df, r$chi2, r$p))
}

cat("\n=== FAMILY 1 Mare Nostrum / Triton-Sophia split inside P1 ===\n")
print(etable(f1_mn_iom_nb, f1_mn_iom_nb_ctl, f1_mn_utd_nb, f1_mn_utd_nb_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM (no ctl)", "IOM + ctl",
                         "UNITED (no ctl)", "UNITED + ctl")))
cat("\nTest 1a (Mare Nostrum) vs 1b (Triton/Sophia) slope equality:\n")
wald_1a_vs_1b(f1_mn_iom_nb,     "IOM NB (no ctl)")
wald_1a_vs_1b(f1_mn_iom_nb_ctl, "IOM NB + ctl")
wald_1a_vs_1b(f1_mn_utd_nb,     "UNITED NB (no ctl)")
wald_1a_vs_1b(f1_mn_utd_nb_ctl, "UNITED NB + ctl")

cat("\n=== FAMILY 1 boundary leakage (drop +/-30d around cuts) ===\n")
cat(sprintf("Sample: N=%d days (drops %d around 2017-02-02 and 2020-10-21).\n\n",
            nrow(d_bl), nrow(d) - nrow(d_bl)))
print(etable(f1_3p_iom_nb_bl, f1_3p_iom_nb_bl_ctl, f1_3p_utd_nb_bl, f1_3p_utd_nb_bl_ctl,
             vcov = NW(14), se.below = TRUE,
             headers = c("IOM (no ctl)", "IOM + ctl",
                         "UNITED (no ctl)", "UNITED + ctl")))
cat("\nWald test (3 period gradients equal) under boundary leakage:\n")
for (i in seq_len(nrow(wald_tbl_bl))) {
  r <- wald_tbl_bl[i, ]
  cat(sprintf("  %-25s chi2(%d) = %.3f, p = %.4f\n",
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
cat("\n--- 3-period LPM, NW(14) ---\n")
print(etable(f3_3p_iom_lpm, f3_3p_utd_lpm, vcov = NW(14), se.below = TRUE,
             headers = c("IOM LPM", "UNITED LPM")))
cat("\n--- 3-period logit, cluster(month_year) ---\n")
print(etable(f3_3p_iom_lgt, f3_3p_utd_lgt, vcov = ~month_year_fac,
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

# ── 7. LaTeX table (\input'd by paper/thesis.qmd) ─────────────
cat("\n--- 7. Writing LaTeX table ---\n")

sig_stars <- function(p) {
  ifelse(p < 0.001, "^{***}",
  ifelse(p < 0.01,  "^{**}",
  ifelse(p < 0.05,  "^{*}", "")))
}
fcoef <- function(b, p) sprintf("$%+.3f%s$", b, sig_stars(p))
fse   <- function(se)   sprintf("(%.3f)", se)
fint  <- function(x)    formatC(round(x), format = "d", big.mark = ",")
fp    <- function(p)    ifelse(p < 0.001, "$< 0.001$", sprintf("$%.3f$", p))

# Extract period-specific coefficients from each fit.
get_period <- function(m, p_label) {
  ct <- coeftable(m, vcov = NW(14))
  rn <- grep(paste0("swh_prev5days:period", p_label), rownames(ct), fixed = TRUE)
  list(b = ct[rn, 1], se = ct[rn, 2],
       p = 2 * pnorm(-abs(ct[rn, 1] / ct[rn, 2])))
}
periods_labels <- c("1. SAR + border control",
                    "2. MoU + NGO containment",
                    "3. Lamorgese rollback")
u <- lapply(periods_labels, function(p) get_period(f1_3p_utd_nb, p))
i <- lapply(periods_labels, function(p) get_period(f1_3p_iom_nb, p))

# Wald test (3 gradients equal) for each source from wald_tbl (NegBin rows).
w_u <- wald_tbl[wald_tbl$spec == "UNITED NB", ]
w_i <- wald_tbl[wald_tbl$spec == "IOM NB", ]

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[h!]")
add("\\centering")
add("\\small")
add("\\caption{Period-specific SWH--mortality gradient. NegBin with month--year FE and NW(14) SEs. UNITED is fit from 2013-01-01 to take advantage of its longer span; IOM is fit from 2014-01-01, when MMP coverage begins. Both samples end on the eve of the January 2023 Piantedosi decree.}")
add("\\label{tab:period}")
add("\\begin{tabular}{lcc}")
add("\\hline")
add(sprintf("SAR/policy regime                                & UNITED (%d--%d) & IOM (%d--%d) \\\\",
            year(min(d_utd_3p$date)), year(max(d_utd_3p$date)),
            year(min(d_iom_3p$date)), year(max(d_iom_3p$date))))
add("\\hline")
add("1. SAR + border control (through Feb 2017)       & ",
    fcoef(u[[1]]$b, u[[1]]$p), " & ", fcoef(i[[1]]$b, i[[1]]$p), " \\\\")
add("                                                 & ",
    fse(u[[1]]$se), " & ", fse(i[[1]]$se), " \\\\")
add("2. MoU + NGO containment (Feb 2017--Oct 2020)    & ",
    fcoef(u[[2]]$b, u[[2]]$p), " & ", fcoef(i[[2]]$b, i[[2]]$p), " \\\\")
add("                                                 & ",
    fse(u[[2]]$se), " & ", fse(i[[2]]$se), " \\\\")
add("3. Lamorgese rollback (Oct 2020--Jan 2023)       & ",
    fcoef(u[[3]]$b, u[[3]]$p), " & ", fcoef(i[[3]]$b, i[[3]]$p), " \\\\")
add("                                                 & ",
    fse(u[[3]]$se), " & ", fse(i[[3]]$se), " \\\\")
add("\\hline")
add("Wald $H_{0}$: all equal, $\\chi^{2}(2)$           & ",
    sprintf("$%.2f$", w_u$chi2), " & ", sprintf("$%.2f$", w_i$chi2), " \\\\")
add("$p$-value                                        & ",
    fp(w_u$p), " & ", fp(w_i$p), " \\\\")
add("Observations (days)                              & ",
    fint(nrow(d_utd_3p)), " & ", fint(nrow(d_iom_3p)), " \\\\")
add("\\hline")
add("\\multicolumn{3}{l}{\\footnotesize Stars denote two-sided $p$-values: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\multicolumn{3}{l}{\\footnotesize Newey--West standard errors in parentheses.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
out_period <- tbl_path("05_analysis", "04_period_sar_gradient.tex")
writeLines(L, out_period)
cat(sprintf("  Saved: %s\n", out_period))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
