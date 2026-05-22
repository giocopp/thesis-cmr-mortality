# 26_push_factor_decomposition.R
# ==============================
# Channel decomposition for the SWH x post_MoU shift on CMR deaths,
# using ACLED Libya conflict (composite: battles + explosive violence +
# violence against civilians) as a political-risk moderator.
#
# Logic:
#   The count-model SWH slope nets the DETERRENCE channel (SWH suppresses
#   departures) with the MORTALITY channel (SWH raises per-crossing risk).
#   ACLED conflict in Libya measures broader political instability. The
#   triple interaction SWH x post_MoU x ACLED tests whether the post-MoU
#   SWH-deaths shift varies with conflict intensity:
#     b_triple = 0   -> no detected heterogeneity by ACLED intensity
#     b_triple < 0   -> shift smaller on high-conflict days
#     b_triple > 0   -> shift bigger on high-conflict days
#
# Three layers:
#   L1: deaths ~ SWH + SWH:post_mou + ACLED_z + ACLED_z:post_mou
#                + SWH:ACLED_z + SWH:post_mou:ACLED_z | month_year_fac
#       UNITED primary + IOM comparison; NB + Poisson.
#   L2: crossing_attempts ~ SWH + SWH:post_mou + ACLED_z + ACLED_z:post_mou
#                            | month_year_fac
#       Deterrence model on the common denominator used by the rate
#       model in 20: crossing_attempts = frx_persons + lcg_tcg_pushbacks
#       + n_dead_missing. UNITED deaths are NOT in the denominator, so
#       the L1-L2 decomposition with UNITED on the death side is clean
#       of source-specific circularity.
#   L3: count spec stratified by ACLED_lag1w median (low vs high).
#       UNITED + IOM, NB.
#
# Plus:
#   L1b: L1 spec + symmetric SAR-share triple, on the SAR-restricted
#        sample. Tests whether the ACLED triple is just proxy-ing SAR
#        capacity heterogeneity.
#
# Note: libya_conflict (composite in the panel) == libya_battles +
# libya_expvio + libya_violciv (verified by row-wise sum; correlation =
# 1.000). Protests and riots are NOT in the composite. ACLED is broadcast
# weekly in acled_daily.RDS; ~4 distinct values per month, so within
# month_year FE identifying variation is week-to-week.
#
# In:  analysis/data/daily_panel_complete.RDS (has libya_conflict)
#      data/processed/{iom_mmp_incidents,united_incidents,core_corridor}.RDS
# Out: output/tables/26_push_factor_decomposition.txt
#      output/figures/26_push_factor_decomposition.png

library(tidyverse)
library(lubridate)
library(fixest)
library(sf)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-02-02")

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("26  PUSH-FACTOR DECOMPOSITION (libya_conflict x SWH x post_MoU)\n")
cat("============================================================\n\n")

# -- 1. Load data + build lagged ACLED moderator ------------
cat("--- 1. Loading data ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

iom_daily    <- build_iom_daily()
united_daily <- build_united_daily()

panel <- panel |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  left_join(united_daily,                                       by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) |>
  arrange(date) |>
  mutate(
    libya_conflict_lag1w = dplyr::lag(libya_conflict, 7),

    # SAR-capacity moderator for Layer 1b: weekly-lagged SAR share
    sar_events_pw = dplyr::lag(zoo::rollsumr(frx_n_sar,     k = 7, fill = NA), 1),
    incidents_pw  = dplyr::lag(zoo::rollsumr(frx_incidents, k = 7, fill = NA), 1),
    sar_share_pw  = ifelse(incidents_pw > 0,
                            sar_events_pw / incidents_pw, NA_real_),

    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    log1p_lc_lag14 = log1p(lc_lag14),

    unit           = 1L,
    year_fac       = factor(year),
    month_year_fac = factor(month_year)
  )

# Sample: same as 20 plus require ACLED lag non-NA.
d <- panel |>
  filter(!is.na(lc_lag14), !is.na(swh_prev5days),
         !is.na(libya_conflict_lag1w))

d <- d |>
  mutate(libya_conflict_z = (libya_conflict_lag1w - mean(libya_conflict_lag1w)) /
                              sd(libya_conflict_lag1w))

d_sar <- d |>
  filter(!is.na(sar_share_pw)) |>
  mutate(sar_share_pw_z = (sar_share_pw - mean(sar_share_pw)) /
                           sd(sar_share_pw))

cat(sprintf("  Main sample (d):        N = %d days (%s to %s)\n",
            nrow(d), min(d$date), max(d$date)))
cat(sprintf("  SAR-restricted (d_sar): N = %d (drop %d for NA sar_share)\n",
            nrow(d_sar), nrow(d) - nrow(d_sar)))
cat(sprintf("  UNITED deaths in d: %.0f over %d death-days\n",
            sum(d$n_dead_united), sum(d$n_dead_united > 0)))
cat(sprintf("  IOM    deaths in d: %.0f over %d death-days\n",
            sum(d$n_dead_iom),    sum(d$n_dead_iom    > 0)))

cat(sprintf("\n  libya_conflict_lag1w: mean=%.2f sd=%.2f median=%.0f range=[%.0f, %.0f]\n",
            mean(d$libya_conflict_lag1w), sd(d$libya_conflict_lag1w),
            median(d$libya_conflict_lag1w),
            min(d$libya_conflict_lag1w),  max(d$libya_conflict_lag1w)))
cat("\n  Key correlations:\n")
cat(sprintf("    cor(ACLED_z, swh_prev5days)     = %+.3f\n",
            cor(d$libya_conflict_z, d$swh_prev5days)))
cat(sprintf("    cor(ACLED_z, crossing_attempts) = %+.3f\n",
            cor(d$libya_conflict_z, d$crossing_attempts)))
cat(sprintf("    cor(ACLED_z, post_mou)          = %+.3f\n",
            cor(d$libya_conflict_z, d$post_mou)))

# -- 2. Helper functions ---------------------------------
fit_one <- function(dep, rhs, family, dat) {
  fn  <- if (family == "NegBin") fenegbin else fepois
  fml <- as.formula(sprintf("%s ~ %s", dep, rhs))
  tryCatch(
    fn(fml, data = dat, vcov = NW(14), panel.id = ~unit + date),
    error = function(e) {
      message(sprintf("[FIT FAILED] dep=%s family=%s: %s",
                      dep, family, conditionMessage(e)))
      NULL
    }
  )
}

extract_coef <- function(m, term) {
  if (is.null(m)) return(c(coef = NA_real_, se = NA_real_, p = NA_real_))
  ct <- coeftable(m, vcov = NW(14))
  rn <- rownames(ct)
  r  <- which(rn == term)
  if (length(r) == 0) {
    parts <- strsplit(term, ":")[[1]]
    if (length(parts) == 2) {
      r <- which(rn == paste(parts[2], parts[1], sep = ":"))
    }
  }
  if (length(r) == 0) {
    return(c(coef = NA_real_, se = NA_real_, p = NA_real_))
  }
  c(coef = ct[r, 1], se = ct[r, 2],
    p    = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

# -- 3. LAYER 1: triple interaction on deaths ---------------
cat("\n--- 2. LAYER 1: deaths ~ SWH x post_MoU x libya_conflict_z ---\n")

mod_z <- "libya_conflict_z"
triple_term <- sprintf("swh_prev5days:post_mou:%s", mod_z)
het_term    <- sprintf("swh_prev5days:%s", mod_z)
acled_term  <- mod_z
acledpost   <- sprintf("%s:post_mou", mod_z)
triple_sar  <- "swh_prev5days:post_mou:sar_share_pw_z"

l1_terms <- c(
  "swh_prev5days",
  "swh_prev5days:post_mou",
  mod_z,
  sprintf("%s:post_mou",                mod_z),
  sprintf("swh_prev5days:%s",           mod_z),
  sprintf("swh_prev5days:post_mou:%s",  mod_z)
)
rhs_l1 <- paste(paste(l1_terms, collapse = " + "), "| month_year_fac")

l1 <- list(
  united_nb = fit_one("n_dead_united", rhs_l1, "NegBin",  d),
  united_po = fit_one("n_dead_united", rhs_l1, "Poisson", d),
  iom_nb    = fit_one("n_dead_iom",    rhs_l1, "NegBin",  d),
  iom_po    = fit_one("n_dead_iom",    rhs_l1, "Poisson", d)
)

cat("\n  Layer 1 full coefficients, NW(14):\n")
print(etable(l1$united_nb, l1$united_po, l1$iom_nb, l1$iom_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED NB", "UNITED Pois",
                         "IOM NB",    "IOM Pois")))

cat("\n  L1 key coefficients:\n")
cat(sprintf("    %-12s %-8s  %+10s  %10s  %10s  %s\n",
            "source", "family", "coef", "SE", "p", "term"))
for (src in c("UNITED", "IOM")) {
  for (fam in c("NegBin", "Poisson")) {
    key  <- paste0(tolower(src), "_", if (fam == "NegBin") "nb" else "po")
    mfit <- l1[[key]]
    for (tm in c("swh_prev5days", "swh_prev5days:post_mou",
                  acled_term, acledpost, het_term, triple_term)) {
      r <- extract_coef(mfit, tm)
      star <- if (!is.na(r["p"]) && r["p"] < 0.05) " *" else ""
      cat(sprintf("    %-12s %-8s  %+10.3f  %10.3f  %10.4f  %s%s\n",
                  src, fam, r["coef"], r["se"], r["p"], tm, star))
    }
  }
}

# -- 4. LAYER 1b: + SAR-share triple (channel-confounding check) ----
cat("\n--- 3. LAYER 1b: + symmetric SAR-share triple (SAR-restricted) ---\n")

l1_sarsmp <- list(
  united_nb = fit_one("n_dead_united", rhs_l1, "NegBin",  d_sar),
  united_po = fit_one("n_dead_united", rhs_l1, "Poisson", d_sar),
  iom_nb    = fit_one("n_dead_iom",    rhs_l1, "NegBin",  d_sar),
  iom_po    = fit_one("n_dead_iom",    rhs_l1, "Poisson", d_sar)
)

l1b_terms <- c(
  l1_terms,
  "sar_share_pw_z",
  "sar_share_pw_z:post_mou",
  "swh_prev5days:sar_share_pw_z",
  "swh_prev5days:post_mou:sar_share_pw_z"
)
rhs_l1b <- paste(paste(l1b_terms, collapse = " + "), "| month_year_fac")

l1b <- list(
  united_nb = fit_one("n_dead_united", rhs_l1b, "NegBin",  d_sar),
  united_po = fit_one("n_dead_united", rhs_l1b, "Poisson", d_sar),
  iom_nb    = fit_one("n_dead_iom",    rhs_l1b, "NegBin",  d_sar),
  iom_po    = fit_one("n_dead_iom",    rhs_l1b, "Poisson", d_sar)
)

cat("\n  Layer 1b coefficients, NW(14):\n")
print(etable(l1b$united_nb, l1b$united_po, l1b$iom_nb, l1b$iom_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED NB", "UNITED Pois", "IOM NB", "IOM Pois")))

cat("\n  L1 -> L1b: ACLED triple stability check\n")
cat(sprintf("  %-12s %-8s %-14s  %+10s  %10s  %10s  %s\n",
            "source", "family", "spec", "ACLED triple", "SE", "p", "N"))
for (src in c("UNITED", "IOM")) {
  for (fam in c("NegBin", "Poisson")) {
    key  <- paste0(tolower(src), "_", if (fam == "NegBin") "nb" else "po")
    r_b   <- extract_coef(l1[[key]],         triple_term)
    r_s   <- extract_coef(l1_sarsmp[[key]],  triple_term)
    r_c   <- extract_coef(l1b[[key]],        triple_term)
    n_b <- if (is.null(l1[[key]]))        NA_integer_ else nobs(l1[[key]])
    n_s <- if (is.null(l1_sarsmp[[key]])) NA_integer_ else nobs(l1_sarsmp[[key]])
    n_c <- if (is.null(l1b[[key]]))       NA_integer_ else nobs(l1b[[key]])
    st_b <- if (!is.na(r_b["p"]) && r_b["p"] < 0.05) " *" else ""
    st_s <- if (!is.na(r_s["p"]) && r_s["p"] < 0.05) " *" else ""
    st_c <- if (!is.na(r_c["p"]) && r_c["p"] < 0.05) " *" else ""
    cat(sprintf("  %-12s %-8s %-14s  %+10.3f  %10.3f  %10.4f%s  %d\n",
                src, fam, "base (full)",    r_b["coef"], r_b["se"], r_b["p"], st_b, n_b))
    cat(sprintf("  %-12s %-8s %-14s  %+10.3f  %10.3f  %10.4f%s  %d\n",
                src, fam, "base (SAR smp)", r_s["coef"], r_s["se"], r_s["p"], st_s, n_s))
    cat(sprintf("  %-12s %-8s %-14s  %+10.3f  %10.3f  %10.4f%s  %d\n",
                src, fam, "+SAR triple",    r_c["coef"], r_c["se"], r_c["p"], st_c, n_c))
  }
}

cat("\n  Layer 1b: SAR-share triple (own channel diagnostic):\n")
for (src in c("UNITED", "IOM")) {
  for (fam in c("NegBin", "Poisson")) {
    key <- paste0(tolower(src), "_", if (fam == "NegBin") "nb" else "po")
    r   <- extract_coef(l1b[[key]], triple_sar)
    star <- if (!is.na(r["p"]) && r["p"] < 0.05) " *" else ""
    cat(sprintf("    %-12s %-8s  %+10.3f  %10.3f  %10.4f%s\n",
                src, fam, r["coef"], r["se"], r["p"], star))
  }
}

# -- 5. LAYER 2: deterrence model on crossing_attempts ----
cat("\n--- 4. LAYER 2: crossing_attempts ~ SWH x post_mou + ACLED + ACLED:post_mou ---\n")

l2_terms <- c(
  "swh_prev5days",
  "swh_prev5days:post_mou",
  mod_z,
  sprintf("%s:post_mou", mod_z)
)
rhs_l2 <- paste(paste(l2_terms, collapse = " + "), "| month_year_fac")

l2 <- list(
  ca_nb = fit_one("crossing_attempts", rhs_l2, "NegBin",  d),
  ca_po = fit_one("crossing_attempts", rhs_l2, "Poisson", d)
)

cat("\n  Layer 2 coefficients, NW(14):\n")
print(etable(l2$ca_nb, l2$ca_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("CrossAttempts NB", "CrossAttempts Pois")))

cat("\n  L2 key coefficients:\n")
cat(sprintf("    %-20s %-8s  %+10s  %10s  %10s  %s\n",
            "outcome", "family", "coef", "SE", "p", "term"))
for (fam in c("NegBin", "Poisson")) {
  key  <- paste0("ca_", if (fam == "NegBin") "nb" else "po")
  mfit <- l2[[key]]
  for (tm in c("swh_prev5days", "swh_prev5days:post_mou",
                acled_term, acledpost)) {
    r <- extract_coef(mfit, tm)
    star <- if (!is.na(r["p"]) && r["p"] < 0.05) " *" else ""
    cat(sprintf("    %-20s %-8s  %+10.3f  %10.3f  %10.4f  %s%s\n",
                "crossing_attempts", fam, r["coef"], r["se"], r["p"], tm, star))
  }
}

# -- 6. LAYER 3: stratified pre/post by ACLED median ------
cat("\n--- 5. LAYER 3: median-split stratification by ACLED ---\n")

raw_var <- "libya_conflict_lag1w"
med <- median(d[[raw_var]])
d_low  <- d |> filter(.data[[raw_var]] <= med)
d_high <- d |> filter(.data[[raw_var]] >  med)

cat(sprintf("  Median ACLED (lag1w): %.0f. Low N=%d, High N=%d.\n",
            med, nrow(d_low), nrow(d_high)))
cat(sprintf("  Low:  UNITED=%.0f, IOM=%.0f\n",
            sum(d_low$n_dead_united),  sum(d_low$n_dead_iom)))
cat(sprintf("  High: UNITED=%.0f, IOM=%.0f\n",
            sum(d_high$n_dead_united), sum(d_high$n_dead_iom)))

rhs_l3 <- "swh_prev5days + swh_prev5days:post_mou | month_year_fac"
l3 <- list(
  utd_low  = fit_one("n_dead_united", rhs_l3, "NegBin", d_low),
  utd_high = fit_one("n_dead_united", rhs_l3, "NegBin", d_high),
  iom_low  = fit_one("n_dead_iom",    rhs_l3, "NegBin", d_low),
  iom_high = fit_one("n_dead_iom",    rhs_l3, "NegBin", d_high)
)

slope_decomp <- function(mfit, src, stratum) {
  if (is.null(mfit)) {
    return(tibble(source = src, stratum = stratum, n_obs = NA_integer_,
                  b_pre = NA_real_,   se_pre   = NA_real_, p_pre   = NA_real_,
                  b_shift = NA_real_, se_shift = NA_real_, p_shift = NA_real_,
                  b_post = NA_real_,  se_post  = NA_real_, p_post  = NA_real_))
  }
  co <- coef(mfit); V <- vcov(mfit, vcov = NW(14))
  bp <- unname(co["swh_prev5days"])
  bs <- unname(co["swh_prev5days:post_mou"])
  sp <- sqrt(V["swh_prev5days", "swh_prev5days"])
  ss <- sqrt(V["swh_prev5days:post_mou", "swh_prev5days:post_mou"])
  bpost <- bp + bs
  vpost <- V["swh_prev5days", "swh_prev5days"] +
    V["swh_prev5days:post_mou", "swh_prev5days:post_mou"] +
    2 * V["swh_prev5days", "swh_prev5days:post_mou"]
  spost <- sqrt(vpost)
  tibble(source = src, stratum = stratum, n_obs = nobs(mfit),
         b_pre   = bp,    se_pre   = sp,    p_pre   = 2 * pnorm(-abs(bp / sp)),
         b_shift = bs,    se_shift = ss,    p_shift = 2 * pnorm(-abs(bs / ss)),
         b_post  = bpost, se_post  = spost, p_post  = 2 * pnorm(-abs(bpost / spost)))
}

l3_tbl <- bind_rows(
  slope_decomp(l3$utd_low,  "UNITED", "low push"),
  slope_decomp(l3$utd_high, "UNITED", "high push"),
  slope_decomp(l3$iom_low,  "IOM",    "low push"),
  slope_decomp(l3$iom_high, "IOM",    "high push")
)

cat("\n  L3 slopes (NegBin, NW(14), month_year FE):\n")
cat(sprintf("  %-7s %-9s %6s  %+10s %10s %10s  %+10s %10s %10s  %+10s %10s %10s\n",
            "source", "stratum", "N",
            "b_pre", "SE", "p", "b_shift", "SE", "p", "b_post", "SE", "p"))
for (i in seq_len(nrow(l3_tbl))) {
  r <- l3_tbl[i, ]
  cat(sprintf("  %-7s %-9s %6d  %+10.3f %10.3f %10.4f  %+10.3f %10.3f %10.4f  %+10.3f %10.3f %10.4f\n",
              r$source, r$stratum, r$n_obs,
              r$b_pre,   r$se_pre,   r$p_pre,
              r$b_shift, r$se_shift, r$p_shift,
              r$b_post,  r$se_post,  r$p_post))
}

# -- 7. Coefficient plot ----------------------------------
cat("\n--- 6. Plot ---\n")

build_triple_row <- function(mfit, src, fam, spec) {
  r <- extract_coef(mfit, triple_term)
  tibble(source = src, family = fam, spec = spec,
         coef = r["coef"], se = r["se"])
}

plot_df <- bind_rows(
  build_triple_row(l1$united_nb,  "UNITED", "NegBin",  "L1 base"),
  build_triple_row(l1$united_po,  "UNITED", "Poisson", "L1 base"),
  build_triple_row(l1$iom_nb,     "IOM",    "NegBin",  "L1 base"),
  build_triple_row(l1$iom_po,     "IOM",    "Poisson", "L1 base"),
  build_triple_row(l1b$united_nb, "UNITED", "NegBin",  "L1b +SAR ctrl"),
  build_triple_row(l1b$united_po, "UNITED", "Poisson", "L1b +SAR ctrl"),
  build_triple_row(l1b$iom_nb,    "IOM",    "NegBin",  "L1b +SAR ctrl"),
  build_triple_row(l1b$iom_po,    "IOM",    "Poisson", "L1b +SAR ctrl")
) |>
  mutate(ci_lo  = coef - 1.96 * se,
         ci_hi  = coef + 1.96 * se,
         source = factor(source, levels = c("UNITED", "IOM")),
         spec   = factor(spec,   levels = c("L1 base", "L1b +SAR ctrl")))

p <- ggplot(plot_df,
            aes(coef, family, colour = source, shape = spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.55)) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                width = 0.2, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c("UNITED" = "#B2182B", "IOM" = "#2166AC")) +
  scale_shape_manual(values  = c("L1 base" = 16, "L1b +SAR ctrl" = 17)) +
  facet_wrap(~ source, ncol = 1) +
  labs(
    title    = "ACLED Libya x SWH x post_MoU triple interaction on deaths",
    subtitle = "L1 base vs L1b (+ SAR-share triple). Stable across L1->L1b => not just SAR-capacity proxy.",
    x        = "Coefficient on SWH x post_MoU x ACLED_z (per 1m SWH x 1-SD ACLED)",
    y        = NULL, colour = "Source", shape = "Spec"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

push_fig <- fig_path("06_robustness", "04_push_factor_decomposition.png")
ggsave(push_fig, p, width = 11, height = 6, dpi = 200)
cat(sprintf("  Saved: %s\n", push_fig))

# -- 8. Save text output -----------------------------------
cat("\n--- 7. Saving results ---\n")

sink_file <- tbl_path("06_robustness", "04_push_factor_decomposition.txt")
sink(sink_file)

cat("26  PUSH-FACTOR DECOMPOSITION (libya_conflict x SWH x post_MoU)\n")
cat("============================================================\n\n")
cat(sprintf("Main sample (d):        N = %d days (%s to %s)\n",
            nrow(d), min(d$date), max(d$date)))
cat(sprintf("SAR-restricted (d_sar): N = %d\n", nrow(d_sar)))
cat(sprintf("UNITED deaths: %.0f over %d death-days\n",
            sum(d$n_dead_united), sum(d$n_dead_united > 0)))
cat(sprintf("IOM    deaths: %.0f over %d death-days\n",
            sum(d$n_dead_iom),    sum(d$n_dead_iom    > 0)))
cat(sprintf("libya_conflict_lag1w: mean=%.2f sd=%.2f median=%.0f range=[%.0f, %.0f]\n",
            mean(d$libya_conflict_lag1w), sd(d$libya_conflict_lag1w),
            median(d$libya_conflict_lag1w),
            min(d$libya_conflict_lag1w),  max(d$libya_conflict_lag1w)))
cat(sprintf("cor(ACLED_z, swh_prev5days)     = %+.3f\n",
            cor(d$libya_conflict_z, d$swh_prev5days)))
cat(sprintf("cor(ACLED_z, crossing_attempts) = %+.3f\n",
            cor(d$libya_conflict_z, d$crossing_attempts)))
cat(sprintf("cor(ACLED_z, post_mou)          = %+.3f\n\n",
            cor(d$libya_conflict_z, d$post_mou)))

cat("Identification logic:\n")
cat("  L2 outcome is crossing_attempts (= frx_persons + lcg_tcg_pushbacks +\n")
cat("  n_dead_missing). UNITED deaths are NOT in the denominator, so the\n")
cat("  L1-L2 decomposition with UNITED on the death side is clean of\n")
cat("  source-specific circularity.\n\n")

cat("=== LAYER 1: deaths ~ SWH + SWH:post_mou + ACLED_z\n")
cat("            + ACLED_z:post_mou + SWH:ACLED_z\n")
cat("            + SWH:post_mou:ACLED_z | month_year_fac ===\n\n")
print(etable(l1$united_nb, l1$united_po, l1$iom_nb, l1$iom_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED NB", "UNITED Pois", "IOM NB", "IOM Pois")))

cat("\n--- L1 key coefficients ---\n")
cat(sprintf("  %-12s %-8s  %+10s  %10s  %10s  %s\n",
            "source", "family", "coef", "SE", "p", "term"))
for (src in c("UNITED", "IOM")) {
  for (fam in c("NegBin", "Poisson")) {
    key  <- paste0(tolower(src), "_", if (fam == "NegBin") "nb" else "po")
    mfit <- l1[[key]]
    for (tm in c("swh_prev5days", "swh_prev5days:post_mou",
                  acled_term, acledpost, het_term, triple_term)) {
      r <- extract_coef(mfit, tm)
      star <- if (!is.na(r["p"]) && r["p"] < 0.05) " *" else ""
      cat(sprintf("  %-12s %-8s  %+10.3f  %10.3f  %10.4f  %s%s\n",
                  src, fam, r["coef"], r["se"], r["p"], tm, star))
    }
  }
}

cat(sprintf("\n\n=== LAYER 1b: ACLED triple + SAR-share triple (SAR-restricted N=%d) ===\n\n",
            nrow(d_sar)))
print(etable(l1b$united_nb, l1b$united_po, l1b$iom_nb, l1b$iom_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("UNITED NB", "UNITED Pois", "IOM NB", "IOM Pois")))

cat("\n--- L1 -> L1b: ACLED triple coefficient across specs ---\n")
cat(sprintf("  %-12s %-8s %-14s  %+10s  %10s  %10s  %s\n",
            "source", "family", "spec", "ACLED triple", "SE", "p", "N"))
for (src in c("UNITED", "IOM")) {
  for (fam in c("NegBin", "Poisson")) {
    key  <- paste0(tolower(src), "_", if (fam == "NegBin") "nb" else "po")
    r_b   <- extract_coef(l1[[key]],         triple_term)
    r_s   <- extract_coef(l1_sarsmp[[key]],  triple_term)
    r_c   <- extract_coef(l1b[[key]],        triple_term)
    n_b <- if (is.null(l1[[key]]))        NA_integer_ else nobs(l1[[key]])
    n_s <- if (is.null(l1_sarsmp[[key]])) NA_integer_ else nobs(l1_sarsmp[[key]])
    n_c <- if (is.null(l1b[[key]]))       NA_integer_ else nobs(l1b[[key]])
    st_b <- if (!is.na(r_b["p"]) && r_b["p"] < 0.05) " *" else ""
    st_s <- if (!is.na(r_s["p"]) && r_s["p"] < 0.05) " *" else ""
    st_c <- if (!is.na(r_c["p"]) && r_c["p"] < 0.05) " *" else ""
    cat(sprintf("  %-12s %-8s %-14s  %+10.3f  %10.3f  %10.4f%s  %d\n",
                src, fam, "base (full)",    r_b["coef"], r_b["se"], r_b["p"], st_b, n_b))
    cat(sprintf("  %-12s %-8s %-14s  %+10.3f  %10.3f  %10.4f%s  %d\n",
                src, fam, "base (SAR smp)", r_s["coef"], r_s["se"], r_s["p"], st_s, n_s))
    cat(sprintf("  %-12s %-8s %-14s  %+10.3f  %10.3f  %10.4f%s  %d\n",
                src, fam, "+SAR triple",    r_c["coef"], r_c["se"], r_c["p"], st_c, n_c))
  }
}

cat("\n--- Layer 1b: SAR-share triple (own channel diagnostic) ---\n")
cat(sprintf("  %-12s %-8s  %+10s  %10s  %10s\n",
            "source", "family", "coef", "SE", "p"))
for (src in c("UNITED", "IOM")) {
  for (fam in c("NegBin", "Poisson")) {
    key <- paste0(tolower(src), "_", if (fam == "NegBin") "nb" else "po")
    r   <- extract_coef(l1b[[key]], triple_sar)
    star <- if (!is.na(r["p"]) && r["p"] < 0.05) " *" else ""
    cat(sprintf("  %-12s %-8s  %+10.3f  %10.3f  %10.4f%s\n",
                src, fam, r["coef"], r["se"], r["p"], star))
  }
}

cat("\n\n=== LAYER 2: crossing_attempts ~ SWH + SWH:post_mou + ACLED_z\n")
cat("            + ACLED_z:post_mou | month_year_fac ===\n\n")
print(etable(l2$ca_nb, l2$ca_po,
             vcov = NW(14), se.below = TRUE,
             headers = c("CrossAttempts NB", "CrossAttempts Pois")))

cat("\n--- L2 key coefficients ---\n")
cat(sprintf("  %-20s %-8s  %+10s  %10s  %10s  %s\n",
            "outcome", "family", "coef", "SE", "p", "term"))
for (fam in c("NegBin", "Poisson")) {
  key  <- paste0("ca_", if (fam == "NegBin") "nb" else "po")
  mfit <- l2[[key]]
  for (tm in c("swh_prev5days", "swh_prev5days:post_mou",
                acled_term, acledpost)) {
    r <- extract_coef(mfit, tm)
    star <- if (!is.na(r["p"]) && r["p"] < 0.05) " *" else ""
    cat(sprintf("  %-20s %-8s  %+10.3f  %10.3f  %10.4f  %s%s\n",
                "crossing_attempts", fam, r["coef"], r["se"], r["p"], tm, star))
  }
}

cat(sprintf("\n\n=== LAYER 3: stratified by libya_conflict_lag1w median = %.0f ===\n",
            med))
cat(sprintf("  Low:  N=%d  UNITED=%.0f  IOM=%.0f\n",
            nrow(d_low),  sum(d_low$n_dead_united),  sum(d_low$n_dead_iom)))
cat(sprintf("  High: N=%d  UNITED=%.0f  IOM=%.0f\n\n",
            nrow(d_high), sum(d_high$n_dead_united), sum(d_high$n_dead_iom)))

cat(sprintf("  %-7s %-9s %6s  %+10s %10s %10s  %+10s %10s %10s  %+10s %10s %10s\n",
            "source", "stratum", "N",
            "b_pre", "SE", "p", "b_shift", "SE", "p", "b_post", "SE", "p"))
for (i in seq_len(nrow(l3_tbl))) {
  r <- l3_tbl[i, ]
  cat(sprintf("  %-7s %-9s %6d  %+10.3f %10.3f %10.4f  %+10.3f %10.3f %10.4f  %+10.3f %10.3f %10.4f\n",
              r$source, r$stratum, r$n_obs,
              r$b_pre,   r$se_pre,   r$p_pre,
              r$b_shift, r$se_shift, r$p_shift,
              r$b_post,  r$se_post,  r$p_post))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
