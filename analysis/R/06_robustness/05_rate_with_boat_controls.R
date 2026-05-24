# ── Rate model with boat-composition controls ──────────────────────────────
# Tests whether the SWH:post_mou shift survives inflatable/wooden share
# controls (Deiana-style composition probe). V1: rate alone; V2: + shares;
# V3: + SWH × inflatable. Poisson, NW(14), IOM + UNITED.

library(tidyverse)
library(lubridate)
library(fixest)
library(sf)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("27  RATE MODEL WITH BOAT COMPOSITION CONTROLS\n")
cat("    Extension of m_rate / m_rate_u in 20_primary_model.R\n")
cat("============================================================\n\n")

# ── 1. Load data (mirror of 20's data prep) ──────────────────
cat("--- 1. Loading data + 20-style prep ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

# IOM daily (same helper / default filter as 20)
iom_daily <- build_iom_daily()

# UNITED daily via the shared builder. Defaults (corridor spatial join;
# country in CMR+Med; manner drowned/other_unknown) replicate the previous
# inline filter exactly — single source of truth, see _helpers.R.
united_daily <- build_united_daily()

panel <- panel |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  left_join(united_daily, by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) |>
  add_crossing_exposure() |>   # living_crossings used for lag controls
  mutate(
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    log1p_lc_lag14 = log1p(lc_lag14),
    unit           = 1L,
    year_fac       = factor(year),
    month_year_fac = factor(month_year)
  )

# 20's rate-model sample: crossing_attempts > 0 + lc_lag14/swh_prev5days
# non-NA. Common denominator -> one sample for both sources.
d_rate_full <- panel |>
  filter(!is.na(lc_lag14), !is.na(swh_prev5days),
         crossing_attempts > 0) |>
  mutate(log_crossing_attempts = log(crossing_attempts))

# 27 boat-observable sub-sample: additionally require frx_incidents > 0
# so that boat-composition shares are defined.
d_rate_boat <- d_rate_full |>
  filter(frx_incidents > 0)

cat(sprintf("  Rate sample (20 baseline):     N = %d days\n",
            nrow(d_rate_full)))
cat(sprintf("  Boat-observable sub-sample:    N = %d days (%.1f%% of full)\n",
            nrow(d_rate_boat),
            100 * nrow(d_rate_boat) / nrow(d_rate_full)))
cat(sprintf("  Days dropped (no Frontex events): %d\n",
            nrow(d_rate_full) - nrow(d_rate_boat)))

cat("\n  Boat composition summary (IOM boat-observable sample):\n")
cat(sprintf("    Inflatable share: mean = %.3f  sd = %.3f  range = [%.2f, %.2f]\n",
            mean(d_rate_boat$frx_inflatable_share),
            sd(d_rate_boat$frx_inflatable_share),
            min(d_rate_boat$frx_inflatable_share),
            max(d_rate_boat$frx_inflatable_share)))
cat(sprintf("    Wooden share:     mean = %.3f  sd = %.3f  range = [%.2f, %.2f]\n",
            mean(d_rate_boat$frx_wooden_share),
            sd(d_rate_boat$frx_wooden_share),
            min(d_rate_boat$frx_wooden_share),
            max(d_rate_boat$frx_wooden_share)))

cat("\n  Pre/post-MoU mean inflatable share (descriptive, matches Deiana Fig 9):\n")
print(d_rate_boat |>
        mutate(period = ifelse(date >= MOU_DATE, "Post-MoU", "Pre-MoU")) |>
        group_by(period) |>
        summarise(n_days = n(),
                  mean_inflatable = mean(frx_inflatable_share),
                  mean_wooden     = mean(frx_wooden_share),
                  .groups = "drop"))

# ── 2. Three specifications -- IOM ───────────────────────────
cat("\n--- 2. IOM rate model: V1, V2, V3 ---\n")

# V1: baseline rate model (matches 20 m_rate spec)
v1_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

# V2: + boat composition (additive)
v2_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               frx_inflatable_share + frx_wooden_share +
               log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

# V3: + boat × SWH interaction (Deiana-style mediator probe)
v3_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               frx_inflatable_share + frx_wooden_share +
               swh_prev5days:frx_inflatable_share +
               log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

cat("\n  IOM, NW(14):\n")
print(etable(v1_iom, v2_iom, v3_iom,
             vcov = NW(14), se.below = TRUE,
             headers = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH")))

# ── 3. Three specifications -- UNITED ────────────────────────
cat("\n--- 3. UNITED rate model: V1, V2, V3 ---\n")

v1_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

v2_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  frx_inflatable_share + frx_wooden_share +
                  log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

v3_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  frx_inflatable_share + frx_wooden_share +
                  swh_prev5days:frx_inflatable_share +
                  log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

cat("\n  UNITED, NW(14):\n")
print(etable(v1_united, v2_united, v3_united,
             vcov = NW(14), se.below = TRUE,
             headers = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH")))

# ── 4. Comparison summary ────────────────────────────────────
cat("\n--- 4. Summary: SWH x post_MoU across V1, V2, V3 ---\n")

extract_b3 <- function(m, label, source) {
  ct <- coeftable(m, vcov = NW(14))
  r <- which(rownames(ct) == "swh_prev5days:post_mou")
  if (length(r) == 0) {
    return(tibble(spec = label, source = source,
                  coef = NA_real_, se = NA_real_, p = NA_real_))
  }
  tibble(spec = label, source = source,
         coef = ct[r, 1], se = ct[r, 2],
         p = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

extract_swh_x_inflatable <- function(m, label, source) {
  ct <- coeftable(m, vcov = NW(14))
  r <- which(rownames(ct) == "swh_prev5days:frx_inflatable_share")
  if (length(r) == 0) {
    return(tibble(spec = label, source = source,
                  coef = NA_real_, se = NA_real_, p = NA_real_))
  }
  tibble(spec = label, source = source,
         coef = ct[r, 1], se = ct[r, 2],
         p = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

summary_b3 <- bind_rows(
  extract_b3(v1_iom,    "V1: Baseline",     "IOM"),
  extract_b3(v2_iom,    "V2: +boat ctrl",   "IOM"),
  extract_b3(v3_iom,    "V3: +boat x SWH",  "IOM"),
  extract_b3(v1_united, "V1: Baseline",     "UNITED"),
  extract_b3(v2_united, "V2: +boat ctrl",   "UNITED"),
  extract_b3(v3_united, "V3: +boat x SWH",  "UNITED")
)

cat("\n  SWH x post_MoU (recorded-death rate, NW(14)):\n\n")
for (i in seq_len(nrow(summary_b3))) {
  r <- summary_b3[i, ]
  star <- if (!is.na(r$p) && r$p < 0.05) " *" else ""
  cat(sprintf("    %-8s %-18s  b = %+.3f (SE=%.3f)  p = %.4f%s\n",
              r$source, r$spec, r$coef, r$se, r$p, star))
}

v3_inflatable <- bind_rows(
  extract_swh_x_inflatable(v3_iom,    "V3", "IOM"),
  extract_swh_x_inflatable(v3_united, "V3", "UNITED")
)

cat("\n  SWH x inflatable_share (V3 only, Deiana-style mediator probe):\n\n")
for (i in seq_len(nrow(v3_inflatable))) {
  r <- v3_inflatable[i, ]
  star <- if (!is.na(r$p) && r$p < 0.05) " *" else ""
  cat(sprintf("    %-8s %-18s  b = %+.3f (SE=%.3f)  p = %.4f%s\n",
              r$source, r$spec, r$coef, r$se, r$p, star))
}

# ── 5. Coefficient plot ──────────────────────────────────────
cat("\n--- 5. Plot ---\n")

plot_df <- summary_b3 |>
  mutate(
    ci_lo = coef - 1.96 * se,
    ci_hi = coef + 1.96 * se,
    spec_f = factor(spec,
                    levels = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH"))
  )

p <- ggplot(plot_df, aes(coef, spec_f, colour = source)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                width = 0.2, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  labs(
    title    = "Recorded-death rate: SWH x post_MoU across boat-control specifications",
    subtitle = "Poisson with log(crossing_attempts) free covariate. Month-year FE. NW(14) SEs.",
    x        = "SWH x post_MoU coefficient (per 1m SWH, log rate)",
    y        = NULL,
    colour   = "Source"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

boat_fig <- fig_path("06_robustness", "05_rate_with_boat_controls_coefplot.png")
ggsave(boat_fig, p, width = 11, height = 5, dpi = 200)
cat(sprintf("  Saved: %s\n", boat_fig))

# ── 6. Save text output ──────────────────────────────────────
cat("\n--- 6. Saving results ---\n")

sink_file <- tbl_path("06_robustness", "05_rate_with_boat_controls.txt")
sink(sink_file)

cat("27  RATE MODEL WITH BOAT COMPOSITION CONTROLS\n")
cat("Extension of m_rate / m_rate_u in 20_primary_model.R\n")
cat("=====================================================\n\n")

cat(sprintf("Sample (boat-observable):  %s to %s\n",
            min(d_rate_boat$date), max(d_rate_boat$date)))
cat(sprintf("                           N = %d days (single sample, both sources)\n",
            nrow(d_rate_boat)))
cat(sprintf("                           (vs. 20 rate sample: N = %d days)\n",
            nrow(d_rate_full)))
cat("                           Sample is restricted to frx_incidents > 0 so boat shares are defined.\n")
cat("Outcome: recorded deaths (IOM or UNITED) on a common-denominator rate.\n")
cat("Standard errors: Newey-West (lag = 14)\n")
cat("Common denom: crossing_attempts = frx_persons + lcg_tcg_pushbacks + n_dead_missing\n")
cat("log(crossing_attempts) enters as a free covariate (not a forced offset).\n")
cat("N = estimation N after fixest drops all-zero FE cells.\n\n")

cat("Three specifications:\n")
cat("  V1: Baseline (matches 20 m_rate spec, on the boat-observable sample)\n")
cat("       deaths ~ swh + swh:post_mou + log(crossing_attempts) | month_year_fac\n")
cat("  V2: + boat composition (additive)\n")
cat("       deaths ~ swh + swh:post_mou + inflatable_share + wooden_share\n")
cat("              + log(crossing_attempts) | month_year_fac\n")
cat("  V3: + boat x SWH interaction (Deiana-style mediator probe)\n")
cat("       deaths ~ swh + swh:post_mou + swh:inflatable_share +\n")
cat("              + inflatable_share + wooden_share + log(crossing_attempts) | month_year_fac\n\n")

cat("=== Boat composition: pre/post-MoU descriptive ===\n")
print(d_rate_boat |>
        mutate(period = ifelse(date >= MOU_DATE, "Post-MoU", "Pre-MoU")) |>
        group_by(period) |>
        summarise(n_days = n(),
                  mean_inflatable = mean(frx_inflatable_share),
                  mean_wooden     = mean(frx_wooden_share),
                  .groups = "drop"))

cat("\n=== IOM rate model -- V1, V2, V3 (NW(14)) ===\n")
print(etable(v1_iom, v2_iom, v3_iom,
             vcov = NW(14), se.below = TRUE,
             headers = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH")))

cat("\n=== UNITED rate model -- V1, V2, V3 (NW(14)) ===\n")
print(etable(v1_united, v2_united, v3_united,
             vcov = NW(14), se.below = TRUE,
             headers = c("V1: Baseline", "V2: +boat ctrl", "V3: +boat x SWH")))

cat("\n=== SUMMARY: SWH x post_MoU coefficient across versions ===\n\n")
for (i in seq_len(nrow(summary_b3))) {
  r <- summary_b3[i, ]
  star <- if (!is.na(r$p) && r$p < 0.05) " *" else ""
  cat(sprintf("  %-8s %-18s  b = %+.3f (SE=%.3f)  p = %.4f%s\n",
              r$source, r$spec, r$coef, r$se, r$p, star))
}

cat("\n=== SWH x inflatable_share interaction (V3 only, Deiana-style) ===\n\n")
for (i in seq_len(nrow(v3_inflatable))) {
  r <- v3_inflatable[i, ]
  star <- if (!is.na(r$p) && r$p < 0.05) " *" else ""
  cat(sprintf("  %-8s %-18s  b = %+.3f (SE=%.3f)  p = %.4f%s\n",
              r$source, r$spec, r$coef, r$se, r$p, star))
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

get_coef <- function(m, name) {
  ct <- coeftable(m, vcov = NW(14))
  if (!(name %in% rownames(ct))) return(list(present = FALSE))
  r <- ct[name, ]
  list(present = TRUE, b = r[1], se = r[2],
       p = 2 * pnorm(-abs(r[1] / r[2])))
}

# Sample sizes after fixest drops all-zero FE cells.
N_utd <- nobs(v1_united)
N_iom <- nobs(v1_iom)

# Helper to render one (b, p)+(se) row pair for a triple V1/V2/V3.
row_triple <- function(label, models, name) {
  cs <- lapply(models, get_coef, name = name)
  bs <- sapply(cs, function(c) if (c$present) fcoef(c$b, c$p) else "")
  ss <- sapply(cs, function(c) if (c$present) fse(c$se)     else "")
  list(
    sprintf("%-31s & %-14s & %-14s & %-14s \\\\", label, bs[1], bs[2], bs[3]),
    sprintf("%-31s & %-14s & %-14s & %-14s \\\\", "",    ss[1], ss[2], ss[3])
  )
}

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[H]")
add("\\centering")
add("\\small")
add("\\caption{Volume-controlled Poisson with boat-composition controls. Boat-observable")
add(sprintf("sample: %s days (UNITED), %s days (IOM).}", fint(N_utd), fint(N_iom)))
add("\\label{tab:appx-boat}")
add("\\begin{tabular}{lccc}")
add("\\hline")
add("                                & V1: baseline & V2: + boat shares & V3: + boat $\\times$ SWH \\\\")
add("\\hline")
add("\\multicolumn{4}{l}{\\textit{UNITED}} \\\\")
models_u <- list(v1_united, v2_united, v3_united)
for (line in row_triple("SWH$_{t-1:t-5}$",                models_u, "swh_prev5days"))            add(line)
for (line in row_triple("SWH $\\times$ Post-MoU",         models_u, "swh_prev5days:post_mou"))   add(line)
for (line in row_triple("$\\log C_t$",                     models_u, "log_crossing_attempts"))    add(line)
for (line in row_triple("SWH $\\times$ inflatable share", models_u, "swh_prev5days:frx_inflatable_share")) add(line)
add("\\multicolumn{4}{l}{\\textit{IOM (comparison)}} \\\\")
models_i <- list(v1_iom, v2_iom, v3_iom)
for (line in row_triple("SWH$_{t-1:t-5}$",                models_i, "swh_prev5days"))            add(line)
for (line in row_triple("SWH $\\times$ Post-MoU",         models_i, "swh_prev5days:post_mou"))   add(line)
for (line in row_triple("$\\log C_t$",                     models_i, "log_crossing_attempts"))    add(line)
for (line in row_triple("SWH $\\times$ inflatable share", models_i, "swh_prev5days:frx_inflatable_share")) add(line)
add("\\hline")
add("Month-year FE                   & Yes            & Yes            & Yes            \\\\")
add("Newey-West SEs (lag 14)         & Yes            & Yes            & Yes            \\\\")
add("\\hline")
add("\\multicolumn{4}{l}{\\footnotesize Newey--West standard errors in parentheses. Stars: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
out_boat <- tbl_path("06_robustness", "05_rate_with_boat_controls.tex")
writeLines(L, out_boat)
cat(sprintf("  Saved: %s\n", out_boat))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
