# Primary spec from 20_primary_model.R run across SWH lag windows
# (1, 2, 3, 4, 5, 7 days) and future placebo windows (1, 3, 5, 7 days),
# on two death sources: UNITED primary, IOM comparison.
#
# For each (source x window x family) combination fits:
#   n_dead ~ swh_win<k> + swh_win<k>:post_mou | month_year_fac
# with both NegBin (fenegbin) and Poisson QMLE (fepois), NW(14) SEs,
# sample = primary filter plus longest lag/lead windows
# (!is.na(lc_lag14) & !is.na(swh_win7) & !is.na(swh_lead7)).
#
# UNITED is filtered to match IOM's primary as closely as possible:
#   - country_of_death in {Algeria, Italy, Libya, Malta, Tunisia, Mediterranean}
#   - manner_of_death in {"drowned", "other_unknown"}     [~ IOM Drowning + Mixed]
#   - spatial join to core corridor polygon (same polygon as build_iom_daily)
# IOM uses the default build_iom_daily() (incident only — split EXCLUDED — drown+mixed, central).
#
# Several lag windows exist in the panel (swh_lag1, swh_prev3days,
# swh_prev5days, swh_prevweek); 2-day, 4-day, and future windows are computed
# on the fly.
#
# Out: output/tables/32_lag_iom_vs_united.txt
#      output/tables/32_lag_iom_vs_united.csv

library(tidyverse)
library(lubridate)
library(fixest)
library(sf)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-02-02")

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("32  SWH WINDOW GRID x UNITED primary + IOM comparison\n")
cat("============================================================\n\n")

# ── 1. Load panel + compute all 6 lag-window SWH measures ────
cat("--- 1. Loading daily panel + building windows ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS")) |>
  select(-n_dead_missing) |>
  arrange(date) |>
  mutate(
    swh_win1 = swh_lag1,        # existing: 1-day lag
    swh_win2 = zoo::rollmeanr(dplyr::lag(swh, 1), k = 2, fill = NA),
    swh_win3 = swh_prev3days,   # existing
    swh_win4 = zoo::rollmeanr(dplyr::lag(swh, 1), k = 4, fill = NA),
    swh_win5 = swh_prev5days,   # existing
    swh_win7 = swh_prevweek,    # existing
    swh_lead1 = dplyr::lead(swh, 1),
    swh_lead3 = zoo::rollmean(dplyr::lead(swh, 1), k = 3,
                              fill = NA, align = "left"),
    swh_lead5 = zoo::rollmean(dplyr::lead(swh, 1), k = 5,
                              fill = NA, align = "left"),
    swh_lead7 = zoo::rollmean(dplyr::lead(swh, 1), k = 7,
                              fill = NA, align = "left"),
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag14 = dplyr::lag(zoo::rollmeanr(living_crossings, k = 7,
                                           fill = NA, align = "right"), 8),
    unit           = 1L,
    year           = year(date),
    month_year_fac = factor(month_year)
  )

cat(sprintf("  Panel: %d days (%s to %s)\n",
            nrow(panel), min(panel$date), max(panel$date)))

# ── 2. IOM deaths via shared helper (default primary filter) ──
cat("\n--- 2. IOM deaths (build_iom_daily defaults) ---\n")
iom_daily <- build_iom_daily() |> rename(n_dead_iom = n_dead_missing)
cat(sprintf("  IOM: %d death-days, %.0f deaths total\n",
            nrow(iom_daily), sum(iom_daily$n_dead_iom)))

# ── 3. UNITED deaths, matched to IOM primary ──────────────────
cat("\n--- 3. UNITED deaths (matched filter) ---\n")

# UNITED daily via the shared builder. Defaults (corridor spatial join;
# country in CMR+Med; manner drowned/other_unknown) replicate the previous
# inline filter exactly — single source of truth, see _helpers.R.
united_daily <- build_united_daily()
cat(sprintf("  UNITED (build_united_daily): %d death-days, %.0f deaths\n",
            nrow(united_daily), sum(united_daily$n_dead_united)))

# ── 4. Build modelling frame, fix sample for all windows ──────
cat("\n--- 4. Modelling frame ---\n")

d <- panel |>
  left_join(iom_daily,    by = "date") |>
  left_join(united_daily, by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) |>
  filter(!is.na(lc_lag14), !is.na(swh_win7), !is.na(swh_lead7))  # shared sample

cat(sprintf("  N = %d days (%s to %s)\n",
            nrow(d), min(d$date), max(d$date)))
cat(sprintf("  IOM total deaths in sample:    %.0f\n", sum(d$n_dead_iom)))
cat(sprintf("  UNITED total deaths in sample: %.0f\n", sum(d$n_dead_united)))

# ── 5. Fit one (source, window, family) combination ───────────
fit_combo <- function(outcome, window, family) {
  fn <- if (family == "NegBin") fenegbin else fepois
  fml <- as.formula(sprintf("%s ~ %s + %s:post_mou | month_year_fac",
                              outcome, window, window))
  m <- tryCatch(
    fn(fml, data = d, vcov = NW(14), panel.id = ~unit + date),
    error = function(e) NULL
  )
  if (is.null(m)) {
    return(tibble(source = outcome, family = family, window = window,
                   b1 = NA_real_, se1 = NA_real_, p1 = NA_real_,
                   b3 = NA_real_, se3 = NA_real_, p3 = NA_real_))
  }
  ct <- coeftable(m, vcov = NW(14))
  r1 <- which(rownames(ct) == window)
  r3 <- grep(":post_mou$", rownames(ct))
  r3 <- r3[grepl(window, rownames(ct)[r3], fixed = TRUE)]
  tibble(
    source = outcome,
    family = family,
    window = window,
    b1  = ct[r1, 1], se1 = ct[r1, 2],
    p1  = 2 * pnorm(-abs(ct[r1, 1] / ct[r1, 2])),
    b3  = ct[r3, 1], se3 = ct[r3, 2],
    p3  = 2 * pnorm(-abs(ct[r3, 1] / ct[r3, 2]))
  )
}

# ── 6. Run the full grid ──────────────────────────────────────
cat("\n--- 5. Fitting 10 windows x 2 sources x 2 families = 40 models ---\n")

window_defs <- tibble(
  window = c("swh_win1", "swh_win2", "swh_win3", "swh_win4", "swh_win5",
             "swh_win7", "swh_lead1", "swh_lead3", "swh_lead5", "swh_lead7"),
  timing = c(rep("past", 6), rep("future_placebo", 4)),
  window_label = c("lag 1d", "lag 1-2d", "lag 1-3d", "lag 1-4d", "lag 1-5d",
                   "lag 1-7d", "lead 1d", "lead 1-3d", "lead 1-5d",
                   "lead 1-7d")
)
outcomes <- c("n_dead_united", "n_dead_iom")
families <- c("NegBin", "Poisson")

grid <- expand_grid(outcome = outcomes, window_defs, family = families)

results <- pmap_dfr(grid, \(outcome, window, timing, window_label, family)
  fit_combo(outcome, window, family)) |>
  bind_cols(grid |> select(timing, window_label)) |>
  mutate(
    source = if_else(source == "n_dead_iom", "IOM", "UNITED"),
    source = factor(source, levels = c("UNITED", "IOM")),
    timing = factor(timing, levels = c("past", "future_placebo")),
    window_order = match(window, window_defs$window)
  ) |>
  arrange(source, timing, window_order, family)

# ── 7. Print results ──────────────────────────────────────────
cat("\n--- 6. Results (b3 = SWH window:post_mou) ---\n\n")

cat(sprintf("  %-7s %-14s %-10s %-8s  %+10s  %10s  %10s\n",
            "source", "timing", "window", "family", "b3", "SE", "p"))
cat(sprintf("  %-7s %-14s %-10s %-8s  %10s  %10s  %10s\n",
            "------", "------", "------", "------", "----------", "----------",
            "----------"))
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  star <- if (!is.na(r$p3) && r$p3 < 0.05) " *" else ""
  cat(sprintf("  %-7s %-14s %-10s %-8s  %+10.3f  %10.3f  %10.4f%s\n",
              as.character(r$source), as.character(r$timing), r$window_label,
              r$family,
              r$b3, r$se3, r$p3, star))
}

# ── 8. Save output ────────────────────────────────────────────
cat("\n--- 7. Saving ---\n")

csv_out <- tbl_path("06_robustness", "06_lag_iom_vs_united.csv")
write.csv(results, csv_out, row.names = FALSE)
cat(sprintf("Saved: %s\n", csv_out))

sink_file <- tbl_path("06_robustness", "06_lag_iom_vs_united.txt")
sink(sink_file)

cat("32  SWH WINDOW GRID x UNITED primary + IOM comparison (primary spec)\n")
cat("Aligned with 20_primary_model.R\n")
cat("============================================\n\n")
cat("Spec: n_dead ~ SWH_window + SWH_window:post_mou | month_year_fac\n")
cat("NegBin (fenegbin) + Poisson QMLE (fepois), NW(14) SEs.\n")
cat("Past windows are lagged SWH; future windows are placebo leads.\n")
cat("Shared sample (!is.na(lc_lag14) & !is.na(swh_win7) & !is.na(swh_lead7)):\n")
cat(sprintf("  N = %d days, %s to %s.\n",
            nrow(d), min(d$date), max(d$date)))
cat(sprintf("  UNITED deaths: %.0f\n", sum(d$n_dead_united)))
cat(sprintf("  IOM deaths:    %.0f\n", sum(d$n_dead_iom)))
cat("\nUNITED filter: country in {Algeria, Italy, Libya, Malta, Tunisia,\n")
cat("   Mediterranean}, manner_of_death in {drowned, other_unknown},\n")
cat("   spatial join to core corridor polygon (same as build_iom_daily).\n\n")

cat("=== b1 (SWH pre-MoU slope) ===\n")
cat(sprintf("  %-7s %-14s %-10s %-8s  %+10s  %10s  %10s\n",
            "source", "timing", "window", "family", "b1", "SE", "p"))
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  star <- if (!is.na(r$p1) && r$p1 < 0.05) " *" else ""
  cat(sprintf("  %-7s %-14s %-10s %-8s  %+10.3f  %10.3f  %10.4f%s\n",
              as.character(r$source), as.character(r$timing), r$window_label,
              r$family,
              r$b1, r$se1, r$p1, star))
}

cat("\n=== b3 (SWH x post_MoU interaction) ===\n")
cat(sprintf("  %-7s %-14s %-10s %-8s  %+10s  %10s  %10s\n",
            "source", "timing", "window", "family", "b3", "SE", "p"))
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  star <- if (!is.na(r$p3) && r$p3 < 0.05) " *" else ""
  cat(sprintf("  %-7s %-14s %-10s %-8s  %+10.3f  %10.3f  %10.4f%s\n",
              as.character(r$source), as.character(r$timing), r$window_label,
              r$family,
              r$b3, r$se3, r$p3, star))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
