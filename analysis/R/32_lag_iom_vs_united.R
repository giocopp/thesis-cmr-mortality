# 32_lag_iom_vs_united.R
# ======================
# Primary spec from 20_primary_model.R run across SWH lag windows
# (1, 2, 3, 4, 5, 7 days) on two death sources: IOM MMP vs UNITED.
#
# For each (source x window x family) combination fits:
#   n_dead ~ swh_win<k> + swh_win<k>:post_mou | month_year_fac
# with both NegBin (fenegbin) and Poisson QMLE (fepois), NW(14) SEs,
# sample = 05d primary filter (!is.na(lc_lag14) & !is.na(swh_win7)).
#
# UNITED is filtered to match IOM's primary as closely as possible:
#   - country_of_death in {Algeria, Italy, Libya, Malta, Tunisia, Mediterranean}
#   - manner_of_death in {"drowned", "other_unknown"}     [~ IOM Drowning + Mixed]
#   - spatial join to core corridor polygon (same polygon as build_iom_daily)
# IOM uses the default build_iom_daily() (incident only — split EXCLUDED — drown+mixed, central).
#
# Both windows exist in the panel (swh_lag1, swh_prev3days, swh_prev5days,
# swh_prevweek); 2-day and 4-day windows are computed on the fly.
#
# Out: output/tables/32_lag_iom_vs_united.txt
#      output/tables/32_lag_iom_vs_united.csv

library(tidyverse)
library(lubridate)
library(fixest)
library(sf)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("32  LAG GRID x IOM vs UNITED (primary spec from 20)\n")
cat("============================================================\n\n")

# ── 1. Load panel + compute all 6 lag-window SWH measures ────
cat("--- 1. Loading daily panel + building windows ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS")) %>%
  select(-n_dead_missing) %>%
  arrange(date) %>%
  mutate(
    swh_win1 = swh_lag1,        # existing: 1-day lag
    swh_win2 = zoo::rollmeanr(dplyr::lag(swh, 1), k = 2, fill = NA),
    swh_win3 = swh_prev3days,   # existing
    swh_win4 = zoo::rollmeanr(dplyr::lag(swh, 1), k = 4, fill = NA),
    swh_win5 = swh_prev5days,   # existing
    swh_win7 = swh_prevweek,    # existing
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
iom_daily <- build_iom_daily() %>% rename(n_dead_iom = n_dead_missing)
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

d <- panel %>%
  left_join(iom_daily,    by = "date") %>%
  left_join(united_daily, by = "date") %>%
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) %>%
  filter(!is.na(lc_lag14), !is.na(swh_win7))  # longest window => shared sample

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
cat("\n--- 5. Fitting 6 windows x 2 sources x 2 families = 24 models ---\n")

windows  <- c("swh_win1", "swh_win2", "swh_win3", "swh_win4",
              "swh_win5", "swh_win7")
outcomes <- c("n_dead_iom", "n_dead_united")
families <- c("NegBin", "Poisson")

grid <- expand_grid(outcome = outcomes, window = windows, family = families)

results <- pmap_dfr(grid, \(outcome, window, family)
  fit_combo(outcome, window, family)) %>%
  mutate(
    source = if_else(source == "n_dead_iom", "IOM", "UNITED"),
    lag_d  = as.integer(str_remove(window, "swh_win"))
  ) %>%
  arrange(source, family, lag_d)

# ── 7. Print results ──────────────────────────────────────────
cat("\n--- 6. Results (b3 = swh_win:post_mou) ---\n\n")

cat(sprintf("  %-7s %-8s %-4s  %+10s  %10s  %10s\n",
            "source", "family", "lag", "b3", "SE", "p"))
cat(sprintf("  %-7s %-8s %-4s  %10s  %10s  %10s\n",
            "------", "------", "---", "----------", "----------", "----------"))
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  star <- if (!is.na(r$p3) && r$p3 < 0.05) " *" else ""
  cat(sprintf("  %-7s %-8s %-4s  %+10.3f  %10.3f  %10.4f%s\n",
              r$source, r$family, sprintf("%dd", r$lag_d),
              r$b3, r$se3, r$p3, star))
}

# ── 8. Save output ────────────────────────────────────────────
cat("\n--- 7. Saving ---\n")

csv_path <- file.path(BASE_DIR, "output", "tables",
                        "32_lag_iom_vs_united.csv")
write.csv(results %>% select(-window), csv_path, row.names = FALSE)
cat(sprintf("Saved: %s\n", csv_path))

sink_file <- file.path(BASE_DIR, "output", "tables",
                        "32_lag_iom_vs_united.txt")
sink(sink_file)

cat("32  LAG GRID x IOM vs UNITED (primary spec)\n")
cat("Aligned with 20_primary_model.R\n")
cat("============================================\n\n")
cat("Spec: n_dead ~ swh_win<k> + swh_win<k>:post_mou | month_year_fac\n")
cat("NegBin (fenegbin) + Poisson QMLE (fepois), NW(14) SEs.\n")
cat("Shared sample (!is.na(lc_lag14) & !is.na(swh_win7)):\n")
cat(sprintf("  N = %d days, %s to %s.\n",
            nrow(d), min(d$date), max(d$date)))
cat(sprintf("  IOM deaths:    %.0f\n", sum(d$n_dead_iom)))
cat(sprintf("  UNITED deaths: %.0f\n", sum(d$n_dead_united)))
cat("\nUNITED filter: country in {Algeria, Italy, Libya, Malta, Tunisia,\n")
cat("   Mediterranean}, manner_of_death in {drowned, other_unknown},\n")
cat("   spatial join to core corridor polygon (same as build_iom_daily).\n\n")

cat("=== b1 (SWH pre-MoU slope) ===\n")
cat(sprintf("  %-7s %-8s %-4s  %+10s  %10s  %10s\n",
            "source", "family", "lag", "b1", "SE", "p"))
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  star <- if (!is.na(r$p1) && r$p1 < 0.05) " *" else ""
  cat(sprintf("  %-7s %-8s %-4s  %+10.3f  %10.3f  %10.4f%s\n",
              r$source, r$family, sprintf("%dd", r$lag_d),
              r$b1, r$se1, r$p1, star))
}

cat("\n=== b3 (SWH x post_MoU interaction) ===\n")
cat(sprintf("  %-7s %-8s %-4s  %+10s  %10s  %10s\n",
            "source", "family", "lag", "b3", "SE", "p"))
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  star <- if (!is.na(r$p3) && r$p3 < 0.05) " *" else ""
  cat(sprintf("  %-7s %-8s %-4s  %+10.3f  %10.3f  %10.4f%s\n",
              r$source, r$family, sprintf("%dd", r$lag_d),
              r$b3, r$se3, r$p3, star))
}

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
