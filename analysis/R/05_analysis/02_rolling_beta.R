# ── Rolling-window beta(SWH) ────────────────────────────────────────────────
# Poisson QMLE n_dead_iom ~ swh_prev5days | month_year_fac per window, NW(14).
# 2-year window stepped weekly; corridor-wide + AFR/EU bloc panels.

library(tidyverse)
library(fixest)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

STEP_DAYS      <- 7
MIN_DEATH_DAYS <- 30
MIN_ROWS       <- 200
WINDOW_PRIMARY <- 730

cat("============================================================\n")
cat("053  ROLLING-WINDOW BETA(SWH)\n")
cat("============================================================\n\n")

# ── 1. Load data + 05d data prep ───────────────────────────
cat("--- 1. Loading panels + 05d data prep ---\n")

panel_daily <- readRDS(file.path(BASE_DIR, "analysis", "data",
                                   "daily_panel_complete.RDS"))

iom_daily <- build_iom_daily()

panel_daily <- panel_daily |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  replace_na(list(n_dead_iom = 0)) |>
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag7  = dplyr::lag(zoo::rollmeanr(living_crossings, k = 7,
                                           fill = NA, align = "right"), 1),
    lc_lag14 = dplyr::lag(zoo::rollmeanr(living_crossings, k = 7,
                                           fill = NA, align = "right"), 8),
    unit           = 1L,
    year_fac       = factor(year),
    month_year_fac = factor(month_year)
  )

# 05d primary sample
sample_dates <- panel_daily |>
  filter(!is.na(lc_lag14), !is.na(swh_prev5days)) |>
  pull(date)

da <- panel_daily |>
  filter(date %in% sample_dates) |>
  arrange(date)

zp <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_zone.RDS")) |>
  filter(date %in% sample_dates) |>
  rename(n_dead_iom = n_dead_missing)

# Collapse zone -> 2 blocs: sum deaths, first SWH (SWH constant across
# zones under A2).
bloc <- zp |>
  group_by(date, sar_bloc) |>
  summarise(n_dead_iom   = sum(n_dead_iom),
            swh_prev5days = first(swh_prev5days),
            .groups = "drop") |>
  mutate(year = year(date))
dim(bloc$date) <- NULL

cat(sprintf("  daily-agg range:  %s to %s (N = %d)\n",
    min(da$date), max(da$date), nrow(da)))
cat(sprintf("  2-bloc panel:     %d rows (%d dates x 2 blocs)\n",
    nrow(bloc), length(unique(bloc$date))))

# ── 2. Rolling estimator ────────────────────────────────────

run_rolling <- function(df, wind = WINDOW_PRIMARY, step = STEP_DAYS) {
  fit_one <- function(sub, mid_date) {
    if (nrow(sub) < MIN_ROWS) return(NULL)
    if (sum(sub$n_dead_iom > 0) < MIN_DEATH_DAYS) return(NULL)

    sub$unit <- 1L

    if (length(unique(sub$month_year_fac)) < 2) return(NULL)

    m <- tryCatch(
      fepois(n_dead_iom ~ swh_prev5days | month_year_fac,
             data = sub, vcov = NW(14), panel.id = ~unit + date),
      error = function(e) NULL
    )
    if (is.null(m) || !"swh_prev5days" %in% names(coef(m))) return(NULL)

    ct <- tryCatch(
      coeftable(m, vcov = NW(14)),
      error = function(e) NULL
    )
    if (is.null(ct)) return(NULL)
    row <- which(rownames(ct) == "swh_prev5days")
    tibble(date_mid = mid_date,
           beta     = ct[row, 1],
           se       = ct[row, 2],
           n_obs    = nrow(sub),
           n_death_days = sum(sub$n_dead_iom > 0))
  }

  dates_all <- sort(unique(df$date))
  start_pool <- seq(min(dates_all) + wind / 2,
                     max(dates_all) - wind / 2,
                     by = step)

  map_dfr(start_pool, function(mid) {
    lo <- mid - wind / 2
    hi <- mid + wind / 2
    sub <- df |> filter(date >= lo, date <= hi) |>
      mutate(month_year_fac = factor(format(date, "%Y-%m")))
    fit_one(sub, mid)
  })
}

# ── 3. Run rolling by flavor ────────────────────────────────
cat("\n--- 2. Running rolling estimators ---\n")

# (1) Daily-agg: 2-year window
cat("\n  [1] daily-agg, 2-year window ...\n")
t0 <- Sys.time()
res_da_2y <- run_rolling(da, wind = WINDOW_PRIMARY) |>
  mutate(flavor = "daily-agg", label = "daily-agg",
         window = "2-year")
cat(sprintf("    %d windows in %.1f s\n",
    nrow(res_da_2y), as.numeric(Sys.time() - t0, units = "secs")))

# (2) 2-bloc: per-bloc rolling, 2-year window
cat("\n  [2] 2-bloc (AFR, EU), 2-year window ...\n")
res_bloc_list <- list()
for (bl in c("AFR", "EU")) {
  t0 <- Sys.time()
  sub_bl <- bloc |> filter(sar_bloc == bl)
  r <- run_rolling(sub_bl, wind = WINDOW_PRIMARY) |>
    mutate(flavor = "2-bloc", label = bl, window = "2-year")
  cat(sprintf("    %s: %d windows in %.1f s\n",
      bl, nrow(r), as.numeric(Sys.time() - t0, units = "secs")))
  res_bloc_list[[bl]] <- r
}
res_bloc <- bind_rows(res_bloc_list)

# ── 4. Combine and save ─────────────────────────────────────
cat("\n--- 3. Save tables ---\n")

all_res <- bind_rows(res_da_2y, res_bloc)
saveRDS(all_res, tbl_path("05_analysis", "02_rolling_beta.RDS"))

# Expected (flavor, label, window) template so any series with zero windows
# shows up explicitly with n_windows = 0 instead of being silently dropped.
expected_series <- bind_rows(
  tibble(flavor = "daily-agg", label = "daily-agg", window = "2-year"),
  tibble(flavor = "2-bloc",    label = c("AFR", "EU"), window = "2-year")
)

# Per-series summary with pre-MoU and post-MoU means.
#   "Fully pre"  = window's RIGHT edge at or before MoU (date_mid <= MoU - W/2)
#   "Fully post" = window's LEFT  edge at or after  MoU (date_mid >= MoU + W/2)
# Windows that straddle MoU are excluded from both means. The means are
# unweighted across windows; adjacent windows overlap heavily (week-stepped
# 2-year windows share 723 of 730 days), so beta_diff is descriptive only —
# do NOT compute a naive t-test on it.
HALF_WIND <- WINDOW_PRIMARY / 2

per_series <- all_res |>
  group_by(flavor, label, window) |>
  group_modify(~ {
    fully_pre  <- .x |> filter(date_mid <= MOU_DATE - HALF_WIND)
    fully_post <- .x |> filter(date_mid >= MOU_DATE + HALF_WIND)
    tibble(
      n_windows      = nrow(.x),
      n_fully_pre    = nrow(fully_pre),
      n_fully_post   = nrow(fully_post),
      date_range     = sprintf("%s to %s",
                                min(.x$date_mid), max(.x$date_mid)),
      beta_min       = round(min(.x$beta), 3),
      beta_max       = round(max(.x$beta), 3),
      beta_pre_mean  = if (nrow(fully_pre)  > 0) round(mean(fully_pre$beta),  3) else NA_real_,
      beta_post_mean = if (nrow(fully_post) > 0) round(mean(fully_post$beta), 3) else NA_real_
    )
  }) |>
  ungroup()

window_summary <- expected_series |>
  left_join(per_series, by = c("flavor", "label", "window")) |>
  mutate(
    n_windows    = coalesce(n_windows, 0L),
    n_fully_pre  = coalesce(n_fully_pre, 0L),
    n_fully_post = coalesce(n_fully_post, 0L),
    beta_diff    = round(beta_post_mean - beta_pre_mean, 3)
  )

# Loud warning if any series produced zero windows.
zero_series <- window_summary |> filter(n_windows == 0)
if (nrow(zero_series) > 0) {
  cat("\n  WARNING: the following series produced 0 rolling windows ",
      "(skip rules: <", MIN_DEATH_DAYS, " death-days OR <", MIN_ROWS, " rows):\n", sep = "")
  for (i in seq_len(nrow(zero_series))) {
    z <- zero_series[i, ]
    cat(sprintf("    %-12s %-16s %-6s\n", z$flavor, z$label, z$window))
  }
}

sink_file <- tbl_path("05_analysis", "02_rolling_beta.txt")
sink(sink_file)
old_opts <- options(tibble.width = Inf, tibble.print_max = Inf)
on.exit(options(old_opts), add = TRUE)
cat("053  ROLLING-WINDOW BETA(SWH)\n")
cat("Aligned with 20_primary_model.R (Poisson arm)\n")
cat("=============================================\n")
cat(sprintf("Window: %d days (2y), step %d days, anchor = window midpoint.\n",
    WINDOW_PRIMARY, STEP_DAYS))
cat("Model: fepois(n_dead_iom ~ swh_prev5days | month_year_fac), NW(14) SEs.\n")
cat("Sample: 05d primary (!is.na(lc_lag14) & !is.na(swh_prev5days)).\n")
cat(sprintf("Skip rules: <%d death-days OR <%d rows OR <2 month_year levels.\n",
    MIN_DEATH_DAYS, MIN_ROWS))
cat("\n")
cat("Pre-MoU mean : avg beta over windows whose RIGHT edge is <= MoU\n")
cat("               (date_mid <= ", as.character(MOU_DATE), " - W/2).\n", sep = "")
cat("Post-MoU mean: avg beta over windows whose LEFT edge is >= MoU\n")
cat("               (date_mid >= ", as.character(MOU_DATE), " + W/2).\n", sep = "")
cat("Windows that straddle MoU are excluded from both means.\n")
cat("Means are unweighted across heavily-overlapping windows; beta_diff\n")
cat("is descriptive only — do NOT compute a naive t-test on it.\n\n")

cat("Per-series summary:\n")
print(window_summary, n = Inf)

if (nrow(zero_series) > 0) {
  cat("\nSeries with 0 windows (failed skip rules):\n")
  print(zero_series |> select(flavor, label, window), n = Inf)
}

cat("\nFull rolling estimates are in 22_rolling_beta.RDS\n")
sink()
cat(sprintf("Saved: %s\n", sink_file))

# ── 5. Plot ─────────────────────────────────────────────────
cat("\n--- 4. Plot ---\n")

# Break the line at large time gaps. EU 2-bloc has multi-year gaps where
# no 2-year window had enough death-days; without breaking, geom_line
# interpolates a misleading straight line across them. We assign a
# "segment" id that increments every time the gap between consecutive
# windows exceeds 2 steps (14 days). ggplot's group aesthetic then draws
# each segment separately.
add_segments <- function(df, gap_threshold = 14) {
  df |>
    group_by(flavor, label, window) |>
    arrange(date_mid) |>
    mutate(gap = c(0, as.numeric(diff(date_mid))),
           segment = cumsum(gap > gap_threshold)) |>
    ungroup() |>
    mutate(series_id = paste(flavor, label, window, segment, sep = "_"))
}

plot_df <- all_res |>
  mutate(ci_lo = beta - 1.96 * se,
         ci_hi = beta + 1.96 * se) |>
  add_segments()

# NOTE: post-hoc rolling-mean smoothing was removed (2026-04-13). The 730-day
# (2y) window already smooths heavily; double-smoothing the rolling estimates
# obscured the actual noise level and made the CI ribbons visually misleading.
# Rolling estimates and CIs are now plotted as-is.

# Common x-axis range across both panels for visual comparability.
x_min <- min(plot_df$date_mid, na.rm = TRUE)
x_max <- max(plot_df$date_mid, na.rm = TRUE)
x_lims <- as.Date(c(x_min, x_max))

# Report the gaps in the log
gap_summary <- plot_df |>
  group_by(flavor, label, window) |>
  summarise(n_segments = n_distinct(segment),
            max_gap = max(gap, na.rm = TRUE),
            .groups = "drop") |>
  filter(n_segments > 1)
if (nrow(gap_summary) > 0) {
  cat("  Gaps detected (windows with insufficient death-days):\n")
  for (i in seq_len(nrow(gap_summary))) {
    g <- gap_summary[i, ]
    cat(sprintf("    %-12s %-20s %-6s: %d segments, max gap = %d days\n",
        g$flavor, g$label, g$window, g$n_segments, g$max_gap))
  }
}

# Pre/post mean segments for the daily-agg panel — one per window size.
# These overlay the rolling line as horizontal dashed segments showing the
# average rolling beta in the fully-pre-MoU and fully-post-MoU regions.
# Visually answers "did the average move?" without misleading the reader
# into reading a discrete jump that isn't there.
x_min_da <- min(plot_df$date_mid[plot_df$flavor == "daily-agg"], na.rm = TRUE)
x_max_da <- max(plot_df$date_mid[plot_df$flavor == "daily-agg"], na.rm = TRUE)

da_summary <- window_summary |> filter(flavor == "daily-agg")

mean_segs <- bind_rows(
  da_summary |>
    filter(!is.na(beta_pre_mean)) |>
    transmute(x    = x_min_da,
              xend = MOU_DATE - HALF_WIND,
              y    = beta_pre_mean,
              yend = beta_pre_mean,
              side = "pre"),
  da_summary |>
    filter(!is.na(beta_post_mean)) |>
    transmute(x    = MOU_DATE + HALF_WIND,
              xend = x_max_da,
              y    = beta_post_mean,
              yend = beta_post_mean,
              side = "post")
)

# Panel (1): daily-agg, 2-year rolling beta
p1 <- plot_df |>
  filter(flavor == "daily-agg") |>
  ggplot(aes(x = date_mid, y = beta, group = series_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dotted",
             colour = "#D32F2F", linewidth = 0.6) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.15, fill = "#2166AC", colour = NA) +
  geom_line(linewidth = 0.7, colour = "#2166AC") +
  geom_segment(data = mean_segs,
               aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 1.1, linetype = "longdash", colour = "#2166AC",
               inherit.aes = FALSE) +
  scale_x_date(limits = x_lims, date_breaks = "2 years",
               date_labels = "%Y") +
  labs(
    title = "(1) Daily-aggregate rolling beta (2-year window)",
    subtitle = paste("Poisson QMLE, NW(14).",
                     "Long-dashed segments = mean beta over windows fully pre/post MoU.",
                     "Red dotted = MoU (2017-07-01).", sep = " "),
    x = NULL, y = expression(beta(SWH[prevweek]))
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

# Panel (2): 2-bloc (AFR, EU)
p2 <- plot_df |>
  filter(flavor == "2-bloc") |>
  ggplot(aes(x = date_mid, y = beta, colour = label, fill = label,
             group = series_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dotted",
             colour = "#D32F2F", linewidth = 0.6) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("AFR" = "#D6604D", "EU" = "#2166AC")) +
  scale_fill_manual(values   = c("AFR" = "#D6604D", "EU" = "#2166AC")) +
  scale_x_date(limits = x_lims, date_breaks = "2 years",
               date_labels = "%Y") +
  labs(
    title = "(2) 2-bloc: African SAR vs EU SAR (2y window)",
    subtitle = "AFR = Libya + Tunisia; EU = Italy + Malta. Raw rolling estimates.",
    x = NULL, y = expression(beta(SWH[prevweek])),
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

combined_plot <- p1 / p2

fig_out <- fig_path("05_analysis", "02_rolling_beta.png")
ggsave(fig_out, combined_plot, width = 10, height = 8, dpi = 200)
cat(sprintf("Saved: %s\n", fig_out))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
