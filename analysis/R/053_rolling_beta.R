# 053_rolling_beta.R
# ==================
# Enhancement #4: ROLLING-WINDOW beta(SWH).
#
# Motivation: the user explicitly rejected sharp placebo-date tests. The
# assumption that MoU produced a discrete June/August 2017 break is
# implausible — any real structural change should be GRADUAL. A rolling
# window visualises the smoothness directly.
#
# Specification:
#   - Windows: 730 days (2 years) primary + 365 days (1 year) robustness
#   - Step:    7 days (weekly)
#   - Model:   fepois (Poisson QMLE) — faster, more stable in short windows
#              than fenegbin, and fine asymptotically
#   - Formula: n_dead_missing ~ swh_prevweek | month_year
#              (main effect only, NO post_mou, NO interaction)
#   - SE:      NW(14) per window
#   - Skip:    windows with <30 death-days or <200 total days
#
# Flavors (4 panels in the output figure):
#   (1) Daily-agg — one overall series, with 2-year AND 1-year windows
#       as a robustness check on window-size sensitivity
#   (2) 2-bloc    — AFR vs EU SAR (2-year window, two lines)
#   (3) 4-country — Libya, Tunisia, Italy, Malta (2-year window, four lines)
#
# Purpose of each panel:
#   (1) shows window-size robustness (same data, different smoothing)
#   (2) shows bloc-level heterogeneity
#   (3) shows country-level heterogeneity — which country drives the shift?
#
# Output: output/tables/053_rolling_beta.txt
#         output/tables/053_rolling_beta.RDS  (full long-format data)
#         output/figures/053_rolling_beta.png

library(tidyverse)
library(fixest)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
STEP_DAYS <- 7
MIN_DEATH_DAYS <- 30
MIN_ROWS <- 200
WINDOW_PRIMARY <- 730
WINDOW_ROBUST  <- 365

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("053  ROLLING-WINDOW BETA(SWH)\n")
cat("============================================================\n\n")

# ── 1. Load data ────────────────────────────────────────────
cat("--- 1. Loading panels ---\n")

# Drop the panel's broad n_dead_missing and replace with the analytical
# series via the shared helper. Default = incident-only, core corridor,
# all causes. Change the call to test sensitivity variants.
da <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_complete.RDS")) %>%
  select(-n_dead_missing) %>%
  left_join(build_iom_daily(), by = "date") %>%
  replace_na(list(n_dead_missing = 0)) %>%
  arrange(date)

zp <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_zone.RDS"))

# Collapse zone -> 2 blocs: sum deaths, mean SWH (SWH constant across
# zones under A2).
bloc <- zp %>%
  group_by(date, sar_bloc) %>%
  summarise(n_dead_missing = sum(n_dead_missing),
            swh_prevweek   = mean(swh_prevweek, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(year = year(date))
dim(bloc$date) <- NULL

cat(sprintf("  daily-agg range: %s to %s (N = %d)\n",
    min(da$date), max(da$date), nrow(da)))
cat(sprintf("  zone 4-country range: %s to %s (N = %d, 4 zones)\n",
    min(zp$date), max(zp$date), nrow(zp)))

# ── 2. Rolling estimator ────────────────────────────────────

run_rolling <- function(df, wind = WINDOW_PRIMARY, step = STEP_DAYS) {
  fit_one <- function(sub, mid_date) {
    if (nrow(sub) < MIN_ROWS) return(NULL)
    if (sum(sub$n_dead_missing > 0) < MIN_DEATH_DAYS) return(NULL)

    sub$unit <- 1L

    if (length(unique(sub$month_year)) < 2) return(NULL)

    m <- tryCatch(
      fepois(n_dead_missing ~ swh_prevweek | month_year,
             data = sub, vcov = NW(14), panel.id = ~unit + date),
      error = function(e) NULL
    )
    if (is.null(m) || !"swh_prevweek" %in% names(coef(m))) return(NULL)

    ct <- tryCatch(
      coeftable(m, vcov = NW(14)),
      error = function(e) NULL
    )
    if (is.null(ct)) return(NULL)
    row <- which(rownames(ct) == "swh_prevweek")
    tibble(date_mid = mid_date,
           beta     = ct[row, 1],
           se       = ct[row, 2],
           n_obs    = nrow(sub),
           n_death_days = sum(sub$n_dead_missing > 0))
  }

  dates_all <- sort(unique(df$date))
  start_pool <- seq(min(dates_all) + wind / 2,
                     max(dates_all) - wind / 2,
                     by = step)

  map_dfr(start_pool, function(mid) {
    lo <- mid - wind / 2
    hi <- mid + wind / 2
    sub <- df %>% filter(date >= lo, date <= hi) %>%
      mutate(month_year = factor(format(date, "%Y-%m")))
    fit_one(sub, mid)
  })
}

# ── 3. Run rolling by flavor ────────────────────────────────
cat("\n--- 2. Running rolling estimators ---\n")

# (1) Daily-agg: 2-year and 1-year windows
cat("\n  [1] daily-agg, 2-year window ...\n")
t0 <- Sys.time()
res_da_2y <- run_rolling(da, wind = WINDOW_PRIMARY) %>%
  mutate(flavor = "daily-agg", label = "daily-agg (2y)",
         window = "2-year")
cat(sprintf("    %d windows in %.1f s\n",
    nrow(res_da_2y), as.numeric(Sys.time() - t0, units = "secs")))

cat("  [1] daily-agg, 1-year window (robustness) ...\n")
t0 <- Sys.time()
res_da_1y <- run_rolling(da, wind = WINDOW_ROBUST) %>%
  mutate(flavor = "daily-agg", label = "daily-agg (1y)",
         window = "1-year")
cat(sprintf("    %d windows in %.1f s\n",
    nrow(res_da_1y), as.numeric(Sys.time() - t0, units = "secs")))

# (2) 2-bloc: per-bloc rolling, 2-year window
cat("\n  [2] 2-bloc (AFR, EU), 2-year window ...\n")
res_bloc_list <- list()
for (bl in c("AFR", "EU")) {
  t0 <- Sys.time()
  sub_bl <- bloc %>% filter(sar_bloc == bl)
  r <- run_rolling(sub_bl, wind = WINDOW_PRIMARY) %>%
    mutate(flavor = "2-bloc", label = bl, window = "2-year")
  cat(sprintf("    %s: %d windows in %.1f s\n",
      bl, nrow(r), as.numeric(Sys.time() - t0, units = "secs")))
  res_bloc_list[[bl]] <- r
}
res_bloc <- bind_rows(res_bloc_list)

# (3) 4-country: per-country rolling, 2-year window
cat("\n  [3] 4-country (Libya, Tunisia, Italy, Malta), 2-year window ...\n")
res_country_list <- list()
for (ct in c("Libya", "Tunisia", "Italy", "Malta")) {
  t0 <- Sys.time()
  sub_ct <- zp %>% filter(country == ct)
  r <- run_rolling(sub_ct, wind = WINDOW_PRIMARY) %>%
    mutate(flavor = "4-country", label = ct, window = "2-year")
  cat(sprintf("    %s: %d windows in %.1f s\n",
      ct, nrow(r), as.numeric(Sys.time() - t0, units = "secs")))
  res_country_list[[ct]] <- r
}
res_country <- bind_rows(res_country_list)

# ── 4. Combine and save ─────────────────────────────────────
cat("\n--- 3. Save tables ---\n")

all_res <- bind_rows(res_da_2y, res_da_1y, res_bloc, res_country)
saveRDS(all_res, file.path(BASE_DIR, "output", "tables", "053_rolling_beta.RDS"))

# Expected (flavor, label, window) template so any series with zero windows
# (e.g., Italy when its zone has too few death-days for the threshold) shows
# up explicitly with n_windows = 0 instead of being silently dropped.
expected_series <- bind_rows(
  tibble(flavor = "daily-agg", label = "daily-agg (2y)", window = "2-year"),
  tibble(flavor = "daily-agg", label = "daily-agg (1y)", window = "1-year"),
  tibble(flavor = "2-bloc",    label = c("AFR", "EU"),   window = "2-year"),
  tibble(flavor = "4-country",
         label = c("Libya", "Tunisia", "Italy", "Malta"),
         window = "2-year")
)

# Per-series summary with pre-MoU and post-MoU means.
#   "Fully pre"  = window's RIGHT edge at or before MoU (date_mid <= MoU - W/2)
#   "Fully post" = window's LEFT  edge at or after  MoU (date_mid >= MoU + W/2)
# Windows that straddle MoU are excluded from both means. The means are
# unweighted across windows; adjacent windows overlap heavily (week-stepped
# 2-year windows share 723 of 730 days), so beta_diff is descriptive only —
# do NOT compute a naive t-test on it.
half_days <- function(w) if (w == "2-year") WINDOW_PRIMARY / 2 else WINDOW_ROBUST / 2

per_series <- all_res %>%
  group_by(flavor, label, window) %>%
  group_modify(~ {
    half <- half_days(.y$window)
    fully_pre  <- .x %>% filter(date_mid <= MOU_DATE - half)
    fully_post <- .x %>% filter(date_mid >= MOU_DATE + half)
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
  }) %>%
  ungroup()

window_summary <- expected_series %>%
  left_join(per_series, by = c("flavor", "label", "window")) %>%
  mutate(
    n_windows    = coalesce(n_windows, 0L),
    n_fully_pre  = coalesce(n_fully_pre, 0L),
    n_fully_post = coalesce(n_fully_post, 0L),
    beta_diff    = round(beta_post_mean - beta_pre_mean, 3)
  )

# Loud warning if any series produced zero windows (Italy is the recurring
# culprit on the 4-country panel).
zero_series <- window_summary %>% filter(n_windows == 0)
if (nrow(zero_series) > 0) {
  cat("\n  WARNING: the following series produced 0 rolling windows ",
      "(skip rules: <", MIN_DEATH_DAYS, " death-days OR <", MIN_ROWS, " rows):\n", sep = "")
  for (i in seq_len(nrow(zero_series))) {
    z <- zero_series[i, ]
    cat(sprintf("    %-12s %-16s %-6s\n", z$flavor, z$label, z$window))
  }
}

sink_file <- file.path(BASE_DIR, "output", "tables", "053_rolling_beta.txt")
sink(sink_file)
old_opts <- options(tibble.width = Inf, tibble.print_max = Inf)
on.exit(options(old_opts), add = TRUE)
cat("053  ROLLING-WINDOW BETA(SWH)\n")
cat("=============================\n")
cat(sprintf("Primary window: %d days (2y). Robustness window: %d days (1y).\n",
    WINDOW_PRIMARY, WINDOW_ROBUST))
cat(sprintf("Step: %d days. Anchor: window midpoint.\n", STEP_DAYS))
cat("Model: fepois(n_dead_missing ~ swh_prevweek | month_year), NW(14) SEs.\n")
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
  print(zero_series %>% select(flavor, label, window), n = Inf)
}

cat("\nFull rolling estimates are in 053_rolling_beta.RDS\n")
sink()
cat(sprintf("Saved: %s\n", sink_file))

# ── 5. Plot ─────────────────────────────────────────────────
cat("\n--- 4. Plot ---\n")

# Break the line at large time gaps. Malta 4-country and EU 2-bloc have
# multi-year gaps where no 2-year window had enough death-days; without
# breaking, geom_line interpolates a misleading straight line across them.
# We assign a "segment" id that increments every time the gap between
# consecutive windows exceeds 2 steps (14 days). ggplot's group aesthetic
# then draws each segment separately.
add_segments <- function(df, gap_threshold = 14) {
  df %>%
    group_by(flavor, label, window) %>%
    arrange(date_mid) %>%
    mutate(gap = c(0, as.numeric(diff(date_mid))),
           segment = cumsum(gap > gap_threshold)) %>%
    ungroup() %>%
    mutate(series_id = paste(flavor, label, window, segment, sep = "_"))
}

plot_df <- all_res %>%
  mutate(ci_lo = beta - 1.96 * se,
         ci_hi = beta + 1.96 * se) %>%
  add_segments()

# NOTE: post-hoc rolling-mean smoothing was removed (2026-04-13). The 730-day
# (2y) window already smooths heavily; double-smoothing the rolling estimates
# obscured the actual noise level and made the CI ribbons visually misleading.
# Rolling estimates and CIs are now plotted as-is.

# Common x-axis range across all three panels for visual comparability.
x_min <- min(plot_df$date_mid, na.rm = TRUE)
x_max <- max(plot_df$date_mid, na.rm = TRUE)
x_lims <- as.Date(c(x_min, x_max))

# Report the gaps in the log
gap_summary <- plot_df %>%
  group_by(flavor, label, window) %>%
  summarise(n_segments = n_distinct(segment),
            max_gap = max(gap, na.rm = TRUE),
            .groups = "drop") %>%
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

da_summary <- window_summary %>%
  filter(flavor == "daily-agg") %>%
  mutate(half_d = if_else(window == "2-year", WINDOW_PRIMARY / 2, WINDOW_ROBUST / 2))

mean_segs <- bind_rows(
  da_summary %>%
    filter(!is.na(beta_pre_mean)) %>%
    transmute(window,
              x    = x_min_da,
              xend = MOU_DATE - half_d,
              y    = beta_pre_mean,
              yend = beta_pre_mean,
              side = "pre"),
  da_summary %>%
    filter(!is.na(beta_post_mean)) %>%
    transmute(window,
              x    = MOU_DATE + half_d,
              xend = x_max_da,
              y    = beta_post_mean,
              yend = beta_post_mean,
              side = "post")
)

# Panel (1): daily-agg with 2y and 1y windows
p1 <- plot_df %>%
  filter(flavor == "daily-agg") %>%
  ggplot(aes(x = date_mid, y = beta, colour = window, fill = window,
             group = series_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dotted",
             colour = "#D32F2F", linewidth = 0.6) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_segment(data = mean_segs,
               aes(x = x, xend = xend, y = y, yend = yend, colour = window),
               linewidth = 1.1, linetype = "longdash",
               inherit.aes = FALSE) +
  scale_colour_manual(values = c("2-year" = "#2166AC", "1-year" = "#D6604D")) +
  scale_fill_manual(values   = c("2-year" = "#2166AC", "1-year" = "#D6604D")) +
  scale_x_date(limits = x_lims, date_breaks = "2 years",
               date_labels = "%Y") +
  labs(
    title = "(1) Daily-aggregate rolling beta: window-size robustness",
    subtitle = paste("Poisson QMLE, NW(14). Raw rolling estimates (no post-hoc smoothing).",
                     "Long-dashed segments = mean beta over windows fully pre/post MoU.",
                     "Red dotted = MoU (2017-07-01).", sep = " "),
    x = NULL, y = expression(beta(SWH[prevweek])),
    colour = "Window", fill = "Window"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

# Panel (2): 2-bloc (AFR, EU)
p2 <- plot_df %>%
  filter(flavor == "2-bloc") %>%
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

# Panel (3): 4-country (Libya, Tunisia, Italy, Malta) with CI ribbons
country_colours <- c(
  "Libya"   = "#B2182B",
  "Tunisia" = "#F4A582",
  "Italy"   = "#2166AC",
  "Malta"   = "#4393C3"
)

p3 <- plot_df %>%
  filter(flavor == "4-country") %>%
  ggplot(aes(x = date_mid, y = beta, colour = label, fill = label,
             group = series_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dotted",
             colour = "#D32F2F", linewidth = 0.6) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = country_colours) +
  scale_fill_manual(values = country_colours) +
  scale_x_date(limits = x_lims, date_breaks = "2 years",
               date_labels = "%Y") +
  labs(
    title = "(3) 4-country: per-SAR rolling beta (2y window)",
    subtitle = "Libya, Tunisia, Italy, Malta estimated separately. 95% CI ribbons, raw rolling estimates.",
    x = NULL, y = expression(beta(SWH[prevweek])),
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

combined_plot <- p1 / p2 / p3

fig_path <- file.path(BASE_DIR, "output", "figures", "053_rolling_beta.png")
ggsave(fig_path, combined_plot, width = 10, height = 12, dpi = 200)
cat(sprintf("Saved: %s\n", fig_path))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
