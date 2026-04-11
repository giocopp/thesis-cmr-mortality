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
#   - Formula: n_dead_missing ~ swh_prevweek_z | month_year
#              (main effect only, NO post_mou, NO interaction)
#   - SE:      NW(28) per window
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
library(zoo)

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

    sub$swh_prevweek_z <- as.numeric(scale(sub$swh_prevweek))
    sub$unit <- 1L

    if (length(unique(sub$month_year)) < 2) return(NULL)

    m <- tryCatch(
      fepois(n_dead_missing ~ swh_prevweek_z | month_year,
             data = sub, vcov = NW(28), panel.id = ~unit + date),
      error = function(e) NULL
    )
    if (is.null(m) || !"swh_prevweek_z" %in% names(coef(m))) return(NULL)

    ct <- tryCatch(
      coeftable(m, vcov = NW(28)),
      error = function(e) NULL
    )
    if (is.null(ct)) return(NULL)
    row <- which(rownames(ct) == "swh_prevweek_z")
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

sink_file <- file.path(BASE_DIR, "output", "tables", "053_rolling_beta.txt")
sink(sink_file)
cat("053  ROLLING-WINDOW BETA(SWH)\n")
cat("=============================\n")
cat(sprintf("Primary window: %d days (2y). Robustness window: %d days (1y).\n",
    WINDOW_PRIMARY, WINDOW_ROBUST))
cat(sprintf("Step: %d days. Anchor: window midpoint.\n", STEP_DAYS))
cat("Model: fepois(n_dead_missing ~ swh_prevweek_z | month_year), NW(28) SEs.\n\n")

cat("Summary by flavor / label / window:\n")
print(all_res %>%
        group_by(flavor, label, window) %>%
        summarise(n_windows = n(),
                  date_range = sprintf("%s to %s",
                                        min(date_mid), max(date_mid)),
                  beta_min = round(min(beta), 3),
                  beta_max = round(max(beta), 3),
                  beta_at_mou = {
                    idx <- which.min(abs(date_mid - MOU_DATE))
                    round(beta[idx], 3)
                  },
                  .groups = "drop"),
      n = 20)

cat("\nFull table is in 053_rolling_beta.RDS\n")
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

# Light rolling-mean smoothing (k = 5 points = ~5 weeks of weekly-stepped
# estimates). Applied WITHIN each series_id so the smoothing does not leak
# across the multi-year gaps in Malta / EU. Short segments (< k points)
# are left unsmoothed.
smooth_cols <- function(df, k = 5) {
  df %>%
    group_by(series_id) %>%
    arrange(date_mid) %>%
    mutate(
      beta  = if (dplyr::n() >= k) zoo::rollmean(beta,  k = k, fill = NA, align = "center") else beta,
      ci_lo = if (dplyr::n() >= k) zoo::rollmean(ci_lo, k = k, fill = NA, align = "center") else ci_lo,
      ci_hi = if (dplyr::n() >= k) zoo::rollmean(ci_hi, k = k, fill = NA, align = "center") else ci_hi
    ) %>%
    ungroup()
}
plot_df <- plot_df %>% smooth_cols(k = 5)

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
  scale_colour_manual(values = c("2-year" = "#2166AC", "1-year" = "#D6604D")) +
  scale_fill_manual(values   = c("2-year" = "#2166AC", "1-year" = "#D6604D")) +
  scale_x_date(limits = x_lims, date_breaks = "2 years",
               date_labels = "%Y") +
  labs(
    title = "(1) Daily-aggregate rolling beta: window-size robustness",
    subtitle = "Poisson QMLE, NW(28). Lines are 5-week rolling-mean smoothed. Red dotted = MoU (2017-07-01).",
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
    subtitle = "AFR = Libya + Tunisia; EU = Italy + Malta. 5-week rolling-mean smoothed.",
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
    subtitle = "Libya, Tunisia, Italy, Malta estimated separately. 95% CI ribbons, 5-week rolling-mean smoothed.",
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
