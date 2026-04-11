# 02_fatality_rate_timeseries.R
# =============================
# Fatality rate time series: weekly and monthly.
# Compares: deaths/(deaths+arrivals) vs deaths/(deaths+arrivals+interceptions)
# Monthly also shows IOM official rate.
#
# Input:
#   analysis/data/daily_panel.RDS
#   data/processed/iom_med_crossings_monthly.RDS  (for IOM official rate)
#
# Output:
#   output/figures/fatality_rate_timeseries.png

Sys.setlocale("LC_TIME", "en_US.UTF-8")

library(tidyverse)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")

cat("============================================================\n")
cat("FATALITY RATE TIME SERIES\n")
cat("============================================================\n\n")

# ============================================================
# 1. Load daily panel
# ============================================================
cat("--- Loading daily panel ---\n")
daily <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS")) %>%
  filter(!is.na(arrivals))

cat(sprintf("  %d days with arrivals data\n", nrow(daily)))

# ============================================================
# 2. Weekly aggregation
# ============================================================
weekly <- daily %>%
  group_by(iso_week) %>%
  summarise(
    date_start           = min(date),
    deaths               = sum(deaths),
    arrivals             = sum(arrivals),
    interceptions_approx = sum(intercept_per_day),
    n_days               = n(),
    .groups              = "drop"
  ) %>%
  filter(n_days >= 5, date_start <= as.Date("2025-10-01")) %>%
  mutate(
    fr    = if_else((deaths + arrivals) > 0,
                     deaths / (deaths + arrivals), NA_real_),
    fr_ic = if_else((deaths + arrivals + interceptions_approx) > 0,
                     deaths / (deaths + arrivals + interceptions_approx), NA_real_)
  )

cat(sprintf("Weekly panel: %d weeks\n", nrow(weekly)))

# ============================================================
# 3. Monthly aggregation (with IOM official rate)
# ============================================================
official_rate <- readRDS(file.path(BASE_DIR, "data", "processed",
                                     "iom_med_crossings_monthly.RDS")) %>%
  transmute(ym = as.Date(date), official_rate = cmr_rate_of_death)

monthly <- daily %>%
  group_by(ym) %>%
  summarise(deaths = sum(deaths), arrivals = sum(arrivals),
            interceptions = sum(intercept_per_day), .groups = "drop") %>%
  left_join(official_rate, by = "ym") %>%
  filter(!is.na(official_rate), ym <= as.Date("2025-10-01")) %>%
  mutate(
    fr    = if_else((deaths + arrivals) > 0,
                     deaths / (deaths + arrivals), NA_real_),
    fr_ic = if_else((deaths + arrivals + interceptions) > 0,
                     deaths / (deaths + arrivals + interceptions), NA_real_)
  )

cat(sprintf("Monthly panel: %d months\n\n", nrow(monthly)))

# ============================================================
# 4. Plots
# ============================================================
cat("--- Generating plots ---\n")

x_limits <- c(as.Date("2014-01-01"), as.Date("2025-11-01"))
shared_x <- scale_x_date(limits = x_limits, date_breaks = "1 year", date_labels = "%Y")
mou_line  <- geom_vline(xintercept = MOU_DATE, linetype = "dashed",
                          colour = "grey40", linewidth = 0.4)
mou_label <- annotate("text", x = MOU_DATE, y = Inf, label = "MoU",
                        vjust = 1.5, hjust = -0.1, size = 3, colour = "grey40")

# --- Weekly ---
weekly_long <- weekly %>%
  select(date_start, fr, fr_ic) %>%
  pivot_longer(-date_start, names_to = "source", values_to = "fr_val") %>%
  mutate(source = recode(source,
    fr    = "Deaths / (deaths + arrivals)",
    fr_ic = "Deaths / (deaths + arrivals + interceptions*)"
  ))

p_weekly <- ggplot(weekly_long %>% filter(!is.na(fr_val)),
                   aes(x = date_start, y = fr_val * 100, colour = source)) +
  geom_line(alpha = 0.4, linewidth = 0.3) +
  geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 0.8) +
  mou_line + mou_label + shared_x +
  scale_colour_manual(values = c(
    "Deaths / (deaths + arrivals)"                    = "steelblue",
    "Deaths / (deaths + arrivals + interceptions*)"   = "#2CA02C"
  )) +
  labs(title = "Weekly fatality rate",
       subtitle = "All CMR deaths. *Interceptions approximated from monthly data.",
       y = "Fatality rate (%)", x = NULL, colour = NULL) +
  coord_cartesian(ylim = c(0, 30)) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

# --- Monthly ---
monthly_long <- monthly %>%
  select(ym, fr, fr_ic, official_rate) %>%
  pivot_longer(-ym, names_to = "source", values_to = "fr_val") %>%
  mutate(source = recode(source,
    fr            = "Our method",
    fr_ic         = "Our method + interceptions",
    official_rate = "IOM official"
  )) %>%
  filter(!is.na(fr_val))

p_monthly <- ggplot(monthly_long, aes(x = ym, y = fr_val * 100, colour = source)) +
  geom_line(alpha = 0.4, linewidth = 0.3) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, linewidth = 0.8) +
  mou_line + mou_label + shared_x +
  scale_colour_manual(values = c(
    "Our method"                 = "steelblue",
    "Our method + interceptions" = "#2CA02C",
    "IOM official"               = "#D4820E"
  )) +
  labs(title = "Monthly fatality rate",
       subtitle = "Blue = deaths/(deaths+arrivals) | Green = with interceptions | Orange = IOM official",
       y = "Fatality rate (%)", x = NULL, colour = NULL) +
  coord_cartesian(ylim = c(0, 15)) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

p_out <- p_weekly / p_monthly
ggsave(file.path(BASE_DIR, "output", "figures", "fatality_rate_timeseries.png"),
       p_out, width = 12, height = 10, dpi = 200)
cat("Saved: output/figures/fatality_rate_timeseries.png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
