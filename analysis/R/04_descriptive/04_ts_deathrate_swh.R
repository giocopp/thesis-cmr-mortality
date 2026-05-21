# Extends Camarena et al. (2024, Fig. 2) to the full CMR analysis period
# (2014-01-01 – 2023-05-31).
#
# Red (Y1): LOESS smooth of weekly fatality rate (deaths / crossings),
#   where weekly sums exclude days with zero crossings (death-attribution
#   timing on no-crossing days is noise). LOESS is weighted by weekly
#   crossings so low-volume weeks (e.g. a single fatal incident with
#   few survivors) don't yank the smooth — same logic as a ratio-of-sums
#   rolling mean, but with smoother curvature and no edge truncation.
#   span = LOESS_SPAN (default 0.10, i.e. each fitted point uses the
#   nearest ~10% of weeks ≈ 24-week effective bandwidth — comparable to
#   a 24-week centered rolling mean but smoother). deaths = UNITED sea
#   deaths (CMR); crossings = deaths + frx_persons + lcg_tcg_pushbacks.
#   Raw weekly rates plotted as small points beneath.
#
# Blue (Y2): raw prior-week SWH (swh_prevweek), i.e. 7-day trailing mean.
#
# Dashed line: Italy-Libya MoU, 2017-02-02.
#
# In:  analysis/data/daily_panel_complete.RDS
#      data/processed/united_incidents.RDS
# Out: output/figures/13_ts_deathrate_swh.png

library(tidyverse)
library(lubridate)
library(scales)

Sys.setlocale("LC_TIME", "C")

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

# ── 1. UNITED daily sea deaths, CMR ──
united_daily <- readRDS(file.path(BASE_DIR, "data", "processed",
                                   "united_incidents.RDS")) |>
  filter(country_of_death %in% CMR_INCIDENT_COUNTRIES,
         (manner_of_death == "drowned" & !is.na(manner_of_death)) |
         (transport_means == "boat_ship_ferry" & !is.na(transport_means))) |>
  group_by(date = incident_date_clean) |>
  summarise(deaths = sum(n_deaths, na.rm = TRUE), .groups = "drop")

# ── 2. Daily panel (kept daily for SWH line) ──
daily <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS")) |>
  left_join(united_daily, by = "date") |>
  mutate(
    deaths      = replace_na(deaths, 0),
    crossings_u = frx_persons + lcg_tcg_pushbacks + deaths
  ) |>
  arrange(date) |>
  filter(!is.na(swh_prevweek))

# ── 3. Weekly aggregates (crossing-days only) ──
weekly <- daily |>
  filter(crossings_u > 0) |>
  mutate(week = floor_date(date, "week", week_start = 1)) |>
  group_by(week) |>
  summarise(
    deaths    = sum(deaths),
    crossings = sum(crossings_u),
    n_days    = n(),
    .groups   = "drop"
  ) |>
  arrange(week) |>
  mutate(fatality_rate = deaths / crossings)

# Weighted LOESS smooth (weights = crossings) so low-volume weeks
# don't dominate the curve.
loess_fit <- loess(fatality_rate ~ as.numeric(week),
                   data    = weekly,
                   weights = weekly$crossings,
                   span    = LOESS_SPAN)
weekly$loess_rate <- predict(loess_fit)

cat(sprintf("Daily obs (SWH line): %d (%s to %s)\n",
            nrow(daily), min(daily$date), max(daily$date)))
cat(sprintf("Daily obs dropped (crossings == 0): %d (%.1f%%)\n",
            sum(daily$crossings_u == 0),
            100 * mean(daily$crossings_u == 0)))
cat(sprintf("Weekly obs after filter: %d; UNITED deaths summed=%.0f\n",
            nrow(weekly), sum(weekly$deaths)))
cat(sprintf("LOESS-smoothed rate (span=%.2f): mean=%.4f, median=%.4f, max=%.3f\n",
            LOESS_SPAN,
            mean(weekly$loess_rate),
            median(weekly$loess_rate),
            max(weekly$loess_rate)))

# ── 4. Dual-axis scaling (use loess_rate range, ignore raw point spikes) ──
dr_max  <- max(weekly$loess_rate,   na.rm = TRUE)
swh_max <- max(daily$swh_prevweek,  na.rm = TRUE)
scl     <- dr_max / swh_max

# ── 5. Plot ──
p <- ggplot() +
  geom_line(data = daily,
            aes(x = date, y = swh_prevweek * scl, colour = "Wave Height"),
            linewidth = 0.35, alpha = 0.85) +
  geom_point(data = weekly,
             aes(x = week, y = fatality_rate, colour = "Death Rate"),
             size = 0.5, alpha = 0.25) +
  geom_line(data = weekly,
            aes(x = week, y = loess_rate, colour = "Death Rate"),
            linewidth = 0.7) +
  geom_vline(xintercept = MOU_SIGN_DATE,
             linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  annotate("text", x = MOU_SIGN_DATE, y = dr_max * 1.02,
           label = "Italy-Libya MoU (Feb 2017)",
           hjust = -0.05, vjust = 1, size = 3, colour = "grey30") +
  scale_colour_manual(values = c("Death Rate"  = "red",
                                 "Wave Height" = "blue")) +
  scale_x_date(breaks = seq(as.Date("2014-01-01"),
                            as.Date("2023-06-01"), by = "1 year"),
               date_labels = "%Y",
               expand = c(0.01, 0.01)) +
  scale_y_continuous(
    name     = sprintf("Death Rate (LOESS, span=%.2f)", LOESS_SPAN),
    labels   = label_percent(accuracy = 1),
    sec.axis = sec_axis(~ . / scl, name = "Wave Height, prior week (m)")
  ) +
  coord_cartesian(ylim = c(0, dr_max * 1.05)) +
  labs(x = "Date", colour = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position    = "top",
    legend.direction   = "horizontal",
    legend.key.width   = unit(1.2, "cm"),
    axis.title.y.left  = element_text(colour = "red"),
    axis.title.y.right = element_text(colour = "blue"),
    axis.text.y.left   = element_text(colour = "red"),
    axis.text.y.right  = element_text(colour = "blue"),
    panel.grid         = element_blank()
  )

ggsave(fig_path("04_descriptive", "04_ts_deathrate_swh.png"),
       plot = p, width = 11, height = 5.5, dpi = 300)
