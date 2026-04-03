# 04_swh_vs_fatality_rate.R
# =========================
# Dual-axis time series: LOWESS-smoothed death rate + raw SWH.
# Camarena et al. (2020) Fig 2 style, using our data.
#
# Death rate = deaths / crossings (deaths + arrivals + interceptions)
# LOWESS: bw = 0.025, iter = 0 (matches Stata twoway lowess default)
# SWH: raw prior-week average (7-day rolling mean of lag-1)
#
# Input:
#   analysis/data/daily_panel.RDS
#
# Output:
#   output/figures/swh_vs_fatality_rate.png       (2016-2021)
#   output/figures/swh_vs_fatality_rate_2016_2018.png (zoomed)

Sys.setlocale("LC_TIME", "en_US.UTF-8")

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")

cat("============================================================\n")
cat("SWH vs FATALITY RATE\n")
cat("============================================================\n\n")

# ============================================================
# 1. Load daily panel
# ============================================================
cat("--- Loading daily panel ---\n")
daily <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS"))
cat(sprintf("  %d days\n", nrow(daily)))

# ============================================================
# 2. Camarena-style plot helper
# ============================================================
# Death rate uses crossings (no interceptions) to match Camarena's definition.
# SWH is raw prior-week average (not smoothed).
# Death rate is LOWESS-smoothed (bw=0.025, iter=0).

camarena_plot <- function(df, date_lim, x_breaks, x_fmt, fr_lim, wh_lim,
                          legend_pos = c(0.15, 0.92)) {

  scale_f  <- diff(fr_lim) / diff(wh_lim)
  offset_f <- fr_lim[1] - wh_lim[1] * scale_f

  plot_df <- df %>%
    filter(date >= date_lim[1], date <= date_lim[2]) %>%
    mutate(deathrate = if_else(crossings > 0, deaths / crossings, NA_real_))

  dr_df <- plot_df %>% filter(!is.na(deathrate))
  lo <- lowess(as.numeric(dr_df$date), dr_df$deathrate, f = 0.025, iter = 0)
  dr_df$dr_smooth <- lo$y

  ggplot() +
    geom_line(data = plot_df %>% filter(!is.na(swh_prevweek)),
              aes(x = date, y = swh_prevweek * scale_f + offset_f),
              colour = "blue", linewidth = 0.5) +
    geom_line(data = dr_df, aes(x = date, y = dr_smooth),
              colour = "red", linewidth = 0.5) +
    # Legend entries
    geom_line(data = dr_df[1:2, ], aes(x = date, y = dr_smooth, colour = "Death Rate"),
              linewidth = 0.5, show.legend = TRUE) +
    geom_line(data = plot_df[1:2, ] %>% filter(!is.na(swh_prevweek)),
              aes(x = date, y = swh_prevweek * scale_f + offset_f, colour = "Wave Height"),
              linewidth = 0.5, show.legend = TRUE) +
    scale_colour_manual(values = c("Death Rate" = "red", "Wave Height" = "blue")) +
    scale_x_date(breaks = x_breaks, date_labels = x_fmt,
                 limits = date_lim) +
    scale_y_continuous(
      name = "Death Rate",
      limits = fr_lim,
      sec.axis = sec_axis(~ (. - offset_f) / scale_f,
                          name = "Wave Height, prior week")
    ) +
    labs(x = "Date", colour = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.ticks       = element_line(colour = "black", linewidth = 0.3),
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      legend.position  = legend_pos,
      legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
      legend.key.size  = unit(0.8, "lines"),
      legend.text      = element_text(size = 8)
    )
}

# ============================================================
# 3. Full-period plot (2016-2021)
# ============================================================
cat("--- Generating plots ---\n")

daily_plot <- daily %>% filter(!is.na(arrivals))

p_full <- camarena_plot(
  daily_plot,
  date_lim = as.Date(c("2016-01-01", "2021-12-31")),
  x_breaks = seq(as.Date("2016-01-01"), as.Date("2022-01-01"), by = "1 year"),
  x_fmt    = "%Y",
  fr_lim   = c(0, 0.10),
  wh_lim   = c(0.3, 2.5),
  legend_pos = c(0.10, 0.92)
) +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x = MOU_DATE, y = 0.95, label = "MoU",
           hjust = -0.1, size = 3, colour = "grey40")

ggsave(file.path(BASE_DIR, "output", "figures", "swh_vs_fatality_rate.png"),
       p_full, width = 12, height = 5, dpi = 200)
cat("Saved: output/figures/swh_vs_fatality_rate.png\n")

# ============================================================
# 4. Zoomed plot (Jan 2016 - Mar 2018)
# ============================================================
x_breaks_zoom <- as.Date(seq(20461, 21272, by = 60), origin = "1960-01-01")

p_zoom <- camarena_plot(
  daily_plot,
  date_lim = as.Date(c("2016-01-08", "2018-03-29")),
  x_breaks = x_breaks_zoom,
  x_fmt    = "%b-%Y",
  fr_lim   = c(0, 0.15),
  wh_lim   = c(0.3, 2.5)
)

ggsave(file.path(BASE_DIR, "output", "figures", "swh_vs_fatality_rate_2016_2018.png"),
       p_zoom, width = 10, height = 5, dpi = 200)
cat("Saved: output/figures/swh_vs_fatality_rate_2016_2018.png\n")

# ============================================================
# 5. Comparison: three panels stacked (2016-2018)
# ============================================================
# Each panel = dual-axis Camarena-style (LOWESS death rate + raw SWH).
# Panel A: Camarena et al. data
# Panel B: Our data, deaths / (deaths + arrivals)
# Panel C: Our data, deaths / (deaths + arrivals + interceptions)

library(patchwork)

cam <- haven::read_dta(file.path(BASE_DIR, "replication", "camarena-et-al",
                                  "MiM_Replication", "DATA", "time_series.dta")) %>%
  filter(sample == 1) %>%
  mutate(date = as.Date(edate, origin = "1960-01-01"))

date_lim <- as.Date(c("2016-01-08", "2018-03-29"))
x_brk    <- as.Date(seq(20461, 21272, by = 60), origin = "1960-01-01")
wh_lim   <- c(0.3, 2.5)

# Shared theme for all three panels
shared_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(colour = "black", linewidth = 0.3),
    legend.position  = "none"
  )

# Helper: build one Camarena-style panel
make_panel <- function(dr_df, swh_df, fr_lim, title_text) {
  scale_f  <- diff(fr_lim) / diff(wh_lim)
  offset_f <- fr_lim[1] - wh_lim[1] * scale_f

  lo <- lowess(as.numeric(dr_df$date), dr_df$dr, f = 0.025, iter = 0)
  dr_df$smooth <- lo$y

  ggplot() +
    geom_line(data = swh_df, aes(x = date, y = swh * scale_f + offset_f),
              colour = "blue", linewidth = 0.4) +
    geom_line(data = dr_df, aes(x = date, y = smooth),
              colour = "red", linewidth = 0.5) +
    scale_x_date(breaks = x_brk, date_labels = "%b-%Y", limits = date_lim) +
    scale_y_continuous(
      name = "Death Rate",
      limits = fr_lim,
      sec.axis = sec_axis(~ (. - offset_f) / scale_f,
                          name = "Wave Height, prior week")
    ) +
    labs(title = title_text, x = NULL) +
    shared_theme
}

# --- Panel A: Camarena data ---
cam_dr <- cam %>% filter(!is.na(deathrate)) %>%
  transmute(date, dr = deathrate)
cam_swh <- cam %>% transmute(date, swh = wave_height_prevweek)

p_a <- make_panel(cam_dr, cam_swh, fr_lim = c(0, 1),
                  "A. Camarena et al. data")

# --- Panel B: Our data, no interceptions ---
own <- daily_plot %>%
  filter(date >= date_lim[1], date <= date_lim[2])

own_dr_no_ic <- own %>%
  mutate(dr = if_else(crossings_no_ic > 0, deaths / crossings_no_ic, NA_real_)) %>%
  filter(!is.na(dr)) %>% select(date, dr)
own_swh <- own %>% filter(!is.na(swh_prevweek)) %>%
  transmute(date, swh = swh_prevweek)

p_b <- make_panel(own_dr_no_ic, own_swh, fr_lim = c(0, 1),
                  "B. Our data — deaths / (deaths + arrivals)")

# --- Panel C: Our data, with interceptions ---
own_dr_ic <- own %>%
  mutate(dr = if_else(crossings > 0, deaths / crossings, NA_real_)) %>%
  filter(!is.na(dr)) %>% select(date, dr)

p_c <- make_panel(own_dr_ic, own_swh, fr_lim = c(0, 0.15),
                  "C. Our data — deaths / (deaths + arrivals + interceptions)")

p_comp <- p_a / p_b / p_c
ggsave(file.path(BASE_DIR, "output", "figures", "death_rate_comparison.png"),
       p_comp, width = 10, height = 12, dpi = 200)
cat("Saved: output/figures/death_rate_comparison.png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
