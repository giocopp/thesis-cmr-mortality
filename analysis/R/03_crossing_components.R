# 03_crossing_components.R
# ========================
# Stacked area plots of crossing components (arrivals, interceptions, deaths)
# with fatality rate overlaid as black line.
#
# Weekly: UNHCR daily arrivals (from Oct 2015), interceptions from monthly.
# Monthly: IOM official arrivals (from 2014), interceptions, deaths.
#
# Input:
#   analysis/data/daily_panel.RDS
#
# Output:
#   output/figures/crossing_components.png

library(tidyverse)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")

cat("============================================================\n")
cat("CROSSING COMPONENTS\n")
cat("============================================================\n\n")

# ============================================================
# 1. Load daily panel
# ============================================================
cat("--- Loading daily panel ---\n")
daily <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS"))
cat(sprintf("  %d days\n", nrow(daily)))

# ============================================================
# 2. Weekly aggregation (arrivals NA before Oct 2015)
# ============================================================
weekly <- daily %>%
  group_by(iso_week) %>%
  summarise(
    date_start    = min(date),
    arrivals      = if (all(is.na(arrivals))) NA_real_ else sum(arrivals, na.rm = TRUE),
    deaths        = sum(deaths),
    interceptions = sum(intercept_per_day),
    n_days        = n(),
    .groups       = "drop"
  ) %>%
  filter(n_days >= 5, date_start <= as.Date("2025-10-01")) %>%
  mutate(
    total = if_else(!is.na(arrivals), arrivals + deaths + interceptions, NA_real_),
    fr    = if_else(!is.na(total) & total > 0, deaths / total, NA_real_)
  )

weekly_stack <- weekly %>%
  mutate(arrivals_plot = replace_na(arrivals, 0)) %>%
  select(date_start, arrivals_plot, interceptions, deaths) %>%
  pivot_longer(-date_start, names_to = "component", values_to = "count") %>%
  mutate(component = factor(component,
    levels = c("arrivals_plot", "interceptions", "deaths"),
    labels = c("Arrivals", "Interceptions (approx.)", "Deaths + missing")))

cat(sprintf("Weekly panel: %d weeks\n", nrow(weekly)))

# ============================================================
# 3. Monthly aggregation (official arrivals from 2014)
# ============================================================
monthly <- daily %>%
  group_by(ym) %>%
  summarise(deaths = sum(deaths), .groups = "drop") %>%
  left_join(
    daily %>%
      distinct(ym, official_arrivals, interceptions) %>%
      filter(!is.na(official_arrivals)),
    by = "ym"
  ) %>%
  filter(!is.na(official_arrivals), ym <= as.Date("2025-10-01")) %>%
  replace_na(list(interceptions = 0)) %>%
  mutate(
    total = official_arrivals + deaths + interceptions,
    fr    = if_else(total > 0, deaths / total, NA_real_)
  )

monthly_stack <- monthly %>%
  select(ym, official_arrivals, interceptions, deaths) %>%
  pivot_longer(-ym, names_to = "component", values_to = "count") %>%
  mutate(component = factor(component,
    levels = c("official_arrivals", "interceptions", "deaths"),
    labels = c("Arrivals", "Interceptions", "Deaths + missing")))

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

fill_colours <- c("Arrivals" = "steelblue",
                   "Interceptions (approx.)" = "#D4820E",
                   "Interceptions" = "#D4820E",
                   "Deaths + missing" = "#B2182B")

# --- Weekly ---
fr_scale_w <- max(weekly$total, na.rm = TRUE) / 0.30

p_weekly <- ggplot() +
  geom_area(data = weekly_stack, aes(x = date_start, y = count, fill = component),
            alpha = 0.7, position = "stack") +
  geom_line(data = weekly %>% filter(!is.na(fr)),
            aes(x = date_start, y = fr * fr_scale_w),
            colour = "black", linewidth = 0.3, alpha = 0.4) +
  geom_smooth(data = weekly %>% filter(!is.na(fr)),
              aes(x = date_start, y = fr * fr_scale_w),
              method = "loess", span = 0.15, se = FALSE,
              colour = "black", linewidth = 0.8) +
  mou_line + mou_label + shared_x +
  scale_fill_manual(values = fill_colours) +
  scale_y_continuous(
    name = "People per week",
    sec.axis = sec_axis(~ . / fr_scale_w * 100, name = "Fatality rate (%)")
  ) +
  labs(
    title = "Weekly crossing components and fatality rate",
    subtitle = "Stacked area = people counts | Black line = fatality rate (deaths/total).\nArrivals from Oct 2015. Interceptions approximated from monthly data.",
    x = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.title.y.right = element_text(colour = "black"))

# --- Monthly ---
fr_scale_m <- max(monthly$total, na.rm = TRUE) / 0.15

p_monthly <- ggplot() +
  geom_area(data = monthly_stack, aes(x = ym, y = count, fill = component),
            alpha = 0.7, position = "stack") +
  geom_line(data = monthly %>% filter(!is.na(fr)),
            aes(x = ym, y = fr * fr_scale_m),
            colour = "black", linewidth = 0.3, alpha = 0.4) +
  geom_smooth(data = monthly %>% filter(!is.na(fr)),
              aes(x = ym, y = fr * fr_scale_m),
              method = "loess", span = 0.2, se = FALSE,
              colour = "black", linewidth = 0.8) +
  mou_line + mou_label + shared_x +
  scale_fill_manual(values = fill_colours) +
  scale_y_continuous(
    name = "People per month",
    sec.axis = sec_axis(~ . / fr_scale_m * 100, name = "Fatality rate (%)")
  ) +
  labs(
    title = "Monthly crossing components and fatality rate (2014-2025)",
    subtitle = "Stacked area = people counts | Black line = fatality rate (deaths/total).\nArrivals from IOM official data. Interceptions = LCG + TCG.",
    x = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.title.y.right = element_text(colour = "black"))

p_out <- p_weekly / p_monthly
ggsave(file.path(BASE_DIR, "output", "figures", "crossing_components.png"),
       p_out, width = 12, height = 10, dpi = 200)
cat("Saved: output/figures/crossing_components.png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
