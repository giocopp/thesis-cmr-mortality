# Paper figure: 2×2 panel of crossings, deaths/fatality, detailed
# composition, and aggregate shares for the Central Mediterranean Route.

library(tidyverse)
library(lubridate)
library(scales)
library(patchwork)
library(cowplot)
library(grid)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))
MOU_DATE <- MOU_SIGN_DATE

# ── 1. Load data ────────────────────────────────────────────────────────────
panel <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel_complete.RDS"))
stopifnot(all(c("date", "frx_persons", "frx_incidents", "n_dead_missing",
                "crossing_attempts", "post_mou",
                "lcg_tcg_pushbacks", "arrivals") %in% names(panel)))

frx <- readRDS(file.path(BASE_DIR, "data", "processed",
                          "frontex_incidents_coords.RDS")) |>
  filter(country_of_departure %in% CMR_DEPARTURES)

PANEL_END <- max(panel$date)
PLOT_XLIM <- c(as.Date("2014-01-01"), PANEL_END)

iom_monthly <- readRDS(file.path(BASE_DIR, "data", "processed",
                                  "iom_med_crossings_monthly.RDS")) |>
  transmute(
    ym = as.Date(date),
    arrivals = as.numeric(sea_arrivals_in_italy),
    lcg = as.numeric(interceptions_by_libyan_coast_guard),
    tcg = as.numeric(interceptions_by_tunisian_coast_guard)
  ) |>
  filter(ym <= floor_date(PANEL_END, "month"))

united_raw <- readRDS(file.path(BASE_DIR, "data", "processed", "united_incidents.RDS"))
UNITED_CMR_COUNTRIES <- c(CMR_INCIDENT_COUNTRIES, "Mediterranean")
united_daily <- united_raw |>
  filter(country_of_death %in% UNITED_CMR_COUNTRIES,
         incident_year >= 2014L,
         incident_date_clean <= PANEL_END,
         (manner_of_death == "drowned" & !is.na(manner_of_death)) |
         (transport_means == "boat_ship_ferry" & !is.na(transport_means))) |>
  group_by(date = incident_date_clean) |>
  summarise(n_dead_united = sum(n_deaths, na.rm = TRUE), .groups = "drop")

# ── 2. Monthly aggregates ───────────────────────────────────────────────────
monthly_united <- panel |>
  select(date, frx_persons, frx_incidents, lcg_tcg_pushbacks, arrivals) |>
  left_join(united_daily, by = "date") |>
  replace_na(list(n_dead_united = 0)) |>
  mutate(ym = floor_date(date, "month")) |>
  group_by(ym) |>
  summarise(
    frx_persons       = sum(frx_persons),
    frx_incidents     = sum(frx_incidents),
    n_dead_missing    = sum(n_dead_united),
    lcg_tcg_pushbacks = sum(lcg_tcg_pushbacks),
    .groups = "drop"
  ) |>
  mutate(crossing_attempts = frx_persons + lcg_tcg_pushbacks + n_dead_missing) |>
  left_join(iom_monthly, by = "ym") |>
  replace_na(list(lcg = 0, tcg = 0))

# ── 3. Event-type classification on Frontex incidents ───────────────────────
frx <- frx |>
  mutate(
    sar_true = !is.na(sar_ops) & sar_ops,
    event_type = case_when(
      sar_true & interceptor_type == "EU_ops"                  ~ "SAR: EU operations",
      sar_true & interceptor_type == "Ita_ops"                 ~ "SAR: Italian authorities",
      sar_true & interceptor_type == "NGO"                     ~ "SAR: NGO operations",
      sar_true & interceptor_type == "EU_Coast_Guard"          ~ "SAR: EU members CG",
      sar_true                                                 ~ "SAR: Mare Nostrum and Others",
      interceptor_type == "EU_Coast_Guard"                     ~ "Non SAR: EU CG operations",
      interceptor_type == "Land_patrol"                        ~ "Not SAR: Land patrol",
      interceptor_type %in% c("NA", "No_intercept")            ~ "Not SAR: No intercep. (detection only)",
      TRUE                                                     ~ "Not SAR: Mare Nostrum and Others"
    )
  ) |>
  select(-sar_true)

etype_order <- c(
  "SAR: Italian authorities",
  "SAR: NGO operations",
  "SAR: EU members CG",
  "SAR: EU operations",
  "SAR: Mare Nostrum and Others",
  " ",  # phantom legend spacer
  "Not SAR: No intercep. (detection only)",
  "Non SAR: EU CG operations",
  "Not SAR: Land patrol",
  "Not SAR: Mare Nostrum and Others",
  "Libyan CG operations",
  "Tunisian CG operations"
)

etype_colours <- c(
  "SAR: Italian authorities"               = "#08519C",
  "SAR: NGO operations"                    = "#3182BD",
  "SAR: EU members CG"                     = "#6BAED6",
  "SAR: EU operations"                     = "#9ECAE1",
  "SAR: Mare Nostrum and Others"           = "#C6DBEF",
  " "                                      = "#FFFFFF00",
  "Non SAR: EU CG operations"              = "#E6550D",
  "Not SAR: Land patrol"                   = "#FD8D3C",
  "Not SAR: No intercep. (detection only)" = "#A63603",
  "Not SAR: Mare Nostrum and Others"       = "#FDBE85",
  "Libyan CG operations"                   = "#252525",
  "Tunisian CG operations"                 = "#969696"
)

frx_etype_p <- frx |>
  mutate(ym = floor_date(date, "month"), etype = event_type) |>
  group_by(ym, etype) |>
  summarise(persons = sum(num_persons, na.rm = TRUE), .groups = "drop")

ic_long <- iom_monthly |>
  filter(!is.na(lcg) | !is.na(tcg)) |>
  pivot_longer(cols = c(lcg, tcg), names_to = "type", values_to = "persons") |>
  filter(!is.na(persons)) |>
  mutate(etype = ifelse(type == "lcg", "Libyan CG operations", "Tunisian CG operations")) |>
  select(ym, etype, persons)

persons_all <- bind_rows(frx_etype_p, ic_long) |>
  filter(ym <= PANEL_END) |>
  mutate(etype = factor(etype, levels = etype_order))

# ── 4. Panel-builder helpers ────────────────────────────────────────────────
comp_order <- c("frontex_persons", "lcg_persons", "tcg_persons", "deaths_missing")
comp_colours <- c(
  "frontex_persons"    = "#1F78B4",
  "lcg_persons"        = "#333333",
  "tcg_persons"        = "#AAAAAA",
  "deaths_missing"     = "#D32F2F"
)

build_crossing_long <- function(m, death_label) {
  comp_labels <- c(
    "frontex_persons"    = "Intercepted persons",
    "lcg_persons"        = "Libyan CG operations",
    "tcg_persons"        = "Tunisian CG operations",
    "deaths_missing"     = death_label
  )
  m |>
    filter(ym <= PANEL_END) |>
    mutate(
      frontex_persons    = replace_na(frx_persons, 0),
      lcg_persons        = replace_na(lcg, 0),
      tcg_persons        = replace_na(tcg, 0),
      deaths_missing     = n_dead_missing
    ) |>
    select(ym, frontex_persons, lcg_persons, tcg_persons, deaths_missing) |>
    pivot_longer(-ym, names_to = "component", values_to = "persons") |>
    mutate(component = factor(component, levels = comp_order,
                              labels = comp_labels[comp_order]))
}

build_cross_panel <- function(m, death_label) {
  cross_long <- build_crossing_long(m, death_label)
  comp_labels_vec <- c(
    "frontex_persons"    = "Intercepted persons",
    "lcg_persons"        = "Libyan CG operations",
    "tcg_persons"        = "Tunisian CG operations",
    "deaths_missing"     = death_label
  )
  max_bar <- cross_long |>
    group_by(ym) |>
    summarise(total = sum(persons, na.rm = TRUE), .groups = "drop") |>
    pull(total) |>
    max(na.rm = TRUE)
  step <- 10^floor(log10(max_bar / 4))
  cross_step <- ceiling((max_bar / 4) / step) * step
  cross_max <- cross_step * 4
  cross_breaks <- seq(cross_step, cross_max, by = cross_step)

  ggplot(cross_long, aes(x = ym, y = persons, fill = component)) +
    geom_col(position = "stack", width = 25) +
    scale_fill_manual(
      values = setNames(comp_colours, comp_labels_vec[names(comp_colours)]),
      breaks = c(
        "Intercepted persons",
        death_label,
        "Libyan CG operations",
        "Tunisian CG operations"
      ),
      name = NULL
    ) +
    geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red", linewidth = 0.5) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(breaks = cross_breaks) +
    coord_cartesian(xlim = PLOT_XLIM, ylim = c(0, cross_max), clip = "off",
                    expand = FALSE) +
    labs(
      title = "(a) Number of persons attempting the crossing, by crossing outcome",
      subtitle = "Red dashed line marks the signature of the 2017 Italy-Libya MoU",
      y = "Number of persons", x = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "red3"),
      plot.margin = margin(t = 16, r = 5, b = 5, l = 5)
    )
}

DEATHS_Y2_MAX <- 20
build_deaths_panel <- function(m, death_label, shared_sf = NULL) {
  dr <- m |>
    filter(ym <= PANEL_END) |>
    mutate(death_rate = ifelse(crossing_attempts > 0,
                                n_dead_missing / crossing_attempts * 100, NA_real_))
  sf <- if (!is.null(shared_sf)) shared_sf else max(dr$n_dead_missing, na.rm = TRUE) / DEATHS_Y2_MAX

  primary_max     <- sf * DEATHS_Y2_MAX
  primary_step    <- primary_max / 4
  primary_breaks  <- seq(primary_step, primary_max, by = primary_step)
  secondary_step  <- DEATHS_Y2_MAX / 4
  secondary_breaks <- seq(secondary_step, DEATHS_Y2_MAX, by = secondary_step)

  ggplot(dr, aes(x = ym)) +
    geom_col(aes(y = n_dead_missing, fill = death_label), width = 25) +
    geom_line(aes(y = death_rate * sf, colour = "Fatality rate (%)"),
              linewidth = 0.7) +
    scale_fill_manual(values = setNames("#D32F2F", death_label), name = NULL) +
    scale_colour_manual(values = c("Fatality rate (%)" = "#333333"), name = NULL) +
    scale_y_continuous(
      name = "Recorded Dead and Missing Migrants",
      breaks = primary_breaks,
      sec.axis = sec_axis(~ . / sf, name = "Fatality rate (%)",
                          breaks = secondary_breaks)
    ) +
    geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red", linewidth = 0.5) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y",
                 expand = expansion(mult = c(0.005, 0.005))) +
    coord_cartesian(xlim = PLOT_XLIM,
                    ylim = if (!is.null(shared_sf)) c(0, shared_sf * DEATHS_Y2_MAX) else NULL,
                    clip = "off", expand = FALSE) +
    labs(title = "(b) Number of persons dead and missing while attempting the crossing", x = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      plot.margin = margin(t = 16, r = 5, b = 5, l = 5)
    )
}

# Shared deaths y-scale so panels (a/b) align with composition panels (c/d).
deaths_max_raw <- max(monthly_united$n_dead_missing[monthly_united$ym <= PANEL_END],
                      na.rm = TRUE)
deaths_step_base <- 10^floor(log10(deaths_max_raw / 4))
deaths_y_step    <- ceiling((deaths_max_raw / 4) / deaths_step_base) * deaths_step_base
deaths_y_max     <- deaths_y_step * 4
deaths_sf        <- deaths_y_max / DEATHS_Y2_MAX

# ── 5. Composition data for panels (c) and (d) ──────────────────────────────
LOESS_SPAN <- 0.15
all_months <- seq(min(persons_all$ym, na.rm = TRUE),
                  floor_date(PANEL_END, "month"),
                  by = "month")

persons_complete <- expand_grid(
  ym = all_months,
  etype = factor(etype_order, levels = etype_order)
) |>
  left_join(persons_all, by = c("ym", "etype")) |>
  replace_na(list(persons = 0)) |>
  group_by(ym) |>
  mutate(
    total = sum(persons),
    pct = ifelse(total > 0, persons / total * 100, 0)
  ) |>
  ungroup()

persons_smoothed <- persons_complete |>
  group_by(etype) |>
  mutate(
    ym_num = as.numeric(ym),
    pct_smooth = predict(loess(pct ~ ym_num, span = LOESS_SPAN))
  ) |>
  ungroup() |>
  mutate(pct_smooth = pmax(pct_smooth, 0)) |>
  group_by(ym) |>
  mutate(pct_smooth = pct_smooth / sum(pct_smooth) * 100) |>
  ungroup()

mou_vline <- geom_vline(
  xintercept = MOU_DATE,
  linetype = "dashed",
  colour = "red",
  linewidth = 0.5
)

p_etype_pct <- ggplot(persons_smoothed, aes(x = ym, y = pct_smooth, fill = etype)) +
  geom_area(position = "stack") +
  scale_fill_manual(values = etype_colours, name = NULL) +
  scale_y_continuous(name = "% of persons intercepted",
                     breaks = c(25, 50, 75, 100)) +
  mou_vline +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  coord_cartesian(xlim = PLOT_XLIM, ylim = c(0, 100), clip = "off",
                  expand = FALSE) +
  labs(
    title = "(b) Detailed composition of interceptions",
    y = "% of persons intercepted",
    x = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(t = -5),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(t = 16, r = 5, b = 5, l = 5)
  ) +
  guides(fill = guide_legend(nrow = 3))

share_colours <- c(
  "SAR operations recorded by Frontex"     = "#1F78B4",
  "Non SAR operations recorded by Frontex" = "#F16913",
  "Libyan CG operations"                   = "#252525",
  "Tunisian CG operations"                 = "#969696"
)

share_lines <- persons_complete |>
  filter(etype != " ") |>
  mutate(group = case_when(
    grepl("^SAR:", etype)               ~ "SAR operations recorded by Frontex",
    grepl("^(Not SAR|Non SAR):", etype) ~ "Non SAR operations recorded by Frontex",
    TRUE                                ~ as.character(etype)
  )) |>
  group_by(ym, group) |>
  summarise(persons = sum(persons), .groups = "drop") |>
  group_by(ym) |>
  mutate(
    total = sum(persons),
    share = ifelse(total > 0, persons / total * 100, 0)
  ) |>
  ungroup() |>
  select(-total) |>
  mutate(group = factor(group, levels = names(share_colours)))

p_share_lines <- ggplot(share_lines, aes(x = ym, y = share, colour = group)) +
  geom_smooth(method = "loess", span = LOESS_SPAN, se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = share_colours, name = NULL) +
  mou_vline +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(breaks = c(25, 50, 75, 100)) +
  coord_cartesian(xlim = PLOT_XLIM, ylim = c(0, 100), clip = "off",
                  expand = FALSE) +
  labs(
    title = "(c) Percentage of persons intercepted",
    y = "% of persons intercepted",
    x = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.margin = margin(t = -5),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(t = 16, r = 5, b = 5, l = 5)
  )

# ── 6. Assemble paper figure (2 × 2) ────────────────────────────────────────
legend_bottom_theme <- theme(
  legend.position = "bottom",
  legend.background = element_blank(),
  legend.box.background = element_rect(fill = "grey97", colour = "grey80", linewidth = 0.5),
  legend.box.margin = margin(3, 5, 3, 5),
  legend.box = "vertical",
  legend.key = element_blank(),
  legend.margin = margin(2, 0, 2, 0),
  legend.text = element_text(size = 14),
  legend.key.size = unit(0.6, "cm"),
  legend.title = element_blank(),
  legend.spacing = unit(0, "cm"),
  legend.spacing.x = unit(0.15, "cm"),
  legend.spacing.y = unit(0.05, "cm")
)

p_cross_united  <- build_cross_panel(monthly_united, "Recorded Dead and Missing Migrants")
p_deaths_united <- build_deaths_panel(monthly_united,
                                      "Recorded Dead and Missing Migrants",
                                      shared_sf = deaths_sf)

panel_2x2_theme <- theme(
  plot.title = element_text(face = "bold", size = 22, lineheight = 0.98),
  plot.subtitle = element_text(size = 16, colour = "red3"),
  axis.title = element_text(size = 18),
  axis.text = element_text(size = 16),
  plot.margin = margin(t = 9, r = 4, b = 2, l = 4)
)

p_cross_united_2x2 <- p_cross_united +
  labs(title = "(a) Totals: persons attempting crossing, by outcome") +
  panel_2x2_theme

p_deaths_united_2x2 <- p_deaths_united +
  labs(title = "(b) Totals: recorded dead/missing and fatality rate") +
  panel_2x2_theme +
  theme(axis.title.y.right = element_blank())

p_etype_pct_2x2 <- p_etype_pct +
  labs(
    title = "(c) Shares: detailed composition of interceptions",
    subtitle = "Red dashed line marks the signature of the 2017 Italy-Libya MoU"
  ) +
  panel_2x2_theme

p_share_lines_2x2 <- p_share_lines +
  labs(title = "(d) Shares: percentage of persons intercepted") +
  panel_2x2_theme

p_cross_united_display <- p_cross_united_2x2 +
  legend_bottom_theme +
  guides(fill = guide_legend(nrow = 2, byrow = FALSE, title = NULL))

p_deaths_united_display <- p_deaths_united_2x2 +
  legend_bottom_theme +
  guides(
    fill   = guide_legend(nrow = 1, order = 1, title = NULL),
    colour = guide_legend(nrow = 1, order = 2, title = NULL)
  )

p_etype_pct_display <- p_etype_pct_2x2 +
  legend_bottom_theme +
  guides(fill = guide_legend(ncol = 2, byrow = FALSE, title = NULL))

p_share_lines_display <- p_share_lines_2x2 +
  legend_bottom_theme +
  guides(colour = guide_legend(nrow = 2, byrow = FALSE, title = NULL))

# Build rows separately so legend-height padding stays within each row: the tall
# (c) legend would otherwise pad (b)/(d) and leave a white band at the figure foot.
panel_top <- cowplot::plot_grid(
  p_cross_united_display, p_etype_pct_display,
  ncol = 2, align = "h", axis = "tb"
)
panel_bottom <- cowplot::plot_grid(
  p_deaths_united_display, p_share_lines_display,
  ncol = 2, align = "h", axis = "tb"
)
panel_event_type <- cowplot::plot_grid(
  panel_top, panel_bottom,
  ncol = 1, align = "v", axis = "lr"
)

# Small even inset so the content does not butt up against the frame.
frame_pad <- 0.012
panel_event_type_framed <- cowplot::ggdraw() +
  cowplot::draw_plot(panel_event_type,
                     x = frame_pad, y = frame_pad,
                     width = 1 - 2 * frame_pad, height = 1 - 2 * frame_pad) +
  cowplot::draw_grob(grid::rectGrob(
    gp = grid::gpar(col = "black", fill = NA, lwd = 2)
  ))

ggsave(
  fig_path("04_descriptive", "fig-panel-event-type.png"),
  panel_event_type_framed,
  width = 18.5, height = 16, dpi = 300
)
