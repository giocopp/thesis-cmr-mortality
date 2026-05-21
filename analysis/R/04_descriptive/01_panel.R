# 10_descriptives_panel.R
# ===========================
# Comprehensive descriptive statistics for the Central Mediterranean Route.
# Monthly stacked bar panels showing crossings, deaths, composition over time.
# Annual summary table with key numbers.
#
# Uses the canonical IOM filter (incident + split incident, all causes,
# CMR countries) — see 01_build_daily_panel.R for the rationale.
#
# Input:
#   data/processed/frontex_incidents_coords.RDS
#   data/processed/iom_mmp_incidents.RDS
#   data/processed/iom_med_crossings_monthly.RDS
#   analysis/data/daily_panel_complete.RDS
#
# Output:
#   output/figures/desc_panel_crossings_iom.png     (crossings + deaths, IOM)
#   output/figures/desc_panel_crossings_united.png  (crossings + deaths, UNITED)
#   output/figures/desc_panel_event_type.png        (events + persons + composition)
#   output/figures/desc_boat_type.png               (line plot)
#   output/tables/desc_annual_summary.html
#   output/tables/desc_annual_summary.tex

library(tidyverse)
library(lubridate)
library(gt)
library(scales)
library(patchwork)
library(cowplot)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))
MOU_DATE <- MOU_SIGN_DATE  # plot annotation uses signing date

cat("============================================================\n")
cat("DESCRIPTIVE STATISTICS — CENTRAL MEDITERRANEAN ROUTE\n")
cat("============================================================\n\n")

# ── 1. Load & validate ──────────────────────────────────
cat("--- 1. Loading data ---\n")

# 1a. Daily panel
panel <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel_complete.RDS"))
stopifnot(all(c("date", "frx_persons", "frx_incidents", "n_dead_missing",
                "crossing_attempts", "post_mou",
                "lcg_tcg_pushbacks", "arrivals") %in% names(panel)))
cat(sprintf("  Daily panel: %d days (%s to %s)\n",
    nrow(panel), min(panel$date), max(panel$date)))

# 1b. Frontex incidents (CMR departures)
frx <- readRDS(file.path(BASE_DIR, "data", "processed", "frontex_incidents_coords.RDS")) |>
  filter(country_of_departure %in% CMR_DEPARTURES)
stopifnot(all(c("date", "boat_category", "interceptor_type", "detector_type",
                "detected_by", "intercepted_by", "num_persons", "sar_ops")
              %in% names(frx)))
FRX_END <- max(frx$date)
cat(sprintf("  Frontex CMR: %d incidents (%s to %s)\n",
    nrow(frx), min(frx$date), FRX_END))

# 1c. IOM MMP incidents — broad descriptive filter via iom_incidents().
iom_cmr <- iom_incidents(
    incident_types = c("incident", "split incident"),
    spatial        = "all_cmr",
    causes         = "all"
  ) |>
  filter(date <= max(panel$date))
cat(sprintf("  IOM CMR: %d incidents\n", nrow(iom_cmr)))

# Panel end (canonical end date for descriptives) — read from the rebuilt
# daily panel, which is bounded by interceptions_daily_disagg.RDS coverage.
PANEL_END <- max(panel$date)

# 1d. IOM monthly crossings (cap at the panel end date)
iom_monthly <- readRDS(file.path(BASE_DIR, "data", "processed",
                                  "iom_med_crossings_monthly.RDS")) |>
  transmute(
    ym = as.Date(date),
    arrivals = as.numeric(sea_arrivals_in_italy),
    lcg = as.numeric(interceptions_by_libyan_coast_guard),
    tcg = as.numeric(interceptions_by_tunisian_coast_guard)
  ) |>
  filter(ym <= floor_date(PANEL_END, "month"))
cat(sprintf("  IOM monthly: %d months\n", nrow(iom_monthly)))


# 1e. UNITED deaths — sea deaths only, restricted to the same 5 CMR countries
#     as IOM plus "Mediterranean" (open-sea), aggregated to daily.
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
cat(sprintf("  UNITED CMR sea deaths: %d days, %.0f total deaths\n",
            nrow(united_daily), sum(united_daily$n_dead_united)))

PLOT_XLIM <- c(as.Date("2014-01-01"), PANEL_END)

# ── 2. Monthly aggregations ─────────────────────────────
cat("\n--- 2. Monthly aggregations ---\n")

# From daily panel
monthly_panel <- panel |>
  mutate(ym = floor_date(date, "month")) |>
  group_by(ym) |>
  summarise(
    frx_persons       = sum(frx_persons),
    frx_incidents     = sum(frx_incidents),
    n_dead_missing    = sum(n_dead_missing),
    lcg_tcg_pushbacks = sum(lcg_tcg_pushbacks),
    unhcr_arrivals    = sum(arrivals, na.rm = TRUE),
    unhcr_days        = sum(!is.na(arrivals)),
    crossing_attempts = sum(crossing_attempts),
    .groups = "drop"
  )

# IOM death incidents aggregated monthly
iom_monthly_events <- iom_cmr |>
  mutate(ym = floor_date(date, "month")) |>
  group_by(ym) |>
  summarise(
    iom_incidents = n(),
    iom_lethal    = sum(dead_missing > 0),
    .groups = "drop"
  )

# Merge everything
monthly <- monthly_panel |>
  left_join(iom_monthly, by = "ym") |>
  left_join(iom_monthly_events, by = "ym") |>
  replace_na(list(iom_incidents = 0, iom_lethal = 0))

cat(sprintf("  Monthly tibble: %d months\n", nrow(monthly)))

# UNITED monthly: same structure but with UNITED deaths
monthly_united_daily <- panel |>
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
    unhcr_arrivals    = sum(arrivals, na.rm = TRUE),
    unhcr_days        = sum(!is.na(arrivals)),
    .groups = "drop"
  ) |>
  mutate(crossing_attempts = frx_persons + lcg_tcg_pushbacks + n_dead_missing)

monthly_united <- monthly_united_daily |>
  left_join(iom_monthly, by = "ym") |>
  replace_na(list(lcg = 0, tcg = 0))

cat(sprintf("  UNITED monthly tibble: %d months\n", nrow(monthly_united)))

# ── 3. Crossings + Deaths panels (separate IOM and UNITED figures) ──
cat("\n--- 3. Crossings + Deaths panels ---\n")

library(gridtext)
library(grid)
library(patchwork)

# Shared aesthetics
comp_order <- c("frontex_persons", "undetected_monthly",
                "lcg_persons", "tcg_persons", "deaths_missing")
comp_colours <- c(
  "frontex_persons"    = "#1F78B4",
  "undetected_monthly" = "#A6CEE3",
  "lcg_persons"        = "#333333",
  "tcg_persons"        = "#AAAAAA",
  "deaths_missing"     = "#D32F2F"
)

# Helper: build crossing-persons long data from a monthly tibble
build_crossing_long <- function(m, death_label) {
  comp_labels <- c(
    "frontex_persons"    = "Intercepted persons",
    "undetected_monthly" = "Undetected arrivals",
    "lcg_persons"        = "LCG pullbacks",
    "tcg_persons"        = "TCG pullbacks",
    "deaths_missing"     = death_label
  )
  m |>
    filter(ym <= PANEL_END) |>
    mutate(
      frontex_persons    = replace_na(frx_persons, 0),
      undetected_monthly = ifelse(unhcr_days > 0,
                                   pmax(unhcr_arrivals - frx_persons, 0), 0),
      lcg_persons        = replace_na(lcg, 0),
      tcg_persons        = replace_na(tcg, 0),
      deaths_missing     = n_dead_missing
    ) |>
    select(ym, frontex_persons, undetected_monthly,
           lcg_persons, tcg_persons, deaths_missing) |>
    pivot_longer(-ym, names_to = "component", values_to = "persons") |>
    mutate(component = factor(component, levels = comp_order,
                              labels = comp_labels[comp_order]))
}

# Shared boxed legend theme — outer-box style so single-guide and multi-guide
# legends render with identical geometry (needed for left-alignment).
cross_legend_box_theme <- theme(
  legend.position = "right",
  legend.background = element_blank(),
  legend.box.background = element_rect(fill = "grey97", colour = "grey80", linewidth = 0.5),
  legend.box.margin = margin(6, 8, 6, 8),
  legend.box = "vertical",
  legend.key = element_blank(),
  legend.margin = margin(0, 0, 0, 0),
  legend.text = element_text(size = 9),
  legend.key.size = unit(0.4, "cm"),
  legend.title = element_text(face = "bold", size = 9, margin = margin(b = 4)),
  legend.spacing = unit(0, "cm"),
  legend.spacing.y = unit(-0.1, "cm")
)

build_note_grob <- function(caption_text) {
  grid::grobTree(
    grid::rectGrob(gp = grid::gpar(fill = "grey97", col = "grey55", lwd = 0.6)),
    gridtext::textbox_grob(
      caption_text,
      x = grid::unit(0.012, "npc"), y = grid::unit(0.99, "npc"),
      width = grid::unit(0.976, "npc"), height = grid::unit(0.98, "npc"),
      hjust = 0, vjust = 1, halign = 0, valign = 1,
      gp = grid::gpar(fontsize = 10, col = "grey20", lineheight = 1.35),
      box_gp = grid::gpar(col = NA, fill = NA),
      padding = grid::unit(c(6, 6, 6, 14), "pt"),
      margin = grid::unit(c(0, 0, 0, 0), "pt")
    )
  )
}

# Panel (a): crossings stacked bars
build_cross_panel <- function(m, death_label) {
  cross_long <- build_crossing_long(m, death_label)
  comp_labels_vec <- c(
    "frontex_persons"    = "Intercepted persons",
    "undetected_monthly" = "Undetected arrivals",
    "lcg_persons"        = "LCG pullbacks",
    "tcg_persons"        = "TCG pullbacks",
    "deaths_missing"     = death_label
  )
  ggplot(cross_long, aes(x = ym, y = persons, fill = component)) +
    geom_col(position = "stack", width = 25) +
    scale_fill_manual(values = setNames(comp_colours, comp_labels_vec[names(comp_colours)]),
                      name = NULL) +
    geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red", linewidth = 0.5) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    coord_cartesian(xlim = PLOT_XLIM, clip = "off") +
    labs(
      title = "(a) Number of persons attempting the crossing, by crossing outcome",
      subtitle = "Red dashed line marks the signature of the 2017 Italy-Libya Memorandum of Understanding",
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

# Panel (b): deaths bars + fatality-rate line (dual axis)
build_deaths_panel <- function(m, death_label, shared_sf = NULL) {
  dr <- m |>
    filter(ym <= PANEL_END) |>
    mutate(death_rate = ifelse(crossing_attempts > 0,
                                n_dead_missing / crossing_attempts * 100, NA_real_))
  sf <- if (!is.null(shared_sf)) shared_sf else max(dr$n_dead_missing, na.rm = TRUE) / 15

  ggplot(dr, aes(x = ym)) +
    geom_col(aes(y = n_dead_missing, fill = death_label), width = 25) +
    geom_line(aes(y = death_rate * sf, colour = "Fatality rate (%)"),
              linewidth = 0.7) +
    scale_fill_manual(values = setNames("#D32F2F", death_label), name = NULL) +
    scale_colour_manual(values = c("Fatality rate (%)" = "#333333"), name = NULL) +
    scale_y_continuous(
      name = "Dead and Missing Migrants",
      sec.axis = sec_axis(~ . / sf, name = "Fatality rate (%)")
    ) +
    geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red", linewidth = 0.5) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    coord_cartesian(xlim = PLOT_XLIM,
                    ylim = if (!is.null(shared_sf)) c(0, shared_sf * 15) else NULL,
                    clip = "off") +
    labs(title = "(b) Number of persons dead and missing while attempting the crossing", x = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      plot.margin = margin(t = 16, r = 5, b = 5, l = 5)
    )
}

# Main figure builder: 2 panels + note + frame + overall title
build_crossings_figure <- function(monthly_data, death_label, deaths_source_label,
                                   main_title, output_path, shared_sf = NULL) {
  p_cross  <- build_cross_panel(monthly_data, death_label)
  p_deaths <- build_deaths_panel(monthly_data, death_label, shared_sf = shared_sf)

  legend_cross <- cowplot::get_legend(
    p_cross +
      theme(legend.position = "right") +
      cross_legend_box_theme +
      guides(fill = guide_legend(ncol = 1, title = "Crossing Outcomes"))
  )
  legend_deaths <- cowplot::get_legend(
    p_deaths +
      theme(legend.position = "right") +
      cross_legend_box_theme +
      guides(fill   = guide_legend(ncol = 1, order = 1, title = "Mortality"),
             colour = guide_legend(ncol = 1, order = 2, title = NULL))
  )

  # Small left offset on both legends so they don't collide with the plot's y-axis
  legend_cross_wrap <- cowplot::ggdraw() +
    cowplot::draw_grob(legend_cross,  x = 0.05, y = 0.5, hjust = 0, vjust = 0.5)
  legend_deaths_wrap <- cowplot::ggdraw() +
    cowplot::draw_grob(legend_deaths, x = 0.05, y = 0.5, hjust = 0, vjust = 0.5)

  caption_text <- paste0(
    "<b>Source:</b> Frontex JORA data (2014-2023), UNHCR sea arrivals to Italy (2014-2023), ",
    "IOM LCG/TCG pullback figures (2014-2023), ", deaths_source_label, " (2014-2023).<br/>",
    "<b>Note:</b> 'Persons attempting crossing' = Frontex persons + undetected arrivals (UNHCR sea arrivals minus ",
    "Frontex persons, floored at zero) + LCG/TCG pullbacks + dead/missing. 'Fatality rate' = dead/missing ÷ ",
    "persons attempting crossing. Undetected arrivals are available only for months with UNHCR data. The Frontex ",
    "'Not SAR: Coast Guard' events and IOM LCG/TCG pullbacks likely partially overlap; the overlap cannot be ",
    "bounded from public data but is reported to be concentrated from 2020 onwards."
  )
  note_grob <- build_note_grob(caption_text)

  layout_design <- "
AAAAAE
BBBBBF
CCCCCC
"
  fig <- (
    p_cross +
      p_deaths +
      wrap_elements(full = note_grob) +
      legend_cross_wrap +
      legend_deaths_wrap
  ) +
    plot_layout(design = layout_design, heights = c(1, 1, 0.22)) +
    plot_annotation(
      title = main_title,
      theme = theme(
        plot.title = element_text(face = "bold", size = 17, hjust = 0.5,
                                  margin = margin(b = 6)),
        plot.margin = margin(14, 14, 14, 14)
      )
    )

  fig_framed <- cowplot::ggdraw(fig) +
    cowplot::draw_grob(grid::rectGrob(
      gp = grid::gpar(col = "black", fill = NA, lwd = 2)
    ))

  ggsave(output_path, fig_framed, width = 14, height = 11, dpi = 300)
  cat("Saved:", output_path, "\n")
}

# Shared y-scale so the IOM and UNITED deaths panels are visually comparable
deaths_sf <- max(
  max(monthly$n_dead_missing[monthly$ym <= PANEL_END], na.rm = TRUE),
  max(monthly_united$n_dead_missing[monthly_united$ym <= PANEL_END], na.rm = TRUE)
) / 15

# IOM figure
build_crossings_figure(
  monthly_data = monthly,
  death_label = "Dead and Missing Migrants",
  deaths_source_label = "IOM Missing Migrants Project",
  main_title = "Crossing the Central Mediterranean: Persons Attempting Crossing and Fatal Outcomes per Month",
  output_path = fig_path("04_descriptive", "01_panel_crossings_iom.png"),
  shared_sf = deaths_sf
)

# UNITED figure
build_crossings_figure(
  monthly_data = monthly_united,
  death_label = "Dead and Missing Migrants",
  deaths_source_label = "UNITED List of Refugee Deaths",
  main_title = "Crossing the Central Mediterranean: Persons Attempting Crossing and Fatal Outcomes",
  output_path = fig_path("04_descriptive", "01_panel_crossings_united.png"),
  shared_sf = deaths_sf
)

# ── 4. Panel 2: Event type composition (detailed) ──────
cat("\n--- 4. Panel: Event type composition ---\n")

library(gridtext)
library(grid)
library(patchwork)

# Build event_type label from the cleaned dataset columns.
# Display spec (11 cats total):
#   SAR: EU operations       <- interceptor_type == "EU_ops"
#   SAR: Italian operations  <- interceptor_type == "Ita_ops" (Marina Militare,
#                               Mare Sicuro, Mare Nostrum)
#   SAR: NGO operations      <- interceptor_type == "NGO"
#   SAR: EU Coast Guard      <- interceptor_type == "EU_Coast_Guard"
#                               (CPV/CPB/OPV — member-state coast-guard vessels)
#   SAR: Other               <- Commercial + Land_patrol + No_intercept + Other + NA
#   Not SAR: EU Coast Guard  <- interceptor_type == "EU_Coast_Guard"
#   Not SAR: Land patrol     <- interceptor_type == "Land_patrol"
#   Not SAR: Self-arrived    <- interceptor_type == "No_intercept"
#   Not SAR: Other           <- EU_ops + Ita_ops + NGO + Commercial + Other + NA
#   LCG pullbacks / TCG pullbacks  (IOM monthly data)
# sar_ops == NA (339 Unknown-SAR rows) is treated as Not-SAR by design.
frx <- frx |>
  mutate(
    sar_true = !is.na(sar_ops) & sar_ops,
    event_type = case_when(
      sar_true & interceptor_type == "EU_ops"                  ~ "SAR: EU operations",
      sar_true & interceptor_type == "Ita_ops"                 ~ "SAR: Italian operations",
      sar_true & interceptor_type == "NGO"                     ~ "SAR: NGO operations",
      sar_true & interceptor_type == "EU_Coast_Guard"          ~ "SAR: EU Coast Guard",
      sar_true                                                 ~ "SAR: Other",
      interceptor_type == "EU_Coast_Guard"                     ~ "Not SAR: EU Coast Guard",
      interceptor_type == "Land_patrol"                        ~ "Not SAR: Land patrol",
      interceptor_type %in% c("NA", "No_intercept")            ~ "Not SAR: No intercep. (detection only)",
      TRUE                                                     ~ "Not SAR: Other"
    )
  ) |>
  select(-sar_true)

etype_order <- c(
  # SAR — darker blue to lighter blue
  "SAR: Italian operations",       # #08519C (darkest)
  "SAR: NGO operations",           # #3182BD
  "SAR: EU Coast Guard",           # #6BAED6
  "SAR: EU operations",            # #9ECAE1
  "SAR: Other",                    # #C6DBEF (lightest)
  # Not SAR — darker brown/red to lighter orange
  "Not SAR: No intercep. (detection only)", # #A63603 (darkest)
  "Not SAR: EU Coast Guard",       # #E6550D
  "Not SAR: Land patrol",          # #FD8D3C
  "Not SAR: Other",                # #FDBE85 (lightest)
  # IOM — darker to lighter grey
  "LCG pullbacks",                 # #252525
  "TCG pullbacks"                  # #969696
)

etype_colours <- c(
  # SAR — original palette. EU Coast Guard takes the slot formerly held by
  # "SAR: Commercial vessels" (#6BAED6) so it is clearly distinct from NGO.
  "SAR: EU operations"            = "#9ECAE1",  # was "SAR: EU operations (IRINI)"
  "SAR: Italian operations"       = "#08519C",  # was "SAR: Italian authorities"
  "SAR: NGO operations"           = "#3182BD",  # was "SAR: NGO"
  "SAR: EU Coast Guard"           = "#6BAED6",  # was "SAR: Commercial vessels" slot
  "SAR: Other"                    = "#C6DBEF",  # was "SAR: Other"
  # Not SAR — original palette
  "Not SAR: EU Coast Guard"       = "#E6550D",  # was "Not SAR: Coast Guard"
  "Not SAR: Land patrol"          = "#FD8D3C",  # was "Not SAR: Land patrol"
  "Not SAR: No intercep. (detection only)" = "#A63603",  # was "Not SAR: Self-arrived" slot (darkest)
  "Not SAR: Other"                = "#FDBE85",  # was "Not SAR: Other"
  # IOM pullbacks — greys
  "LCG pullbacks"                 = "#252525",
  "TCG pullbacks"                 = "#969696"
)

# Frontex interceptions/rescues by event type (counts)
frx_etype_n <- frx |>
  mutate(ym = floor_date(date, "month"), etype = event_type) |>
  count(ym, etype, name = "n") |>
  mutate(etype = factor(etype, levels = etype_order))

# Frontex persons by event type
frx_etype_p <- frx |>
  mutate(ym = floor_date(date, "month"), etype = event_type) |>
  group_by(ym, etype) |>
  summarise(persons = sum(num_persons, na.rm = TRUE), .groups = "drop")

# LCG/TCG pushbacks (persons)
ic_long <- iom_monthly |>
  filter(!is.na(lcg) | !is.na(tcg)) |>
  pivot_longer(cols = c(lcg, tcg), names_to = "type", values_to = "persons") |>
  filter(!is.na(persons)) |>
  mutate(etype = ifelse(type == "lcg", "LCG pullbacks", "TCG pullbacks")) |>
  select(ym, etype, persons)

persons_all <- bind_rows(frx_etype_p, ic_long) |>
  filter(ym <= PANEL_END) |>
  mutate(etype = factor(etype, levels = etype_order))

# Shared line and label for MoU
mou_vline <- geom_vline(
  xintercept = MOU_DATE,
  linetype = "dashed",
  colour = "red",
  linewidth = 0.5
)

mou_label <- annotate(
  "text",
  x = MOU_DATE,
  y = Inf,
  label = "2017 Italy-Libya MoU",
  colour = "red3",
  fontface = "bold",
  size = 3.2,
  vjust = 0.8,
  hjust = -0.02
)

# Shared theme for panels 4a-4c (no individual legends)
etype_theme_noleg <- theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(t = 16, r = 5, b = 5, l = 5)
  )

# 4a. Interceptions/rescues by full event type
p_etype_n <- ggplot(frx_etype_n, aes(x = ym, y = n, fill = etype)) +
  geom_col(position = "stack", width = 25) +
  scale_fill_manual(values = etype_colours, drop = FALSE, name = NULL) +
  mou_vline +
  mou_label +
  coord_cartesian(xlim = PLOT_XLIM, clip = "off") +
  labs(
    title = "Frontex interceptions by event type",
    subtitle = "Monthly count of interception events by SAR and non-SAR categories.",
    y = "Interceptions",
    x = NULL
  ) +
  etype_theme_noleg

# 4b. Persons by full event type
p_etype_p <- ggplot(persons_all, aes(x = ym, y = persons, fill = etype)) +
  geom_col(position = "stack", width = 25) +
  scale_fill_manual(values = etype_colours, name = NULL) +
  mou_vline +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  coord_cartesian(xlim = PLOT_XLIM, clip = "off") +
  labs(
    title = "(a) Number of persons intercepted",
    subtitle = "Red dashed line marks the signature of the 2017 Italy-Libya Memorandum of Understanding",
    y = "Number of persons",
    x = NULL
  ) +
  etype_theme_noleg +
  theme(plot.subtitle = element_text(size = 9, colour = "red3"))

# 4c. 100% stacked area (composition)
all_months <- seq(
  min(persons_all$ym, na.rm = TRUE),
  floor_date(PANEL_END, "month"),
  by = "month"
)

LOESS_SPAN <- 0.15

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

# Smooth shares per category, then re-normalise to 100%
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

p_etype_pct <- ggplot(persons_smoothed, aes(x = ym, y = pct_smooth, fill = etype)) +
  geom_area(position = "stack") +
  scale_fill_manual(values = etype_colours, name = NULL) +
  scale_y_continuous(name = "% of persons intercepted") +
  mou_vline +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  coord_cartesian(xlim = PLOT_XLIM, clip = "off") +
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

# 4d. Share line plot (aggregate: SAR / Not SAR / LCG / TCG)
share_colours <- c(
  "SAR"           = "#08519C",
  "Not SAR"       = "#E6550D",
  "LCG pullbacks" = "#252525",
  "TCG pullbacks" = "#969696"
)

share_lines <- persons_complete |>
  mutate(group = case_when(
    grepl("^SAR:", etype)     ~ "SAR",
    grepl("^Not SAR:", etype) ~ "Not SAR",
    TRUE                      ~ as.character(etype)
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
  coord_cartesian(xlim = PLOT_XLIM, ylim = c(0, 100), clip = "off") +
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

# Note box
caption_text <- paste0(
  "<b>Source:</b> Frontex JORA data (2014-2023), IOM MMP data (2014-2023).<br/>",
  "<b>Note:</b> The numbers under the Frontex category 'Not SAR: Coast Guard' and the IOM category ",
  "'LCG and TCG pullback' likely partially overlap. Frontex aerial assets routinely patrol the central ",
  "Mediterranean and sometimes relay sightings to the Libyan and Tunisian Coast Guards, who then carry out the ",
  "interception (Border Forensics, 2022; Human Rights Watch, 2022; Lighthouse Reports, 2024). ",
  "Such operations are recorded both in the Frontex incident database and in the IOM monthly pullback figures, ",
  "producing some double-counting of persons. The overlap cannot be bounded precisely from public data, but it ",
  "is unlikely to alter the broader picture. The phenomenon is reported to be concentrated from 2020 onwards."
)

note_grob <- grid::grobTree(
  grid::rectGrob(
    gp = grid::gpar(fill = "grey97", col = "grey55", lwd = 0.6)
  ),
  gridtext::textbox_grob(
    caption_text,
    x = grid::unit(0.012, "npc"),
    y = grid::unit(0.99, "npc"),
    width = grid::unit(0.976, "npc"),
    height = grid::unit(0.98, "npc"),
    hjust = 0,
    vjust = 1,
    halign = 0,
    valign = 1,
    gp = grid::gpar(fontsize = 10, col = "grey20", lineheight = 1.35),
    box_gp = grid::gpar(col = NA, fill = NA),
    padding = grid::unit(c(6, 6, 6, 14), "pt"),
    margin = grid::unit(c(0, 0, 0, 0), "pt")
  )
)

# ── Extract boxed side legends ──
legend_box_theme <- theme(
  legend.position = "right",
  legend.background = element_rect(fill = "grey97", colour = "grey80", linewidth = 0.5),
  legend.key = element_blank(),
  legend.margin = margin(10, 12, 10, 12),
  legend.text = element_text(size = 10),
  legend.key.size = unit(0.5, "cm"),
  legend.title = element_text(face = "bold", size = 10, margin = margin(b = 4))
)

legend_ab <- cowplot::get_legend(
  p_etype_pct + legend_box_theme +
    guides(fill = guide_legend(
      ncol = 1,
      title = "Actors involved in operations (detailed)"
    ))
)
legend_c <- cowplot::get_legend(
  p_share_lines + legend_box_theme +
    guides(colour = guide_legend(
      ncol = 1,
      title = "Actor categories involved in operations"
    ))
)

# Left-align both legends against the cell's left edge; vertically centre
legend_ab_wrap <- cowplot::ggdraw() +
  cowplot::draw_grob(legend_ab, x = 0, y = 0.5, hjust = 0, vjust = 0.5)
legend_c_wrap <- cowplot::ggdraw() +
  cowplot::draw_grob(legend_c,  x = 0, y = 0.5, hjust = 0, vjust = 0.5)

# Display versions (no inline legends)
p_etype_pct_noleg    <- p_etype_pct    + theme(legend.position = "none")
p_share_lines_noleg  <- p_share_lines  + theme(legend.position = "none")

# Assemble panel: plots left, legends right (legend_ab spans rows 1-2)
layout_design <- "
AAAAE
BBBBE
CCCCF
DDDDD
"

panel_event_type <- (
  p_etype_p +
    p_etype_pct_noleg +
    p_share_lines_noleg +
    wrap_elements(full = note_grob) +
    legend_ab_wrap +
    legend_c_wrap
) +
  plot_layout(
    design = layout_design,
    heights = c(1, 1, 1, 0.27)
  ) +
  plot_annotation(
    title = "Crossing the Central Mediterranean: Interception of Migrant Boats per Month, by Actor Type",
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, hjust = 0.5,
                                margin = margin(b = 6)),
      plot.margin = margin(14, 14, 14, 14)
    )
  )

# Wrap in ggdraw and add an outer black frame around everything (title included)
panel_event_type_framed <- cowplot::ggdraw(panel_event_type) +
  cowplot::draw_grob(grid::rectGrob(
    gp = grid::gpar(col = "black", fill = NA, lwd = 2)
  ))

ggsave(
  fig_path("04_descriptive", "01_panel_event_type.png"),
  panel_event_type_framed,
  width = 14,
  height = 14,
  dpi = 300
)
# ── 5. Boat type panel (share + avg persons) ───────────
cat("\n--- 5. Boat type ---\n")

BOAT_SPAN <- 0.15
boat_levels  <- c("Inflatable", "Wooden", "Other")
boat_colours <- c("Inflatable" = "#D32F2F", "Wooden" = "#1F78B4", "Other" = "#8B8386")

frx_boat <- frx |>
  mutate(ym = floor_date(date, "month"),
         boat = case_when(
           boat_category == "Inflatable" ~ "Inflatable",
           boat_category == "Wooden"     ~ "Wooden",
           TRUE                          ~ "Other"
         ),
         boat = factor(boat, levels = boat_levels)) |>
  group_by(ym, boat) |>
  summarise(n = n(),
            persons = sum(num_persons, na.rm = TRUE),
            .groups = "drop") |>
  group_by(ym) |>
  mutate(share = n / sum(n) * 100,
         avg_persons = ifelse(n > 0, persons / n, NA_real_)) |>
  ungroup()

# 5a. Share by boat type (LOESS-smoothed)
p_boat_share <- ggplot(frx_boat, aes(x = ym, y = share, colour = boat)) +
  geom_smooth(method = "loess", span = BOAT_SPAN, se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = boat_colours, name = NULL) +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  coord_cartesian(xlim = PLOT_XLIM, ylim = c(0, 100)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Share of Frontex interceptions in the Central Mediterranean, by boat type",
       subtitle = paste0("LOESS-smoothed (span = ", BOAT_SPAN, ")."),
       y = "Share of interception events", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# 5b. Average persons per boat (proxy for boat size, LOESS-smoothed)
p_boat_size <- ggplot(frx_boat, aes(x = ym, y = avg_persons, colour = boat)) +
  geom_smooth(method = "loess", span = BOAT_SPAN, se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = boat_colours, name = NULL) +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  coord_cartesian(xlim = PLOT_XLIM) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Average persons per boat, by boat type",
       subtitle = paste0("Proxy for vessel capacity. LOESS-smoothed (span = ", BOAT_SPAN, ")."),
       y = "Number of persons per boat", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.margin = margin(t = -5))

panel_boat <- p_boat_share / p_boat_size
ggsave(fig_path("04_descriptive", "01_panel_boat_type.png"),
       panel_boat, width = 12, height = 8, dpi = 300)

# ── 6. Annual summary table ──────────────────────────────
cat("\n--- 6. Annual summary table ---\n")

# IOM lethal events by year (incident count + average size)
iom_annual <- iom_cmr |>
  mutate(year = year(date)) |>
  group_by(year) |>
  summarise(
    iom_lethal = sum(dead_missing > 0),
    avg_deaths_per_lethal = ifelse(sum(dead_missing > 0) > 0,
                                    sum(dead_missing) / sum(dead_missing > 0), NA_real_),
    .groups = "drop"
  )

# Panel annual: matches the daily formula crossing_attempts = frx_persons +
# lcg_tcg_pushbacks + n_dead_missing (the single inclusive death count from 02).
panel_annual <- panel |>
  mutate(year = year(date)) |>
  group_by(year) |>
  summarise(
    crossing_attempts = sum(crossing_attempts),
    frx_incidents     = sum(frx_incidents),
    frx_persons       = sum(frx_persons),
    n_dead_missing    = sum(n_dead_missing),
    .groups = "drop"
  ) |>
  mutate(
    death_rate_pct    = ifelse(crossing_attempts > 0,
                                n_dead_missing / crossing_attempts * 100, NA_real_),
    avg_persons_event = ifelse(frx_incidents > 0,
                                frx_persons / frx_incidents, NA_real_)
  )

# IOM monthly arrivals (Italian sea arrivals) by year — used for arrival share
iom_arr_annual <- iom_monthly |>
  mutate(year = year(ym)) |>
  group_by(year) |>
  summarise(arrivals = sum(arrivals, na.rm = TRUE), .groups = "drop")

annual_summary <- panel_annual |>
  left_join(iom_annual,    by = "year") |>
  left_join(iom_arr_annual, by = "year") |>
  rename(total_crossings = crossing_attempts) |>
  mutate(arrival_share_pct = ifelse(total_crossings > 0,
                                     arrivals / total_crossings * 100, NA_real_))

cat("\n  Annual summary:\n")
print(annual_summary, n = 12)

# Build gt table (transposed: variables as rows, years as columns)
var_labels <- c(
  total_crossings       = "Total crossing attempts",
  frx_incidents         = "Frontex interceptions",
  n_dead_missing        = "Dead and missing (IOM)",
  iom_lethal            = "Lethal events",
  death_rate_pct        = "Death rate (%)",
  avg_persons_event     = "Persons / event",
  avg_deaths_per_lethal = "Dead and missing / lethal event",
  arrival_share_pct     = "Arrival share (%)"
)

# Format numbers before pivoting
annual_transposed <- annual_summary |>
  select(year, all_of(names(var_labels))) |>
  mutate(across(c(total_crossings, frx_incidents, n_dead_missing, iom_lethal),
                ~ formatC(round(.x), format = "d", big.mark = ",")),
         across(c(death_rate_pct, avg_persons_event, avg_deaths_per_lethal, arrival_share_pct),
                ~ ifelse(is.na(.x), "—", formatC(round(.x, 1), format = "f", digits = 1)))) |>
  pivot_longer(-year, names_to = "variable", values_to = "value") |>
  pivot_wider(names_from = year, values_from = value) |>
  mutate(variable = var_labels[variable],
         variable = factor(variable, levels = var_labels)) |>
  arrange(variable)

tbl <- annual_transposed |>
  gt(rowname_col = "variable") |>
  tab_header(
    title = "Annual Summary Statistics — Central Mediterranean Route",
    subtitle = "Data: Frontex Themis, IOM MMP, UNHCR (2014–2023)"
  ) |>
  tab_source_note("Total crossing attempts = Frontex persons + LCG/TCG pushbacks + dead/missing (IOM). Lower bound: excludes ~10% undetected arrivals.") |>
  tab_source_note("Death rate = dead+missing / total crossing attempts.") |>
  tab_source_note("Arrival share = Italian sea arrivals (IOM monthly) / total crossing attempts.")

# Save
gtsave(tbl, tbl_path("04_descriptive", "01_annual_summary.html"))

# LaTeX export
latex_out <- as_latex(tbl) |> as.character()
writeLines(latex_out, tbl_path("04_descriptive", "01_annual_summary.tex"))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
