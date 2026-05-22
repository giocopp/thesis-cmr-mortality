# 05_policy_timeline.R
# ===========================
# Timeline of Central Mediterranean policy moments since 2013.
# Top half: specific policies/programs as dated callouts on the timeline spine.
# Bottom half: five phase bands showing the era each policy belongs to.
#
# Output:
#   output/figures/04_descriptive/05_policy_timeline.png

library(tidyverse)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

# The right edge matches the Frontex-bounded analytical panel; the left edge
# includes Mare Nostrum's pre-panel launch date for policy accuracy.
panel_dates <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel_complete.RDS")) |>
  summarise(
    start = min(date, na.rm = TRUE),
    end = max(date, na.rm = TRUE)
  )

timeline_start <- panel_dates$start
timeline_end   <- panel_dates$end
ongoing_end    <- timeline_end
regime_break   <- as.Date("2017-02-02")
mare_nostrum_start <- as.Date("2013-10-18")
phase_band_step <- 0.64
phase_band_half <- 0.25
phase_band_top <- -1.28

# ── 1. Phases ────────────────────────────────────────────
phases <- tibble(
  phase = 1:5,
  label = c(
    "1. State-led\nhumanitarian rescue",
    "2. Border control +\nactive NGO SAR",
    "3. Italy-Libya\nMoU period",
    "4. Closed ports +\nNGO containment",
    "5. Lamorgese\npartial rollback"
  ),
  start = as.Date(c(
    as.character(mare_nostrum_start),
    "2014-11-01", "2017-02-02", "2017-07-31", "2020-10-21"
  )),
  end = as.Date(c(
    "2014-10-31", "2018-06-10",
    as.character(ongoing_end), as.character(ongoing_end), "2023-01-01"
  )),
  family = c(
    "Rescue", "Border control", "Externalisation",
    "NGO containment", "Rollback"
  )
 ) |>
  mutate(
    band_y = phase_band_top - (phase - 1) * phase_band_step
  )

phase_y_bottom <- min(phases$band_y) - phase_band_half - 0.3

phase_colours <- c(
  "1" = "#08519C",  # pre-2017 rescue
  "2" = "#3182BD",  # pre-2017 border control + NGO SAR
  "3" = "#A63603",  # Italy-Libya MoU period
  "4" = "#E6550D",  # ongoing closed ports / NGO containment
  "5" = "#F2C94C"   # Lamorgese partial rollback
)

event_colours <- c(
  "Rescue" = "#08519C",
  "Border control" = "#3182BD",
  "Externalisation" = "#A63603",
  "NGO containment" = "#E6550D",
  "Rollback" = "#B7791F"
)

# ── 2. Policies & programs ───────────────────────────────
events <- tribble(
  ~date,         ~date_label,    ~label,                                  ~family,             ~label_x,       ~label_y,
  "2013-10-18",  "18 Oct 2013",  "Mare Nostrum\nstate-led rescue",         "Rescue",            "2014-05-01",   2.65,
  "2014-11-01",  "01 Nov 2014",  "Operation Triton\nborder control",       "Border control",    "2014-12-15",   1.35,
  "2015-06-22",  "22 Jun 2015",  "Operation Sophia\nanti-smuggling",       "Border control",    "2015-08-15",   2.35,
  "2017-02-02",  "02 Feb 2017",  "Italy-Libya MoU\nLibyan interception",   "Externalisation",   "2016-09-15",   1.35,
  "2017-07-31",  "Jul 2017",     "Minniti Code\nNGO port access\nconditional", "NGO containment", "2017-10-01",   2.75,
  "2018-06-10",  "10 Jun 2018",  "Salvini closed ports\ndelayed disembarkation", "NGO containment", "2018-07-15", 1.35,
  "2019-06-14",  "14 Jun 2019",  "Decreto Sicurezza Bis\nbans + fines",    "NGO containment",   "2019-06-01",   2.55,
  "2020-10-21",  "21 Oct 2020",  "Lamorgese reforms\npartial rollback",    "Rollback",          "2020-09-01",   1.35,
  "2023-01-02",  "02 Jan 2023",  "Piantedosi decree\ndistant ports,\none rescue", "NGO containment", "2022-07-01", 2.45
) |>
  mutate(
    across(c(date, label_x), as.Date),
    callout_body = paste0("\n", label),
    callout_lines = stringr::str_count(label, "\n") + 2L,
    date_y = label_y + (callout_lines - 1L) * 0.118
  )

link_gap_upper <- -0.24
link_gap_lower <- -0.80

event_phase_base <- events |>
  mutate(
    phase_ref = case_when(
      family == "Rescue" ~ 1L,
      family == "Border control" ~ 2L,
      family == "Externalisation" ~ 3L,
      family == "NGO containment" ~ 4L,
      family == "Rollback" ~ 5L
    )
  ) |>
  left_join(
    phases |> select(phase_ref = phase, band_y),
    by = "phase_ref"
  )

event_phase_links <- bind_rows(
  event_phase_base |>
    transmute(
      date,
      family,
      y_start = -0.18,
      y_end = link_gap_upper
    ),
  event_phase_base |>
    transmute(
      date,
      family,
      y_start = link_gap_lower,
      y_end = band_y + phase_band_half
    )
)

# Year tick marks
years <- as.Date(paste0(year(timeline_start):year(timeline_end), "-01-01"))

# ── 3. Plot ─────────────────────────────────────────────
timeline_left <- mare_nostrum_start - days(10)
label_col_left <- timeline_start - days(720)
label_col_right <- timeline_start - days(185)
label_col_marker_x <- label_col_left + days(42)
label_col_x <- label_col_left + days(78)
xlim_dates <- c(label_col_left, timeline_end + days(20))
pre_regime_label_x <- timeline_left + days(round(as.numeric(regime_break - timeline_left) / 2))
post_regime_label_x <- regime_break + days(round(as.numeric(xlim_dates[2] - regime_break) / 2))

p_timeline <- ggplot() +
  annotate(
    "rect",
    xmin = timeline_left, xmax = regime_break,
    ymin = phase_y_bottom, ymax = 3.62,
    fill = "#EAF4F8", alpha = 0.55
  ) +
  annotate(
    "rect",
    xmin = regime_break, xmax = xlim_dates[2],
    ymin = phase_y_bottom, ymax = 3.62,
    fill = "#FFF1E6", alpha = 0.55
  ) +
  annotate(
    "rect",
    xmin = label_col_left, xmax = label_col_right,
    ymin = phase_y_bottom, ymax = phase_band_top + phase_band_half + 0.55,
    fill = "white", colour = "grey80", linewidth = 0.35
  ) +
  annotate(
    "text",
    x = label_col_x,
    y = phase_band_top + phase_band_half + 0.25,
    label = "Policy periods",
    hjust = 0, vjust = 0.5,
    colour = "grey20", size = 3.6, fontface = "bold",
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = pre_regime_label_x, y = 3.33,
    label = "From institutional SAR\nto border control",
    colour = "#266577", size = 3.2, fontface = "bold",
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = post_regime_label_x, y = 3.33,
    label = "Border control externalization\nand deterrence",
    colour = "#8C3A11", size = 3.2, fontface = "bold",
    lineheight = 0.9
  ) +
  # Subtle calendar guides
  geom_segment(
    data = tibble(d = years),
    aes(x = d, xend = d),
    y = phase_y_bottom + 0.15, yend = 3.45,
    colour = "grey92", linewidth = 0.25
  ) +
  # Phase bars
  geom_rect(
    data = phases,
    aes(xmin = start, xmax = end,
        ymin = band_y - phase_band_half,
        ymax = band_y + phase_band_half,
        fill = factor(phase)),
    colour = "white", linewidth = 0.5, alpha = 0.92
  ) +
  # Dotted links: each dated policy moment points to its reference period.
  geom_segment(
    data = event_phase_links,
    aes(
      x = date, xend = date,
      y = y_start, yend = y_end,
      colour = family
    ),
    linetype = "dotted", linewidth = 1.1, alpha = 0.98
  ) +
  # Timeline spine
  annotate(
    "segment",
    x = timeline_left, xend = xlim_dates[2] - days(6),
    y = 0, yend = 0,
    linewidth = 1.35, colour = "grey8",
    arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed")
  ) +
  # Year tick marks
  geom_segment(
    data = tibble(d = years),
    aes(x = d, xend = d, y = -0.18, yend = 0.18),
    colour = "grey8", linewidth = 0.85
  ) +
  geom_text(
    data = tibble(d = years),
    aes(x = d, y = -0.52, label = format(d, "%Y")),
    size = 3.65, colour = "grey12", fontface = "bold"
  ) +
  geom_text(
    data = phases,
    aes(x = label_col_x, y = band_y, label = label),
    hjust = 0, vjust = 0.5,
    colour = "grey20", size = 3.25, fontface = "bold",
    lineheight = 0.95
  ) +
  geom_segment(
    data = phases,
    aes(
      x = label_col_marker_x, xend = label_col_marker_x,
      y = band_y - 0.18, yend = band_y + 0.18,
      colour = family
    ),
    linewidth = 2.8, lineend = "round"
  ) +
  # Leader lines from policy moments to callouts
  geom_segment(
    data = events,
    aes(x = date, xend = label_x, y = 0.17, yend = label_y - 0.25,
        colour = family),
    linewidth = 0.35, alpha = 0.75, lineend = "round"
  ) +
  # Event dots on the spine
  geom_point(
    data = events,
    aes(x = date, y = 0, colour = family),
    size = 3.55, shape = 21,
    fill = "white", stroke = 1.15
  ) +
  # Dated event labels above the spine
  geom_label(
    data = events,
    aes(x = label_x, y = label_y, label = callout_body, colour = family),
    fill = "white", linewidth = 0.45,
    label.padding = grid::unit(0.28, "lines"),
    label.r = grid::unit(0.10, "lines"),
    size = 3.25, lineheight = 0.92
  ) +
  geom_text(
    data = events,
    aes(x = label_x, y = date_y, label = date_label, colour = family),
    size = 3.25, fontface = "bold", lineheight = 0.92
  ) +
  scale_fill_manual(values = phase_colours, guide = "none") +
  scale_colour_manual(values = event_colours, guide = "none") +
  scale_x_date(limits = xlim_dates, expand = c(0.005, 0.005)) +
  scale_y_continuous(
    breaks = NULL,
    labels = NULL,
    limits = c(phase_y_bottom, 3.68),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    plot.margin = margin(8, 16, 8, 10)
  )

ggsave(
  fig_path("04_descriptive", "05_policy_timeline.png"),
  p_timeline,
  width = 11.8, height = 5.35, dpi = 300
)

cat("Saved:", fig_path("04_descriptive", "05_policy_timeline.png"), "\n")
