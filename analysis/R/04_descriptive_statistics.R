# 04_descriptive_statistics.R
# ===========================
# Comprehensive descriptive statistics for the Central Mediterranean Route.
# Monthly stacked bar panels showing crossings, deaths, composition over time.
# Annual summary table with key numbers.
#
# Uses the canonical IOM filter (incident + split incident, all causes,
# CMR countries) — see 02_build_daily_panel.R for the rationale.
#
# Input:
#   data/processed/frontex_incidents.RDS
#   data/processed/iom_mmp_incidents.RDS
#   data/processed/iom_med_crossings_monthly.RDS
#   analysis/data/daily_panel_complete.RDS
#
# Output:
#   output/figures/desc_panel_crossings.png   (persons + events + deaths/rate)
#   output/figures/desc_panel_event_type.png  (events + persons + composition)
#   output/figures/desc_boat_type.png         (line plot)
#   output/tables/desc_annual_summary.html
#   output/tables/desc_annual_summary.tex

library(tidyverse)
library(lubridate)
library(gt)
library(scales)
library(patchwork)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
CMR_DEPARTURES <- c("Libya", "Tunisia", "Algeria")
CMR_INCIDENT_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia")

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
frx <- readRDS(file.path(BASE_DIR, "data", "processed", "frontex_incidents.RDS")) %>%
  filter(country_of_departure %in% CMR_DEPARTURES)
stopifnot(all(c("date", "boat_category", "event_type", "detected_by",
                "num_persons", "sar_flag") %in% names(frx)))
FRX_END <- max(frx$date)
cat(sprintf("  Frontex CMR: %d incidents (%s to %s)\n",
    nrow(frx), min(frx$date), FRX_END))

# 1c. IOM MMP incidents — canonical filter (incident + split, no cause filter,
#     CMR countries; see 02_build_daily_panel.R header for rationale).
iom_raw <- readRDS(file.path(BASE_DIR, "data", "processed", "iom_mmp_incidents.RDS"))
iom_cmr <- iom_raw %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) %in% c("incident", "split incident"),
         `Country of Incident` %in% CMR_INCIDENT_COUNTRIES) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date), date <= max(panel$date))
cat(sprintf("  IOM CMR: %d incidents\n", nrow(iom_cmr)))

# Panel end (canonical end date for descriptives) — read from the rebuilt
# daily panel, which is bounded by interceptions_daily_disagg.RDS coverage.
PANEL_END <- max(panel$date)

# 1d. IOM monthly crossings (cap at the panel end date)
iom_monthly <- readRDS(file.path(BASE_DIR, "data", "processed",
                                  "iom_med_crossings_monthly.RDS")) %>%
  transmute(
    ym = as.Date(date),
    arrivals = as.numeric(sea_arrivals_in_italy),
    lcg = as.numeric(interceptions_by_libyan_coast_guard),
    tcg = as.numeric(interceptions_by_tunisian_coast_guard)
  ) %>%
  filter(ym <= floor_date(PANEL_END, "month"))
cat(sprintf("  IOM monthly: %d months\n", nrow(iom_monthly)))


PLOT_XLIM <- c(as.Date("2014-01-01"), PANEL_END)

# ── 2. Monthly aggregations ─────────────────────────────
cat("\n--- 2. Monthly aggregations ---\n")

# From daily panel
monthly_panel <- panel %>%
  mutate(ym = floor_date(date, "month")) %>%
  group_by(ym) %>%
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
iom_monthly_events <- iom_cmr %>%
  mutate(ym = floor_date(date, "month")) %>%
  group_by(ym) %>%
  summarise(
    iom_incidents = n(),
    iom_lethal    = sum(dead_missing > 0),
    .groups = "drop"
  )

# Merge everything
monthly <- monthly_panel %>%
  left_join(iom_monthly, by = "ym") %>%
  left_join(iom_monthly_events, by = "ym") %>%
  replace_na(list(iom_incidents = 0, iom_lethal = 0))

cat(sprintf("  Monthly tibble: %d months\n", nrow(monthly)))

# ── 3. Panel 1: Crossings + Deaths ───────────────────────
cat("\n--- 3. Panel: Crossings + Deaths ---\n")

# 3a. Crossing attempts by persons
# The daily panel uses crossing_attempts = Frontex + LCG/TCG pushbacks + deaths
# (a lower bound). For the monthly plot we can add undetected arrivals
# computed at the monthly level, where UNHCR-Frontex timing mismatches
# wash out. The daily subtraction inflates the gap 3-27x, but at the
# monthly level the gap is genuine (~10% of UNHCR arrivals).
crossing_persons <- monthly %>%
  filter(ym <= PANEL_END) %>%
  mutate(
    frontex_persons    = replace_na(frx_persons, 0),
    undetected_monthly = ifelse(unhcr_days > 0,
                                 pmax(unhcr_arrivals - frx_persons, 0), 0),
    lcg_persons        = replace_na(lcg, 0),
    tcg_persons        = replace_na(tcg, 0),
    deaths_missing     = n_dead_missing
  ) %>%
  select(ym, frontex_persons, undetected_monthly,
         lcg_persons, tcg_persons, deaths_missing) %>%
  pivot_longer(-ym, names_to = "component", values_to = "persons")

comp_order <- c("frontex_persons", "undetected_monthly",
                "lcg_persons", "tcg_persons", "deaths_missing")
comp_labels <- c(
  "frontex_persons"    = "Persons intercepted/rescued during operations (Frontex)",
  "undetected_monthly" = "Undetected arrivals (monthly UNHCR - Frontex)",
  "lcg_persons"        = "LCG pushbacks (IOM)",
  "tcg_persons"        = "TCG pushbacks (IOM)",
  "deaths_missing"     = "Deaths and missing persons (IOM)"
)
comp_colours <- c(
  "frontex_persons"    = "#1F78B4",
  "undetected_monthly" = "#A6CEE3",
  "lcg_persons"        = "#333333",
  "tcg_persons"        = "#AAAAAA",
  "deaths_missing"     = "#D32F2F"
)

crossing_persons <- crossing_persons %>%
  mutate(component = factor(component, levels = comp_order,
                            labels = comp_labels[comp_order]))

legend_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        legend.margin = margin(t = -5))

p_cross_persons <- ggplot(crossing_persons, aes(x = ym, y = persons, fill = component)) +
  geom_col(position = "stack", width = 25) +
  scale_fill_manual(values = setNames(comp_colours, comp_labels[names(comp_colours)]),
                    name = NULL) +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  coord_cartesian(xlim = PLOT_XLIM) +
  labs(title = "Number of persons attempting crossing in the Central Mediterranean, by crossing outcome",
       subtitle = "Monthly totals. Undetected arrivals computed as max(UNHCR daily arrivals in Italy - Frontex interceptions/rescues, 0) at monthly level.",
       y = "Number of persons", x = NULL) +
  legend_theme +
  guides(fill = guide_legend(nrow = 2))

# 3b. Deaths + death rate (dual-axis, single plot)
deaths_rate <- monthly %>%
  filter(ym <= PANEL_END) %>%
  mutate(death_rate = ifelse(crossing_attempts > 0,
                              n_dead_missing / crossing_attempts * 100, NA_real_))

scale_factor <- max(deaths_rate$n_dead_missing, na.rm = TRUE) / 15

p_deaths <- ggplot(deaths_rate, aes(x = ym)) +
  geom_col(aes(y = n_dead_missing, fill = "Deaths and missing persons (IOM)"), width = 25) +
  geom_line(aes(y = death_rate * scale_factor, colour = "Fatality rate (%)"),
            linewidth = 0.7) +
  scale_fill_manual(values = c("Deaths and missing persons (IOM)" = "#D32F2F"), name = NULL) +
  scale_colour_manual(values = c("Fatality rate (%)" = "#333333"), name = NULL) +
  scale_y_continuous(
    name = "Number of dead and missing persons",
    sec.axis = sec_axis(~ . / scale_factor, name = "Fatality rate (%)")
  ) +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  coord_cartesian(xlim = PLOT_XLIM) +
  labs(title = "Number of dead and missing persons and fatality rate",
       subtitle = "Monthly totals. Fatality rate = deaths + missing / crossing attempts (Frontex events + LCG/TCG pushbacks + deaths + missing).",
       x = NULL) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  legend_theme

panel_crossings <- p_cross_persons / p_deaths
ggsave(file.path(BASE_DIR, "output", "figures", "desc_panel_crossings.png"),
       panel_crossings, width = 12, height = 8, dpi = 300)
cat("Saved: output/figures/desc_panel_crossings.png\n")

# ── 4. Panel 2: Event type composition (detailed) ──────
cat("\n--- 4. Panel: Event type composition ---\n")

# Full breakdown: SAR subtypes / Not SAR subtypes / LCG / TCG
etype_order <- c(
  "SAR: EU operations (IRINI)", "SAR: Commercial vessels",
  "SAR: NGO", "SAR: Italian authorities", "SAR: Other",
  "Not SAR: Coast Guard", "Not SAR: Land patrol",
  "Not SAR: Self-arrived", "Not SAR: Other",
  "LCG pushbacks", "TCG pushbacks"
)
etype_colours <- c(
  "SAR: EU operations (IRINI)" = "#9ECAE1",
  "SAR: Commercial vessels"    = "#6BAED6",
  "SAR: NGO"                   = "#3182BD",
  "SAR: Italian authorities"   = "#08519C",
  "SAR: Other"                 = "#C6DBEF",
  "Not SAR: Coast Guard"       = "#E6550D",
  "Not SAR: Land patrol"       = "#FD8D3C",
  "Not SAR: Self-arrived"      = "#A63603",
  "Not SAR: Other"             = "#FDBE85",
  "LCG pushbacks"          = "#252525",
  "TCG pushbacks"          = "#969696"
)

# Frontex interceptions/rescues by full event type (counts)
frx_etype_n <- frx %>%
  mutate(ym = floor_date(date, "month"),
         etype = event_type) %>%
  count(ym, etype, name = "n") %>%
  mutate(etype = factor(etype, levels = etype_order))

# Frontex persons by full event type
frx_etype_p <- frx %>%
  mutate(ym = floor_date(date, "month"),
         etype = event_type) %>%
  group_by(ym, etype) %>%
  summarise(persons = sum(num_persons, na.rm = TRUE), .groups = "drop")

# LCG/TCG pushbacks (persons)
ic_long <- iom_monthly %>%
  filter(!is.na(lcg) | !is.na(tcg)) %>%
  pivot_longer(cols = c(lcg, tcg), names_to = "type", values_to = "persons") %>%
  filter(!is.na(persons)) %>%
  mutate(etype = ifelse(type == "lcg", "LCG pushbacks", "TCG pushbacks")) %>%
  select(ym, etype, persons)

persons_all <- bind_rows(frx_etype_p, ic_long) %>%
  filter(ym <= PANEL_END) %>%
  mutate(etype = factor(etype, levels = etype_order))

# Shared theme for panels 4a-4c (no individual legends)
etype_theme_noleg <- theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# 4a. Interceptions/rescues by full event type
p_etype_n <- ggplot(frx_etype_n, aes(x = ym, y = n, fill = etype)) +
  geom_col(position = "stack", width = 25) +
  scale_fill_manual(values = etype_colours, drop = FALSE, name = NULL) +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red", linewidth = 0.5) +
  coord_cartesian(xlim = PLOT_XLIM) +
  labs(title = "Frontex interceptions/rescues by event type",
       subtitle = "Monthly count of interception/rescue events by SAR and non-SAR categories.",
       y = "Interceptions/rescues", x = NULL) +
  etype_theme_noleg

# 4b. Persons by full event type
p_etype_p <- ggplot(persons_all, aes(x = ym, y = persons, fill = etype)) +
  geom_col(position = "stack", width = 25) +
  scale_fill_manual(values = etype_colours, name = NULL) +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  coord_cartesian(xlim = PLOT_XLIM) +
  labs(title = "Number of persons intercepted/rescued in the Central Mediterranean, by interception/rescue type",
       subtitle = "Monthly totals. Frontex persons by SAR/non-SAR category; LCG and TCG pushbacks from IOM monthly data.",
       y = "Number of persons", x = NULL) +
  etype_theme_noleg

# 4c. 100% stacked area (composition)
all_months <- seq(min(persons_all$ym, na.rm = TRUE),
                  floor_date(PANEL_END, "month"), by = "month")

LOESS_SPAN <- 0.15

persons_complete <- expand_grid(
  ym = all_months,
  etype = factor(etype_order, levels = etype_order)
) %>%
  left_join(persons_all, by = c("ym", "etype")) %>%
  replace_na(list(persons = 0)) %>%
  group_by(ym) %>%
  mutate(total = sum(persons),
         pct = ifelse(total > 0, persons / total * 100, 0)) %>%
  ungroup()

# Smooth shares per category, then re-normalise to 100%
persons_smoothed <- persons_complete %>%
  group_by(etype) %>%
  mutate(ym_num = as.numeric(ym),
         pct_smooth = predict(loess(pct ~ ym_num, span = LOESS_SPAN))) %>%
  ungroup() %>%
  mutate(pct_smooth = pmax(pct_smooth, 0)) %>%
  group_by(ym) %>%
  mutate(pct_smooth = pct_smooth / sum(pct_smooth) * 100) %>%
  ungroup()

p_etype_pct <- ggplot(persons_smoothed, aes(x = ym, y = pct_smooth, fill = etype)) +
  geom_area(position = "stack") +
  scale_fill_manual(values = etype_colours, name = NULL) +
  scale_y_continuous(name = "Share of persons intercepted/rescued (%)") +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  coord_cartesian(xlim = PLOT_XLIM) +
  labs(title = "Composition of interception/rescue events, in percentage",
       subtitle = paste0("LOESS-smoothed (span = ", LOESS_SPAN, ")."),
       x = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(t = -5)) +
  guides(fill = guide_legend(nrow = 3))

# 4d. Share line plot (aggregate: SAR / Not SAR / LCG / TCG)
share_colours <- c(
  "SAR"                = "#08519C",
  "Not SAR"            = "#E6550D",
  "LCG pushbacks"  = "#252525",
  "TCG pushbacks"  = "#969696"
)

share_lines <- persons_complete %>%
  mutate(group = case_when(
    grepl("^SAR:", etype)     ~ "SAR",
    grepl("^Not SAR:", etype) ~ "Not SAR",
    TRUE                      ~ as.character(etype)
  )) %>%
  group_by(ym, group) %>%
  summarise(persons = sum(persons), .groups = "drop") %>%
  group_by(ym) %>%
  mutate(total = sum(persons),
         share = ifelse(total > 0, persons / total * 100, 0)) %>%
  ungroup() %>%
  select(-total) %>%
  mutate(group = factor(group, levels = names(share_colours)))

p_share_lines <- ggplot(share_lines, aes(x = ym, y = share, colour = group)) +
  geom_smooth(method = "loess", span = LOESS_SPAN, se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = share_colours, name = NULL) +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  coord_cartesian(xlim = PLOT_XLIM, ylim = c(0, 100)) +
  labs(title = "Aggregate composition of interception/rescue events, in percentage",
       subtitle = paste0("SAR, non-SAR, LCG and TCG pushbacks as share of total persons. LOESS-smoothed (span = ", LOESS_SPAN, ")."),
       y = "Share of persons intercepted/rescued (%)", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.margin = margin(t = -5))

# Single column: persons bar, composition area (with shared legend), aggregate lines
panel_event_type <- p_etype_p / p_etype_pct / p_share_lines

ggsave(file.path(BASE_DIR, "output", "figures", "desc_panel_event_type.png"),
       panel_event_type, width = 12, height = 14, dpi = 300)
cat("Saved: output/figures/desc_panel_event_type.png\n")

# ── 5. Boat type panel (share + avg persons) ───────────
cat("\n--- 5. Boat type ---\n")

BOAT_SPAN <- 0.15
boat_levels  <- c("Inflatable", "Wooden", "Other")
boat_colours <- c("Inflatable" = "#D32F2F", "Wooden" = "#1F78B4", "Other" = "#8B8386")

frx_boat <- frx %>%
  mutate(ym = floor_date(date, "month"),
         boat = case_when(
           boat_category == "Inflatable" ~ "Inflatable",
           boat_category == "Wooden"     ~ "Wooden",
           TRUE                          ~ "Other"
         ),
         boat = factor(boat, levels = boat_levels)) %>%
  group_by(ym, boat) %>%
  summarise(n = n(),
            persons = sum(num_persons, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(ym) %>%
  mutate(share = n / sum(n) * 100,
         avg_persons = ifelse(n > 0, persons / n, NA_real_)) %>%
  ungroup()

# 5a. Share by boat type (LOESS-smoothed)
p_boat_share <- ggplot(frx_boat, aes(x = ym, y = share, colour = boat)) +
  geom_smooth(method = "loess", span = BOAT_SPAN, se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = boat_colours, name = NULL) +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  coord_cartesian(xlim = PLOT_XLIM, ylim = c(0, 100)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Share of Frontex interceptions/rescues in the Central Mediterranean, by boat type",
       subtitle = paste0("LOESS-smoothed (span = ", BOAT_SPAN, ")."),
       y = "Share of interception/rescue events (%)", x = NULL) +
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
ggsave(file.path(BASE_DIR, "output", "figures", "desc_boat_type.png"),
       panel_boat, width = 12, height = 8, dpi = 300)
cat("Saved: output/figures/desc_boat_type.png\n")

# ── 6. Annual summary table ──────────────────────────────
cat("\n--- 6. Annual summary table ---\n")

# IOM lethal events by year (incident count + average size)
iom_annual <- iom_cmr %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    iom_lethal = sum(dead_missing > 0),
    avg_deaths_per_lethal = ifelse(sum(dead_missing > 0) > 0,
                                    sum(dead_missing) / sum(dead_missing > 0), NA_real_),
    .groups = "drop"
  )

# Panel annual: matches the daily formula crossing_attempts = frx_persons +
# lcg_tcg_pushbacks + n_dead_missing (the single inclusive death count from 02).
panel_annual <- panel %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    crossing_attempts = sum(crossing_attempts),
    frx_incidents     = sum(frx_incidents),
    frx_persons       = sum(frx_persons),
    n_dead_missing    = sum(n_dead_missing),
    .groups = "drop"
  ) %>%
  mutate(
    death_rate_pct    = ifelse(crossing_attempts > 0,
                                n_dead_missing / crossing_attempts * 100, NA_real_),
    avg_persons_event = ifelse(frx_incidents > 0,
                                frx_persons / frx_incidents, NA_real_)
  )

# IOM monthly arrivals (Italian sea arrivals) by year — used for arrival share
iom_arr_annual <- iom_monthly %>%
  mutate(year = year(ym)) %>%
  group_by(year) %>%
  summarise(arrivals = sum(arrivals, na.rm = TRUE), .groups = "drop")

annual_summary <- panel_annual %>%
  left_join(iom_annual,    by = "year") %>%
  left_join(iom_arr_annual, by = "year") %>%
  rename(total_crossings = crossing_attempts) %>%
  mutate(arrival_share_pct = ifelse(total_crossings > 0,
                                     arrivals / total_crossings * 100, NA_real_))

cat("\n  Annual summary:\n")
print(annual_summary, n = 12)

# Build gt table (transposed: variables as rows, years as columns)
var_labels <- c(
  total_crossings       = "Total crossing attempts",
  frx_incidents         = "Frontex interceptions/rescues",
  n_dead_missing        = "Dead and missing (IOM)",
  iom_lethal            = "Lethal events",
  death_rate_pct        = "Death rate (%)",
  avg_persons_event     = "Persons / event",
  avg_deaths_per_lethal = "Dead and missing / lethal event",
  arrival_share_pct     = "Arrival share (%)"
)

# Format numbers before pivoting
annual_transposed <- annual_summary %>%
  select(year, all_of(names(var_labels))) %>%
  mutate(across(c(total_crossings, frx_incidents, n_dead_missing, iom_lethal),
                ~ formatC(round(.x), format = "d", big.mark = ",")),
         across(c(death_rate_pct, avg_persons_event, avg_deaths_per_lethal, arrival_share_pct),
                ~ ifelse(is.na(.x), "—", formatC(round(.x, 1), format = "f", digits = 1)))) %>%
  pivot_longer(-year, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = year, values_from = value) %>%
  mutate(variable = var_labels[variable],
         variable = factor(variable, levels = var_labels)) %>%
  arrange(variable)

tbl <- annual_transposed %>%
  gt(rowname_col = "variable") %>%
  tab_header(
    title = "Annual Summary Statistics — Central Mediterranean Route",
    subtitle = "Data: Frontex Themis, IOM MMP, UNHCR (2014–2023)"
  ) %>%
  tab_source_note("Total crossing attempts = Frontex persons + LCG/TCG pushbacks + dead/missing (IOM). Lower bound: excludes ~10% undetected arrivals.") %>%
  tab_source_note("Death rate = dead+missing / total crossing attempts.") %>%
  tab_source_note("Arrival share = Italian sea arrivals (IOM monthly) / total crossing attempts.")

# Save
gtsave(tbl, file.path(BASE_DIR, "output", "tables", "desc_annual_summary.html"))
cat("Saved: output/tables/desc_annual_summary.html\n")

# LaTeX export
latex_out <- as_latex(tbl) %>% as.character()
writeLines(latex_out, file.path(BASE_DIR, "output", "tables", "desc_annual_summary.tex"))
cat("Saved: output/tables/desc_annual_summary.tex\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
