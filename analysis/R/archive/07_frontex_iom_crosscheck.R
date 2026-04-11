# 07_frontex_iom_crosscheck.R
# ===========================
# Cross-check Frontex Themis data with IOM MMP and daily panel.
# Also validate Rodriguez-Sanchez NGO SAR timeline against Frontex SAR flag.
#
# Input:
#   data/processed/frontex_incidents.RDS  (Frontex Themis incidents)
#   data/processed/iom_mmp_incidents.RDS (IOM MMP)
#   analysis/data/daily_panel.RDS (existing daily TS)
#   data/processed/archive/sar_ngo_ops_daily_RS.RDS (Rodriguez-Sanchez NGO timeline)
# Output:
#   output/figures/frontex_iom_crosscheck.png
#   output/tables/frontex_iom_crosscheck.txt

library(tidyverse)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
CMR_INCIDENT_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia")
SEA_CAUSES <- c("Drowning", "Mixed or unknown")

cat("============================================================\n")
cat("FRONTEX / IOM / SAR CROSS-CHECK\n")
cat("============================================================\n\n")

# ── 1. Load all datasets ──────────────────────────────────
cat("--- 1. Loading data ---\n")

# 1a. Frontex Themis
frx <- readRDS(file.path(BASE_DIR, "data", "processed", "frontex_incidents.RDS"))
cat(sprintf("  Frontex Themis: %d incidents, %s to %s\n",
    nrow(frx), min(frx$date), max(frx$date)))

# Aggregate Frontex to daily (CMR departures only)
frx_daily <- frx %>%
  filter(country_of_departure %in% c("Libya", "Tunisia", "Algeria")) %>%
  group_by(date) %>%
  summarise(
    frx_incidents  = n(),
    frx_persons    = sum(num_persons, na.rm = TRUE),
    frx_deaths     = sum(num_deaths, na.rm = TRUE),
    frx_migrants   = sum(num_migrants, na.rm = TRUE),
    frx_pct_sar    = mean(sar_flag == 1, na.rm = TRUE) * 100,
    frx_n_sar      = sum(sar_flag == 1, na.rm = TRUE),
    frx_n_ngo      = sum(ngo_involved == TRUE, na.rm = TRUE),
    frx_ngo_persons = sum(num_persons[ngo_involved == TRUE], na.rm = TRUE),
    frx_pct_ngo    = mean(ngo_involved == TRUE, na.rm = TRUE) * 100,
    frx_n_inflatable = sum(grepl("inflatable|rubber|zodiac|dinghy",
                                  transport_type, ignore.case = TRUE), na.rm = TRUE),
    .groups = "drop"
  )

# 1b. IOM MMP — CMR sea deaths (core corridor)
iom_raw <- readRDS(file.path(BASE_DIR, "data", "processed",
                               "iom_mmp_incidents.RDS"))
iom_cmr <- iom_raw %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) == "incident",
         `Cause of death (category)` %in% SEA_CAUSES) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE),
         in_cmr_countries = `Country of Incident` %in% CMR_INCIDENT_COUNTRIES) %>%
  filter(!is.na(date))

iom_daily <- iom_cmr %>%
  group_by(date) %>%
  summarise(
    iom_incidents    = n(),
    iom_dead_missing = sum(dead_missing),
    iom_core_deaths  = sum(dead_missing * in_cmr_countries),
    .groups = "drop"
  )

# 1c. Existing daily panel
panel <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS")) %>%
  select(date, deaths, n_incidents, arrivals, crossings)

# 1d. Rodriguez-Sanchez SAR
sar <- readRDS(file.path(BASE_DIR, "data", "processed", "archive", "sar_ngo_ops_daily_RS.RDS"))

cat(sprintf("  IOM CMR sea deaths: %d incidents\n", nrow(iom_cmr)))
cat(sprintf("  Daily panel: %d days\n", nrow(panel)))
cat(sprintf("  SAR daily: %d days\n", nrow(sar)))

# ── 2. Merge all to daily ─────────────────────────────────
cat("\n--- 2. Merging to daily panel ---\n")

# Coverage boundaries
FRX_START <- min(frx$date)
FRX_END   <- max(frx$date)
RS_NGO_END <- as.Date("2021-10-31")  # last date with R-S NGO vessel data

combined <- tibble(date = seq(as.Date("2014-01-01"), as.Date("2023-12-31"), by = "day")) %>%
  left_join(frx_daily, by = "date") %>%
  left_join(iom_daily, by = "date") %>%
  left_join(panel, by = "date") %>%
  left_join(sar, by = "date") %>%
  mutate(
    # Zero-fill Frontex only within its coverage period
    across(c(frx_incidents, frx_persons, frx_deaths, frx_migrants,
             frx_pct_sar, frx_n_sar, frx_n_ngo, frx_ngo_persons, frx_pct_ngo,
             frx_n_inflatable),
           ~ ifelse(date >= FRX_START & date <= FRX_END & is.na(.x), 0, .x)),
    # Zero-fill IOM and deaths (IOM covers full period)
    across(c(iom_incidents, iom_dead_missing, iom_core_deaths,
             deaths, n_incidents),
           ~ replace_na(.x, 0)),
    # R-S NGO data: set to NA after coverage ends (zeros are not real)
    n_ngo_vessels = ifelse(date > RS_NGO_END, NA_real_, n_ngo_vessels),
    year = year(date), month = month(date)
  )

cat(sprintf("  Combined panel: %d days\n", nrow(combined)))

# ── 3. Compare deaths: Frontex vs IOM ─────────────────────
cat("\n--- 3. Deaths comparison: Frontex vs IOM ---\n")

yearly_deaths <- combined %>%
  group_by(year) %>%
  summarise(
    frx_deaths       = sum(frx_deaths),
    iom_all_cmr      = sum(iom_dead_missing),
    iom_core         = sum(iom_core_deaths),
    panel_deaths     = sum(deaths),
    .groups = "drop"
  )
cat("\n  Annual deaths by source:\n")
print(yearly_deaths, n = 12)

cat("\n  Daily correlation (days with any deaths in either source):\n")
death_days <- combined %>% filter(frx_deaths > 0 | iom_dead_missing > 0)
cat(sprintf("    N days: %d\n", nrow(death_days)))
if (nrow(death_days) > 10) {
  cat(sprintf("    Cor(frx_deaths, iom_all_cmr): %.3f\n",
      cor(death_days$frx_deaths, death_days$iom_dead_missing, use = "complete.obs")))
}

# ── 4. Compare volume: Frontex persons vs panel arrivals ──
cat("\n--- 4. Volume comparison: Frontex persons vs panel ---\n")

yearly_volume <- combined %>%
  group_by(year) %>%
  summarise(
    frx_persons  = sum(frx_persons),
    frx_migrants = sum(frx_migrants),
    panel_arrivals  = sum(arrivals, na.rm = TRUE),
    panel_crossings = sum(crossings, na.rm = TRUE),
    .groups = "drop"
  )
cat("\n  Annual volume by source:\n")
print(yearly_volume, n = 12)

# ── 5. NGO SAR validation: Frontex NGO flag vs R-S NGO count
cat("\n--- 5. NGO SAR validation: Frontex vs Rodriguez-Sanchez ---\n")

# Weekly aggregation for comparison
weekly_sar <- combined %>%
  filter(year >= 2014, year <= 2021) %>%
  mutate(week = floor_date(date, "week")) %>%
  group_by(week) %>%
  summarise(
    frx_ngo_incidents = sum(frx_n_ngo),
    frx_total_incidents = sum(frx_incidents),
    frx_ngo_pct = ifelse(sum(frx_incidents) > 0,
                          sum(frx_n_ngo) / sum(frx_incidents) * 100, NA),
    rs_ngo_vessels = mean(n_ngo_vessels, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(frx_ngo_pct), frx_total_incidents > 0)

cat(sprintf("  Weekly obs with Frontex data: %d\n", nrow(weekly_sar)))
if (nrow(weekly_sar) > 10) {
  cat(sprintf("  Cor(frx_ngo_pct, rs_ngo_vessels): %.3f\n",
      cor(weekly_sar$frx_ngo_pct, weekly_sar$rs_ngo_vessels, use = "complete.obs")))
}

# Monthly for a cleaner view
monthly_sar <- combined %>%
  filter(year >= 2014, year <= 2021) %>%
  mutate(ym = floor_date(date, "month")) %>%
  group_by(ym) %>%
  summarise(
    frx_ngo_pct = ifelse(sum(frx_incidents) > 0,
                          sum(frx_n_ngo) / sum(frx_incidents) * 100, NA),
    rs_ngo_vessels = mean(n_ngo_vessels, na.rm = TRUE),
    frx_incidents = sum(frx_incidents),
    .groups = "drop"
  ) %>%
  filter(!is.na(frx_ngo_pct), frx_incidents > 0)

cat(sprintf("  Monthly obs: %d\n", nrow(monthly_sar)))
if (nrow(monthly_sar) > 10) {
  cat(sprintf("  Cor(frx_ngo_pct, rs_ngo_vessels) monthly: %.3f\n",
      cor(monthly_sar$frx_ngo_pct, monthly_sar$rs_ngo_vessels, use = "complete.obs")))
}

cat("\n  Monthly SAR comparison (sample):\n")
monthly_sar %>%
  mutate(ym = format(ym, "%Y-%m")) %>%
  filter(ym %in% c("2015-06","2016-06","2017-06","2017-12","2018-06","2019-06","2020-06","2021-06")) %>%
  print()

# ── 6. Frontex-unique info: boat type breakdown ───────────
cat("\n--- 6. What Frontex adds: boat type over time ---\n")

boat_yearly <- frx %>%
  mutate(year = year(date),
         boat_cat = case_when(
           grepl("inflatable|rubber|zodiac|dinghy", transport_type, ignore.case = TRUE) ~ "Inflatable",
           grepl("wooden|wood", transport_type, ignore.case = TRUE) ~ "Wooden",
           grepl("metal|makeshift", transport_type, ignore.case = TRUE) ~ "Metal/makeshift",
           grepl("fibre|fiber|glass", transport_type, ignore.case = TRUE) ~ "Fibreglass",
           grepl("fishing", transport_type, ignore.case = TRUE) ~ "Fishing",
           grepl("sailing", transport_type, ignore.case = TRUE) ~ "Sailing",
           TRUE ~ "Other"
         )) %>%
  group_by(year, boat_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(year) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup()

cat("\n  Inflatable share by year:\n")
boat_yearly %>% filter(boat_cat == "Inflatable") %>%
  select(year, n, pct) %>% print(n = 12)

# ── 7. Plots ──────────────────────────────────────────────
cat("\n--- 7. Plots ---\n")

# Monthly aggregation for plots
monthly <- combined %>%
  mutate(ym = floor_date(date, "month")) %>%
  group_by(ym) %>%
  summarise(
    frx_deaths = if (all(is.na(frx_deaths))) NA_real_ else sum(frx_deaths, na.rm = TRUE),
    iom_deaths = sum(iom_dead_missing),
    frx_persons = if (all(is.na(frx_persons))) NA_real_ else sum(frx_persons, na.rm = TRUE),
    panel_arrivals = if (all(is.na(arrivals))) NA_real_ else sum(arrivals, na.rm = TRUE),
    frx_ngo_count = if (all(is.na(frx_n_ngo))) NA_real_ else sum(frx_n_ngo, na.rm = TRUE),
    frx_ngo_persons = if (all(is.na(frx_ngo_persons))) NA_real_ else sum(frx_ngo_persons, na.rm = TRUE),
    frx_total = if (all(is.na(frx_incidents))) NA_real_ else sum(frx_incidents, na.rm = TRUE),
    rs_ngo = mean(n_ngo_vessels, na.rm = TRUE),
    .groups = "drop"
  )

# Total crossings per month = arrivals + LCG/TCG interceptions + deaths
ic_monthly <- readRDS(file.path(BASE_DIR, "data", "processed",
                                 "iom_med_crossings_monthly.RDS")) %>%
  transmute(ym = as.Date(date),
            lcg = as.numeric(interceptions_by_libyan_coast_guard),
            tcg = as.numeric(interceptions_by_tunisian_coast_guard))

monthly <- monthly %>%
  left_join(ic_monthly, by = "ym") %>%
  mutate(
    # Arrivals: use UNHCR when available, Frontex persons as fallback
    arrivals_best = ifelse(!is.na(panel_arrivals), panel_arrivals, frx_persons),
    total_crossings = replace_na(arrivals_best, 0) +
                      replace_na(lcg, 0) + replace_na(tcg, 0) + iom_deaths,
    # NGO persons as share of total crossings (all in persons)
    frx_ngo_share = ifelse(!is.na(frx_ngo_persons) & total_crossings > 0,
                            frx_ngo_persons / total_crossings * 100, NA_real_)
  )

p1 <- ggplot(monthly, aes(x = ym)) +
  geom_line(aes(y = frx_deaths, colour = "deaths recorded by Frontex during interceptions"), linewidth = 0.6) +
  geom_line(aes(y = iom_deaths, colour = "deaths recorded by IOM MMP"), linewidth = 0.6) +
  scale_colour_manual(values = c("deaths recorded by Frontex during interceptions" = "#D32F2F", 
                                 "deaths recorded by IOM MMP" = "#1565C0")) +
  labs(title = "Monthly death count in Central Mediterranean",
       y = "Number of deaths and missing persons", x = NULL, colour = "Source") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

p2 <- ggplot(monthly, aes(x = ym)) +
  geom_line(aes(y = frx_persons, colour = "number of persons detected (Frontex)"), linewidth = 0.6) +
  geom_line(aes(y = panel_arrivals, colour = "number of arrivals to Italy (UNHCR)"), linewidth = 0.6) +
  scale_colour_manual(values = c("number of persons detected (Frontex)" = "#D32F2F", 
                                 "number of arrivals to Italy (UNHCR)" = "#E08A00")) +
  labs(title = "Monthly volume of crossings in Central Medierranean",
       y = "Number of persons", x = NULL, colour = "Source") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# p3a: NGO operation counts — Frontex NGO incidents vs R-S vessel count
p3a <- ggplot(monthly, aes(x = ym)) +
  geom_line(data = monthly %>% filter(!is.na(frx_ngo_count)),
            aes(y = frx_ngo_count, colour = "Number of NGO operations (Frontex)"), linewidth = 0.6) +
  geom_line(data = monthly %>% filter(!is.na(rs_ngo)),
            aes(y = rs_ngo * 8, colour = "Active NGO vessels (Rodriguez-Sanchez et al., 2023)"), linewidth = 0.6) +
  scale_colour_manual(values = c("Number of NGO operations (Frontex)" = "#D32F2F",
                                  "Active NGO vessels (Rodriguez-Sanchez et al., 2023)" = "#2E7D32")) +
  scale_y_continuous(
    name = "Number of NGO operations",
    limits = c(0, 104),
    sec.axis = sec_axis(~ . / 8, name = "Number of active NGO vessels")
  ) +
  labs(title = "NGO operations: Frontex incident count vs active vessel count",
       x = NULL, colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# p3b: NGO rescue share of total crossings vs R-S vessel count
p3b <- ggplot(monthly, aes(x = ym)) +
  geom_line(data = monthly %>% filter(!is.na(frx_ngo_share)),
            aes(y = frx_ngo_share, colour = "Share of people crossing rescued by NGOs"), linewidth = 0.6) +
  geom_line(data = monthly %>% filter(!is.na(rs_ngo)),
            aes(y = rs_ngo * 8, colour = "Active NGO vessels (Rodriguez-Sanchez et al., 2023)"), linewidth = 0.6) +
  scale_colour_manual(values = c("Share of people crossing rescued by NGOs" = "#D32F2F",
                                  "Active NGO vessels (Rodriguez-Sanchez et al., 2023)" = "#2E7D32")) +
  scale_y_continuous(
    name = "Percentage of people crossing rescued by NGOs",
    limits = c(0, 104),
    sec.axis = sec_axis(~ . / 8, name = "Number of active NGO vessels")
  ) +
  labs(title = "NGO rescue share and vessel presence in the Central Mediterranean",
       x = NULL, colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# p4: Stacked area — event type composition over time
# Monthly counts by event_type (CMR departures only)
frx_cmr <- frx %>%
  filter(country_of_departure %in% c("Libya", "Tunisia", "Algeria"))

monthly_type <- frx_cmr %>%
  mutate(ym = floor_date(date, "month")) %>%
  group_by(ym, event_type) %>%
  summarise(n = n(), .groups = "drop")

# Order: SAR categories first (bottom), then Not SAR (top)
type_order <- c("SAR: EU operations (IRINI)",
                "SAR: Commercial vessels", "SAR: NGO",
                "SAR: Italian authorities", "SAR: Other",
                "Not SAR: Coast Guard", "Not SAR: Land patrol",
                "Not SAR: Self-arrived", "Not SAR: Other")

type_colours <- c(
  "SAR: EU operations (IRINI)" = "#7570B3",
  "SAR: Commercial vessels"    = "#E7298A",
  "SAR: NGO"                   = "#D95F02",
  "SAR: Italian authorities"   = "#1F78B4",
  "SAR: Other"                 = "#A6CEE3",
  "Not SAR: Coast Guard"       = "#B2DF8A",
  "Not SAR: Land patrol"       = "#FDBF6F",
  "Not SAR: Self-arrived"      = "#FB9A99",
  "Not SAR: Other"             = "#CAB2D6"
)

type_order_b <- c(type_order,
                  "LCG interceptions", "TCG interceptions")
type_colours_b <- c(type_colours,
                    "LCG interceptions" = "#333333",
                    "TCG interceptions" = "#AAAAAA")

monthly_type <- monthly_type %>%
  mutate(event_type = factor(event_type, levels = type_order_b))

# Monthly persons by event_type (CMR departures)
monthly_persons <- frx_cmr %>%
  mutate(ym = floor_date(date, "month")) %>%
  group_by(ym, event_type) %>%
  summarise(persons = sum(num_persons, na.rm = TRUE), .groups = "drop") %>%
  mutate(event_type = factor(event_type, levels = type_order))

# LCG/TCG monthly interceptions (persons) from clean IOM crossings data
ic <- readRDS(file.path(BASE_DIR, "data", "processed",
                         "iom_med_crossings_monthly.RDS")) %>%
  transmute(ym = as.Date(date),
            lcg = as.numeric(interceptions_by_libyan_coast_guard),
            tcg = as.numeric(interceptions_by_tunisian_coast_guard))

PLOT_XLIM <- c(as.Date("2014-01-01"), FRX_END)

# Panel A: incidents by event type
p4a <- ggplot(monthly_type, aes(x = ym, y = n, fill = event_type)) +
  geom_col(position = "stack", width = 25) +
  scale_fill_manual(values = type_colours_b, drop = FALSE, name = NULL) +
  coord_cartesian(xlim = PLOT_XLIM) +
  labs(title = "Number of Frontex detections by event type",
       y = "Number of detections", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Panel B: persons by event type + LCG/TCG as stacked bars
# Add interceptions as extra event types
ic_long <- ic %>%
  filter(!is.na(lcg) | !is.na(tcg)) %>%
  tidyr::pivot_longer(cols = c(lcg, tcg), names_to = "type", values_to = "persons") %>%
  filter(!is.na(persons)) %>%
  mutate(event_type = case_when(
    type == "lcg" ~ "LCG interceptions",
    type == "tcg" ~ "TCG interceptions"
  )) %>%
  select(ym, event_type, persons)

# Combine Frontex persons + interceptions
persons_all <- bind_rows(
  monthly_persons %>% select(ym, event_type, persons),
  ic_long
)

persons_all <- persons_all %>%
  filter(ym <= FRX_END) %>%
  mutate(event_type = factor(event_type, levels = type_order_b))

# Total crossings = arrivals + deaths + interceptions (monthly)
crossings_monthly <- readRDS(file.path(BASE_DIR, "data", "processed",
                                        "iom_med_crossings_monthly.RDS")) %>%
  transmute(
    ym = as.Date(date),
    arrivals = as.numeric(sea_arrivals_in_italy),
    interceptions = replace_na(as.numeric(interceptions_by_libyan_coast_guard), 0) +
                    replace_na(as.numeric(interceptions_by_tunisian_coast_guard), 0)
  )

iom_monthly_deaths <- readRDS(file.path(BASE_DIR, "data", "processed",
                                         "iom_mmp_incidents.RDS")) %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) != "sub-incident") %>%
  mutate(ym = floor_date(as.Date(incident_date_clean), "month"),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(ym)) %>%
  group_by(ym) %>%
  summarise(deaths = sum(dead_missing), .groups = "drop")

crossings_monthly <- crossings_monthly %>%
  left_join(iom_monthly_deaths, by = "ym") %>%
  replace_na(list(deaths = 0)) %>%
  mutate(total_crossings = replace_na(arrivals, 0) + deaths + interceptions) %>%
  filter(ym <= FRX_END)


p4b <- ggplot(persons_all, aes(x = ym, y = persons, fill = event_type)) +
  geom_col(position = "stack", width = 25) +
  geom_line(data = crossings_monthly,
            aes(x = ym, y = total_crossings, fill = NULL),
            colour = "black", linewidth = 1) +
  scale_fill_manual(values = type_colours_b, name = NULL) +
  coord_cartesian(xlim = PLOT_XLIM) +
  labs(title = "Number of persons crossing by event type",
       y = "Number of persons", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Panel C: 100% stacked area (composition)
# Complete grid: every month × every event type, fill 0s
all_months <- seq(min(persons_all$ym, na.rm = TRUE),
                  floor_date(FRX_END, "month"), by = "month")

complete_grid <- tidyr::expand_grid(
  ym = all_months,
  event_type = factor(type_order_b, levels = type_order_b)
)

persons_complete <- complete_grid %>%
  left_join(persons_all %>% mutate(event_type = factor(event_type, levels = type_order_b)),
            by = c("ym", "event_type")) %>%
  replace_na(list(persons = 0)) %>%
  group_by(ym) %>%
  mutate(total = sum(persons),
         pct = ifelse(total > 0, persons / total * 100, 0)) %>%
  ungroup()

p4c <- ggplot(persons_complete, aes(x = ym, y = pct, fill = event_type)) +
  geom_area(position = "stack") +
  scale_fill_manual(values = type_colours_b, name = NULL) +
  scale_y_continuous(name = "Share of detected persons (%)") +
  coord_cartesian(xlim = PLOT_XLIM) +
  labs(title = "CMR crossing composition by event type (%)",
       x = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8))

p4b <- p4b + theme(legend.position = "right",
                    legend.justification = "center",
                    legend.text = element_text(size = 9),
                    legend.key.size = unit(0.4, "cm"))

p4_combined <- (p4a + theme(legend.position = "none")) /
               p4b /
               (p4c + theme(legend.position = "none"))
ggsave(file.path(BASE_DIR, "output", "figures", "frontex_event_types.png"),
       p4_combined, width = 12, height = 14, dpi = 200)
cat("Saved: output/figures/frontex_event_types.png\n")

# Figure 1: p1 + p2 (deaths and volume comparison)
fig_crosscheck <- p1 / p2
ggsave(file.path(BASE_DIR, "output", "figures", "frontex_iom_crosscheck.png"),
       fig_crosscheck, width = 11, height = 8, dpi = 200)
cat("Saved: output/figures/frontex_iom_crosscheck.png\n")

# Figure 2: p3a + p3b (NGO operations and rescue share)
fig_ngo <- p3a / p3b
ggsave(file.path(BASE_DIR, "output", "figures", "frontex_ngo_validation.png"),
       fig_ngo, width = 11, height = 8, dpi = 200)
cat("Saved: output/figures/frontex_ngo_validation.png\n")

# Figure 4: Crossing attempt components
# Build monthly components: Frontex persons, UNHCR extra, LCG, TCG, deaths
crossing_components <- monthly %>%
  filter(ym <= FRX_END) %>%
  mutate(
    frontex_persons = replace_na(frx_persons, 0),
    unhcr_extra     = pmax(replace_na(panel_arrivals, 0) - replace_na(frx_persons, 0), 0),
    lcg_persons     = replace_na(lcg, 0),
    tcg_persons     = replace_na(tcg, 0),
    deaths_missing  = iom_deaths
  ) %>%
  select(ym, frontex_persons, unhcr_extra, lcg_persons, tcg_persons, deaths_missing) %>%
  tidyr::pivot_longer(-ym, names_to = "component", values_to = "persons")

comp_order <- c("frontex_persons", "unhcr_extra",
                "lcg_persons", "tcg_persons", "deaths_missing")
comp_labels <- c(
  "frontex_persons" = "Persons detected during operations (Frontex)",
  "unhcr_extra"     = "Arrivals not detected during operations (UNHCR - Frontex)",
  "lcg_persons"     = "LCG interceptions (IOM)",
  "tcg_persons"     = "TCG interceptions (IOM)",
  "deaths_missing"  = "Deaths and missing persons (IOM)"
)
comp_colours <- c(
  "frontex_persons" = "#1F78B4",
  "unhcr_extra"     = "#A6CEE3",
  "lcg_persons"     = "#333333",
  "tcg_persons"     = "#AAAAAA",
  "deaths_missing"  = "#D32F2F"
)

crossing_components <- crossing_components %>%
  mutate(component = factor(component, levels = comp_order, labels = comp_labels[comp_order]))

p5 <- ggplot(crossing_components, aes(x = ym, y = persons, fill = component)) +
  geom_col(position = "stack", width = 25) +
  scale_fill_manual(values = setNames(comp_colours, comp_labels[names(comp_colours)]),
                    name = NULL) +
  coord_cartesian(xlim = PLOT_XLIM) +
  labs(title = "Estimated total crossing attempts by component (monthly)",
       y = "Number of persons", x = NULL) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9))

ggsave(file.path(BASE_DIR, "output", "figures", "crossing_attempts_components.png"),
       p5, width = 12, height = 6, dpi = 200)
cat("Saved: output/figures/crossing_attempts_components.png\n")

# Figure 6: Detection method composition (bar plot)
frx_det <- frx_cmr %>%
  filter(!is.na(detected_by)) %>%
  mutate(
    ym = floor_date(date, "month"),
    det_method = case_when(
      grepl("NGO vessel", detected_by, ignore.case = TRUE) ~ "NGO vessel",
      grepl("FWA|RPAS|HELO", detected_by, ignore.case = TRUE) ~ "Aerial surveillance",
      grepl("CPV|CPB|Land Patrol|OPV", detected_by, ignore.case = TRUE) ~ "Patrol vessel",
      grepl("Marina|MAS|Mare Sicuro|Mare Nostrum|EUNAVFOR",
            detected_by, ignore.case = TRUE) ~ "Navy (Italian/EU)",
      grepl("Call-Migrant|Call-Civilian|Call-Pleasure", detected_by, ignore.case = TRUE) ~ "Distress call",
      grepl("Commercial|fishing|Merchant", detected_by, ignore.case = TRUE) ~ "Commercial vessel",
      TRUE ~ "Other"
    )
  )

# Order: most common at bottom
det_levels <- c("Patrol vessel", "Navy (Italian/EU)", "Aerial surveillance",
                "NGO vessel", "Commercial vessel", "Distress call", "Other")
frx_det$det_method <- factor(frx_det$det_method, levels = rev(det_levels))

det_monthly <- frx_det %>% count(ym, det_method)

det_colours <- c(
  "Patrol vessel"       = "#2166AC",
  "Navy (Italian/EU)"   = "#4393C3",
  "Aerial surveillance" = "#E69F00",
  "NGO vessel"          = "#4DAF4A",
  "Commercial vessel"   = "#A6761D",
  "Distress call"       = "#D32F2F",
  "Other"               = "#BDBDBD"
)

p6 <- ggplot(det_monthly, aes(x = ym, y = n, fill = det_method)) +
  geom_col(width = 25) +
  scale_fill_manual(values = det_colours, name = "Detection method") +
  geom_vline(xintercept = as.Date("2017-07-01"), linetype = "dashed",
             colour = "red", linewidth = 0.5) +
  annotate("text", x = as.Date("2017-07-01"), y = Inf, label = "MoU",
           vjust = 2, hjust = -0.2, colour = "red", size = 3) +
  labs(
    title = "Detection Method Composition — Central Mediterranean",
    subtitle = "Monthly Frontex detections by method (CMR departures)",
    x = NULL, y = "Detections"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey50"),
    axis.text = element_text(size = 8, colour = "grey50"),
    axis.title = element_text(size = 9, colour = "grey50")
  )

ggsave(file.path(BASE_DIR, "output", "figures", "frx_detection_composition_time.png"),
       p6, width = 12, height = 6, dpi = 300)
cat("Saved: output/figures/frx_detection_composition_time.png\n")

# ── 8. Save text summary ──────────────────────────────────
cat("\n--- 8. Saving summary ---\n")

sink(file.path(BASE_DIR, "output", "tables", "frontex_iom_crosscheck.txt"))
cat("FRONTEX / IOM / SAR CROSS-CHECK\n")
cat("================================\n\n")

cat("DEATHS BY YEAR AND SOURCE:\n")
print(yearly_deaths, n = 12)

cat("\nVOLUME BY YEAR AND SOURCE:\n")
print(yearly_volume, n = 12)

cat("\nFRONTEX BOAT TYPE (inflatable share by year):\n")
boat_yearly %>% filter(boat_cat == "Inflatable") %>%
  select(year, n, pct) %>% print(n = 12)

cat("\nNGO SAR VALIDATION (monthly):\n")
if (nrow(monthly_sar) > 0) {
  cat(sprintf("Cor(Frontex NGO%%, R-S NGO vessels): %.3f\n",
      cor(monthly_sar$frx_ngo_pct, monthly_sar$rs_ngo_vessels, use = "complete.obs")))
}
sink()
cat("Saved: output/tables/frontex_iom_crosscheck.txt\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
