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
CORE <- list(lon_min = 10.0, lon_max = 15.1, lat_min = 32.4, lat_max = 37.8)
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

# Aggregate Frontex to daily
frx_daily <- frx %>%
  group_by(date) %>%
  summarise(
    frx_incidents  = n(),
    frx_persons    = sum(num_persons, na.rm = TRUE),
    frx_deaths     = sum(num_deaths, na.rm = TRUE),
    frx_migrants   = sum(num_migrants, na.rm = TRUE),
    frx_pct_sar    = mean(sar_flag == 1, na.rm = TRUE) * 100,
    frx_n_sar      = sum(sar_flag == 1, na.rm = TRUE),
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
         in_core = Longitude >= CORE$lon_min & Longitude <= CORE$lon_max &
                   Latitude >= CORE$lat_min & Latitude <= CORE$lat_max) %>%
  filter(!is.na(date))

iom_daily <- iom_cmr %>%
  group_by(date) %>%
  summarise(
    iom_incidents    = n(),
    iom_dead_missing = sum(dead_missing),
    iom_core_deaths  = sum(dead_missing * in_core),
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

combined <- tibble(date = seq(as.Date("2014-01-01"), as.Date("2023-12-31"), by = "day")) %>%
  left_join(frx_daily, by = "date") %>%
  left_join(iom_daily, by = "date") %>%
  left_join(panel, by = "date") %>%
  left_join(sar, by = "date") %>%
  replace_na(list(
    frx_incidents = 0, frx_persons = 0, frx_deaths = 0, frx_migrants = 0,
    frx_pct_sar = 0, frx_n_sar = 0, frx_n_inflatable = 0,
    iom_incidents = 0, iom_dead_missing = 0, iom_core_deaths = 0,
    deaths = 0, n_incidents = 0
  )) %>%
  mutate(year = year(date), month = month(date))

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

# ── 5. SAR validation: Frontex SAR flag vs R-S NGO count ──
cat("\n--- 5. SAR validation: Frontex vs Rodriguez-Sanchez ---\n")

# Weekly aggregation for comparison
weekly_sar <- combined %>%
  filter(year >= 2014, year <= 2021) %>%
  mutate(week = floor_date(date, "week")) %>%
  group_by(week) %>%
  summarise(
    frx_sar_incidents = sum(frx_n_sar),
    frx_total_incidents = sum(frx_incidents),
    frx_sar_pct = ifelse(sum(frx_incidents) > 0,
                          sum(frx_n_sar) / sum(frx_incidents) * 100, NA),
    rs_ngo_vessels = mean(n_ngo_vessels, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(frx_sar_pct), frx_total_incidents > 0)

cat(sprintf("  Weekly obs with Frontex data: %d\n", nrow(weekly_sar)))
if (nrow(weekly_sar) > 10) {
  cat(sprintf("  Cor(frx_sar_pct, rs_ngo_vessels): %.3f\n",
      cor(weekly_sar$frx_sar_pct, weekly_sar$rs_ngo_vessels, use = "complete.obs")))
}

# Monthly for a cleaner view
monthly_sar <- combined %>%
  filter(year >= 2014, year <= 2021) %>%
  mutate(ym = floor_date(date, "month")) %>%
  group_by(ym) %>%
  summarise(
    frx_sar_pct = ifelse(sum(frx_incidents) > 0,
                          sum(frx_n_sar) / sum(frx_incidents) * 100, NA),
    rs_ngo_vessels = mean(n_ngo_vessels, na.rm = TRUE),
    frx_incidents = sum(frx_incidents),
    .groups = "drop"
  ) %>%
  filter(!is.na(frx_sar_pct), frx_incidents > 0)

cat(sprintf("  Monthly obs: %d\n", nrow(monthly_sar)))
if (nrow(monthly_sar) > 10) {
  cat(sprintf("  Cor(frx_sar_pct, rs_ngo_vessels) monthly: %.3f\n",
      cor(monthly_sar$frx_sar_pct, monthly_sar$rs_ngo_vessels, use = "complete.obs")))
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
    frx_deaths = sum(frx_deaths),
    iom_deaths = sum(iom_dead_missing),
    frx_persons = sum(frx_persons),
    panel_arrivals = sum(arrivals, na.rm = TRUE),
    frx_sar_pct = ifelse(sum(frx_incidents) > 0,
                          sum(frx_n_sar) / sum(frx_incidents) * 100, NA),
    rs_ngo = mean(n_ngo_vessels, na.rm = TRUE),
    .groups = "drop"
  )

p1 <- ggplot(monthly, aes(x = ym)) +
  geom_line(aes(y = frx_deaths, colour = "Frontex"), linewidth = 0.6) +
  geom_line(aes(y = iom_deaths, colour = "IOM MMP"), linewidth = 0.6) +
  scale_colour_manual(values = c("Frontex" = "#D32F2F", "IOM MMP" = "#1565C0")) +
  labs(title = "Monthly deaths: Frontex vs IOM MMP",
       y = "Deaths", x = NULL, colour = "Source") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

p2 <- ggplot(monthly, aes(x = ym)) +
  geom_line(aes(y = frx_persons, colour = "Frontex persons"), linewidth = 0.6) +
  geom_line(aes(y = panel_arrivals, colour = "Panel arrivals"), linewidth = 0.6) +
  scale_colour_manual(values = c("Frontex persons" = "#D32F2F", "Panel arrivals" = "#1565C0")) +
  labs(title = "Monthly volume: Frontex total persons vs daily panel arrivals",
       y = "Persons", x = NULL, colour = "Source") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

p3 <- monthly %>%
  filter(!is.na(frx_sar_pct), !is.na(rs_ngo)) %>%
  ggplot(aes(x = ym)) +
  geom_line(aes(y = frx_sar_pct, colour = "Frontex SAR %"), linewidth = 0.6) +
  geom_line(aes(y = rs_ngo * 10, colour = "R-S NGO vessels (x10)"), linewidth = 0.6) +
  scale_colour_manual(values = c("Frontex SAR %" = "#D32F2F", "R-S NGO vessels (x10)" = "#2E7D32")) +
  scale_y_continuous(
    name = "Frontex SAR %",
    sec.axis = sec_axis(~ . / 10, name = "NGO vessels (R-S)")
  ) +
  labs(title = "SAR validation: Frontex SAR flag vs Rodriguez-Sanchez NGO count",
       x = NULL, colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

p_combined <- p1 / p2 / p3
ggsave(file.path(BASE_DIR, "output", "figures", "frontex_iom_crosscheck.png"),
       p_combined, width = 11, height = 12, dpi = 200)
cat("Saved: output/figures/frontex_iom_crosscheck.png\n")

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

cat("\nSAR VALIDATION (monthly):\n")
if (nrow(monthly_sar) > 0) {
  cat(sprintf("Cor(Frontex SAR%%, R-S NGO vessels): %.3f\n",
      cor(monthly_sar$frx_sar_pct, monthly_sar$rs_ngo_vessels, use = "complete.obs")))
}
sink()
cat("Saved: output/tables/frontex_iom_crosscheck.txt\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
