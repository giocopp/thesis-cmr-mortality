# 00b_build_daily_panel.R
# =======================
# Build the daily panel dataset used by analysis scripts.
# Merges clean datasets from data/processed/ into a single daily panel.
#
# Sources (all from data/processed/):
#   - era5_swh_daily.RDS         (ERA5 daily SWH, spatial mean over core corridor)
#   - iom_mmp_incidents.RDS      (IOM MMP incident-level, 2014+)
#   - iom_med_crossings_monthly.RDS (monthly interceptions + official arrivals)
#   - unhcr_daily_arrivals.RDS   (UNHCR daily arrivals to Italy)
#
# Optional:
#   - data/processed/archive/migrant_files_cmr_pre_iom.RDS (2008-2013, pre-IOM)
#
# Deaths merge strategy:
#   2008-2013: Migrant Files (researcher-compiled, only source available)
#   2014+:     IOM MMP (systematic institutional monitoring)
#   No overlap — IOM takes precedence for all years it covers.
#
# Output:
#   analysis/data/daily_panel.RDS

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()

cat("============================================================\n")
cat("BUILD DAILY PANEL\n")
cat("============================================================\n\n")

# ============================================================
# 1. ERA5 daily SWH
# ============================================================
cat("--- 1. ERA5 daily SWH ---\n")

weather <- readRDS(file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS"))
cat(sprintf("  %d weather days (%s to %s)\n",
    nrow(weather), min(weather$date), max(weather$date)))

# ============================================================
# 2. IOM deaths (CMR, all incident types except sub-incidents)
# ============================================================
cat("--- 2. IOM deaths ---\n")

iom <- readRDS(file.path(BASE_DIR, "data", "processed", "iom_mmp_incidents.RDS")) %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) != "sub-incident") %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date))

daily_iom <- iom %>%
  group_by(date) %>%
  summarise(deaths = sum(dead_missing),
            n_incidents = n(),
            .groups = "drop")

cat(sprintf("  %d incidents, %.0f deaths\n", nrow(iom), sum(iom$dead_missing)))

# ============================================================
# 2b. Migrant Files deaths (2008-2013, pre-IOM period)
# ============================================================
cat("--- 2b. Migrant Files deaths (2008-2013) ---\n")

mf_path <- file.path(BASE_DIR, "data", "processed", "archive",
                      "migrant_files_cmr_pre_iom.RDS")
if (file.exists(mf_path)) {
  mf <- readRDS(mf_path)
  daily_mf <- mf %>%
    group_by(date) %>%
    summarise(deaths = sum(dead_missing),
              n_incidents = n(),
              .groups = "drop")
  cat(sprintf("  %d incidents, %.0f deaths (2008-2013)\n",
      nrow(mf), sum(mf$dead_missing)))
} else {
  daily_mf <- tibble(date = as.Date(character()),
                     deaths = numeric(), n_incidents = integer())
  cat("  migrant_files_cmr_pre_iom.RDS not found — skipping pre-IOM period\n")
}

# Combine: Migrant Files for 2008-2013, IOM for 2014+
daily_deaths_combined <- bind_rows(
  daily_mf %>% filter(date < as.Date("2014-01-01")),
  daily_iom
)
cat(sprintf("  Combined: %d days with incidents\n", nrow(daily_deaths_combined)))

# ============================================================
# 3. UNHCR daily arrivals
# ============================================================
cat("--- 3. UNHCR arrivals ---\n")

arrivals <- readRDS(file.path(BASE_DIR, "data", "processed",
                              "unhcr_daily_arrivals.RDS"))
cat(sprintf("  %s to %s\n", min(arrivals$date), max(arrivals$date)))

# ============================================================
# 4. Monthly interceptions and official arrivals
# ============================================================
cat("--- 4. Monthly interceptions ---\n")

monthly_official <- readRDS(file.path(BASE_DIR, "data", "processed",
                                       "iom_med_crossings_monthly.RDS")) %>%
  transmute(
    ym = as.Date(date),
    official_arrivals = as.numeric(sea_arrivals_in_italy),
    interceptions = replace_na(as.numeric(interceptions_by_libyan_coast_guard), 0) +
                    replace_na(as.numeric(interceptions_by_tunisian_coast_guard), 0),
    intercept_per_day = interceptions / days_in_month(as.Date(date))
  )

cat(sprintf("  %d months\n", nrow(monthly_official)))

# ============================================================
# 5. Assemble daily panel
# ============================================================
cat("--- 5. Assembling daily panel ---\n")

daily <- tibble(date = seq(as.Date("2008-01-01"),
                           as.Date("2025-12-31"), by = "day")) %>%
  left_join(weather, by = "date") %>%
  left_join(daily_deaths_combined, by = "date") %>%
  left_join(arrivals, by = "date") %>%
  replace_na(list(deaths = 0, n_incidents = 0L)) %>%
  mutate(ym = floor_date(date, "month"),
         iso_week = paste0(isoyear(date), "_w",
                           sprintf("%02d", isoweek(date)))) %>%
  left_join(monthly_official, by = "ym") %>%
  replace_na(list(intercept_per_day = 0, interceptions = 0)) %>%
  mutate(crossings       = deaths + arrivals + intercept_per_day,
         crossings_no_ic = deaths + arrivals)

cat(sprintf("  %d days\n", nrow(daily)))
cat(sprintf("  Arrivals available: %d days\n", sum(!is.na(daily$arrivals))))
cat(sprintf("  SWH available: %d days\n", sum(!is.na(daily$swh))))

# ============================================================
# 6. Save
# ============================================================
out_dir <- file.path(BASE_DIR, "analysis", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(daily, file.path(out_dir, "daily_panel.RDS"))
cat("\nSaved: analysis/data/daily_panel.RDS\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
