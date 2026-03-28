# 02_weather_danger_analysis.R
# ============================
# Descriptive analysis: weather conditions and crossing danger.
#
# Questions:
#   Q1) How dangerous is being at sea during different weather conditions?
#   Q2) Is being at sea during rough conditions more dangerous post-MoU?
#   Q3) Have people at sea during rough conditions decreased or increased post-MoU?
#
# SWH measure: spatial mean over Camarena-style sea zone, avg lag 1-3
# Sea state cutoffs: Calm < 0.5m | Medium 0.5-1.25m | Rough > 1.25m
# Fatality rate: deaths / (deaths + arrivals)
# Temporal aggregation: daily and weekly (compared)
#
# Input:
#   data/raw/era5/era5_daily_cmr_wave_YYYY.nc
#   data/processed/iom_mmp_incidents_2014_2025_reg.csv
#   data/raw/unhcr/unhcr_daily_arrivals_italy.csv
#
# Output:
#   output/figures/q1_q2_q3_main.png
#   output/figures/q1_q2_q3_yearly.png
#   output/tables/fatality_rate_by_sea_state.csv

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ncdf4)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
ERA5_DIR <- file.path(BASE_DIR, "data", "raw", "era5")
MOU_DATE <- as.Date("2017-07-01")

cat("============================================================\n")
cat("WEATHER & CROSSING DANGER ANALYSIS\n")
cat("============================================================\n\n")

# ============================================================
# 1. Sea zone polygon (Camarena-style)
# ============================================================
cat("--- 1. Sea zone ---\n")

outer_coords <- matrix(c(
  15.5, 31.0, 15.5, 37.0, 12.4, 37.8, 11.0, 37.1,
  9.0, 34.0, 9.0, 31.0, 15.5, 31.0
), ncol = 2, byrow = TRUE)

sea_zone <- st_sfc(st_polygon(list(outer_coords)), crs = 4326) %>%
  st_difference(st_union(ne_countries(scale = "medium", returnclass = "sf")))

# ============================================================
# 2. Daily SWH from ERA5 (spatial mean over sea zone)
# ============================================================
cat("--- 2. ERA5 daily SWH ---\n")

# Build ocean mask
nc0 <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
wave_lon <- ncvar_get(nc0, "longitude")
wave_lat <- ncvar_get(nc0, "latitude")
nc_close(nc0)

grid_pts <- expand.grid(lon = wave_lon, lat = wave_lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
in_zone <- st_intersects(grid_pts, sea_zone, sparse = FALSE)[, 1]

nc_chk <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
swh_day1 <- ncvar_get(nc_chk, "swh", start = c(1, 1, 1), count = c(-1, -1, 1))
nc_close(nc_chk)
mask <- matrix(!is.na(as.vector(swh_day1)) & in_zone, nrow = length(wave_lon))
cat(sprintf("  Ocean cells in zone: %d\n", sum(mask)))

# Helper: extract dates from netcdf
get_nc_dates <- function(nc) {
  tn <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  tv <- ncvar_get(nc, tn)
  tu <- ncatt_get(nc, tn, "units")$value
  if (grepl("seconds since 1970", tu)) {
    as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC"))
  } else if (grepl("hours since 1900", tu)) {
    as.Date(as.POSIXct("1900-01-01", tz = "UTC") + tv * 3600)
  } else {
    ref <- sub("(hours|seconds|days) since ", "", tu)
    mult <- ifelse(grepl("hours", tu), 3600, ifelse(grepl("days", tu), 86400, 1))
    as.Date(as.POSIXct(ref, tz = "UTC") + tv * mult)
  }
}

# Extract daily spatial-mean SWH
weather <- map_dfr(2014:2025, function(yr) {
  f <- file.path(ERA5_DIR, paste0("era5_daily_cmr_wave_", yr, ".nc"))
  if (!file.exists(f)) return(tibble())
  nc <- nc_open(f)
  dates <- get_nc_dates(nc)
  swh_3d <- ncvar_get(nc, "swh")
  nc_close(nc)
  tibble(
    date = dates,
    swh  = map_dbl(seq_along(dates), ~ mean(swh_3d[, , .x][mask], na.rm = TRUE))
  )
}) %>%
  arrange(date) %>%
  mutate(swh_avg13 = (lag(swh, 1) + lag(swh, 2) + lag(swh, 3)) / 3)

cat(sprintf("  Weather days: %d\n", nrow(weather)))

# ============================================================
# 3. IOM incidents (filtered to sea zone)
# ============================================================
cat("--- 3. IOM incidents ---\n")

iom_raw <- read_csv(file.path(BASE_DIR, "data", "processed",
                                "iom_mmp_incidents_2014_2025_reg.csv"),
                     show_col_types = FALSE) %>%
  filter(Route == "Central Mediterranean") %>%
  mutate(
    lat = as.numeric(Latitude),
    lon = as.numeric(Longitude),
    date = as.Date(incident_date_clean),
    dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)
  ) %>%
  filter(!is.na(lat), !is.na(lon), !is.na(date))

pts <- st_as_sf(iom_raw, coords = c("lon", "lat"), crs = 4326)
iom_raw$in_zone <- st_intersects(pts, sea_zone, sparse = FALSE)[, 1]

daily_iom <- iom_raw %>%
  filter(in_zone) %>%
  group_by(date) %>%
  summarise(n_incidents = n(), deaths = sum(dead_missing), .groups = "drop")

cat(sprintf("  Incidents in zone: %d\n", sum(iom_raw$in_zone)))

# ============================================================
# 4. UNHCR arrivals
# ============================================================
cat("--- 4. UNHCR arrivals ---\n")

arrivals <- read_csv(file.path(BASE_DIR, "data", "raw", "unhcr",
                                 "unhcr_daily_arrivals_italy.csv"),
                      show_col_types = FALSE) %>%
  transmute(date = as.Date(data_date), arrivals = individuals)

cat(sprintf("  Arrivals: %s to %s\n", min(arrivals$date), max(arrivals$date)))

# ============================================================
# 5. Assemble daily panel
# ============================================================
cat("--- 5. Daily panel ---\n")

daily <- tibble(date = seq(as.Date("2014-01-01"), as.Date("2025-12-31"), by = "day")) %>%
  left_join(weather, by = "date") %>%
  left_join(daily_iom, by = "date") %>%
  left_join(arrivals, by = "date") %>%
  replace_na(list(n_incidents = 0L, deaths = 0)) %>%
  filter(!is.na(arrivals), !is.na(swh_avg13), year(date) <= 2024) %>%
  mutate(
    crossings = deaths + arrivals,
    sea       = cut(swh_avg13,
                    breaks = c(-Inf, 0.5, 1.25, Inf),
                    labels = c("Calm", "Medium", "Rough"),
                    right  = FALSE),
    post_mou  = if_else(date >= MOU_DATE, "Post-MoU", "Pre-MoU") %>%
                  factor(levels = c("Pre-MoU", "Post-MoU")),
    year      = year(date),
    iso_week  = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))
  )

cat(sprintf("  Days: %d | Sea state: %s\n", nrow(daily),
    paste(table(daily$sea), collapse = " / ")))

# ============================================================
# 6. Assemble weekly panel
# ============================================================
cat("--- 6. Weekly panel ---\n")

weekly <- daily %>%
  group_by(iso_week) %>%
  summarise(
    deaths   = sum(deaths),
    arrivals = sum(arrivals),
    n_incidents = sum(n_incidents),
    swh_avg13 = mean(swh_avg13),
    date_start = min(date),
    n_days   = n(),
    .groups  = "drop"
  ) %>%
  filter(n_days >= 5) %>%
  mutate(
    crossings = deaths + arrivals,
    sea       = cut(swh_avg13,
                    breaks = c(-Inf, 0.5, 1.25, Inf),
                    labels = c("Calm", "Medium", "Rough"),
                    right  = FALSE),
    post_mou  = if_else(date_start >= MOU_DATE, "Post-MoU", "Pre-MoU") %>%
                  factor(levels = c("Pre-MoU", "Post-MoU")),
    year      = isoyear(date_start)
  )

cat(sprintf("  Weeks: %d | Sea state: %s\n", nrow(weekly),
    paste(table(weekly$sea), collapse = " / ")))

# ============================================================
# 7. Compute tables
# ============================================================

# Helper: fatality rate summary
fr_summary <- function(df, ...) {
  df %>%
    filter(crossings > 0) %>%
    group_by(...) %>%
    summarise(
      obs        = n(),
      deaths     = sum(deaths),
      crossings  = sum(crossings),
      fat_rate   = round(100 * sum(deaths) / sum(crossings), 2),
      .groups    = "drop"
    )
}

# Q1: overall
q1_daily  <- fr_summary(daily, sea) %>% mutate(level = "Daily")
q1_weekly <- fr_summary(weekly, sea) %>% mutate(level = "Weekly")
q1 <- bind_rows(q1_daily, q1_weekly)

# Q2: by period
q2_daily  <- fr_summary(daily, post_mou, sea) %>% mutate(level = "Daily")
q2_weekly <- fr_summary(weekly, post_mou, sea) %>% mutate(level = "Weekly")
q2 <- bind_rows(q2_daily, q2_weekly)

# Q3: crossing share
share_fn <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(total_crossings = sum(crossings), obs = n(), .groups = "drop") %>%
    group_by(post_mou) %>%
    mutate(share_pct = round(100 * total_crossings / sum(total_crossings), 1)) %>%
    ungroup()
}
q3_daily  <- share_fn(daily, post_mou, sea) %>% mutate(level = "Daily")
q3_weekly <- share_fn(weekly, post_mou, sea) %>% mutate(level = "Weekly")
q3 <- bind_rows(q3_daily, q3_weekly)

# Year-by-year
yr_fr_daily  <- fr_summary(daily, year, sea) %>% mutate(level = "Daily")
yr_fr_weekly <- fr_summary(weekly, year, sea) %>% mutate(level = "Weekly")
yr_fr <- bind_rows(yr_fr_daily, yr_fr_weekly)

yr_share_daily <- daily %>%
  group_by(year, sea) %>%
  summarise(total = sum(crossings), .groups = "drop") %>%
  group_by(year) %>%
  mutate(share = round(100 * total / sum(total), 1), level = "Daily") %>%
  ungroup()
yr_share_weekly <- weekly %>%
  group_by(year, sea) %>%
  summarise(total = sum(crossings), .groups = "drop") %>%
  group_by(year) %>%
  mutate(share = round(100 * total / sum(total), 1), level = "Weekly") %>%
  ungroup()
yr_share <- bind_rows(yr_share_daily, yr_share_weekly)

# ============================================================
# 8. Print results
# ============================================================
cat("\n============================================================\n")
cat("Q1: FATALITY RATE BY SEA STATE (overall)\n")
cat("============================================================\n")
print(q1 %>% select(level, sea, obs, deaths, crossings, fat_rate))

cat("\n============================================================\n")
cat("Q2: FATALITY RATE BY SEA STATE x PERIOD\n")
cat("============================================================\n")
print(q2 %>% select(level, post_mou, sea, obs, deaths, crossings, fat_rate))

cat("\n============================================================\n")
cat("Q3: CROSSING SHARE BY SEA STATE x PERIOD\n")
cat("============================================================\n")
print(q3 %>% select(level, post_mou, sea, total_crossings, share_pct, obs))

# Save
write_csv(q2, file.path(BASE_DIR, "output", "tables", "fatality_rate_by_sea_state.csv"))
cat("\nSaved: output/tables/fatality_rate_by_sea_state.csv\n")

# ============================================================
# 9. Plots — main (Q1, Q2, Q3)
# ============================================================
cat("\n--- Generating plots ---\n")

fill_period <- scale_fill_manual(values = c("Pre-MoU" = "#D4820E", "Post-MoU" = "#D32F2F"))
fill_sea    <- scale_fill_manual(values = c("Calm" = "#2166AC", "Medium" = "#D4820E", "Rough" = "#B2182B"))

p_q1 <- ggplot(q1, aes(x = sea, y = fat_rate, fill = sea)) +
  geom_col(alpha = 0.8, width = 0.6) +
  facet_wrap(~level) +
  fill_sea + guides(fill = "none") +
  labs(title = "Q1: Fatality rate by sea state",
       subtitle = "SWH = avg lag 1-3 | Calm < 0.5m | Medium 0.5-1.25m | Rough > 1.25m",
       y = "Fatality rate (%)", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_q2 <- ggplot(q2, aes(x = sea, y = fat_rate, fill = post_mou)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  facet_wrap(~level) +
  fill_period +
  labs(title = "Q2: Fatality rate by sea state and period",
       y = "Fatality rate (%)", x = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_q3 <- ggplot(q3, aes(x = sea, y = share_pct, fill = post_mou)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  facet_wrap(~level) +
  fill_period +
  labs(title = "Q3: Share of crossings by sea state and period",
       y = "Share of crossings (%)", x = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_main <- p_q1 / p_q2 / p_q3
ggsave(file.path(BASE_DIR, "output", "figures", "q1_q2_q3_main.png"),
       p_main, width = 12, height = 12, dpi = 200)

# ============================================================
# 10. Plots — year by year
# ============================================================

p_yr_fr <- ggplot(yr_fr, aes(x = factor(year), y = fat_rate, fill = sea)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  facet_wrap(~level) +
  fill_sea +
  labs(title = "Fatality rate by sea state, year by year",
       subtitle = "SWH = avg lag 1-3 | Calm < 0.5m | Medium 0.5-1.25m | Rough > 1.25m",
       y = "Fatality rate (%)", x = NULL, fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_yr_share <- ggplot(yr_share, aes(x = factor(year), y = share, fill = sea)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  facet_wrap(~level) +
  fill_sea +
  labs(title = "Share of crossings by sea state, year by year",
       y = "Share (%)", x = NULL, fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_yearly <- p_yr_fr / p_yr_share
ggsave(file.path(BASE_DIR, "output", "figures", "q1_q2_q3_yearly.png"),
       p_yearly, width = 12, height = 9, dpi = 200)

cat("Saved: output/figures/q1_q2_q3_main.png\n")
cat("Saved: output/figures/q1_q2_q3_yearly.png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
