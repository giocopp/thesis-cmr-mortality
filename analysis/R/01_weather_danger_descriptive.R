# 01_weather_danger_descriptive.R
# ================================
# Descriptive analysis: how does weather relate to crossing danger?
#
# Questions:
#   a0) How dangerous is being at sea during a storm?
#   a1) Did weather-danger change after the MoU?
#   a2) Did crossing composition shift toward/away from rough weather?
#
# Data:
#   - IOM MMP incidents (filtered to new CMR sea zone polygon)
#   - UNHCR daily arrivals to Italy
#   - ERA5 daily SWH (spatial mean over sea zone)
#   - Monthly interceptions (for monthly robustness)
#
# Geography: Camarena-style sea zone polygon (coastline-following)
# Aggregation: weekly (primary), monthly (robustness with interceptions)
#
# Output:
#   output/figures/weather_danger_descriptive.png
#   output/tables/weather_danger_descriptive.csv

library(ncdf4)
library(sf)
library(rnaturalearth)
library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)

BASE_DIR <- here::here()
ERA5_DIR <- file.path(BASE_DIR, "data", "raw", "era5")
MOU_DATE <- as.Date("2017-07-01")

cat("============================================================\n")
cat("WEATHER & CROSSING DANGER: DESCRIPTIVE ANALYSIS\n")
cat("============================================================\n\n")

# ============================================================
# 1. Define sea zone polygon (Camarena-style)
# ============================================================
cat("--- 1. Sea zone polygon ---\n")
outer_coords <- matrix(c(
  15.5, 31.0, 15.5, 37.0, 12.4, 37.8, 11.0, 37.1,
  9.0, 34.0, 9.0, 31.0, 15.5, 31.0
), ncol = 2, byrow = TRUE)
outer_poly <- st_sfc(st_polygon(list(outer_coords)), crs = 4326)
world <- ne_countries(scale = "medium", returnclass = "sf")
sea_zone <- st_difference(outer_poly, st_union(world))
cat("Sea zone built\n")

# ============================================================
# 2. Compute daily SWH over sea zone from ERA5
# ============================================================
cat("\n--- 2. Daily SWH from ERA5 ---\n")

# Build mask
nc0 <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
wave_lon <- ncvar_get(nc0, "longitude")
wave_lat <- ncvar_get(nc0, "latitude")
nc_close(nc0)

grid <- expand.grid(lon = wave_lon, lat = wave_lat)
grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
in_zone <- st_intersects(grid_sf, sea_zone, sparse = FALSE)[, 1]

nc_chk <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
swh_chk <- ncvar_get(nc_chk, "swh", start = c(1,1,1), count = c(-1,-1,1))
nc_close(nc_chk)
grid$is_ocean <- !is.na(as.vector(swh_chk)) & in_zone
mask <- matrix(grid$is_ocean, nrow = length(wave_lon), ncol = length(wave_lat))
cat(sprintf("Ocean cells in sea zone: %d\n", sum(mask)))

get_nc_dates <- function(nc) {
  tn <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  tv <- ncvar_get(nc, tn)
  tu <- ncatt_get(nc, tn, "units")$value
  if (grepl("seconds since 1970", tu)) as.Date(as.POSIXct(tv, origin="1970-01-01", tz="UTC"))
  else if (grepl("hours since 1900", tu)) as.Date(as.POSIXct("1900-01-01", tz="UTC") + tv*3600)
  else { ref <- sub("(hours|seconds|days) since ", "", tu)
         mult <- ifelse(grepl("hours",tu),3600,ifelse(grepl("days",tu),86400,1))
         as.Date(as.POSIXct(ref, tz="UTC") + tv*mult) }
}

all_wx <- list()
for (yr in 2014:2025) {
  f <- file.path(ERA5_DIR, paste0("era5_daily_cmr_wave_", yr, ".nc"))
  if (!file.exists(f)) { cat(yr, "MISSING\n"); next }
  nc <- nc_open(f); dates <- get_nc_dates(nc); swh_3d <- ncvar_get(nc, "swh"); nc_close(nc)
  swh_daily <- sapply(seq_along(dates), function(t) mean(swh_3d[,,t][mask], na.rm=TRUE))
  all_wx[[as.character(yr)]] <- data.table(date = dates, swh = swh_daily)
  cat(yr, ":", length(dates), "days\n")
}
weather <- rbindlist(all_wx)
cat("Total weather days:", nrow(weather), "\n")

# ============================================================
# 3. Load and filter IOM incidents to sea zone
# ============================================================
cat("\n--- 3. IOM incidents (sea zone filter) ---\n")

iom <- fread(file.path(BASE_DIR, "data", "processed", "iom_mmp_incidents_2014_2025_reg.csv"))
cmr <- iom[Route == "Central Mediterranean"]
cmr[, `:=`(lat = as.numeric(Latitude), lon = as.numeric(Longitude),
           date = as.Date(incident_date_clean),
           dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE))]
cmr <- cmr[!is.na(lat) & !is.na(lon) & !is.na(date)]

pts <- st_as_sf(cmr, coords = c("lon", "lat"), crs = 4326)
cmr$in_zone <- st_intersects(pts, sea_zone, sparse = FALSE)[, 1]
cmr_zone <- cmr[in_zone == TRUE]

cat(sprintf("CMR incidents total: %d\n", nrow(cmr)))
cat(sprintf("In sea zone: %d (%.1f%%)\n", nrow(cmr_zone), 100*nrow(cmr_zone)/nrow(cmr)))
cat(sprintf("Deaths in zone: %.0f\n", sum(cmr_zone$dead_missing)))

# Aggregate to daily
daily_iom <- cmr_zone[, .(n_incidents = .N, deaths = sum(dead_missing)), by = date]

# ============================================================
# 4. Load UNHCR arrivals
# ============================================================
cat("\n--- 4. UNHCR arrivals ---\n")
arr <- fread(file.path(BASE_DIR, "data", "raw", "unhcr", "unhcr_daily_arrivals_italy.csv"))
arr[, date := as.Date(data_date)]
arr <- arr[, .(date, arrivals = individuals)]
cat(sprintf("UNHCR arrivals: %s to %s\n", min(arr$date), max(arr$date)))

# ============================================================
# 5. Build weekly panel
# ============================================================
cat("\n--- 5. Weekly panel ---\n")

# Full daily grid
all_dates <- data.table(date = seq(as.Date("2014-01-01"), as.Date("2025-12-31"), by = "day"))
daily <- merge(all_dates, weather, by = "date", all.x = TRUE)
daily <- merge(daily, daily_iom, by = "date", all.x = TRUE)
daily <- merge(daily, arr, by = "date", all.x = TRUE)
daily[is.na(n_incidents), n_incidents := 0L]
daily[is.na(deaths), deaths := 0]

daily[, iso_week := paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))]

weekly <- daily[, .(
  n_incidents = sum(n_incidents),
  deaths = sum(deaths),
  arrivals = sum(arrivals, na.rm = TRUE),
  arr_days = sum(!is.na(arrivals)),
  swh = mean(swh, na.rm = TRUE),
  date_start = min(date),
  year = isoyear(min(date))
), by = iso_week]

# Filter: need arrivals data (>= 5 days in the week)
w <- weekly[arr_days >= 5 & year >= 2016 & year <= 2024]
w[, crossings := deaths + arrivals]
w[, fat_rate := fifelse(crossings > 0, deaths / crossings, NA_real_)]
w[, post_mou := factor(fifelse(date_start >= MOU_DATE, "Post-MoU", "Pre-MoU"),
                        levels = c("Pre-MoU", "Post-MoU"))]

# SWH terciles
brks <- quantile(w$swh, c(0, 1/3, 2/3, 1))
w[, sea := cut(swh, breaks = brks, labels = c("Calm", "Medium", "Rough"), include.lowest = TRUE)]

cat(sprintf("Weekly panel: %d weeks (%s to %s)\n", nrow(w),
    as.character(min(w$date_start)), as.character(max(w$date_start))))
cat(sprintf("  Pre-MoU: %d weeks, Post-MoU: %d weeks\n",
    sum(w$post_mou == "Pre-MoU"), sum(w$post_mou == "Post-MoU")))

# ============================================================
# 6. Answer the questions
# ============================================================

cat("\n============================================================\n")
cat("a0) HOW DANGEROUS IS BEING AT SEA DURING A STORM?\n")
cat("============================================================\n\n")
a0 <- w[crossings > 0, .(
  weeks = .N,
  aggregate_fat_rate = round(100 * sum(deaths) / sum(crossings), 2),
  total_deaths = sum(deaths),
  total_crossings = sum(crossings)
), by = sea]
print(a0)

cat("\n============================================================\n")
cat("a1) DID WEATHER-DANGER CHANGE AFTER THE MOU?\n")
cat("============================================================\n\n")
a1 <- w[crossings > 0, .(
  weeks = .N,
  fat_rate_pct = round(100 * sum(deaths) / sum(crossings), 2),
  total_deaths = sum(deaths),
  total_crossings = sum(crossings)
), by = .(post_mou, sea)]
print(a1[order(post_mou, sea)])

cat("\n============================================================\n")
cat("a2) CROSSING COMPOSITION: WHO IS AT SEA IN ROUGH WEATHER?\n")
cat("============================================================\n\n")
a2 <- w[, .(
  weeks = .N,
  mean_crossings_week = round(mean(crossings)),
  total_crossings = sum(crossings)
), by = .(post_mou, sea)]
a2[, share_pct := round(100 * total_crossings / sum(total_crossings), 1), by = post_mou]
print(a2[order(post_mou, sea)])

# ============================================================
# 7. Plots
# ============================================================
cat("\n--- Generating plots ---\n")

# a0+a1: Fatality rate by weather x period
plot_a1 <- w[crossings > 0, .(fat_rate = 100 * sum(deaths)/sum(crossings)), by = .(post_mou, sea)]
p1 <- ggplot(plot_a1, aes(x = sea, y = fat_rate, fill = post_mou)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Pre-MoU" = "#D4820E", "Post-MoU" = "#D32F2F")) +
  labs(title = "Fatality rate by sea conditions and period",
       subtitle = "Fatality rate = deaths / (deaths + arrivals)",
       y = "Fatality rate (%)", x = "Sea conditions (SWH tercile)", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

# a2: Crossing share by weather x period
p2 <- ggplot(a2, aes(x = sea, y = share_pct, fill = post_mou)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Pre-MoU" = "#D4820E", "Post-MoU" = "#D32F2F")) +
  labs(title = "Share of crossings by sea conditions",
       subtitle = "% of total crossings in each period",
       y = "Share of crossings (%)", x = "Sea conditions (SWH tercile)", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

# Fatality rate vs SWH, by period (LOESS)
p3 <- ggplot(w[crossings > 10], aes(x = swh, y = fat_rate * 100, colour = post_mou)) +
  geom_point(alpha = 0.25, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  scale_colour_manual(values = c("Pre-MoU" = "#D4820E", "Post-MoU" = "#D32F2F")) +
  labs(title = "Fatality rate vs. SWH by period",
       subtitle = "Each dot = one week (crossings > 10). LOESS fit.",
       x = "Weekly mean SWH (m)", y = "Fatality rate (%)", colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# Year-by-year fatality rate by SWH tercile
yearly_fr <- w[crossings > 0, .(fat_rate = 100 * sum(deaths)/sum(crossings)), by = .(year, sea)]
p4 <- ggplot(yearly_fr, aes(x = factor(year), y = fat_rate, fill = sea)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Calm" = "#2166AC", "Medium" = "#D4820E", "Rough" = "#B2182B")) +
  labs(title = "Fatality rate by sea conditions, year by year",
       y = "Fatality rate (%)", x = NULL, fill = "SWH tercile") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_out <- (p1 | p2) / (p3 | p4)
ggsave(file.path(BASE_DIR, "output", "figures", "weather_danger_descriptive.png"),
       p_out, width = 14, height = 11, dpi = 200)
cat("Saved: output/figures/weather_danger_descriptive.png\n")

# Save tables
fwrite(a1, file.path(BASE_DIR, "output", "tables", "weather_danger_by_period.csv"))
fwrite(yearly_fr, file.path(BASE_DIR, "output", "tables", "weather_danger_by_year.csv"))
cat("Saved: output/tables/weather_danger_by_period.csv\n")
cat("Saved: output/tables/weather_danger_by_year.csv\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
