# Reproducible EDA plots for Central Mediterranean (IOM TS + Missing Migrants)
# Uses:
#   - med_crossings_monthlyTS.csv
#   - iom_mmp_incidents_2014_2025_reg.csv
#
# Toggle:
#   - impute_missing_to_zero = TRUE  -> NA -> 0 (Rodriguez/Sanchez style)
#   - impute_missing_to_zero = FALSE -> NA stays NA (strict missingness)

library(tidyverse)
library(lubridate)
library(zoo)
library(sf)
library(ggspatial)

# ----------------------------
# Paths (edit these as needed)
# ----------------------------
data_dir <- "/Users/giocopp/Desktop/Uni/Hertie School/6th Semester/Thesis-MDS/IOM Data/Clean"
out_dir  <- "/Users/giocopp/Desktop/Uni/Hertie School/6th Semester/Thesis-MDS/IOM Data/EDA"

ts_file  <- file.path(data_dir, "med_crossings_monthlyTS.csv")
inc_file <- file.path(data_dir, "iom_mmp_incidents_2014_2025_reg.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Switch: NA -> 0 imputation
# ----------------------------
impute_missing_to_zero <- TRUE  # TRUE = Rodriguez/Sanchez style; FALSE = strict missingness

# Your treatment cutoff (used in plots/markers if desired)
treatment_cutoff <- as.Date("2017-05-01")
mou_date <- as.Date("2017-02-01")

# Helper: sum but keep NA if all values are NA
sum_na <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}

# ----------------------------
# 1) Monthly crossings time series (IOM-compiled)
# ----------------------------
ts <- readr::read_csv(ts_file, show_col_types = FALSE)

# Guard: ensure interceptions columns exist
req_int_cols <- c("interceptions_by_libyan_coast_guard", "interceptions_by_tunisian_coast_guard")
missing_int <- setdiff(req_int_cols, names(ts))
if (length(missing_int) > 0) {
  stop("Missing interceptions column(s) in TS CSV: ", paste(missing_int, collapse = ", "))
}

ts <- ts %>%
  mutate(
    date = as.Date(date),
    
    # Make sure numeric (handles if they came in as character)
    sea_arrivals_in_italy = suppressWarnings(as.numeric(sea_arrivals_in_italy)),
    sea_arrivals_in_malta = suppressWarnings(as.numeric(sea_arrivals_in_malta)),
    cmr                   = suppressWarnings(as.numeric(cmr)),
    interceptions_by_libyan_coast_guard   = suppressWarnings(as.numeric(interceptions_by_libyan_coast_guard)),
    interceptions_by_tunisian_coast_guard = suppressWarnings(as.numeric(interceptions_by_tunisian_coast_guard))
  )

if (impute_missing_to_zero) {
  # ---- NA -> 0 style (produces attempts/mortality even before interception series starts) ----
  ts <- ts %>%
    mutate(
      # arrivals: sum available; if both missing -> 0 (not NA)
      arrivals_cmr = rowSums(across(c(sea_arrivals_in_italy, sea_arrivals_in_malta)), na.rm = TRUE),
      
      # pushbacks/interceptions: NA -> 0
      lcg = replace_na(interceptions_by_libyan_coast_guard, 0),
      tcg = replace_na(interceptions_by_tunisian_coast_guard, 0),
      interceptions = lcg + tcg,
      
      # deaths: NA -> 0
      deaths_cmr = replace_na(cmr, 0),
      
      # attempts and mortality always defined (except attempts==0 -> NA mortality)
      attempts_cmr = arrivals_cmr + interceptions + deaths_cmr,
      mortality_per1000 = if_else(
        attempts_cmr > 0,
        (deaths_cmr / attempts_cmr) * 1000,
        NA_real_
      )
    )
} else {
  # ---- Strict missingness style (your original principle) ----
  ts <- ts %>%
    mutate(
      arrivals_cmr = dplyr::if_else(
        is.na(sea_arrivals_in_italy) & is.na(sea_arrivals_in_malta),
        NA_real_,
        rowSums(across(c(sea_arrivals_in_italy, sea_arrivals_in_malta)), na.rm = TRUE)
      ),
      
      interceptions = dplyr::if_else(
        is.na(interceptions_by_libyan_coast_guard) & is.na(interceptions_by_tunisian_coast_guard),
        NA_real_,
        rowSums(across(c(interceptions_by_libyan_coast_guard, interceptions_by_tunisian_coast_guard)), na.rm = TRUE)
      ),
      
      deaths_cmr = cmr,
      
      attempts_cmr = arrivals_cmr + interceptions + deaths_cmr,
      mortality_per1000 = if_else(
        attempts_cmr > 0,
        (deaths_cmr / attempts_cmr) * 1000,
        NA_real_
      )
    )
}

# Plot A: components of attempted crossings (2016+ when interceptions series typically begins)
p_attempt_components <- ts %>%
  filter(date >= as.Date("2016-01-01")) %>%
  transmute(
    date,
    `Arrivals (Italy+Malta)`        = arrivals_cmr,
    `Interceptions (Libya+Tunisia)` = interceptions,
    `Recorded deaths (CMR)`         = deaths_cmr
  ) %>%
  pivot_longer(-date, names_to = "component", values_to = "count") %>%
  ggplot(aes(date, count, fill = component)) +
  geom_area(alpha = 0.8) +
  labs(
    title = "Central Mediterranean: components of 'attempted crossings' (monthly, 2016–2025)",
    x = NULL, y = "People", fill = NULL
  ) +
  theme_minimal()

# Plot B: mortality per 1,000 attempts + policy markers (WITH labels)
policy_dates <- tibble::tribble(
  ~label, ~d,
  "Mare Nostrum start",           as.Date("2013-10-01"),
  "Mare Nostrum ends",            as.Date("2014-10-31"),
  "Italy–Libya MoU",              mou_date,
  "NGO code-of-conduct",          as.Date("2017-07-01"),
  "Piantedosi decree",            as.Date("2023-01-02"),
  "Flussi decree",                as.Date("2024-12-01")
)

ts_m <- ts %>%
  arrange(date) %>%
  mutate(
    mortality_ma3 = zoo::rollmean(mortality_per1000, 3, fill = NA, align = "right")
  )

y_max <- suppressWarnings(max(ts_m$mortality_per1000, na.rm = TRUE))
y_lab <- if (is.finite(y_max)) y_max * 0.7 else 0
policy_lab <- policy_dates %>% mutate(y = y_lab)

p_mortality <- ggplot(ts_m, aes(date)) +
  geom_line(aes(y = mortality_per1000), alpha = 0.35) +
  geom_line(aes(y = mortality_ma3), linewidth = 0.9) +
  geom_vline(data = policy_dates, aes(xintercept = d), linetype = "dashed") +
  geom_text(
    data = policy_lab,
    aes(x = d, y = y, label = label),
    angle = 90, vjust = 1.2, hjust = 0, size = 3
  ) +
  labs(
    title = "Central Mediterranean mortality rate (deaths per 1,000 attempted crossings)",
    x = NULL, y = "Deaths per 1,000 attempts"
  ) +
  theme_minimal()

# ----------------------------
# 2) Missing Migrants incidents (create `inc`, `inc_year`, and `p_compare`)
# ----------------------------
inc_raw <- readr::read_csv(inc_file, show_col_types = FALSE)

# Column-name compatibility guards
date_col <- dplyr::case_when(
  "Incident date" %in% names(inc_raw) ~ "Incident date",
  "incident_date" %in% names(inc_raw) ~ "incident_date",
  TRUE ~ NA_character_
)

dm_col <- dplyr::case_when(
  "No. dead/missing" %in% names(inc_raw) ~ "No. dead/missing",
  "No. dead/ missing" %in% names(inc_raw) ~ "No. dead/ missing",
  "No. dead / missing" %in% names(inc_raw) ~ "No. dead / missing",
  TRUE ~ NA_character_
)

if (is.na(date_col)) stop("Could not find an incident date column in incidents CSV.")
if (is.na(dm_col))   stop("Could not find a deaths+missing column in incidents CSV.")

inc <- inc_raw %>%
  rename(incident_date = !!sym(date_col)) %>%
  mutate(
    incident_date = as.Date(incident_date),
    deaths_missing = readr::parse_number(as.character(.data[[dm_col]])),
    Latitude  = suppressWarnings(as.numeric(Latitude)),
    Longitude = suppressWarnings(as.numeric(Longitude)),
    `Cause of death (category)` = na_if(trimws(as.character(`Cause of death (category)`)), "")
  ) %>%
  filter(Route == "Central Mediterranean") %>%
  filter(
    incident_date >= as.Date("2014-01-01"),
    incident_date <= as.Date("2025-12-31")
  )

inc_year <- inc %>%
  filter(!is.na(incident_date)) %>%
  mutate(year = lubridate::year(incident_date)) %>%
  group_by(year) %>%
  summarise(
    deaths_missing = sum_na(deaths_missing),  # NA if all missing
    incidents      = n(),
    .groups = "drop"
  ) %>%
  arrange(year)

# ----------------------------
# Plot C: yearly deaths (bars) and incident counts (line; scaled)
# ----------------------------
max_dm  <- suppressWarnings(max(inc_year$deaths_missing, na.rm = TRUE))
max_inc <- suppressWarnings(max(inc_year$incidents, na.rm = TRUE))

scale_fac <- if (is.finite(max_dm) && is.finite(max_inc) && max_inc > 0) max_dm / max_inc else 1
if (!is.finite(scale_fac) || scale_fac == 0) scale_fac <- 1

p_yearly <- ggplot(inc_year, aes(x = year)) +
  geom_col(aes(y = deaths_missing)) +
  geom_line(aes(y = incidents * scale_fac), linewidth = 0.8) +
  geom_point(aes(y = incidents * scale_fac), size = 1.8) +
  scale_y_continuous(
    name = "Deaths + Missing",
    sec.axis = sec_axis(~ . / scale_fac, name = "Number of recorded incidents")
  ) +
  labs(
    title = "Central Mediterranean incidents (IOM MMP): deaths+missing (bars) vs incident counts (line)",
    x = "Year"
  ) +
  theme_minimal()

# ----------------------------
# Plot D: top cause-of-death categories (by number of incidents)
# ----------------------------
p_causes <- inc %>%
  mutate(cause_cat = replace_na(`Cause of death (category)`, "Unknown/unspecified")) %>%
  count(cause_cat, sort = TRUE) %>%
  slice_head(n = 10) %>%
  ggplot(aes(x = reorder(cause_cat, n), y = n)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Central Mediterranean incidents: top reported cause-of-death categories",
    x = NULL, y = "Number of incidents"
  ) +
  theme_minimal()

# ----------------------------
# Plot E: spatial scatter pre/post treatment cutoff (OSM tiles if available)
# ----------------------------
inc_pts <- inc %>%
  filter(!is.na(Latitude), !is.na(Longitude), !is.na(incident_date)) %>%
  mutate(
    period = if_else(incident_date <= treatment_cutoff, "Pre treatment", "Post treatment")
  )

inc_sf <- st_as_sf(
  inc_pts, coords = c("Longitude", "Latitude"),
  crs = 4326, remove = FALSE
) %>%
  st_transform(3857)  # Web Mercator for tiles

bbox_ll <- st_bbox(c(
  xmin = 5.0,
  xmax = 23.5,
  ymin = 30.0,
  ymax = 41.8
), crs = st_crs(4326))

bbox_3857 <- st_as_sfc(bbox_ll) %>% st_transform(3857) %>% st_bbox()

p_spatial_osm <- tryCatch({
  ggplot() +
    annotation_map_tile(type = "osm", zoom = 7) +
    geom_sf(data = inc_sf, aes(color = period), alpha = 0.45, size = 0.9) +
    coord_sf(
      xlim = c(bbox_3857["xmin"], bbox_3857["xmax"]),
      ylim = c(bbox_3857["ymin"], bbox_3857["ymax"]),
      expand = FALSE
    ) +
    labs(title = "Central Mediterranean deadly incidents (pre/post May 2017)", color = NULL) +
    theme_minimal()
}, error = function(e) {
  message("Tile fetch failed (no internet or OSM blocked). Falling back to plain point map.")
  ggplot(inc_pts, aes(Longitude, Latitude, color = period)) +
    geom_point(alpha = 0.45, size = 0.9) +
    coord_cartesian(xlim = c(5.0, 23.5), ylim = c(30.0, 41.8), expand = FALSE) +
    labs(title = "Central Mediterranean deadly incidents (no tiles)", color = NULL) +
    theme_minimal()
})

# ----------------------------
# Plot F: definition check / compare monthly fatalities (TS vs incidents)
# ----------------------------
inc_m <- inc %>%
  filter(!is.na(incident_date)) %>%
  mutate(month = floor_date(incident_date, "month")) %>%
  group_by(month) %>%
  summarise(deaths_missing = sum_na(deaths_missing), .groups = "drop")

compare <- ts %>%
  select(date, deaths_cmr) %>%
  left_join(inc_m, by = c("date" = "month"))

p_compare <- compare %>%
  pivot_longer(c(deaths_cmr, deaths_missing), names_to = "series", values_to = "count") %>%
  mutate(series = recode(series,
                         deaths_cmr = "TS deaths (CMR)",
                         deaths_missing = "Incidents (dead+missing)"
  )) %>%
  ggplot(aes(date, count, color = series)) +
  geom_line() +
  labs(
    title = "Definition check: fatalities in TS vs Missing Migrants incidents (monthly)",
    x = NULL, y = "People", color = NULL
  ) +
  theme_minimal()

# ---- Save plots ----
ggsave(file.path(out_dir, "plot_ts_attempt_components.png"),
       p_attempt_components, width = 11, height = 4.2, dpi = 200)
ggsave(file.path(out_dir, "plot_ts_mortality_policy.png"),
       p_mortality, width = 11, height = 4.2, dpi = 200)
ggsave(file.path(out_dir, "plot_inc_yearly.png"),
       p_yearly, width = 11, height = 4.2, dpi = 200)
ggsave(file.path(out_dir, "plot_inc_causes.png"),
       p_causes, width = 11, height = 4.2, dpi = 200)
ggsave(file.path(out_dir, "plot_inc_spatial_osm.png"),
       p_spatial_osm, width = 7.2, height = 6.2, dpi = 200)
ggsave(file.path(out_dir, "plot_compare.png"),
       p_compare, width = 11, height = 4.2, dpi = 200)

message("Done. Plots saved to: ", out_dir)
