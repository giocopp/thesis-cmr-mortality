# Merge Triton 2014-2017 coordinates into the pad-194 Themis panel
# via group-match on (date, transport, detected_by, intercept_by, sar).

library(dplyr)
library(readxl)
library(tidyr)

BASE_DIR <- here::here()

# ── 1. Load raw ─────────────────────────────────────────────────────────────
tri_raw <- read_excel(
  file.path(BASE_DIR, "data", "raw", "frontex",
            "Triton 2014 2017 Incident Details.xlsx"),
  sheet = 1
)

pad_raw <- read_excel(
  file.path(BASE_DIR, "data", "raw", "frontex",
            "pad-194_themis_2014_2023.xlsx"),
  sheet = "Sheet1"
)

# ── 2. Standardise + fix Triton coord-sum bug ───────────────────────────────
# Triton aggregated rows store SUMS of lat/lon over n_incidents; divide back.
tri <- tri_raw |>
  transmute(
    tri_row      = row_number(),
    date         = as.Date(`Detection date`),
    op_name      = `Operation name`,
    sar          = `Search and rescue involved`,
    n_incidents  = `Incidents on irregular migration`,
    n_migrants   = `Total number of irregular migrants`,
    n_deaths     = `Death cases`,
    transport    = `Transport type`,
    detected_by  = `Type of detected by`,
    intercept_by = `Type of intercepted by`,
    det_lat_raw  = `Detection latitude`,
    det_lon_raw  = `Detection longitude`,
    int_lat_raw  = `Interception latitude`,
    int_lon_raw  = `Interception longitude`
  ) |>
  mutate(divisor = pmax(n_incidents, 1),
         det_lat = det_lat_raw / divisor,
         det_lon = det_lon_raw / divisor,
         int_lat = int_lat_raw / divisor,
         int_lon = int_lon_raw / divisor)

MED_LAT <- c(30, 46)
MED_LON <- c(-10, 40)
tri <- tri |> mutate(
  bad_coord = !is.na(det_lat) & (det_lat < MED_LAT[1] | det_lat > MED_LAT[2] |
                                  det_lon < MED_LON[1] | det_lon > MED_LON[2]),
  det_lat = ifelse(bad_coord, NA_real_, det_lat),
  det_lon = ifelse(bad_coord, NA_real_, det_lon),
  int_lat = ifelse(bad_coord, NA_real_, int_lat),
  int_lon = ifelse(bad_coord, NA_real_, int_lon),
  has_coord = !is.na(det_lat)
) |>
  select(-det_lat_raw, -det_lon_raw, -int_lat_raw, -int_lon_raw,
         -divisor, -bad_coord)

pad <- pad_raw |>
  transmute(
    pad_row      = row_number(),
    incident_id  = IncidentNumber,
    date         = as.Date(DetectionDate),
    op_name      = OperationName,
    sar          = SAR,
    n_migrants   = num_total_irreg_migrants,
    n_deaths     = num_DeathCases,
    transport    = TransportType,
    detected_by  = TypeOfDetectedBy,
    intercept_by = TypeOfInterceptedBy
  )

# ── 3. Window pad-194 to Triton date range ──────────────────────────────────
tri_min <- min(tri$date, na.rm = TRUE)
tri_max <- max(tri$date, na.rm = TRUE)
pad_win <- pad |> filter(date >= tri_min & date <= tri_max)

# ── 4. Group match ──────────────────────────────────────────────────────────
key_cols <- c("date", "transport", "detected_by", "intercept_by", "sar")

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; rd <- pi / 180
  dlat <- (lat2 - lat1) * rd
  dlon <- (lon2 - lon1) * rd
  a <- sin(dlat / 2)^2 + cos(lat1 * rd) * cos(lat2 * rd) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

AMB_SPREAD_MAX_KM <- 28   # one ERA5 0.25 deg cell at Mediterranean latitude

tri_centroid <- tri |>
  group_by(across(all_of(key_cols))) |>
  mutate(
    w = pmax(n_incidents, 1),
    cent_det_lat = weighted.mean(det_lat, w, na.rm = TRUE),
    cent_det_lon = weighted.mean(det_lon, w, na.rm = TRUE),
    cent_int_lat = weighted.mean(int_lat, w, na.rm = TRUE),
    cent_int_lon = weighted.mean(int_lon, w, na.rm = TRUE)
  ) |>
  mutate(across(starts_with("cent_"),
                ~ ifelse(is.nan(.x), NA_real_, .x))) |>
  mutate(dist_km = ifelse(has_coord,
                          haversine_km(det_lat, det_lon,
                                        cent_det_lat, cent_det_lon),
                          NA_real_)) |>
  ungroup()

tri_group <- tri_centroid |>
  group_by(across(all_of(key_cols))) |>
  summarise(
    tri_rows_in_group = n(),
    tri_n_inc     = sum(n_incidents),
    tri_n_mig     = sum(n_migrants),
    tri_n_dea     = sum(n_deaths),
    tri_has_coord = all(has_coord),
    det_lat       = first(cent_det_lat),
    det_lon       = first(cent_det_lon),
    int_lat       = first(cent_int_lat),
    int_lon       = first(cent_int_lon),
    spread_km     = suppressWarnings(max(dist_km, na.rm = TRUE)),
    tri_row       = first(tri_row),
    .groups = "drop"
  ) |>
  mutate(spread_km = ifelse(is.infinite(spread_km), NA_real_, spread_km))

pad_group <- pad_win |>
  group_by(across(all_of(key_cols))) |>
  summarise(
    pad_n       = n(),
    pad_mig_sum = sum(n_migrants),
    pad_dea_sum = sum(n_deaths, na.rm = TRUE),
    pad_rows    = list(pad_row),
    .groups = "drop"
  )

joined <- pad_group |>
  inner_join(tri_group, by = key_cols) |>
  mutate(
    unique_tri = tri_rows_in_group == 1,
    cnt_match  = pad_n == tri_n_inc,
    mig_match  = pad_mig_sum == tri_n_mig,
    dea_match  = pad_dea_sum == tri_n_dea,
    sums_ok    = cnt_match & mig_match & dea_match,
    perfect    = unique_tri & sums_ok,
    amb_usable = (!unique_tri) & sums_ok & tri_has_coord &
                 !is.na(spread_km) & spread_km <= AMB_SPREAD_MAX_KM
  )

# ── 5. Explode verified groups back to pad-194 row level ────────────────────
matched <- joined |>
  filter(perfect | amb_usable) |>
  mutate(match_quality = case_when(
    perfect & !tri_has_coord             ~ "no_coord",
    perfect & tri_n_inc == 1             ~ "1-1_match",
    perfect                              ~ "avg_match",
    amb_usable                           ~ "amb_avg_match"
  )) |>
  select(all_of(key_cols),
         det_lat, det_lon, int_lat, int_lon,
         match_quality, tri_row, tri_n_inc, pad_rows) |>
  unnest(pad_rows) |>
  rename(pad_row = pad_rows) |>
  select(pad_row, det_lat, det_lon, int_lat, int_lon,
         match_quality, tri_row, tri_n_inc)

stopifnot(!any(duplicated(matched$pad_row)))
stopifnot(all(matched$pad_row %in% pad_win$pad_row))

# ── 6. Assemble enriched output ─────────────────────────────────────────────
frx <- readRDS(file.path(BASE_DIR, "data", "processed", "frontex_incidents.RDS"))
stopifnot(nrow(frx) == nrow(pad))

frx <- frx |> mutate(pad_row = row_number())

enriched <- frx |>
  left_join(matched |> select(pad_row, det_lat, det_lon, int_lat, int_lon,
                                match_quality, tri_n_inc),
            by = "pad_row") |>
  mutate(match_quality = ifelse(!is.na(match_quality), match_quality,
                                NA_character_)) |>
  select(-pad_row)

stopifnot(nrow(enriched) == nrow(frx))
stopifnot("incident_id" %in% names(enriched))

saveRDS(enriched,
        file.path(BASE_DIR, "data", "processed", "frontex_incidents_coords.RDS"))
