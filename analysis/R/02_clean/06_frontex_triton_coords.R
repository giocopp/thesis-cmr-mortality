# 04b_merge_triton_coords.R
# =========================
# Merge Triton 2014-2017 coordinates into the pad-194 Themis panel.
#
# The pad-194 Themis extract (used in 04_clean_frontex.R) has no coordinates.
# The Triton 2014-2017 extract has detection/interception lat-lon but no
# IncidentNumber, so the two cannot be joined on a key.
#
# Diagnostic finding 1 (aggregation):
#   - Both files cover exactly the same dates (2014-11-01 to 2017-09-02) and
#     have the same total deaths (1,013) and near-identical total migrants
#     (447,990 vs 448,529).
#   - On 284/696 days pad-194 has MORE rows than Triton, but the daily sum of
#     migrants is identical. Triton aggregates multi-vessel events while
#     pad-194 records each boat individually.
#
# Diagnostic finding 2 (coordinate SUM bug):
#   - Triton lat/lon for aggregated rows store SUMS of the underlying
#     individual coordinates, not means. A row with n_incidents = 20 has
#     lat = sum of 20 individual lats (e.g. lat = 666.67, not 33.33).
#   - Fix: divide lat/lon by pmax(n_incidents, 1). After correction, 2,021 /
#     2,027 Triton rows fall inside the Mediterranean bbox [30,46] x [-10,40].
#     The 6 rows that remain outside are dropped as unrecoverable.
#
# Matching strategy (group-match):
#   For each (date, transport, detected_by, intercept_by, sar) combination we
#   check whether the pad-194 rows in that group are EXACTLY the set of
#   incidents that the Triton aggregate summarises. A group is a verified
#   match if and only if:
#     1) Triton has exactly ONE row for that group (unambiguous)
#     2) pad's row count equals Triton's n_incidents
#     3) pad's sum of migrants equals Triton's n_migrants
#     4) pad's sum of deaths equals Triton's n_deaths
#   When all four hold, we attach the (corrected) Triton coordinate to every
#   pad-194 row in the group. Rows that do not fall in a verified group get
#   NA coordinates — we do NOT fall back to day-wide centroids because those
#   average across unrelated operations and miss the ERA5 cell.
#
# Diagnostic finding 3 (no-coord Triton rows still carry metadata):
#   Of 2,713 Triton rows, only 2,021 have a valid coordinate after the SUM
#   fix. The remaining 692 have NA coords but otherwise valid metadata
#   (date, classification, counts). The group match is run against the FULL
#   Triton set so that pad rows can still be verified against a no-coord
#   Triton row — they just end up tagged "no_coord" instead of "no_match".
#
# match_quality values:
#   "1-1_match"     : unique Triton row with n_incidents = 1 (single vessel),
#                     1:1 match with a pad row; coord is Triton's own
#   "avg_match"     : unique Triton row with n_incidents > 1 (multi-vessel
#                     aggregate); coord is the mean over n_incidents already
#                     stored by Triton (after SUM fix)
#   "amb_avg_match" : MULTIPLE Triton rows share the same classification key;
#                     aggregate-level sums (count + migrants + deaths) still
#                     agree with pad; coord is the n_incidents-weighted
#                     centroid of the Triton rows, and the max distance from
#                     that centroid is <= 28 km (one ERA5 cell)
#   "no_coord"      : unique Triton row found and metadata verified, but
#                     Triton has NA coordinates for that event
#   "no_match"      : inside Triton window but no verified group match,
#                     OR an ambiguous group whose spread > 28 km; coords NA
#   NA              : outside Triton window (2014-11 to 2017-09); coords NA
#
# Input:  data/raw/frontex/pad-194_themis_2014_2023.xlsx
#         data/raw/frontex/Triton 2014 2017 Incident Details.xlsx
# Output: data/processed/frontex_incidents_coords.RDS
#         output/tables/04b_triton_merge_diagnostics.txt

library(dplyr)
library(readxl)
library(tidyr)

BASE_DIR <- here::here()

cat("============================================================\n")
cat("MERGE TRITON COORDINATES INTO PAD-194 (group-match)\n")
cat("============================================================\n\n")

# ── 1. Load raw data ─────────────────────────────────────
cat("--- 1. Loading files ---\n")

tri_raw <- read_excel(
  file.path(BASE_DIR, "data", "raw", "frontex",
            "Triton 2014 2017 Incident Details.xlsx"),
  sheet = 1
)
cat(sprintf("  Triton:  %d rows, %d cols\n", nrow(tri_raw), ncol(tri_raw)))

pad_raw <- read_excel(
  file.path(BASE_DIR, "data", "raw", "frontex",
            "pad-194_themis_2014_2023.xlsx"),
  sheet = "Sheet1"
)
cat(sprintf("  pad-194: %d rows, %d cols\n", nrow(pad_raw), ncol(pad_raw)))

# ── 2. Standardize + fix Triton coord SUM bug ────────────
cat("\n--- 2. Standardising and fixing coordinates ---\n")

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

# Drop coords outside the Mediterranean bbox after the fix (but KEEP the row
# so its metadata can still be used for group matching — see Finding 3 below).
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
)
cat(sprintf("  Triton rows total:           %d\n", nrow(tri)))
cat(sprintf("  Triton rows with coord:      %d\n", sum(tri$has_coord)))
cat(sprintf("  Triton rows without coord:   %d (metadata still usable)\n",
            sum(!tri$has_coord)))
tri <- tri |> select(-det_lat_raw, -det_lon_raw, -int_lat_raw, -int_lon_raw,
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

# ── 3. Window pad-194 to Triton date range ───────────────
tri_min <- min(tri$date, na.rm = TRUE)
tri_max <- max(tri$date, na.rm = TRUE)
cat(sprintf("\n  Triton window: %s to %s\n", format(tri_min), format(tri_max)))

pad_win <- pad |> filter(date >= tri_min & date <= tri_max)
pad_out <- pad |> filter(date < tri_min | date > tri_max)
cat(sprintf("  pad-194 rows in window:       %d\n", nrow(pad_win)))
cat(sprintf("  pad-194 rows outside window:  %d (kept, no match attempted)\n",
            nrow(pad_out)))

# ── 4. Group match ───────────────────────────────────────
cat("\n--- 3. Group match ---\n")

key_cols <- c("date", "transport", "detected_by", "intercept_by", "sar")

# Haversine distance in km (for computing intra-group spread)
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; rd <- pi / 180
  dlat <- (lat2 - lat1) * rd
  dlon <- (lon2 - lon1) * rd
  a <- sin(dlat / 2)^2 + cos(lat1 * rd) * cos(lat2 * rd) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

AMB_SPREAD_MAX_KM <- 28   # one ERA5 0.25 deg cell at Mediterranean latitude

# Use ALL Triton rows for the group match, including those without coords.
# This lets us verify metadata against no-coord Triton rows and mark pad rows
# as "no_coord" (verified event but no lat/lon available) rather than
# "no_match" (genuinely cannot verify).
#
# For groups where Triton has multiple rows sharing the same key, we compute
# a weighted centroid (weighted by n_incidents — so the coordinate reflects
# the true mean over all underlying individual events, not a mean-of-means)
# and the max distance of any row from that centroid. If the spread fits
# inside one ERA5 cell (<= AMB_SPREAD_MAX_KM), the group is usable.
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
    # Unique-Triton groups: a "perfect" match requires all sums agree
    perfect    = unique_tri & sums_ok,
    # Ambiguous groups (Triton has >1 row): verified if sums agree AND all
    # Triton rows have coords AND the intra-group spread fits one ERA5 cell
    amb_usable = (!unique_tri) & sums_ok & tri_has_coord &
                 !is.na(spread_km) & spread_km <= AMB_SPREAD_MAX_KM
  )

n_perfect       <- sum(joined$perfect)
n_amb_usable    <- sum(joined$amb_usable)
n_near_miss     <- sum(joined$unique_tri & !joined$perfect)
n_ambiguous_grp <- sum(!joined$unique_tri)

cat(sprintf("  Joined groups:                     %d\n", nrow(joined)))
cat(sprintf("  Perfect unique matches:            %d\n", n_perfect))
cat(sprintf("  Ambiguous but usable (spread<=%dkm): %d\n",
            AMB_SPREAD_MAX_KM, n_amb_usable))
cat(sprintf("  Near-miss (sums disagree):         %d\n", n_near_miss))
cat(sprintf("  Ambiguous (>1 Triton row, total):  %d\n", n_ambiguous_grp))

# ── 5. Explode perfect groups back to pad-194 row level ──
cat("\n--- 4. Assigning coordinates ---\n")

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

n_11   <- sum(matched$match_quality == "1-1_match")
n_avg  <- sum(matched$match_quality == "avg_match")
n_amb  <- sum(matched$match_quality == "amb_avg_match")
n_ncrd <- sum(matched$match_quality == "no_coord")
cat(sprintf("  pad-194 rows in verified groups: %d (%.1f%% of in-window)\n",
            nrow(matched), 100 * nrow(matched) / nrow(pad_win)))
cat(sprintf("    1-1_match     (unique Triton n_inc=1, with coord): %d\n", n_11))
cat(sprintf("    avg_match     (unique Triton n_inc>1, with coord): %d\n", n_avg))
cat(sprintf("    amb_avg_match (>1 Triton rows, spread<=%dkm):      %d\n",
            AMB_SPREAD_MAX_KM, n_amb))
cat(sprintf("    no_coord      (verified but Triton lat/lon NA):    %d\n", n_ncrd))

# Sanity checks
stopifnot(!any(duplicated(matched$pad_row)))
stopifnot(all(matched$pad_row %in% pad_win$pad_row))

# ── 6. Assemble final enriched dataset ───────────────────
# Load the full frontex_incidents (produced by 04_clean_frontex.R)
# so ALL original variables are preserved in the output.
cat("\n--- 5. Assembling final output ---\n")

frx <- readRDS(file.path(BASE_DIR, "data", "processed", "frontex_incidents.RDS"))
stopifnot(nrow(frx) == nrow(pad))   # same raw source, same row order

frx <- frx |> mutate(pad_row = row_number())

enriched <- frx |>
  left_join(matched |> select(pad_row, det_lat, det_lon, int_lat, int_lon,
                                match_quality, tri_n_inc),
            by = "pad_row") |>
  mutate(match_quality = ifelse(!is.na(match_quality), match_quality,
                                NA_character_)) |>
  select(-pad_row)

# Sanity checks
stopifnot(nrow(enriched) == nrow(frx))
stopifnot("incident_id" %in% names(enriched))

# Breakdown
qb <- enriched |> count(match_quality) |>
  mutate(pct   = round(100 * n / nrow(enriched), 1),
         label = ifelse(is.na(match_quality), "NA (out of window)",
                                                match_quality))
cat("\n  Match-quality breakdown (full pad-194):\n")
for (i in seq_len(nrow(qb))) {
  cat(sprintf("    %-22s %6d (%4.1f%%)\n",
              qb$label[i], qb$n[i], qb$pct[i]))
}

with_coord <- sum(!is.na(enriched$det_lat))
cat(sprintf("\n  Rows with detection coord:   %d (%.1f%%)\n",
            with_coord, 100 * with_coord / nrow(enriched)))

in_window <- enriched |> filter(!is.na(match_quality))
with_coord_in_window <- sum(!is.na(in_window$det_lat))
cat(sprintf("  In-window rows with coord:   %d / %d  (%.1f%%)\n",
            with_coord_in_window, nrow(in_window),
            100 * with_coord_in_window / nrow(in_window)))

# ── 7. Save output ───────────────────────────────────────
cat("\n--- 6. Saving ---\n")

out_rds <- file.path(BASE_DIR, "data", "processed", "frontex_incidents_coords.RDS")
saveRDS(enriched, out_rds)
cat(sprintf("  Saved: %s\n", sub(BASE_DIR, "", out_rds, fixed = TRUE)))

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))
diag_path <- tbl_path("02_clean", "06_frontex_triton_coords_diagnostics.txt")
sink(diag_path)
cat("Triton <-> pad-194 merge diagnostics (group match)\n")
cat("====================================================\n\n")
cat(sprintf("Triton rows:                 %d\n", nrow(tri)))
cat(sprintf("Triton rows with valid coord: %d\n", sum(tri$has_coord)))
cat(sprintf("Triton rows without coord:    %d\n", sum(!tri$has_coord)))
cat(sprintf("pad-194 rows (total):        %d\n", nrow(pad)))
cat(sprintf("pad-194 rows in window:      %d\n", nrow(pad_win)))
cat(sprintf("pad-194 rows out of window:  %d\n", nrow(pad_out)))
cat(sprintf("Triton date range:           %s to %s\n\n",
            format(tri_min), format(tri_max)))
cat("Joined groups summary:\n")
cat(sprintf("  Joined groups:             %d\n", nrow(joined)))
cat(sprintf("  Perfect matches:           %d\n", n_perfect))
cat(sprintf("  Near-miss (sums disagree): %d\n", n_near_miss))
cat(sprintf("  Ambiguous (>1 Triton row): %d\n\n", n_ambiguous_grp))
cat("Match-quality breakdown (full pad-194):\n")
for (i in seq_len(nrow(qb))) {
  cat(sprintf("  %-22s %6d (%4.1f%%)\n",
              qb$label[i], qb$n[i], qb$pct[i]))
}
cat(sprintf("\nIn-window coverage: %d/%d (%.1f%%)\n",
            with_coord_in_window, nrow(in_window),
            100 * with_coord_in_window / nrow(in_window)))
sink()
cat(sprintf("  Saved: %s\n",
            sub(BASE_DIR, "", diag_path, fixed = TRUE)))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
