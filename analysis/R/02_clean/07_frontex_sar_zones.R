# 04c_assign_sar_zones.R
# ======================
# Assign SAR-zone membership to Frontex incidents.
#
# Input:
#   data/processed/frontex_incidents_coords.RDS  (output of 04b, with det_lat/lon)
#   data/processed/sar_zones.RDS              (4 polygons: Italy, Libya, Malta,
#                                              Tunisia; built by
#                                              analysis/R/00_define_sea_zones.R)
#
# Output:
#   data/processed/frontex_with_sar_zones.RDS
#   output/tables/04c_sar_zone_assignment.txt
#
# Overlap handling:
#   In sar_zones.RDS, Tunisia is already deduplicated against the other three.
#   Italy and Malta still overlap (~41,350 km2 — Malta's SRR sits inside
#   Italy's MRCC ROMA region). Libya does not overlap with the others in the
#   processed file. To preserve the overlap information we attach:
#     sar_italy, sar_malta, sar_libya, sar_tunisia   (4 booleans)
#     sar_countries   (string listing all zones the point falls in)
#     sar_n_zones     (count: 0, 1, or 2)
#   We do NOT pick a single "primary" zone here — that is an analysis choice
#   left to downstream scripts.
#
# Rows without coordinates get NA for all sar_* flags and sar_countries.
# Row count is preserved.

library(dplyr)
library(sf)
library(tidyr)

BASE_DIR <- here::here()

cat("============================================================\n")
cat("ASSIGN SAR ZONES TO FRONTEX INCIDENTS\n")
cat("============================================================\n\n")

# ── 1. Load inputs ───────────────────────────────────────
cat("--- 1. Loading inputs ---\n")

frx <- readRDS(file.path(BASE_DIR, "data", "processed",
                         "frontex_incidents_coords.RDS"))
cat(sprintf("  Frontex rows: %d\n", nrow(frx)))

zones <- readRDS(file.path(BASE_DIR, "data", "processed", "sar_zones.RDS"))
cat(sprintf("  SAR zones:    %d (%s)\n",
            nrow(zones), paste(zones$country, collapse = ", ")))

# Disable s2 — the Italian SRR polygon has self-intersecting edges, same as
# in 00_define_sea_zones.R
sf_use_s2(FALSE)

# ── 2. Spatial join ──────────────────────────────────────
cat("\n--- 2. Spatial join on points with coordinates ---\n")

frx_geo <- frx |> filter(!is.na(det_lat) & !is.na(det_lon))
cat(sprintf("  Rows with detection coord: %d\n", nrow(frx_geo)))

# Convert to sf without dropping the original lat/lon columns
frx_sf <- st_as_sf(frx_geo,
                   coords = c("det_lon", "det_lat"),
                   crs = 4326,
                   remove = FALSE)

# st_join keeps one row per matching zone — a point in Italy AND Malta
# produces two rows. We collapse back to one row per incident below.
joined <- st_join(frx_sf, zones[, "country"], left = TRUE)

# ── 3. Collapse to one row per incident ───────────────────────
wide <- joined |>
  st_drop_geometry() |>
  group_by(incident_id) |>
  summarise(
    sar_italy   = any(country == "Italy",   na.rm = TRUE),
    sar_malta   = any(country == "Malta",   na.rm = TRUE),
    sar_libya   = any(country == "Libya",   na.rm = TRUE),
    sar_tunisia = any(country == "Tunisia", na.rm = TRUE),
    sar_countries = {
      cs <- sort(unique(na.omit(country)))
      if (length(cs) == 0) NA_character_ else paste(cs, collapse = ", ")
    },
    .groups = "drop"
  ) |>
  mutate(sar_n_zones = sar_italy + sar_malta + sar_libya + sar_tunisia)

# Sanity: one row per incident_id, every incident_id in frx_geo accounted for
stopifnot(!any(duplicated(wide$incident_id)))
stopifnot(nrow(wide) == nrow(frx_geo))

cat("\n  Coverage among geo rows:\n")
cat(sprintf("    in any SAR zone : %d (%.1f%%)\n",
            sum(wide$sar_n_zones > 0),
            100 * mean(wide$sar_n_zones > 0)))
cat(sprintf("    in 1 zone only  : %d (%.1f%%)\n",
            sum(wide$sar_n_zones == 1),
            100 * mean(wide$sar_n_zones == 1)))
cat(sprintf("    in 2 zones      : %d (%.1f%%)\n",
            sum(wide$sar_n_zones == 2),
            100 * mean(wide$sar_n_zones == 2)))
cat(sprintf("    in 0 zones      : %d (%.1f%%)\n",
            sum(wide$sar_n_zones == 0),
            100 * mean(wide$sar_n_zones == 0)))

cat("\n  Distribution by country (can double-count on overlaps):\n")
for (c in c("sar_italy", "sar_malta", "sar_libya", "sar_tunisia")) {
  cat(sprintf("    %-12s %d\n", c, sum(wide[[c]])))
}

cat("\n  Most common sar_countries strings:\n")
cs_tbl <- wide |> count(sar_countries, sort = TRUE) |> head(10)
for (i in seq_len(nrow(cs_tbl))) {
  lbl <- if (is.na(cs_tbl$sar_countries[i])) "(none)" else cs_tbl$sar_countries[i]
  cat(sprintf("    %-40s %6d\n", lbl, cs_tbl$n[i]))
}

# ── 4. Merge back onto full Frontex panel ────────────────
cat("\n--- 3. Merging back onto full Frontex panel ---\n")

enriched <- frx |>
  left_join(wide, by = "incident_id") |>
  mutate(
    # Rows with coordinates but outside all zones: FALSE (genuinely not in any zone)
    # Rows without coordinates: NA (unknown, not FALSE)
    sar_italy   = ifelse(is.na(det_lat), NA, sar_italy),
    sar_malta   = ifelse(is.na(det_lat), NA, sar_malta),
    sar_libya   = ifelse(is.na(det_lat), NA, sar_libya),
    sar_tunisia = ifelse(is.na(det_lat), NA, sar_tunisia),
    sar_n_zones = ifelse(is.na(det_lat), NA_integer_, sar_n_zones)
  )

# Sanity checks
stopifnot(nrow(enriched) == nrow(frx))
stopifnot(!any(duplicated(enriched$incident_id)))

# Rows without coords should have NA SAR flags
no_coord_rows <- enriched |> filter(is.na(det_lat))
stopifnot(all(is.na(no_coord_rows$sar_italy)))
stopifnot(all(is.na(no_coord_rows$sar_malta)))
stopifnot(all(is.na(no_coord_rows$sar_libya)))
stopifnot(all(is.na(no_coord_rows$sar_tunisia)))
stopifnot(all(is.na(no_coord_rows$sar_n_zones)))

cat(sprintf("  Enriched rows: %d (== input)\n", nrow(enriched)))

# ── 5. Save ──────────────────────────────────────────────
cat("\n--- 4. Saving ---\n")

out_rds <- file.path(BASE_DIR, "data", "processed",
                     "frontex_with_sar_zones.RDS")
saveRDS(enriched, out_rds)
cat(sprintf("  Saved: %s\n", sub(BASE_DIR, "", out_rds, fixed = TRUE)))

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))
diag_path <- tbl_path("02_clean", "07_frontex_sar_zones_diagnostics.txt")
sink(diag_path)
cat("SAR zone assignment diagnostics\n")
cat("===============================\n\n")
cat(sprintf("Input rows (Frontex):        %d\n", nrow(frx)))
cat(sprintf("Rows with coordinates:       %d\n", nrow(frx_geo)))
cat(sprintf("Rows in any SAR zone:        %d\n", sum(wide$sar_n_zones > 0)))
cat(sprintf("Rows in 1 zone only:         %d\n", sum(wide$sar_n_zones == 1)))
cat(sprintf("Rows in 2 zones (overlap):   %d\n", sum(wide$sar_n_zones == 2)))
cat(sprintf("Rows in 0 zones:             %d\n\n",
            sum(wide$sar_n_zones == 0)))
cat("Per-country counts (double-counts overlaps):\n")
for (c in c("sar_italy", "sar_malta", "sar_libya", "sar_tunisia")) {
  cat(sprintf("  %-12s %d\n", c, sum(wide[[c]])))
}
cat("\nMost common sar_countries values:\n")
for (i in seq_len(nrow(cs_tbl))) {
  lbl <- if (is.na(cs_tbl$sar_countries[i])) "(none)" else cs_tbl$sar_countries[i]
  cat(sprintf("  %-40s %6d\n", lbl, cs_tbl$n[i]))
}
sink()
cat(sprintf("  Saved: %s\n",
            sub(BASE_DIR, "", diag_path, fixed = TRUE)))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
