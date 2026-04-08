# 00_define_sea_zones.R
# ====================
# Define geographic zones for the CMR analysis:
#   1. Core Corridor — sea zone used ONLY for wave height calculation
#   2. SAR zones — parsed from IMO GISIS GML boundary files
#
# Output:
#   output/figures/sea_zone_core_map.png
#   output/figures/sar_zones_map.png
#   data/processed/sar_zones.RDS

Sys.setlocale("LC_TIME", "en_US.UTF-8")

library(tidyverse)
library(sf)
library(xml2)
library(rnaturalearth)
library(lubridate)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
GML_DIR  <- file.path(BASE_DIR, "data", "raw", "IMO-SAR-boundaries")

cat("============================================================\n")
cat("DEFINE SEA ZONES (Core Corridor + SAR)\n")
cat("============================================================\n\n")

# Disable s2 — the Italian SRR polygon has self-intersecting edges
sf_use_s2(FALSE)

# ── Shared spatial data ──────────────────────────────────
world <- ne_countries(scale = "medium", returnclass = "sf")
land  <- st_union(world)

# Ocean mask for SAR zone construction
bbox  <- st_as_sfc(st_bbox(c(xmin = -5, ymin = 25, xmax = 40, ymax = 50),
                            crs = 4326))
ocean <- st_difference(bbox, land)

# ── Shared IOM incident data ────────────────────────────
cat("--- Loading IOM incidents ---\n")

df_cmr <- readRDS(file.path(BASE_DIR, "data", "processed",
                             "iom_mmp_incidents.RDS")) %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) != "sub-incident") %>%
  mutate(
    lat = as.numeric(Latitude),
    lon = as.numeric(Longitude),
    date = as.Date(incident_date_clean),
    dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE),
    post_mou = if_else(date >= MOU_DATE, "Post-MoU", "Pre-MoU") %>%
                 factor(levels = c("Pre-MoU", "Post-MoU"))
  ) %>%
  filter(!is.na(lat), !is.na(lon), !is.na(date))

cat(sprintf("  Total CMR incidents: %d\n", nrow(df_cmr)))
cat(sprintf("  Total dead + missing: %.0f\n", sum(df_cmr$dead_missing)))

# ============================================================
# 1. CORE CORRIDOR SEA ZONE
# ============================================================
cat("\n--- 1. Core Corridor zone ---\n")

coords_core <- matrix(c(
  15.2, 30.0,
  15.2, 36.7,
  12.4, 37.8,
  11.0, 37.1,
   9.0, 34.0,
   9.0, 31.0,
  15.2, 30.0
), ncol = 2, byrow = TRUE)

core_poly <- st_sfc(st_polygon(list(coords_core)), crs = 4326)
core_sea  <- st_difference(core_poly, land)
core_proj <- st_transform(core_sea, 3857)
core_buf  <- st_buffer(core_proj, dist = 5000)
core_vis  <- st_transform(core_buf, 4326) %>% st_intersection(core_poly)

cat("  Core Corridor polygon built\n")

# ============================================================
# 2. SAR ZONES FROM IMO GISIS GML
# ============================================================
cat("\n--- 2. Parsing SAR zone GML files ---\n")

parse_gml <- function(file, country) {
  doc <- read_xml(file)
  ns  <- xml_ns(doc)

  coords_text <- xml_text(xml_find_first(doc, ".//gml:posList", ns))
  desc        <- xml_text(xml_find_first(doc, ".//gml:description", ns))

  vals <- as.numeric(strsplit(trimws(coords_text), "\\s+")[[1]])
  lon  <- vals[seq(1, length(vals), 2)]
  lat  <- vals[seq(2, length(vals), 2)]

  # Remove closing vertex (duplicate of first) — re-closed below
  if (lon[1] == lon[length(lon)] && lat[1] == lat[length(lat)]) {
    lon <- lon[-length(lon)]
    lat <- lat[-length(lat)]
  }

  cat(sprintf("  %s (%s): %d boundary vertices\n", country, desc, length(lon)))
  list(lon = lon, lat = lat, country = country, desc = desc)
}

gml_italy   <- parse_gml(file.path(GML_DIR, "ITA-RCCAreas-807.gml"), "Italy")
gml_libya   <- parse_gml(file.path(GML_DIR, "LBY-RCCAreas-2032.gml"), "Libya")
gml_malta   <- parse_gml(file.path(GML_DIR, "MLT-RCCAreas-150.gml"), "Malta")
gml_tunisia <- parse_gml(file.path(GML_DIR, "TUN-RCCAreas-9294.gml"), "Tunisia")

# Inland waypoints for closing coastal zones.
# Routes the closing line deep inland so all bays are captured.
inland_closure <- list(
  Italy   = data.frame(lon = c(10.5), lat = c(46.5)),
  Tunisia = data.frame(lon = c(9.5),  lat = c(34.0))
)

make_gml_polygon <- function(parsed) {
  lon <- parsed$lon
  lat <- parsed$lat

  wp <- inland_closure[[parsed$country]]
  if (!is.null(wp)) {
    lon <- c(lon, wp$lon, lon[1])
    lat <- c(lat, wp$lat, lat[1])
  } else {
    lon <- c(lon, lon[1])
    lat <- c(lat, lat[1])
  }

  st_polygon(list(cbind(lon, lat)))
}

build_sar_zone <- function(parsed) {
  gml_sfc <- st_sfc(make_gml_polygon(parsed), crs = 4326) %>% st_make_valid()
  sea <- st_intersection(gml_sfc, ocean) %>% st_make_valid()

  cat(sprintf("  %s: sea zone built\n", parsed$country))
  st_sf(country = parsed$country, description = parsed$desc, geometry = sea)
}

zones_sea <- do.call(rbind, lapply(
  list(gml_italy, gml_libya, gml_malta, gml_tunisia), build_sar_zone))

# Remove overlaps: Tunisia gets only uniquely assigned area
others <- zones_sea %>% filter(country != "Tunisia") %>% st_union()
tun_idx <- which(zones_sea$country == "Tunisia")
zones_sea$geometry[tun_idx] <- st_difference(
  zones_sea$geometry[tun_idx], others) %>% st_make_valid()

# Display labels for the legend
zones_sea$zone_label <- c(
  "Italy" = "EU: Italy", "Malta" = "EU: Malta",
  "Libya" = "North Africa: Libya", "Tunisia" = "North Africa: Tunisia"
)[zones_sea$country]

cat(sprintf("\n  Total SAR zones: %d\n", nrow(zones_sea)))
print(zones_sea %>% st_drop_geometry())

# SAR zone incident assignment
cat("\n--- SAR zone incident assignment ---\n")
iom_sf <- st_as_sf(df_cmr, coords = c("lon", "lat"), crs = 4326)
zone_assignment <- st_join(iom_sf, zones_sea, left = TRUE)

zone_table <- zone_assignment %>%
  st_drop_geometry() %>%
  group_by(country) %>%
  summarise(n_incidents = n(), deaths = sum(dead_missing), .groups = "drop") %>%
  arrange(desc(deaths))
print(zone_table)

# ============================================================
# 3. FRONTEX OPERATIONAL AREAS
# ============================================================
# Approximate boundaries reconstructed as buffers from the Italian coastline.
# The actual operational plans are classified; boundaries based on publicly
# available distance specifications.
#
# Sources:
#   Triton phase 1 (30 NM from Italian coast):
#     Frontex internal letter, Rosler to Pinto, 25 Nov 2014, quoting
#     Operational Plan Annex 3: "surface patrolling activities will be
#     implemented within the 30 NM limit of the coastal line or islands
#     of Italy and Malta."
#     https://www.statewatch.org/media/documents/news/2017/apr/eu-frontex-rosler-pinto-letter-italy-sar-25-11-14.pdf
#
#   Triton phase 2 (138 NM south of Sicily):
#     Frontex press release, 26 May 2015: "The operational area will be
#     extended to 138 NM south of Sicily."
#     https://www.frontex.europa.eu/media-centre/news/news-release/frontex-expands-its-joint-operation-triton-udpbHP
#     EU Parliament answer E-003910/2015:
#     https://www.europarl.europa.eu/doceo/document/E-8-2015-003910-ASW_EN.html
#
#   Themis (24 NM from Italian coast):
#     No primary Frontex document found. Widely cited in secondary literature:
#     Heinrich Boll Foundation (2021):
#     https://us.boell.org/en/2021/04/16/evolution-eus-naval-operations-central-mediterranean-gradual-shift-away-search-and

cat("\n--- 3. Building Frontex operational areas ---\n")

NM_TO_M <- 1852  # 1 nautical mile = 1852 meters

# Italian coastline (includes Sicily, Sardinia, Lampedusa)
italy_land <- world %>% filter(name == "Italy") %>% st_geometry() %>% st_union()
italy_proj <- st_transform(italy_land, 3857)

op_specs <- list(
  list(name = "Triton (Nov 2014 – May 2015)",  nm = 30),
  list(name = "Triton (May 2015 – Feb 2018)",   nm = 138),
  list(name = "Themis (Feb 2018 – present)",     nm = 24)
)

op_polys <- lapply(op_specs, function(op) {
  buf <- st_buffer(italy_proj, dist = op$nm * NM_TO_M)
  buf_wgs <- st_transform(buf, 4326)
  sea <- st_intersection(buf_wgs, ocean) %>% st_make_valid()
  cat(sprintf("  %s (%d NM): built\n", op$name, op$nm))
  st_sf(operation = op$name, distance_nm = op$nm, geometry = sea)
})

frontex_op <- do.call(rbind, op_polys)

# Order by date for legend
frontex_op$operation <- factor(frontex_op$operation, levels = c(
  "Triton (Nov 2014 – May 2015)",
  "Triton (May 2015 – Feb 2018)",
  "Themis (Feb 2018 – present)"
))

cat(sprintf("  Total operational area polygons: %d\n", nrow(frontex_op)))

# ============================================================
# 4. MAPS (1-row x 3-column panel)
# ============================================================
cat("\n--- 4. Generating maps ---\n")

library(patchwork)

# Shared map extent and theme
MAP_XLIM <- c(5, 22)
MAP_YLIM <- c(29, 42)

map_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 7.5, colour = "grey50"),
    axis.text = element_text(size = 7, colour = "grey50"),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

# IOM incident point aesthetics (shared across both maps)
colour_vals <- c("Pre-MoU" = "#D4820E", "Post-MoU" = "#D32F2F")

country_labels <- data.frame(
  label = c("Sicily", "Tunisia", "Libya", "Malta", "Lampedusa"),
  lon   = c(14.2, 9.5, 15.0, 14.4, 12.6),
  lat   = c(37.5, 35.0, 30.8, 35.7, 35.6)
)

# ── 4a. ERA5 sea weather area ────────────────────────────

# Count incidents inside the SWH polygon
iom_sf <- st_as_sf(df_cmr, coords = c("lon", "lat"), crs = 4326)
n_in_zone <- sum(st_within(iom_sf, core_poly, sparse = FALSE)[, 1])
n_total   <- nrow(df_cmr)

p_core <- ggplot() +
  geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = core_vis, fill = "#2166AC", alpha = 0.15,
          colour = "#2166AC", linewidth = 0.4) +
  geom_point(data = df_cmr, aes(x = lon, y = lat, size = dead_missing,
                                 colour = post_mou), alpha = 0.4, shape = 16) +
  scale_colour_manual(values = colour_vals, name = "Period") +
  scale_size_continuous(range = c(0.3, 4), name = "Dead + missing",
                        breaks = c(1, 10, 50, 200)) +
  guides(colour = guide_legend(order = 1,
           override.aes = list(size = 2.5, alpha = 0.7)),
         size = guide_legend(order = 2)) +
  geom_text(data = country_labels, aes(x = lon, y = lat, label = label),
            colour = "grey40", size = 2.5, fontface = "italic") +
  coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
  labs(title = "ERA5 sea weather area",
       subtitle = sprintf("Polygon used for spatial averaging of SWH. IOM MMP incidents: %s in SWH zone / %s total.",
                           formatC(n_in_zone, big.mark = ","),
                           formatC(n_total, big.mark = ","))) +
  map_theme

# ── 4b. SAR zones ────────────────────────────────────────

zone_colours <- c(
  "EU: Italy" = "#2166AC",
  "EU: Malta" = "#4393C3",
  "North Africa: Libya" = "#D6604D",
  "North Africa: Tunisia" = "#F4A582"
)

p_sar <- ggplot() +
  geom_sf(data = world, fill = "grey92", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = zones_sea, aes(fill = zone_label), alpha = 0.15,
          colour = "#2166AC", linewidth = 0.3) +
  geom_point(data = df_cmr, aes(x = lon, y = lat, size = dead_missing),
             colour = "grey30", alpha = 0.3, shape = 16) +
  scale_fill_manual(values = zone_colours, name = "SAR zone") +
  scale_size_continuous(range = c(0.3, 4), guide = "none") +
  geom_text(data = country_labels, aes(x = lon, y = lat, label = label),
            colour = "#B22222", size = 2.5, fontface = "bold.italic") +
  coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
  labs(title = "IMO search and rescue zones",
       subtitle = sprintf("Zone boundaries from IMO GISIS Global SAR Plan.",
                           formatC(nrow(df_cmr), big.mark = ",")),
       x = NULL, y = NULL) +
  map_theme

# ── Combine into 1x2 panel ──────────────────────────────

panel_maps <- p_core | p_sar
ggsave(file.path(BASE_DIR, "output", "figures", "desc_panel_maps.png"),
       panel_maps, width = 12, height = 6, dpi = 300)
cat("Saved: output/figures/desc_panel_maps.png\n")

# ============================================================
# 5. SAVE
# ============================================================
cat("\n--- 5. Saving ---\n")

saveRDS(zones_sea, file.path(BASE_DIR, "data", "processed", "sar_zones.RDS"))
cat("Saved: data/processed/sar_zones.RDS\n")

saveRDS(frontex_op, file.path(BASE_DIR, "data", "processed", "frontex_op_areas.RDS"))
cat("Saved: data/processed/frontex_op_areas.RDS\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
