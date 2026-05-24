# Define sea zones: area of analysis polygon, SAR zones (IMO GISIS),
# Frontex operational areas. Produces the paper figure
# fig-sea-zones-united-sar-panel.png and the intermediate sf objects used
# downstream (core_corridor, sar_zones, sar_zones_in_corridor, frontex_op_areas).

Sys.setlocale("LC_TIME", "en_US.UTF-8")

library(tidyverse)
library(sf)
library(lwgeom)
library(xml2)
library(rnaturalearth)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

MAP_END_DATE <- as.Date("2023-05-31")
GML_DIR      <- file.path(BASE_DIR, "data", "raw", "IMO-SAR-boundaries")

sf_use_s2(FALSE)

# ── Shared spatial data ─────────────────────────────────────────────────────
world <- ne_countries(scale = "medium", returnclass = "sf")
land  <- st_union(world)

bbox  <- st_as_sfc(st_bbox(c(xmin = -5, ymin = 25, xmax = 40, ymax = 50),
                            crs = 4326))
ocean <- st_difference(bbox, land)

# ── 1. Area-of-analysis polygon ─────────────────────────────────────────────
coords_core <- matrix(c(
  20.5, 32.0,
  15.9, 38.2,
  15.2, 38.1,
  15.2, 37.8,
  12.4, 37.8,
  11.0, 37.1,
   9.0, 34.0,
   9.0, 31.0,
    20, 30.0,
  20.5, 32.0
), ncol = 2, byrow = TRUE)

core_poly <- st_sfc(st_polygon(list(coords_core)), crs = 4326)
core_sea  <- st_difference(core_poly, land) |> st_make_valid()

# ── 2. SAR zones from IMO GISIS GML ─────────────────────────────────────────
parse_gml <- function(file, country) {
  doc <- read_xml(file)
  ns  <- xml_ns(doc)

  coords_text <- xml_text(xml_find_first(doc, ".//gml:posList", ns))
  desc        <- xml_text(xml_find_first(doc, ".//gml:description", ns))

  vals <- as.numeric(strsplit(trimws(coords_text), "\\s+")[[1]])
  lon  <- vals[seq(1, length(vals), 2)]
  lat  <- vals[seq(2, length(vals), 2)]

  if (lon[1] == lon[length(lon)] && lat[1] == lat[length(lat)]) {
    lon <- lon[-length(lon)]
    lat <- lat[-length(lat)]
  }

  list(lon = lon, lat = lat, country = country, desc = desc)
}

gml_italy   <- parse_gml(file.path(GML_DIR, "ITA-RCCAreas-807.gml"), "Italy")
gml_libya   <- parse_gml(file.path(GML_DIR, "LBY-RCCAreas-2032.gml"), "Libya")
gml_malta   <- parse_gml(file.path(GML_DIR, "MLT-RCCAreas-150.gml"), "Malta")
gml_tunisia <- parse_gml(file.path(GML_DIR, "TUN-RCCAreas-9294.gml"), "Tunisia")

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
  gml_sfc <- st_sfc(make_gml_polygon(parsed), crs = 4326) |> st_make_valid()
  sea <- st_intersection(gml_sfc, ocean) |> st_make_valid()
  st_sf(country = parsed$country, description = parsed$desc, geometry = sea)
}

zones_sea <- do.call(rbind, lapply(
  list(gml_italy, gml_libya, gml_malta, gml_tunisia), build_sar_zone))

# Remove overlaps: Tunisia gets only uniquely assigned area
others  <- zones_sea |> filter(country != "Tunisia") |> st_union()
tun_idx <- which(zones_sea$country == "Tunisia")
zones_sea$geometry[tun_idx] <- st_difference(
  zones_sea$geometry[tun_idx], others) |> st_make_valid()

zones_sea$zone_label <- c(
  "Italy" = "EU: Italy", "Malta" = "EU: Malta",
  "Libya" = "North Africa: Libya", "Tunisia" = "North Africa: Tunisia"
)[zones_sea$country]

# ── 2b. SAR ∩ area-of-analysis ──────────────────────────────────────────────
zones_in_corridor <- zones_sea
for (i in seq_len(nrow(zones_sea))) {
  zg <- zones_sea$geometry[i]
  xg <- tryCatch(
    st_intersection(zg, core_sea) |> st_make_valid(),
    error = function(e) st_sfc(st_geometrycollection(), crs = 4326)
  )
  if (length(xg) == 0 || st_is_empty(xg)[1]) {
    zones_in_corridor$geometry[i] <- st_sfc(st_geometrycollection(), crs = 4326)
  } else {
    zones_in_corridor$geometry[i] <- xg
  }
}

zones_in_corridor$full_sar_area_km2 <- as.numeric(st_area(zones_sea)) / 1e6
zones_in_corridor$zone_area_km2     <- as.numeric(st_area(zones_in_corridor)) / 1e6
zones_in_corridor$pct_of_full_sar   <- round(
  100 * zones_in_corridor$zone_area_km2 / zones_in_corridor$full_sar_area_km2, 1)

# ── 3. Frontex operational areas (buffers from Italian coast) ───────────────
NM_TO_M <- 1852

italy_land <- world |> filter(name == "Italy") |> st_geometry() |> st_union()
italy_proj <- st_transform(italy_land, 3857)

op_specs <- list(
  list(name = "Triton (Nov 2014 – May 2015)",  nm = 30),
  list(name = "Triton (May 2015 – Feb 2018)",   nm = 138),
  list(name = "Themis (Feb 2018 – present)",     nm = 24)
)

op_polys <- lapply(op_specs, function(op) {
  buf <- st_buffer(italy_proj, dist = op$nm * NM_TO_M)
  buf_wgs <- st_transform(buf, 4326)
  sea <- st_intersection(buf_wgs, ocean) |> st_make_valid()
  st_sf(operation = op$name, distance_nm = op$nm, geometry = sea)
})

frontex_op <- do.call(rbind, op_polys)
frontex_op$operation <- factor(frontex_op$operation, levels = c(
  "Triton (Nov 2014 – May 2015)",
  "Triton (May 2015 – Feb 2018)",
  "Themis (Feb 2018 – present)"
))

# ── 4. Paper figure: UNITED deaths + SAR zones panel ────────────────────────
MAP_XLIM <- c(5, 22)
MAP_YLIM <- c(29, 42)

df_united <- readRDS(file.path(BASE_DIR, "data", "processed",
                                "united_incidents.RDS")) |>
  filter(is_cmr,
         !is.na(latitude), !is.na(longitude),
         incident_year >= 2014L,
         incident_date_clean <= MAP_END_DATE,
         (manner_of_death == "drowned" & !is.na(manner_of_death)) |
         (transport_means == "boat_ship_ferry" & !is.na(transport_means)))

united_sf     <- st_as_sf(df_united, coords = c("longitude", "latitude"), crs = 4326)
inside_united <- st_within(united_sf, core_poly, sparse = FALSE)[, 1]
df_united$in_corridor <- factor(
  ifelse(inside_united, "Included", "Excluded"),
  levels = c("Included", "Excluded")
)
df_united_plot <- df_united |> arrange(in_corridor == "Included")

country_labels <- data.frame(
  label = c("Sicily", "Tunisia", "Libya", "Malta", "Lampedusa", "Algeria", "Sardinia"),
  lon   = c(14.2,     9.5,       15.0,    14.4,    12.6,        6.0,       9.1),
  lat   = c(37.5,     35.0,      30.8,    35.7,    35.6,        35.5,      40.1)
)

panel_land_fill <- "#F0F0EC"
panel_land_line <- "#B9B9B3"
panel_sea_fill  <- "#F7FBFD"
panel_corridor  <- "#F08A4B"
panel_incident_colours <- c("Included" = "#C9342D",
                            "Excluded" = "#AEB2B5")
panel_zone_colours <- c(
  "EU: Italy" = "#4E79A7",
  "EU: Malta" = "#76B7B2",
  "North Africa: Libya" = "#E15759",
  "North Africa: Tunisia" = "#F28E2B"
)

panel_theme <- theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = panel_sea_fill, colour = NA),
    panel.border = element_rect(fill = NA, colour = "#D6D6D0",
                                linewidth = 0.35),
    plot.title.position = "plot",
    plot.title = element_text(size = 13, face = "bold",
                              colour = "#202020",
                              margin = margin(b = 3)),
    plot.caption.position = "plot",
    plot.caption = element_text(size = 8.5, colour = "#666666",
                                hjust = 0, margin = margin(t = 4)),
    legend.position = "right",
    legend.background = element_blank(),
    legend.box.background = element_rect(
      fill = scales::alpha("white", 0.96), colour = "#D2D2CC",
      linewidth = 0.25),
    legend.box = "vertical",
    legend.box.just = "center",
    legend.box.spacing = unit(1.5, "pt"),
    legend.box.margin = margin(3, 3, 3, 3),
    legend.key = element_rect(fill = scales::alpha("white", 0), colour = NA),
    legend.margin = margin(2, 2, 2, 2),
    legend.text = element_text(size = 11, colour = "#2B2B2B"),
    legend.title = element_text(size = 11, face = "bold",
                                colour = "#2B2B2B"),
    legend.key.size = unit(0.42, "cm"),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_country_labels <- country_labels |>
  mutate(
    lon = case_when(
      label == "Lampedusa" ~ 12.35,
      label == "Malta" ~ 14.85,
      TRUE ~ lon
    ),
    lat = case_when(
      label == "Lampedusa" ~ 35.45,
      label == "Malta" ~ 35.80,
      TRUE ~ lat
    )
  )

label_layer <- shadowtext::geom_shadowtext(
  data = panel_country_labels,
  aes(x = lon, y = lat, label = label),
  colour = "#4F4F4B",
  bg.colour = "white",
  size = 3.2,
  fontface = "bold.italic"
)

p_united_panel <- ggplot() +
  geom_sf(data = world, fill = panel_land_fill, colour = panel_land_line,
          linewidth = 0.18) +
  geom_sf(data = core_sea, fill = panel_corridor, colour = panel_corridor,
          alpha = 0.16, linewidth = 0.35) +
  geom_point(data = df_united_plot,
             aes(x = longitude, y = latitude, size = n_deaths,
                 colour = in_corridor),
             alpha = 0.55, shape = 16) +
  scale_colour_manual(values = panel_incident_colours,
                      name = "Shipwrecks") +
  scale_size_area(max_size = 4.4, name = "Shipwrecks",
                  breaks = c(1, 10, 50, 200)) +
  guides(colour = guide_legend(order = 1, title = "Shipwrecks",
           override.aes = list(size = 3, alpha = 0.9)),
         size   = guide_legend(order = 2, title = NULL,
           override.aes = list(colour = "#2B2B2B", alpha = 0.65))) +
  label_layer +
  coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
  labs(
    title = "(a) Shipwrecks and area of analysis",
    caption = sprintf("Source: UNITED List of Refugee Deaths, 2014 to %s; Natural Earth base map.",
                      format(MAP_END_DATE, "%B %Y")),
    x = NULL, y = NULL
  ) +
  panel_theme

p_sar_panel <- ggplot() +
  geom_sf(data = world, fill = panel_land_fill, colour = panel_land_line,
          linewidth = 0.18) +
  geom_sf(data = zones_sea, aes(fill = zone_label), alpha = 0.22,
          colour = "#2F6FAF", linewidth = 0.32) +
  scale_fill_manual(
    values = panel_zone_colours,
    labels = c("EU: Italy", "EU: Malta",
               "North Africa:\nLibya", "North Africa:\nTunisia"),
    name = "SAR zone"
  ) +
  label_layer +
  coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
  labs(
    title = "(b) Search-and-rescue responsibility zones",
    caption = "Source: IMO GISIS Global SAR Plan boundaries; Natural Earth base map.",
    x = NULL, y = NULL
  ) +
  panel_theme

panel_united_sar <- p_united_panel + p_sar_panel +
  plot_layout(widths = c(1, 1), guides = "keep")

panel_united_sar_boxed <- cowplot::ggdraw(panel_united_sar) +
  cowplot::draw_grob(grid::rectGrob(
    gp = grid::gpar(col = "#202020", fill = NA, lwd = 1.2)
  ))

ggsave(fig_path("03_build", "fig-sea-zones-united-sar-panel.png"),
       panel_united_sar_boxed, width = 11, height = 4.65, dpi = 450,
       device = ragg::agg_png, bg = "white")

# ── 5. Save intermediates ───────────────────────────────────────────────────
saveRDS(core_poly, file.path(BASE_DIR, "data", "processed", "core_corridor.RDS"))
saveRDS(zones_sea, file.path(BASE_DIR, "data", "processed", "sar_zones.RDS"))
saveRDS(zones_in_corridor,
        file.path(BASE_DIR, "data", "processed", "sar_zones_in_corridor.RDS"))
saveRDS(frontex_op,
        file.path(BASE_DIR, "data", "processed", "frontex_op_areas.RDS"))
