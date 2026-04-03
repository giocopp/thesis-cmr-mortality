# 00_define_sea_zone.R
# ====================
# Define the CMR Core Corridor sea zone and generate map visualization.
# Sea zone is used ONLY for wave height calculation (computed in data/scripts/01).
# ALL CMR incidents are included in the analysis regardless of location.
#
# Core Corridor zone:
#   Follows Tunisian/Libyan coastlines, diagonal Cap Bon → W Sicily,
#   Sicily south coast, east edge at 16E.
#
# Output:
#   output/figures/sea_zone_core_map.png

Sys.setlocale("LC_TIME", "en_US.UTF-8")

library(tidyverse)
library(sf)
library(rnaturalearth)
library(lubridate)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")

world <- ne_countries(scale = "medium", returnclass = "sf")
land  <- st_union(world)

# Helper: build sea zone from outer polygon
build_sea_zone <- function(outer_coords) {
  outer_poly <- st_sfc(st_polygon(list(outer_coords)), crs = 4326)
  sea_raw    <- st_difference(outer_poly, land)
  sea_proj   <- st_transform(sea_raw, 3857)
  sea_buf    <- st_buffer(sea_proj, dist = 5000)
  sea_back   <- st_transform(sea_buf, 4326)
  sea_vis    <- st_intersection(sea_back, outer_poly)
  list(outer_poly = outer_poly, sea_vis = sea_vis)
}

# ============================================================
# 1. Zone A: Core Corridor
# ============================================================
cat("--- Zone A: Core Corridor ---\n")

coords_a <- matrix(c(
    17, 30.0,
  15.1, 36.7,
  12.4, 37.8,
  11.0, 37.1,
   9.0, 34.0,
   9.0, 31.0,
    17, 30.0
), ncol = 2, byrow = TRUE)

zone_a <- build_sea_zone(coords_a)
cat("  Sea zone polygon built (not saved — computed inline where needed)\n")

# ============================================================
# 2. Load ALL CMR incidents (no zone filtering)
# ============================================================
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

cat(sprintf("\nTotal CMR incidents: %d\n", nrow(df_cmr)))
cat(sprintf("Total dead + missing: %.0f\n", sum(df_cmr$dead_missing)))

# ============================================================
# 4. Plot helper
# ============================================================
make_zone_plot <- function(df, sea_vis, title_label, xlims, ylims,
                           country_labels) {

  colour_vals <- c("Pre-MoU" = "#D4820E", "Post-MoU" = "#D32F2F")

  p <- ggplot() +
    geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.3) +
    geom_sf(data = sea_vis, fill = "#2166AC", alpha = 0.15,
            colour = "#2166AC", linewidth = 0.5) +
    geom_point(data = df, aes(x = lon, y = lat, size = dead_missing,
                               colour = post_mou), alpha = 0.5, shape = 16) +
    scale_colour_manual(values = colour_vals, name = "Period") +
    scale_size_continuous(range = c(0.5, 5), name = "Dead + missing",
                          breaks = c(1, 10, 50, 200)) +
    guides(colour = guide_legend(order = 1,
             override.aes = list(size = 3, alpha = 0.7))) +
    coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) +
    labs(
      title = sprintf("CMR incidents and sea conditions zone: %s (2014-2025)", title_label),
      subtitle = sprintf("All CMR incidents: %d | Sea zone used for wave height calculation only",
                          nrow(df)),
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")

  for (i in seq_len(nrow(country_labels))) {
    cl <- country_labels[i, ]
    p <- p + annotate("text", x = cl$x, y = cl$y, label = cl$label,
                       size = cl$size, fontface = "italic")
  }

  p
}

# ============================================================
# 5. Generate maps
# ============================================================
cat("\n--- Generating maps ---\n")

labels_a <- tibble(
  x     = c(12.6,  15.0,  9.5, 14.5),
  y     = c(35.6,  32.0, 36.5, 38.5),
  label = c("Lampedusa", "Libya", "Tunisia", "Sicily"),
  size  = c(2.5, 3, 3, 3)
)

p_a <- make_zone_plot(df_cmr, zone_a$sea_vis,
                       "Core Corridor", xlims = c(5, 22), ylims = c(29, 40),
                       country_labels = labels_a)
ggsave(file.path(BASE_DIR, "output", "figures", "sea_zone_core_map.png"),
       p_a, width = 12, height = 7, dpi = 200)
cat("Saved: output/figures/sea_zone_core_map.png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
