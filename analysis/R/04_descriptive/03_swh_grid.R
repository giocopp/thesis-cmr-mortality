# 12_descriptives_swh_grid.R
# ====================
# SWH spatial-homogeneity diagnostic. For each ERA5 ocean cell inside
# the core corridor polygon, we extract the daily SWH time series
# (2014-2025) and compute:
#   (1) pairwise cell-to-cell correlations — how much do different cells
#       within the polygon move together across days?
#   (2) each cell's correlation with the polygon-average daily SWH —
#       how representative is the polygon mean of each individual cell?
#
# Input:  data/processed/core_corridor.RDS
#         data/raw/era5/era5_daily_cmr_wave_[YEAR].nc
# Output: output/tables/12_descriptives_swh_grid.txt
#         output/figures/12_descriptives_swh_grid_hist.png
#         output/figures/12_descriptives_swh_grid_map.png

library(sf)
library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(rnaturalearth)
library(patchwork)
library(fixest)

sf_use_s2(FALSE)

BASE_DIR <- here::here()
ERA5_DIR <- file.path(BASE_DIR, "data", "raw", "era5")
YEARS    <- 2014:2025

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("12_descriptives_swh_grid: spatial homogeneity of SWH in core corridor\n")
cat("============================================================\n\n")

# ── 1. Load polygon ─────────────────────────────────────────
core_poly <- readRDS(file.path(BASE_DIR, "data", "processed",
                                 "core_corridor.RDS"))
cat("Core corridor area:", round(as.numeric(st_area(core_poly)) / 1e6),
    "km^2\n")

# ── 2. Build mask from ERA5 grid ───────────────────────────
nc0      <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
wave_lon <- ncvar_get(nc0, "longitude")
wave_lat <- ncvar_get(nc0, "latitude")
swh1     <- ncvar_get(nc0, "swh", start = c(1, 1, 1), count = c(-1, -1, 1))
nc_close(nc0)

grid_df <- expand.grid(lon = wave_lon, lat = wave_lat)
grid_pts <- st_as_sf(grid_df, coords = c("lon", "lat"), crs = 4326)

ocean   <- !is.na(as.vector(swh1))
inside  <- st_intersects(grid_pts, core_poly, sparse = FALSE)[, 1]
mask_vec <- ocean & inside

cell_idx <- which(mask_vec)                # flat indices
cell_lon <- grid_df$lon[cell_idx]
cell_lat <- grid_df$lat[cell_idx]
n_cells  <- length(cell_idx)

cat("ERA5 grid:", length(wave_lon), "lon x", length(wave_lat), "lat\n")
cat("Cells in mask (ocean + inside polygon):", n_cells, "\n\n")

# ── 3. Extract cell × day SWH matrix for all years ─────────
cat("--- Extracting cell x day SWH matrix ---\n")

get_nc_dates <- function(nc) {
  tn <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  tv <- ncvar_get(nc, tn)
  as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC"))
}

extract_year <- function(yr) {
  f <- file.path(ERA5_DIR, sprintf("era5_daily_cmr_wave_%d.nc", yr))
  if (!file.exists(f)) {
    cat("  ", yr, ": file not found, skipping\n"); return(NULL)
  }
  nc    <- nc_open(f)
  dates <- get_nc_dates(nc)
  swh3  <- ncvar_get(nc, "swh")  # [lon, lat, time]
  nc_close(nc)
  # Flatten each day and keep only the masked cells: result is cells x days
  mat <- vapply(seq_along(dates),
                function(t) as.vector(swh3[, , t])[cell_idx],
                numeric(n_cells))
  cat(sprintf("  %d: %d days x %d cells\n", yr, length(dates), n_cells))
  list(dates = dates, mat = mat)
}

chunks <- map(YEARS, extract_year)
chunks <- chunks[!vapply(chunks, is.null, logical(1))]

all_dates <- do.call(c, lapply(chunks, \(x) x$dates))
all_mat   <- do.call(cbind, lapply(chunks, \(x) x$mat))

cat(sprintf("\nFull matrix: %d cells x %d days\n",
            nrow(all_mat), ncol(all_mat)))

# ── 4. Pairwise cell-to-cell correlations ──────────────────
cat("\n--- Pairwise cell-to-cell correlations ---\n")
cat("(rows of all_mat are cells, so we correlate on the transpose)\n\n")

cor_mat <- cor(t(all_mat), use = "pairwise.complete.obs")
upper   <- cor_mat[upper.tri(cor_mat)]

# ── 5. Each cell's correlation with the polygon-average ────
polygon_avg  <- colMeans(all_mat, na.rm = TRUE)
cor_with_avg <- apply(all_mat, 1,
                       \(cell) cor(cell, polygon_avg,
                                    use = "pairwise.complete.obs"))

# ── 6. Pairwise correlation vs great-circle distance ───────
cell_sf <- st_as_sf(data.frame(lon = cell_lon, lat = cell_lat),
                     coords = c("lon", "lat"), crs = 4326)
dist_m  <- as.matrix(st_distance(cell_sf))
dist_km <- dist_m[upper.tri(dist_m)] / 1000

dist_decay <- cor(dist_km, upper, use = "pairwise.complete.obs")

# ── 6b. Incident-density-weighted SWH (sanity check) ───────
# Build STATIC weights from the overall spatial distribution of CMR
# incidents (not day-specific — that would be endogenous). Compute a
# weighted daily SWH and correlate it against the unweighted polygon
# mean. Correlation is descriptive; the model check in 6c is the direct
# test of whether weighted SWH changes the primary estimates.
cat("\n--- Incident-density-weighted SWH sanity check ---\n")

iom <- iom_incidents(
    incident_types = c("incident", "split incident"),
    spatial        = "all_cmr"
  ) |>
  dplyr::rename(dead = dead_missing) |>
  tidyr::drop_na(lon, lat)

cat(sprintf("  IOM incidents loaded: %d\n", nrow(iom)))

# Snap each incident to its nearest ERA5 cell (regular 0.5 deg grid).
iom$lon_idx <- vapply(iom$lon, \(x) which.min(abs(wave_lon - x)), integer(1))
iom$lat_idx <- vapply(iom$lat, \(x) which.min(abs(wave_lat - x)), integer(1))
iom$flat    <- (iom$lat_idx - 1L) * length(wave_lon) + iom$lon_idx

# Keep only incidents whose nearest cell is inside our analytical mask.
in_mask <- iom$flat %in% cell_idx
cat(sprintf("  Incidents snapping inside the polygon mask: %d (%.1f%%)\n",
            sum(in_mask), 100 * sum(in_mask) / nrow(iom)))
iom_in <- iom[in_mask, ]

# Build per-cell counts and per-cell death totals.
cell_to_row  <- setNames(seq_along(cell_idx), cell_idx)
weight_count <- numeric(n_cells)
weight_death <- numeric(n_cells)
for (i in seq_len(nrow(iom_in))) {
  r <- cell_to_row[[as.character(iom_in$flat[i])]]
  weight_count[r] <- weight_count[r] + 1
  weight_death[r] <- weight_death[r] + iom_in$dead[i]
}

cat(sprintf("  Cells with >= 1 incident: %d / %d\n",
            sum(weight_count > 0), n_cells))
# Concentration: how few cells account for 50/80/90% of deaths?
sorted_d <- sort(weight_death, decreasing = TRUE)
cumsh <- cumsum(sorted_d) / sum(sorted_d)
cat(sprintf("  Cells accounting for 50%% of deaths: %d\n",
            which(cumsh >= 0.5)[1]))
cat(sprintf("  Cells accounting for 80%% of deaths: %d\n",
            which(cumsh >= 0.8)[1]))
cat(sprintf("  Cells accounting for 90%% of deaths: %d\n",
            which(cumsh >= 0.9)[1]))

# Normalise weights to sum to 1 then compute weighted daily SWH series.
w_count_norm <- weight_count / sum(weight_count)
w_death_norm <- weight_death / sum(weight_death)

swh_weighted_count <- vapply(seq_len(ncol(all_mat)), \(t)
  sum(all_mat[, t] * w_count_norm, na.rm = TRUE) /
  sum(w_count_norm[!is.na(all_mat[, t])]),
  numeric(1))
swh_weighted_death <- vapply(seq_len(ncol(all_mat)), \(t)
  sum(all_mat[, t] * w_death_norm, na.rm = TRUE) /
  sum(w_death_norm[!is.na(all_mat[, t])]),
  numeric(1))

cor_polymean_vs_countw <- cor(polygon_avg, swh_weighted_count,
                                use = "pairwise.complete.obs")
cor_polymean_vs_deathw <- cor(polygon_avg, swh_weighted_death,
                                use = "pairwise.complete.obs")
cor_countw_vs_deathw   <- cor(swh_weighted_count, swh_weighted_death,
                                use = "pairwise.complete.obs")

cat(sprintf("\n  cor(polygon-mean SWH, incident-count-weighted SWH): %.4f\n",
            cor_polymean_vs_countw))
cat(sprintf("  cor(polygon-mean SWH, death-count-weighted SWH):    %.4f\n",
            cor_polymean_vs_deathw))
cat(sprintf("  cor(count-weighted, death-weighted):                %.4f\n",
            cor_countw_vs_deathw))

# ── 6c. Re-estimate the primary model with incident-weighted SWH ─────
# The correlation check above tells us how close the weighted and
# unweighted SWH series are. This model check asks whether the paper's
# primary SWH x post-MoU estimates change when the SWH regressor is
# computed as a static incident-count-weighted cell average instead of
# a simple polygon mean. We use static all-period IOM incident counts
# by ERA5 cell as weights; using day-specific weights would condition on
# realized incidents and be endogenous to the outcome process.
cat("\n--- Weighted-SWH model check ---\n")

make_prev5 <- function(x) zoo::rollmeanr(dplyr::lag(x, 1), k = 5, fill = NA)

weighted_weather <- tibble(
  date               = all_dates,
  swh_polygon_grid   = polygon_avg,
  swh_weighted_count = swh_weighted_count,
  swh_weighted_death = swh_weighted_death
) |>
  arrange(date) |>
  mutate(
    swh_polygon_grid_prev5days   = make_prev5(swh_polygon_grid),
    swh_weighted_count_prev5days = make_prev5(swh_weighted_count),
    swh_weighted_death_prev5days = make_prev5(swh_weighted_death)
  )

main_panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                                "daily_panel_complete.RDS"))
iom_daily_primary    <- build_iom_daily()
united_daily_primary <- build_united_daily()

model_panel <- main_panel |>
  left_join(iom_daily_primary |> rename(n_dead_iom = n_dead_missing),
            by = "date") |>
  left_join(united_daily_primary, by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) |>
  add_crossing_exposure() |>
  arrange(date) |>
  mutate(
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    unit = 1L,
    month_year_fac = factor(month_year)
  ) |>
  left_join(weighted_weather, by = "date") |>
  filter(!is.na(lc_lag14),
         !is.na(swh_prev5days),
         !is.na(swh_weighted_count_prev5days),
         !is.na(swh_weighted_death_prev5days))

extract_shift <- function(model, x, source, family, series_label) {
  ct <- coeftable(model, vcov = NW(14))
  co <- coef(model)
  V  <- vcov(model, vcov = NW(14))
  int_name <- grep(paste0("(^", x, ":post_mou$)|(^post_mou:", x, "$)"),
                   names(co), value = TRUE)
  if (length(int_name) != 1L) {
    stop("Could not uniquely identify interaction for ", x, ": ",
         paste(int_name, collapse = ", "))
  }
  b_pre   <- unname(co[x])
  b_shift <- unname(co[int_name])
  b_post  <- b_pre + b_shift
  v_post  <- V[x, x] + V[int_name, int_name] + 2 * V[x, int_name]
  se_post <- sqrt(v_post)
  tibble(
    source = source,
    family = family,
    swh_series = series_label,
    n_obs = nobs(model),
    b_pre = b_pre,
    se_pre = ct[x, 2],
    p_pre = 2 * pnorm(-abs(b_pre / ct[x, 2])),
    b_shift = b_shift,
    se_shift = ct[int_name, 2],
    p_shift = 2 * pnorm(-abs(b_shift / ct[int_name, 2])),
    b_post = b_post,
    se_post = se_post,
    p_post = 2 * pnorm(-abs(b_post / se_post))
  )
}

fit_weighted_set <- function(dep, source, x, series_label) {
  f <- as.formula(sprintf("%s ~ %s + %s:post_mou | month_year_fac",
                          dep, x, x))
  nb <- fenegbin(f, data = model_panel, vcov = NW(14),
                 panel.id = ~unit + date)
  po <- fepois(f, data = model_panel, vcov = NW(14),
               panel.id = ~unit + date)
  bind_rows(
    extract_shift(nb, x, source, "NegBin", series_label),
    extract_shift(po, x, source, "Poisson", series_label)
  )
}

swh_specs <- tibble(
  x = c("swh_prev5days",
        "swh_weighted_count_prev5days",
        "swh_weighted_death_prev5days"),
  label = c("Polygon mean",
            "Incident-count weighted",
            "Death-count weighted")
)

weighted_model_tbl <- bind_rows(
  purrr::pmap_dfr(swh_specs, \(x, label)
    fit_weighted_set("n_dead_united", "UNITED", x, label)),
  purrr::pmap_dfr(swh_specs, \(x, label)
    fit_weighted_set("n_dead_iom", "IOM", x, label))
)

weighted_delta_tbl <- weighted_model_tbl |>
  select(source, family, swh_series, n_obs, b_pre, b_shift, b_post) |>
  group_by(source, family) |>
  mutate(
    d_b_pre   = b_pre   - b_pre[swh_series == "Polygon mean"],
    d_b_shift = b_shift - b_shift[swh_series == "Polygon mean"],
    d_b_post  = b_post  - b_post[swh_series == "Polygon mean"]
  ) |>
  ungroup() |>
  filter(swh_series != "Polygon mean")

weighted_csv <- tbl_path("04_descriptive", "03_swh_weighted_model.csv")
write.csv(weighted_model_tbl, weighted_csv, row.names = FALSE)
cat(sprintf("  Weighted model comparison saved: %s\n", weighted_csv))

# ── 7. Write table ─────────────────────────────────────────
sink_file <- tbl_path("04_descriptive", "03_swh_grid.txt")
sink(sink_file)

cat("SWH spatial-homogeneity diagnostic — core corridor polygon\n")
cat("==========================================================\n\n")
cat(sprintf("Period: %s to %s (%d days)\n",
            min(all_dates), max(all_dates), ncol(all_mat)))
cat(sprintf("Cells inside polygon (ERA5 ocean): %d\n", n_cells))
cat(sprintf("Number of pairwise cell-cell pairs: %d\n\n",
            length(upper)))

cat("== Pairwise cell-to-cell correlation ==\n")
cat(sprintf("  mean   : %.4f\n", mean(upper,   na.rm = TRUE)))
cat(sprintf("  median : %.4f\n", median(upper, na.rm = TRUE)))
cat(sprintf("  sd     : %.4f\n", sd(upper,     na.rm = TRUE)))
cat(sprintf("  min    : %.4f\n", min(upper,    na.rm = TRUE)))
cat(sprintf("  10th pct: %.4f\n", quantile(upper, 0.10, na.rm = TRUE)))
cat(sprintf("  25th pct: %.4f\n", quantile(upper, 0.25, na.rm = TRUE)))
cat(sprintf("  75th pct: %.4f\n", quantile(upper, 0.75, na.rm = TRUE)))
cat(sprintf("  max    : %.4f\n", max(upper,    na.rm = TRUE)))
cat("\n  Fraction of pairs with correlation above threshold:\n")
for (thr in c(0.99, 0.95, 0.90, 0.80, 0.70)) {
  cat(sprintf("    > %.2f : %5.1f%%\n",
              thr, 100 * mean(upper > thr, na.rm = TRUE)))
}

cat("\n== Each cell's correlation with the polygon-average ==\n")
cat(sprintf("  mean   : %.4f\n", mean(cor_with_avg, na.rm = TRUE)))
cat(sprintf("  median : %.4f\n", median(cor_with_avg, na.rm = TRUE)))
cat(sprintf("  min    : %.4f\n", min(cor_with_avg, na.rm = TRUE)))
cat(sprintf("  10th pct: %.4f\n", quantile(cor_with_avg, 0.10, na.rm = TRUE)))

cat("\n== Pairwise correlation vs distance (continuous) ==\n")
cat(sprintf("  Pearson r(distance_km, cell-cell correlation) = %+.4f\n", dist_decay))
cat(sprintf("  Distance range: %.0f - %.0f km (mean %.0f)\n",
            min(dist_km), max(dist_km), mean(dist_km)))

cat("\n== Incident-density-weighted SWH sanity check ==\n")
cat(sprintf("  Incidents snapped to polygon cells: %d\n", sum(in_mask)))
cat(sprintf("  Cells with >= 1 incident: %d / %d\n",
            sum(weight_count > 0), n_cells))
cat(sprintf("  Cells accounting for 50%% of deaths: %d\n",
            which(cumsh >= 0.5)[1]))
cat(sprintf("  Cells accounting for 80%% of deaths: %d\n",
            which(cumsh >= 0.8)[1]))
cat(sprintf("  Cells accounting for 90%% of deaths: %d\n",
            which(cumsh >= 0.9)[1]))
cat("\n  Correlations:\n")
cat(sprintf("    polygon-mean vs incident-count weighted: %.4f\n",
            cor_polymean_vs_countw))
cat(sprintf("    polygon-mean vs death-count weighted:    %.4f\n",
            cor_polymean_vs_deathw))
cat(sprintf("    count-weighted vs death-weighted:        %.4f\n",
            cor_countw_vs_deathw))
cat("\n  Interpretation:\n")
cat("    The weighted series are highly correlated with the polygon mean\n")
cat("    but not identical, so the next check re-estimates the primary\n")
cat("    count model using incident- and death-weighted SWH.\n")

cat("\n== Weighted-SWH primary model check ==\n")
cat("  Static weights are based on all IOM CMR incidents snapped to ERA5 cells.\n")
cat("  The estimation sample and specification match the primary count model:\n")
cat("  deaths ~ SWH_{t-1:t-5} + SWH_{t-1:t-5} x post_MoU | month-year FE,\n")
cat("  estimated by NegBin and Poisson QMLE with NW(14) SEs.\n\n")
cat(sprintf("  Shared model-check rows before FE drops: %d days\n\n",
            nrow(model_panel)))
print(weighted_model_tbl |>
        mutate(across(c(b_pre, se_pre, p_pre,
                        b_shift, se_shift, p_shift,
                        b_post, se_post, p_post), \(x) round(x, 4))),
      n = Inf, width = Inf)
cat("\n  Change relative to polygon mean:\n")
print(weighted_delta_tbl |>
        mutate(across(c(b_pre, b_shift, b_post,
                        d_b_pre, d_b_shift, d_b_post), \(x) round(x, 4))),
      n = Inf, width = Inf)
cat(sprintf("\n  CSV: %s\n", weighted_csv))

sink()
cat(sprintf("Saved: %s\n", sink_file))

# ── 8. Histogram of pairwise correlations ─────────────────
cat("\n--- Plotting histogram ---\n")

p_hist <- ggplot(data.frame(r = upper), aes(x = r)) +
  geom_histogram(bins = 60, fill = "#2166AC", colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = mean(upper, na.rm = TRUE),
             colour = "#B2182B", linetype = "dashed") +
  labs(
    title = "Pairwise SWH correlation between ERA5 cells in core corridor",
    subtitle = sprintf("n cells = %d, n unique pairs = %d, mean r = %.3f",
                        n_cells, length(upper), mean(upper, na.rm = TRUE)),
    x = "Pearson correlation (across days 2014-2025)",
    y = "Number of cell pairs"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

fig_hist <- fig_path("04_descriptive", "03_swh_grid_hist.png")
ggsave(fig_hist, p_hist, width = 9, height = 5, dpi = 200)
cat("Saved:", fig_hist, "\n")

# ── 9. Real-map plots with land basemap + incident overlay ─
cat("\n--- Plotting correlation + incident-density maps ---\n")

# Base data: land polygons, country labels, extent just past the polygon.
world <- ne_countries(scale = "medium", returnclass = "sf")
bb    <- st_bbox(core_poly)
pad   <- 0.6
MAP_XLIM <- c(bb["xmin"] - pad, bb["xmax"] + pad)
MAP_YLIM <- c(bb["ymin"] - pad, bb["ymax"] + pad)

country_labels <- data.frame(
  label = c("Sicily", "Tunisia", "Libya", "Malta", "Calabria", "Algeria"),
  lon   = c(14.2,     9.5,        15.0,    14.4,    16.4,       8.0),
  lat   = c(37.4,     35.0,       30.8,    35.7,    38.8,       35.5)
)

# Incidents for the overlay (only those snapped inside the mask).
iom_overlay <- iom_in |> dplyr::select(lon, lat, dead)

cell_df   <- data.frame(lon = cell_lon, lat = cell_lat,
                        cor_with_avg = cor_with_avg)
weight_df <- data.frame(lon = cell_lon, lat = cell_lat,
                        incidents = weight_count,
                        deaths = weight_death)

base_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(size = 8.5, colour = "grey40"),
    axis.text        = element_text(size = 7, colour = "grey50"),
    legend.position  = "right",
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 8)
  )

base_layer <- function(p) {
  p +
    geom_sf(data = world, fill = "grey93", colour = "grey70", linewidth = 0.2) +
    geom_text(data = country_labels,
              aes(x = lon, y = lat, label = label),
              colour = "grey35", size = 2.7, fontface = "italic") +
    coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
    base_theme
}

# (A) Cell correlation with polygon-average, with incident overlay
p_map <- base_layer(
  ggplot() +
    geom_tile(data = cell_df,
              aes(x = lon, y = lat, fill = cor_with_avg),
              width = 0.5, height = 0.5, alpha = 0.85) +
    geom_point(data = iom_overlay,
               aes(x = lon, y = lat, size = dead),
               shape = 21, fill = NA, colour = "black",
               stroke = 0.25, alpha = 0.35) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(min(cor_with_avg, na.rm = TRUE), 1),
      name   = "r vs\npolygon\nmean"
    ) +
    scale_size_continuous(range = c(0.3, 4.5), name = "Incident\nsize",
                          breaks = c(1, 10, 50, 200, 500),
                          guide = guide_legend(override.aes =
                            list(colour = "black", alpha = 0.7)))
) +
  labs(
    title    = "Per-cell correlation with polygon-mean SWH",
    subtitle = "Incident locations overlaid (hollow circles, sized by deaths)",
    x = NULL, y = NULL
  )

# (B) Incident-density heat map, with raw incident points
p_weights <- base_layer(
  ggplot() +
    geom_tile(data = weight_df,
              aes(x = lon, y = lat, fill = incidents),
              width = 0.5, height = 0.5, alpha = 0.85) +
    geom_point(data = iom_overlay,
               aes(x = lon, y = lat, size = dead),
               shape = 21, fill = NA, colour = "grey10",
               stroke = 0.2, alpha = 0.35) +
    scale_fill_viridis_c(
      option = "inferno",
      trans  = "sqrt",
      name   = "Incidents\nper cell\n(sqrt)",
      na.value = "transparent"
    ) +
    scale_size_continuous(range = c(0.3, 4.5), name = "Incident\nsize",
                          breaks = c(1, 10, 50, 200, 500),
                          guide = guide_legend(override.aes =
                            list(colour = "grey10", alpha = 0.7)))
) +
  labs(
    title    = "Static incident-density weights (incidents per ERA5 cell)",
    subtitle = sprintf("Total incidents snapped inside mask: %d across %d cells",
                        nrow(iom_in), sum(weight_count > 0)),
    x = NULL, y = NULL
  )

# (C) Side-by-side panel for comparison
panel <- p_map + p_weights + plot_layout(ncol = 2)
fig_panel <- fig_path("04_descriptive", "03_swh_grid_panel.png")
ggsave(fig_panel, panel, width = 16, height = 7, dpi = 220)
cat("Saved:", fig_panel, "\n")

# ── 10. LaTeX table for tab:appx-swh-weighted ─────────────────
cat("\n--- 10. Writing LaTeX table ---\n")

fbse <- function(b, se) sprintf("$%+.3f$ (%.3f)", b, se)
fdif <- function(d)     sprintf("$%+.3f$", d)

get_row <- function(src, fam, series) {
  r <- weighted_model_tbl[weighted_model_tbl$source     == src &
                           weighted_model_tbl$family     == fam &
                           weighted_model_tbl$swh_series == series, ]
  list(b = r$b_shift, se = r$se_shift)
}
rows <- list(
  c(src = "UNITED", fam = "NegBin"),
  c(src = "UNITED", fam = "Poisson"),
  c(src = "IOM",    fam = "NegBin"),
  c(src = "IOM",    fam = "Poisson")
)

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[h!]")
add("\\centering")
add("\\small")
add("\\caption{Primary slope-shift estimates using polygon-mean and incident-count-weighted SWH.}")
add("\\label{tab:appx-swh-weighted}")
add("\\begin{tabular}{llccc}")
add("\\hline")
add("Source & Family & Polygon mean $\\beta_3$ & Incident-weighted $\\beta_3$ & Change \\\\")
add("\\hline")
for (r in rows) {
  poly  <- get_row(r["src"], r["fam"], "Polygon mean")
  inc   <- get_row(r["src"], r["fam"], "Incident-count weighted")
  delta <- inc$b - poly$b
  add(sprintf("%-7s & %-8s & %-19s & %-19s & %s \\\\",
              r["src"], r["fam"],
              fbse(poly$b, poly$se),
              fbse(inc$b,  inc$se),
              fdif(delta)))
}
add("\\hline")
add("\\multicolumn{5}{l}{\\footnotesize Standard errors in parentheses. Month-year fixed effects and Newey-West SEs (lag 14).} \\\\")
add("\\multicolumn{5}{l}{\\footnotesize Incident weights are fixed over the full sample and do not vary by day.}")
add("\\end{tabular}")
add("\\end{table}")
out_swhw <- tbl_path("04_descriptive", "03_swh_weighted_model.tex")
writeLines(L, out_swhw)
cat(sprintf("  Saved: %s\n", out_swhw))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
