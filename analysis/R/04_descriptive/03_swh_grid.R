# Appendix A: SWH spatial homogeneity + incident-weighted SWH robustness.
# Produces fig-appx-swh-grid-panel.png (cell correlation + incident weights)
# and tab-appx-swh-weighted.tex (re-estimation with weighted SWH).

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

# ── 1. Polygon + ERA5 grid mask ─────────────────────────────────────────────
core_poly <- readRDS(file.path(BASE_DIR, "data", "processed",
                                 "core_corridor.RDS"))

nc0      <- nc_open(file.path(ERA5_DIR, "era5_daily_cmr_wave_2014.nc"))
wave_lon <- ncvar_get(nc0, "longitude")
wave_lat <- ncvar_get(nc0, "latitude")
swh1     <- ncvar_get(nc0, "swh", start = c(1, 1, 1), count = c(-1, -1, 1))
nc_close(nc0)

grid_df  <- expand.grid(lon = wave_lon, lat = wave_lat)
grid_pts <- st_as_sf(grid_df, coords = c("lon", "lat"), crs = 4326)

ocean    <- !is.na(as.vector(swh1))
inside   <- st_intersects(grid_pts, core_poly, sparse = FALSE)[, 1]
mask_vec <- ocean & inside

cell_idx <- which(mask_vec)
cell_lon <- grid_df$lon[cell_idx]
cell_lat <- grid_df$lat[cell_idx]
n_cells  <- length(cell_idx)

# ── 2. Extract cell × day SWH matrix ────────────────────────────────────────
get_nc_dates <- function(nc) {
  tn <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  tv <- ncvar_get(nc, tn)
  as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC"))
}

extract_year <- function(yr) {
  f <- file.path(ERA5_DIR, sprintf("era5_daily_cmr_wave_%d.nc", yr))
  if (!file.exists(f)) return(NULL)
  nc    <- nc_open(f)
  dates <- get_nc_dates(nc)
  swh3  <- ncvar_get(nc, "swh")
  nc_close(nc)
  mat <- vapply(seq_along(dates),
                function(t) as.vector(swh3[, , t])[cell_idx],
                numeric(n_cells))
  list(dates = dates, mat = mat)
}

chunks    <- map(YEARS, extract_year)
chunks    <- chunks[!vapply(chunks, is.null, logical(1))]
all_dates <- do.call(c, lapply(chunks, \(x) x$dates))
all_mat   <- do.call(cbind, lapply(chunks, \(x) x$mat))

# ── 3. Per-cell correlation with polygon mean ───────────────────────────────
polygon_avg  <- colMeans(all_mat, na.rm = TRUE)
cor_with_avg <- apply(all_mat, 1,
                       \(cell) cor(cell, polygon_avg,
                                    use = "pairwise.complete.obs"))

# ── 4. Incident-density weights (static, all-period IOM) ────────────────────
iom <- iom_incidents(
    incident_types = c("incident", "split incident"),
    spatial        = "all_cmr"
  ) |>
  dplyr::rename(dead = dead_missing) |>
  tidyr::drop_na(lon, lat)

iom$lon_idx <- vapply(iom$lon, \(x) which.min(abs(wave_lon - x)), integer(1))
iom$lat_idx <- vapply(iom$lat, \(x) which.min(abs(wave_lat - x)), integer(1))
iom$flat    <- (iom$lat_idx - 1L) * length(wave_lon) + iom$lon_idx

in_mask <- iom$flat %in% cell_idx
iom_in  <- iom[in_mask, ]

cell_to_row  <- setNames(seq_along(cell_idx), cell_idx)
weight_count <- numeric(n_cells)
weight_death <- numeric(n_cells)
for (i in seq_len(nrow(iom_in))) {
  r <- cell_to_row[[as.character(iom_in$flat[i])]]
  weight_count[r] <- weight_count[r] + 1
  weight_death[r] <- weight_death[r] + iom_in$dead[i]
}

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

# ── 5. Re-estimate primary count model with weighted SWH ────────────────────
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

# ── 6. Paper figure: per-cell correlation + incident weights ────────────────
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

panel <- p_map + p_weights + plot_layout(ncol = 2)
ggsave(fig_path("04_descriptive", "fig-appx-swh-grid-panel.png"),
       panel, width = 16, height = 7, dpi = 220)

# ── 7. Paper table: weighted-SWH primary slope-shift ────────────────────────
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
writeLines(L, tbl_path("04_descriptive", "tab-appx-swh-weighted.tex"))
