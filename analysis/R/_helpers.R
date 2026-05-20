# ── Shared helpers ────────────────────────────────────────────────────────────
# Builders for IOM/UNITED daily series and small path/spatial utilities.
# Sourced by every analytical script.

source(here::here("analysis", "R", "_constants.R"))

# ── Output paths ──────────────────────────────────────────────────────────────
# fig_path("05_analysis", "01_primary.png") returns the full path under
# output/figures/05_analysis/ and creates the directory if needed.
fig_path <- function(section, file) {
  p <- here::here("output", "figures", section, file)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

tbl_path <- function(section, file) {
  p <- here::here("output", "tables", section, file)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

# ── Spatial: restrict to core corridor ────────────────────────────────────────
# Spatial filter to the corridor polygon. Expects non-missing coords; toggles
# off s2 for robust point-in-polygon on lon/lat.
filter_corridor <- function(d, coords = c("lon", "lat"),
                            base_dir = here::here()) {
  core_poly <- readRDS(file.path(base_dir, "data", "processed", "core_corridor.RDS"))
  prev_s2 <- sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(prev_s2), add = TRUE)
  pts    <- sf::st_as_sf(d, coords = coords, crs = 4326, remove = FALSE)
  inside <- lengths(sf::st_within(pts, core_poly)) > 0
  dplyr::filter(d, inside)
}

# ── IOM daily deaths ──────────────────────────────────────────────────────────
# Daily aggregate of dead+missing from raw IOM MMP. Defaults match the primary
# analytical spec: incident-only, corridor spatial join, drowning + mixed cause.
# Override incident_types/spatial/causes for robustness variants.
build_iom_daily <- function(
  incident_types = c("incident"),
  spatial        = c("central", "all_cmr"),
  causes         = c("sea", "all"),
  countries      = CMR_INCIDENT_COUNTRIES,
  base_dir       = here::here()
) {
  spatial <- match.arg(spatial)
  causes  <- match.arg(causes)

  d <- readRDS(file.path(base_dir, "data", "processed", "iom_mmp_incidents.RDS")) |>
    dplyr::filter(
      Route == "Central Mediterranean",
      tolower(`Incident Type`) %in% tolower(incident_types),
      `Country of Incident` %in% countries
    ) |>
    dplyr::transmute(
      date         = as.Date(incident_date_clean),
      dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE),
      cause_cat    = `Cause of death (category)`,
      lon          = as.numeric(Longitude),
      lat          = as.numeric(Latitude)
    ) |>
    tidyr::drop_na(date)

  if (causes == "sea") {
    d <- dplyr::filter(d, cause_cat %in% c("Drowning", "Mixed or unknown"))
  }

  if (spatial == "central") {
    d <- d |>
      tidyr::drop_na(lon, lat) |>
      filter_corridor(coords = c("lon", "lat"), base_dir = base_dir)
  }

  d |>
    dplyr::group_by(date) |>
    dplyr::summarise(n_dead_missing = sum(dead_missing), .groups = "drop") |>
    dplyr::arrange(date)
}

# ── UNITED daily deaths ───────────────────────────────────────────────────────
# Same spatial + cause logic as build_iom_daily() so the two sources are
# directly comparable. UNITED records open-sea deaths under "Mediterranean",
# so that label is added to the country list.
build_united_daily <- function(
  causes    = c("sea", "all"),
  spatial   = c("central", "all_cmr"),
  countries = c(CMR_INCIDENT_COUNTRIES, "Mediterranean"),
  base_dir  = here::here()
) {
  causes  <- match.arg(causes)
  spatial <- match.arg(spatial)

  d <- readRDS(file.path(base_dir, "data", "processed", "united_incidents.RDS")) |>
    dplyr::filter(country_of_death %in% countries) |>
    dplyr::mutate(date = as.Date(incident_date_clean))

  if (causes == "sea") {
    d <- dplyr::filter(d, manner_of_death %in% c("drowned", "other_unknown"))
  }

  if (spatial == "central") {
    d <- d |>
      dplyr::filter(!is.na(latitude), !is.na(longitude)) |>
      filter_corridor(coords = c("longitude", "latitude"), base_dir = base_dir)
  }

  d |>
    dplyr::group_by(date) |>
    dplyr::summarise(n_dead_united = sum(n_deaths, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(date)
}

# ── Crossing-exposure columns ─────────────────────────────────────────────────
# Adds living_crossings (frx_persons + lcg_tcg_pushbacks) and the source-specific
# attempts_iom / attempts_united denominators. Call after the build_*_daily()
# joins and the replace_na() step that fills n_dead_{iom,united}.
add_crossing_exposure <- function(df) {
  stopifnot(all(c("frx_persons", "lcg_tcg_pushbacks",
                  "n_dead_iom", "n_dead_united") %in% names(df)))
  dplyr::mutate(
    df,
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    attempts_iom     = living_crossings + n_dead_iom,
    attempts_united  = living_crossings + n_dead_united
  )
}
