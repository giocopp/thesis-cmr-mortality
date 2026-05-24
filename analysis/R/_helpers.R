# ── Shared helpers: builders + path/spatial utilities ────────────────────────

source(here::here("analysis", "R", "_constants.R"))

# ── Output paths ─────────────────────────────────────────────────────────────
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

# ── Spatial filter: corridor polygon ─────────────────────────────────────────
filter_corridor <- function(d, coords = c("lon", "lat"),
                            base_dir = here::here()) {
  core_poly <- readRDS(file.path(base_dir, "data", "processed", "core_corridor.RDS"))
  prev_s2 <- sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(prev_s2), add = TRUE)
  pts    <- sf::st_as_sf(d, coords = coords, crs = 4326, remove = FALSE)
  inside <- lengths(sf::st_within(pts, core_poly)) > 0
  dplyr::filter(d, inside)
}

# ── IOM incidents ────────────────────────────────────────────────────────────
iom_incidents <- function(
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

  d
}

build_iom_daily <- function(...) {
  iom_incidents(...) |>
    dplyr::group_by(date) |>
    dplyr::summarise(n_dead_missing = sum(dead_missing), .groups = "drop") |>
    dplyr::arrange(date)
}

# ── UNITED incidents ─────────────────────────────────────────────────────────
united_incidents <- function(
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

  d
}

build_united_daily <- function(...) {
  united_incidents(...) |>
    dplyr::group_by(date) |>
    dplyr::summarise(n_dead_united = sum(n_deaths, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(date)
}

# ── Crossing-exposure columns ────────────────────────────────────────────────
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
