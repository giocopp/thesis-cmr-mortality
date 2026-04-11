# _helpers.R
# ==========
# Shared helpers for the analysis pipeline. Not a numbered pipeline step —
# sourced by other scripts.
#
# build_iom_daily()
# -----------------
# Builds a daily aggregate of dead+missing from raw IOM MMP with configurable
# filters. Analytical scripts source this file and call build_iom_daily() with
# the filter combination they want to test — the defaults match the primary
# analytical spec of 05_reduced_form_primary.R.
#
# Switch between sensitivity variants by changing one argument, e.g.:
#   build_iom_daily()                                    # primary
#   build_iom_daily(spatial = "all_cmr")                 # drop spatial filter
#   build_iom_daily(causes  = "sea")                     # sea-only causes
#   build_iom_daily(incident_types = c("incident", "split incident"),
#                   spatial = "all_cmr")                 # broad/descriptive
#
# Arguments
#   incident_types : character vector of IOM `Incident Type` values to keep.
#                    Default c("incident", "split incident"). Split incidents
#                    are IOM's decomposition of one reported event across
#                    multiple dates — dropping them loses ~10% of sample
#                    rows and is more conservative, but the primary spec
#                    keeps them to match the descriptive panel (02) and to
#                    avoid systematically under-sampling large multi-day
#                    events. Pass `"incident"` to exclude them.
#   spatial        : "central" (restrict to points inside the core corridor
#                    polygon via spatial join) or "all_cmr" (no spatial filter).
#                    Default "central".
#   causes         : "sea" or "all". "sea" keeps Drowning + Mixed or unknown —
#                    the cause categories that map most directly to the act of
#                    crossing the sea (where SWH is the relevant exposure).
#                    Other categories (violence, vehicle accident, sickness,
#                    harsh exposure) include events from before/after the
#                    maritime leg of the journey and are dropped. Either choice
#                    introduces some measurement error; drown+mixed are the
#                    overwhelming majority of CMR deaths so they are the
#                    default and "all" is the robustness alternative.
#                    Default "sea".
#   countries      : character vector of allowed `Country of Incident` values.
#                    Default = the 5 CMR countries.
#   base_dir       : project root. Default here::here().
#
# Returns a tibble with columns `date` and `n_dead_missing`, one row per day
# with at least one matching incident. Dates with no matching incidents are
# absent and should be filled with 0 after left-joining onto the panel.
#
# The function avoids loading `sf` at source time; it uses `sf::` prefixed
# calls only when a spatial filter is requested, so callers that do not need
# spatial filtering do not need `library(sf)`.

build_iom_daily <- function(
  incident_types = c("incident", "split incident"),
  spatial        = c("central", "all_cmr"),
  causes         = c("sea", "all"),
  countries      = c("Algeria", "Italy", "Libya", "Malta", "Tunisia"),
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
    d <- tidyr::drop_na(d, lon, lat)
    core_poly <- readRDS(file.path(base_dir, "data", "processed", "core_corridor.RDS"))
    prev_s2 <- sf::sf_use_s2(FALSE)
    on.exit(sf::sf_use_s2(prev_s2), add = TRUE)
    pts <- sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
    d <- d[lengths(sf::st_within(pts, core_poly)) > 0, , drop = FALSE]
  }

  d |>
    dplyr::group_by(date) |>
    dplyr::summarise(n_dead_missing = sum(dead_missing), .groups = "drop") |>
    dplyr::arrange(date)
}
