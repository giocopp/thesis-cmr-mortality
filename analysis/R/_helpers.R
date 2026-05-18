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
#                    Default c("incident") — split incidents EXCLUDED.
#                    "Split incidents" are IOM's decomposition of ONE
#                    reported event across multiple calendar dates. That
#                    date-smearing decouples deaths from the day's recent
#                    SWH, flattens the SWH->mortality gradient, and has no
#                    analogue in UNITED (one event -> one date), so keeping
#                    them breaks IOM/UNITED comparability. Verified in the
#                    P2/P3 decomposition: dropping split moves the IOM
#                    4-period gradient toward UNITED (P2 +0.14 -> +0.39,
#                    P3 +0.04 -> +0.46) and barely touches P1. Costs ~10%
#                    of sample rows. Pass c("incident", "split incident")
#                    to opt back in (e.g. the broad volume-denominator
#                    build in 01_build_daily_panel.R).
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
  incident_types = c("incident"),
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

# build_united_daily()
# --------------------
# UNITED daily death aggregate, constructed with the SAME spatial + cause
# logic as build_iom_daily() so the two sources are directly comparable.
# This is the single source of truth for the UNITED series: scripts must
# call this instead of re-implementing the filter inline (inline copies
# drift — that is exactly how 31 diverged from 20/28).
#
# With the defaults this reproduces, byte-for-byte, the UNITED block in
# 20_primary_model.R / 28_period_sar_gradient.R (country in CMR+Med, manner
# in {drowned, other_unknown}, non-missing lon/lat, spatial join to the
# core_corridor polygon).
#
# Arguments
#   causes    : "sea" (manner_of_death in {drowned, other_unknown}) or
#               "all" (no manner-of-death filter). Default "sea".
#   spatial   : "central" (spatial join to the core corridor polygon;
#               requires non-missing lon/lat) or "all_cmr" (no spatial
#               filter). Default "central".
#   countries : allowed `country_of_death` values. Default = the 5 CMR
#               countries + "Mediterranean" (UNITED records open-sea
#               deaths under "Mediterranean"; the IOM analogue does not
#               need this because build_iom_daily filters Country of
#               Incident on the CMR route instead).
#   base_dir  : project root. Default here::here().
#
# Returns a tibble with columns `date` and `n_dead_united`, one row per
# day with at least one matching incident. Dates with no matching incident
# are absent and should be filled with 0 after left-joining onto the panel.
#
# Like build_iom_daily(), this avoids loading `sf` at source time and only
# uses `sf::` prefixed calls when spatial == "central".

build_united_daily <- function(
  causes    = c("sea", "all"),
  spatial   = c("central", "all_cmr"),
  countries = c("Algeria", "Italy", "Libya", "Malta", "Tunisia", "Mediterranean"),
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
    d <- dplyr::filter(d, !is.na(latitude), !is.na(longitude))
    core_poly <- readRDS(file.path(base_dir, "data", "processed", "core_corridor.RDS"))
    prev_s2 <- sf::sf_use_s2(FALSE)
    on.exit(sf::sf_use_s2(prev_s2), add = TRUE)
    pts <- sf::st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326)
    d <- d[lengths(sf::st_within(pts, core_poly)) > 0, , drop = FALSE]
  }

  d |>
    dplyr::group_by(date) |>
    dplyr::summarise(n_dead_united = sum(n_deaths, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(date)
}

# add_crossing_exposure()
# -----------------------
# Single source of truth for the constructed crossing-exposure denominator
# (shared by the rate-model scripts; previously inline-duplicated in 20/27/28).
#
# Call AFTER the build_*_daily() left-joins and replace_na(): input must
# already have `frx_persons`, `lcg_tcg_pushbacks`, `n_dead_iom`,
# `n_dead_united`.
#
# Adds:
#   living_crossings = frx_persons + lcg_tcg_pushbacks
#   attempts_iom     = living_crossings + n_dead_iom
#   attempts_united  = living_crossings + n_dead_united   (source-specific)
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
