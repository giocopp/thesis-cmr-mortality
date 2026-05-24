# Integrated daily panel: Frontex + IOM + UNITED + UNHCR + LCG/TCG + ERA5 + ACLED.
# Daily spine 2014-01-01 to the end of interceptions_daily_disagg.RDS.

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

# ── 1. Frontex Themis (CMR departures only) ─────────────────────────────────
frx_all <- readRDS(file.path(BASE_DIR, "data", "processed", "frontex_incidents.RDS"))
frx     <- frx_all |> filter(country_of_departure %in% CMR_DEPARTURES)
FRX_END <- max(frx$date)

# ── 2. Aggregate Frontex to daily (interceptor × SAR matrix + detection + boat) ─
frx_daily <- frx |>
  mutate(
    sar_bucket = ifelse(!is.na(sar_ops) & sar_ops, "sar", "notsar"),
    int_lab = dplyr::recode(interceptor_type,
      "NGO" = "ngo", "EU_ops" = "eu", "Ita_ops" = "ita",
      "Commercial" = "comm", "EU_Coast_Guard" = "cg",
      "Land_patrol" = "land", "No_intercept" = "noint",
      "Other" = "other", "NA" = "na")
  ) |>
  group_by(date) |>
  summarise(
    frx_incidents             = n(),
    frx_persons               = sum(num_persons, na.rm = TRUE),
    frx_n_sar                 = sum(sar_bucket == "sar"),

    frx_n_sar_ngo        = sum(sar_bucket == "sar"    & int_lab == "ngo"),
    frx_n_sar_eu         = sum(sar_bucket == "sar"    & int_lab == "eu"),
    frx_n_sar_ita        = sum(sar_bucket == "sar"    & int_lab == "ita"),
    frx_n_sar_comm       = sum(sar_bucket == "sar"    & int_lab == "comm"),
    frx_n_sar_cg         = sum(sar_bucket == "sar"    & int_lab == "cg"),
    frx_n_sar_land       = sum(sar_bucket == "sar"    & int_lab == "land"),
    frx_n_sar_noint      = sum(sar_bucket == "sar"    & int_lab == "noint"),
    frx_n_sar_other      = sum(sar_bucket == "sar"    & int_lab == "other"),
    frx_n_sar_na         = sum(sar_bucket == "sar"    & int_lab == "na"),
    frx_n_notsar_ngo     = sum(sar_bucket == "notsar" & int_lab == "ngo"),
    frx_n_notsar_eu      = sum(sar_bucket == "notsar" & int_lab == "eu"),
    frx_n_notsar_ita     = sum(sar_bucket == "notsar" & int_lab == "ita"),
    frx_n_notsar_comm    = sum(sar_bucket == "notsar" & int_lab == "comm"),
    frx_n_notsar_cg      = sum(sar_bucket == "notsar" & int_lab == "cg"),
    frx_n_notsar_land    = sum(sar_bucket == "notsar" & int_lab == "land"),
    frx_n_notsar_noint   = sum(sar_bucket == "notsar" & int_lab == "noint"),
    frx_n_notsar_other   = sum(sar_bucket == "notsar" & int_lab == "other"),
    frx_n_notsar_na      = sum(sar_bucket == "notsar" & int_lab == "na"),

    frx_persons_sar_ngo      = sum(num_persons[sar_bucket == "sar"    & int_lab == "ngo"],   na.rm = TRUE),
    frx_persons_sar_eu       = sum(num_persons[sar_bucket == "sar"    & int_lab == "eu"],    na.rm = TRUE),
    frx_persons_sar_ita      = sum(num_persons[sar_bucket == "sar"    & int_lab == "ita"],   na.rm = TRUE),
    frx_persons_sar_comm     = sum(num_persons[sar_bucket == "sar"    & int_lab == "comm"],  na.rm = TRUE),
    frx_persons_sar_cg       = sum(num_persons[sar_bucket == "sar"    & int_lab == "cg"],    na.rm = TRUE),
    frx_persons_sar_land     = sum(num_persons[sar_bucket == "sar"    & int_lab == "land"],  na.rm = TRUE),
    frx_persons_sar_noint    = sum(num_persons[sar_bucket == "sar"    & int_lab == "noint"], na.rm = TRUE),
    frx_persons_sar_other    = sum(num_persons[sar_bucket == "sar"    & int_lab == "other"], na.rm = TRUE),
    frx_persons_sar_na       = sum(num_persons[sar_bucket == "sar"    & int_lab == "na"],    na.rm = TRUE),
    frx_persons_notsar_ngo   = sum(num_persons[sar_bucket == "notsar" & int_lab == "ngo"],   na.rm = TRUE),
    frx_persons_notsar_eu    = sum(num_persons[sar_bucket == "notsar" & int_lab == "eu"],    na.rm = TRUE),
    frx_persons_notsar_ita   = sum(num_persons[sar_bucket == "notsar" & int_lab == "ita"],   na.rm = TRUE),
    frx_persons_notsar_comm  = sum(num_persons[sar_bucket == "notsar" & int_lab == "comm"],  na.rm = TRUE),
    frx_persons_notsar_cg    = sum(num_persons[sar_bucket == "notsar" & int_lab == "cg"],    na.rm = TRUE),
    frx_persons_notsar_land  = sum(num_persons[sar_bucket == "notsar" & int_lab == "land"],  na.rm = TRUE),
    frx_persons_notsar_noint = sum(num_persons[sar_bucket == "notsar" & int_lab == "noint"], na.rm = TRUE),
    frx_persons_notsar_other = sum(num_persons[sar_bucket == "notsar" & int_lab == "other"], na.rm = TRUE),
    frx_persons_notsar_na    = sum(num_persons[sar_bucket == "notsar" & int_lab == "na"],    na.rm = TRUE),

    frx_n_det_fwa             = sum(grepl("FWA",  detected_by), na.rm = TRUE),
    frx_n_det_helo            = sum(grepl("HELO", detected_by), na.rm = TRUE),
    frx_n_det_rpas            = sum(grepl("RPAS", detected_by), na.rm = TRUE),
    frx_n_det_mas             = sum(grepl("MAS",  detected_by), na.rm = TRUE),
    frx_n_det_aerial          = sum(grepl("FWA|HELO|RPAS|MAS", detected_by), na.rm = TRUE),
    frx_n_det_ngo             = sum(grepl("NGO vessel", detected_by), na.rm = TRUE),
    frx_n_det_call            = sum(grepl("Call-", detected_by), na.rm = TRUE),
    frx_n_det_cg              = sum(grepl("CPV|CPB|OPV", detected_by), na.rm = TRUE),
    frx_n_det_land            = sum(grepl("Land Patrol", detected_by), na.rm = TRUE),

    frx_n_inflatable          = sum(boat_category == "Inflatable"),
    frx_persons_inflatable    = sum(num_persons[boat_category == "Inflatable"], na.rm = TRUE),
    frx_n_wooden              = sum(boat_category == "Wooden"),
    frx_persons_wooden        = sum(num_persons[boat_category == "Wooden"], na.rm = TRUE),
    frx_n_fibreglass          = sum(boat_category == "Fibre glass"),
    frx_persons_fibreglass    = sum(num_persons[boat_category == "Fibre glass"], na.rm = TRUE),

    frx_dep_libya             = sum(country_of_departure == "Libya"),
    frx_dep_tunisia           = sum(country_of_departure == "Tunisia"),
    frx_dep_algeria           = sum(country_of_departure == "Algeria"),

    frx_n_in_oparea           = sum(in_op_area, na.rm = TRUE),
    frx_n_multi_actors        = sum(multi_actors_inv, na.rm = TRUE),

    .groups = "drop"
  ) |>
  mutate(
    frx_sar_share                 = frx_n_sar / frx_incidents,
    frx_det_aerial_share          = frx_n_det_aerial / frx_incidents,
    frx_det_ngo_share             = frx_n_det_ngo / frx_incidents,
    frx_det_call_share            = frx_n_det_call / frx_incidents,
    frx_det_cg_share              = frx_n_det_cg / frx_incidents,
    frx_det_land_share            = frx_n_det_land / frx_incidents,
    frx_inflatable_share          = frx_n_inflatable / frx_incidents,
    frx_wooden_share              = frx_n_wooden / frx_incidents,
    frx_fibreglass_share          = frx_n_fibreglass / frx_incidents,
    frx_frac_inflatable_persons   = ifelse(frx_persons > 0,
                                            frx_persons_inflatable / frx_persons, NA_real_),
    frx_libya_share               = frx_dep_libya / frx_incidents,
    frx_in_oparea_share           = frx_n_in_oparea / frx_incidents
  )

# ── 3. Daily IOM (broad descriptive) and UNITED (analytical, feeds C_t) ─────
iom_daily <- build_iom_daily(
  incident_types = c("incident", "split incident"),
  spatial        = "all_cmr",
  causes         = "all",
  countries      = CMR_INCIDENT_COUNTRIES,
  base_dir       = BASE_DIR
)

united_for_ct <- build_united_daily() |>
  dplyr::rename(n_dead_united_for_ct = n_dead_united)

# ── 4. Date spine + merge all sources ───────────────────────────────────────
disagg_path <- file.path(BASE_DIR, "analysis", "data", "interceptions_daily_disagg.RDS")
if (!file.exists(disagg_path)) {
  stop("interceptions_daily_disagg.RDS not found — run 02 first.")
}
PANEL_END <- max(readRDS(disagg_path)$date)

spine <- tibble(date = seq(as.Date("2014-01-01"), PANEL_END, by = "day"))

weather <- readRDS(file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS"))

acled <- readRDS(file.path(BASE_DIR, "data", "processed", "acled_daily.RDS")) |>
  select(date, week_date,
         libya_conflict, libya_conflict_fatalities,
         libya_battles, libya_expvio, libya_violciv,
         tunisia_conflict, tunisia_conflict_fatalities,
         tunisia_battles, tunisia_expvio, tunisia_violciv)

lcg_tcg_daily <- readRDS(disagg_path) |>
  select(date, lcg_pushbacks, tcg_pushbacks, lcg_tcg_pushbacks)

unhcr <- readRDS(file.path(BASE_DIR, "data", "processed",
                            "unhcr_daily_arrivals.RDS"))

panel <- spine |>
  mutate(iso_week = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))) |>
  left_join(frx_daily,     by = "date") |>
  left_join(iom_daily,     by = "date") |>
  left_join(united_for_ct, by = "date") |>
  left_join(weather,       by = "date") |>
  left_join(acled,         by = "date") |>
  left_join(lcg_tcg_daily, by = "date") |>
  left_join(unhcr,         by = "date")

int_types <- c("ngo","eu","ita","comm","cg","land","noint","other","na")
frx_matrix_cols <- c(
  paste0("frx_n_sar_",       int_types),
  paste0("frx_n_notsar_",    int_types),
  paste0("frx_persons_sar_",    int_types),
  paste0("frx_persons_notsar_", int_types)
)
count_cols <- c("frx_incidents", "frx_persons", "frx_n_sar",
                frx_matrix_cols,
                "frx_n_det_fwa", "frx_n_det_helo", "frx_n_det_rpas",
                "frx_n_det_mas", "frx_n_det_aerial", "frx_n_det_ngo",
                "frx_n_det_call", "frx_n_det_cg", "frx_n_det_land",
                "frx_n_inflatable", "frx_persons_inflatable",
                "frx_n_wooden", "frx_persons_wooden",
                "frx_n_fibreglass", "frx_persons_fibreglass",
                "frx_dep_libya", "frx_dep_tunisia", "frx_dep_algeria",
                "frx_n_in_oparea", "frx_n_multi_actors",
                "n_dead_missing", "n_dead_united_for_ct",
                "lcg_pushbacks", "tcg_pushbacks", "lcg_tcg_pushbacks")
panel <- panel |>
  mutate(across(all_of(count_cols), ~ replace_na(.x, 0L)))

# ── 5. Derived variables ────────────────────────────────────────────────────
# crossing_attempts is a daily lower-bound (frx_persons + LCG/TCG + UNITED deaths).
panel <- panel |>
  mutate(
    crossing_attempts = frx_persons + lcg_tcg_pushbacks + n_dead_united_for_ct,
    fatality_rate     = ifelse(crossing_attempts > 0,
                                n_dead_united_for_ct / crossing_attempts, NA_real_),
    post_mou          = as.integer(date >= MOU_DATE),
    year              = year(date),
    month_year        = factor(format(date, "%Y-%m"))
  )

# ── 6. Save ─────────────────────────────────────────────────────────────────
out_dir <- file.path(BASE_DIR, "analysis", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(panel, file.path(out_dir, "daily_panel_complete.RDS"))
