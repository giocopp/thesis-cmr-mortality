# Clean IOM Missing Migrants Project incidents and Mediterranean crossings.

library(dplyr)
library(readxl)
library(stringr)
library(lubridate)
library(purrr)
library(tidyr)

BASE_DIR <- here::here()

# Known coordinate corrections
COORDINATE_FIXES <- list(
  "2022.MMP0765" = list(lat = 32.6486, lon = 14.2714,
                        reason = "Khums, Libya — original had Qatar coordinates")
)

# ── Date parser: handles Excel serials, ISO, DMY, ranges, partial dates ─────
parse_iom_date <- function(raw_date) {
  if (is.null(raw_date) || (length(raw_date) == 1 && is.na(raw_date))) {
    return(tibble(incident_date_clean = as.Date(NA),
                  incident_date_raw = "",
                  incident_date_precision = "unknown"))
  }

  if (inherits(raw_date, "POSIXct") || inherits(raw_date, "POSIXlt")) {
    d <- as.Date(raw_date)
    return(tibble(incident_date_clean = d,
                  incident_date_raw = as.character(d),
                  incident_date_precision = "day"))
  }
  if (inherits(raw_date, "Date")) {
    return(tibble(incident_date_clean = raw_date,
                  incident_date_raw = as.character(raw_date),
                  incident_date_precision = "day"))
  }

  s <- str_trim(as.character(raw_date))

  if (tolower(s) %in% c("unknown", "") || is.na(s)) {
    return(tibble(incident_date_clean = as.Date(NA),
                  incident_date_raw = s,
                  incident_date_precision = "unknown"))
  }

  if (grepl("^\\d{5}$", s)) {
    num <- as.numeric(s)
    if (!is.na(num) && num >= 30000 && num <= 60000) {
      d <- as.Date(num, origin = "1899-12-30")
      return(tibble(incident_date_clean = d,
                    incident_date_raw = as.character(d),
                    incident_date_precision = "day"))
    }
  }

  m <- str_match(s, "^(\\d{4}-\\d{2}-\\d{2})\\s+\\d{2}:\\d{2}:\\d{2}")
  if (!is.na(m[1, 1])) {
    d <- as.Date(m[1, 2])
    return(tibble(incident_date_clean = d,
                  incident_date_raw = m[1, 2],
                  incident_date_precision = "day"))
  }

  m <- str_match(s, "^(\\d{4})-(\\d{2})-(\\d{2})$")
  if (!is.na(m[1, 1])) {
    d <- tryCatch(as.Date(s), error = function(e) NA)
    return(tibble(incident_date_clean = d,
                  incident_date_raw = s,
                  incident_date_precision = ifelse(is.na(d), "imprecise", "day")))
  }

  m <- str_match(s, "^(\\d{1,2})\\.(\\d{1,2})\\.(\\d{4})$")
  if (!is.na(m[1, 1])) {
    d <- tryCatch(
      make_date(as.integer(m[1, 4]), as.integer(m[1, 3]), as.integer(m[1, 2])),
      error = function(e) NA
    )
    return(tibble(incident_date_clean = d,
                  incident_date_raw = s,
                  incident_date_precision = ifelse(is.na(d), "imprecise", "day")))
  }

  m <- str_match(s, "^(\\d{1,2})-\\d{1,2}\\.(\\d{1,2})\\.(\\d{4})$")
  if (!is.na(m[1, 1])) {
    d <- tryCatch(
      make_date(as.integer(m[1, 4]), as.integer(m[1, 3]), as.integer(m[1, 2])),
      error = function(e) NA
    )
    return(tibble(incident_date_clean = d,
                  incident_date_raw = s,
                  incident_date_precision = "imprecise"))
  }

  m <- str_match(s, "^(\\d{4})-(\\d{2})-\\?+$")
  if (!is.na(m[1, 1])) {
    d <- tryCatch(
      make_date(as.integer(m[1, 2]), as.integer(m[1, 3]), 1L),
      error = function(e) NA
    )
    return(tibble(incident_date_clean = d,
                  incident_date_raw = s,
                  incident_date_precision = "month"))
  }

  m <- str_match(s, "^\\?\\?[/.](\\d{1,2})[/.](\\d{4})$")
  if (!is.na(m[1, 1])) {
    d <- tryCatch(
      make_date(as.integer(m[1, 3]), as.integer(m[1, 2]), 1L),
      error = function(e) NA
    )
    return(tibble(incident_date_clean = d,
                  incident_date_raw = s,
                  incident_date_precision = "month"))
  }

  d <- suppressWarnings(ymd(s, quiet = TRUE))
  if (!is.na(d)) {
    return(tibble(incident_date_clean = d,
                  incident_date_raw = s,
                  incident_date_precision = "day"))
  }
  d <- suppressWarnings(dmy(s, quiet = TRUE))
  if (!is.na(d)) {
    return(tibble(incident_date_clean = d,
                  incident_date_raw = s,
                  incident_date_precision = "day"))
  }

  tibble(incident_date_clean = as.Date(NA),
         incident_date_raw = s,
         incident_date_precision = "imprecise")
}

# ── PART 1: IOM MMP incidents ───────────────────────────────────────────────
mmp_path <- file.path(BASE_DIR, "data", "raw", "iom",
                       "Missing_Migrants_Global_Figures_allData.xlsx")

df_raw <- read_excel(mmp_path, sheet = "Worksheet", col_types = "text")

df_all <- df_raw |>
  transmute(
    `Main ID`                   = `Main ID`,
    `Incident ID`               = `Incident ID`,
    `Incident Type`             = `Incident Type`,
    `Region of Incident`        = `Region of Incident`,
    `Incident date`             = `Incident Date`,
    `Incident year`             = suppressWarnings(as.integer(`Incident Year`)),
    `Incident month`            = match(Month, month.name),
    `No. dead`                  = suppressWarnings(as.numeric(`Number of Dead`)),
    `No. missing`               = suppressWarnings(as.numeric(`Minimum Estimated Number of Missing`)),
    `No. dead/missing`          = suppressWarnings(as.numeric(`Total Number of Dead and Missing`)),
    `No. survivors`             = suppressWarnings(as.numeric(`Number of Survivors`)),
    `No. Female`                = suppressWarnings(as.numeric(`Number of Females`)),
    `No. Male`                  = suppressWarnings(as.numeric(`Number of Males`)),
    `No. minors`                = suppressWarnings(as.numeric(`Number of Children`)),
    `Country of Origin`         = `Country of Origin`,
    `Region of Origin`          = `Region of Origin`,
    `Cause of death (category)` = `Cause of Death`,
    `Cause of death (reported)` = `Cause of Death`,
    `Route`                     = `Migration Route`,
    `Country of Incident`       = `Country of Incident`,
    `Location of death`         = `Location of Incident`,
    `UNSD region`               = `UNSD Geographical Grouping`,
    `Source`                    = `Information Source`,
    `Link`                      = URL,
    `Source Quality`            = `Source Quality`,
    Latitude  = suppressWarnings(as.numeric(sub(",.*$",           "", Coordinates))),
    Longitude = suppressWarnings(as.numeric(sub("^[^,]*,[[:space:]]*", "", Coordinates)))
  ) |>
  filter(!is.na(`Main ID`))

date_results <- map_dfr(df_all$`Incident date`, parse_iom_date)
df_all <- bind_cols(df_all, date_results)

for (main_id in names(COORDINATE_FIXES)) {
  fix <- COORDINATE_FIXES[[main_id]]
  idx <- which(df_all$`Main ID` == main_id)
  if (length(idx) > 0) {
    df_all$Latitude[idx]  <- fix$lat
    df_all$Longitude[idx] <- fix$lon
  }
}

saveRDS(df_all, file.path(BASE_DIR, "data", "processed", "iom_mmp_incidents.RDS"))

# ── PART 2: Mediterranean crossings (monthly time series) ───────────────────
crossings_path <- file.path(BASE_DIR, "data", "raw", "iom",
                             "ALL MED DATA 2010-2025_12.08.2025.xlsx")

df_raw <- read_excel(crossings_path,
                     sheet = "Crossings 2014-24 ALL ROUTE",
                     col_names = FALSE)

year_boundaries <- list(
  "2014" = c(4, 15),   "2015" = c(17, 28),  "2016" = c(30, 41),
  "2017" = c(43, 54),  "2018" = c(56, 67),  "2019" = c(69, 80),
  "2020" = c(82, 93),  "2021" = c(95, 106), "2022" = c(108, 119),
  "2023" = c(121, 132), "2024" = c(134, 145), "2025" = c(147, 158)
)

month_map <- c(
  "January" = 1, "February" = 2, "March" = 3, "April" = 4,
  "May" = 5, "June" = 6, "July" = 7, "August" = 8,
  "September" = 9, "October" = 10, "November" = 11, "December" = 12
)

col_names <- c(
  "3" = "waar", "4" = "wmr_sea_arrivals", "5" = "wmr_land_arrivals",
  "6" = "total_sea_arrivals_waar_wmr", "7" = "sea_arrivals_in_italy",
  "8" = "sea_arrivals_in_malta", "9" = "land_arrivals_in_greece",
  "10" = "sea_arrivals_in_greece", "11" = "sea_arrivals_in_cyprus",
  "12" = "land_arrivals_in_cyprus", "13" = "total_sea_arrivals_in_europe",
  "14" = "total_arrivals_in_europe",
  "15" = "interceptions_by_turkish_coast_guard",
  "16" = "interceptions_by_libyan_coast_guard",
  "17" = "interceptions_by_tunisian_coast_guard",
  "18" = "interceptions_by_algerian_coast_guard",
  "19" = "total_interceptions",
  "20" = "emr", "21" = "cmr", "22" = "wmr", "23" = "waar_1",
  "24" = "total_deaths_on_maritime_routes_to_europe",
  "25" = "total_attempted_crossings",
  "26" = "mortality_rate_maritime_routes_to_europe",
  "27" = "mortality_rate_1_in",
  "28" = "emr_rate_of_death", "29" = "cmr_rate_of_death",
  "30" = "wmr_proportion_of_deaths_vs_arrivals",
  "31" = "waar_proportion_of_deaths_vs_arrivals"
)

records <- list()

for (year_str in names(year_boundaries)) {
  bounds <- year_boundaries[[year_str]]
  year_data <- df_raw[bounds[1]:bounds[2], ]

  for (i in seq_len(nrow(year_data))) {
    row <- year_data[i, ]
    month_name <- as.character(row[[2]])

    if (is.na(month_name) || grepl("TOTAL", month_name)) next

    month_name <- str_trim(month_name)
    month_num  <- month_map[month_name]
    if (is.na(month_num)) next

    year <- as.integer(year_str)
    record <- tibble(
      date       = sprintf("%d-%02d-01", year, month_num),
      year       = year,
      month      = as.integer(month_num),
      month_name = month_name
    )

    for (col_idx_str in names(col_names)) {
      col_idx  <- as.integer(col_idx_str)
      col_name <- col_names[col_idx_str]
      val <- if (col_idx <= ncol(row)) as.numeric(row[[col_idx]]) else NA_real_
      record[[col_name]] <- val
    }

    records[[length(records) + 1]] <- record
  }
}

df_crossings <- bind_rows(records) |>
  arrange(year, month)

saveRDS(df_crossings, file.path(BASE_DIR, "data", "processed",
                                 "iom_med_crossings_monthly.RDS"))
