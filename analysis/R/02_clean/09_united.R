# Clean UNITED List of Refugee Deaths dataset.

library(dplyr)
library(readxl)
library(stringr)
library(lubridate)
library(purrr)
library(tidyr)

BASE_DIR <- here::here()

MONTH_LOOKUP <- c(
  "jan" = 1L, "january" = 1L,
  "feb" = 2L, "february" = 2L,
  "mar" = 3L, "march" = 3L,
  "apr" = 4L, "april" = 4L,
  "may" = 5L,
  "jun" = 6L, "june" = 6L,
  "jul" = 7L, "july" = 7L,
  "aug" = 8L, "august" = 8L,
  "sep" = 9L, "sept" = 9L, "september" = 9L,
  "oct" = 10L, "october" = 10L,
  "nov" = 11L, "november" = 11L,
  "dec" = 12L, "december" = 12L
)

SEASON_MID <- list(
  spring = list(month = 4L, day = 15L),
  summer = list(month = 7L, day = 15L),
  autumn = list(month = 10L, day = 15L),
  fall   = list(month = 10L, day = 15L),
  winter = list(month = 1L, day = 15L)
)

# ── Date helpers ────────────────────────────────────────────────────────────
resolve_year <- function(yy) {
  yy <- as.integer(yy)
  ifelse(yy <= 30L, 2000L + yy, 1900L + yy)
}

parse_united_date <- function(raw_text, sort_date_posix = NA) {
  result_na <- tibble(incident_date_clean = as.Date(NA),
                      incident_date_raw   = as.character(raw_text),
                      incident_date_precision = "unknown")

  if (is.null(raw_text) || is.na(raw_text) || str_trim(raw_text) == "") {
    result_na$incident_date_raw <- ""
    return(result_na)
  }

  txt <- str_trim(raw_text)

  m <- str_match(txt, "^(\\d{1,2})[/\\-](\\d{1,2})[/\\-](\\d{1,4})$")
  if (!is.na(m[1, 1])) {
    dd <- as.integer(m[1, 2])
    mm <- as.integer(m[1, 3])
    yr_raw <- m[1, 4]
    yr <- if (nchar(yr_raw) == 4) as.integer(yr_raw) else resolve_year(yr_raw)
    d <- tryCatch(make_date(yr, mm, dd), error = function(e) NA)
    if (!is.na(d)) {
      return(tibble(incident_date_clean = d,
                    incident_date_raw   = txt,
                    incident_date_precision = "day"))
    }
  }

  txt_lower <- tolower(txt)

  m <- str_match(txt_lower, "^in\\s+mid\\s+(\\w+)\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    mo <- MONTH_LOOKUP[m[1, 2]]
    if (!is.na(mo)) {
      yr <- if (nchar(m[1, 3]) == 4) as.integer(m[1, 3]) else resolve_year(m[1, 3])
      return(tibble(incident_date_clean = make_date(yr, mo, 15L),
                    incident_date_raw   = txt,
                    incident_date_precision = "month"))
    }
  }

  m <- str_match(txt_lower, "^mid\\s+(\\w+)\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    mo <- MONTH_LOOKUP[m[1, 2]]
    if (!is.na(mo)) {
      yr <- if (nchar(m[1, 3]) == 4) as.integer(m[1, 3]) else resolve_year(m[1, 3])
      return(tibble(incident_date_clean = make_date(yr, mo, 15L),
                    incident_date_raw   = txt,
                    incident_date_precision = "month"))
    }
  }

  m <- str_match(txt_lower, "^beg(?:in)?\\s+(\\w+)\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    mo <- MONTH_LOOKUP[m[1, 2]]
    if (!is.na(mo)) {
      yr <- if (nchar(m[1, 3]) == 4) as.integer(m[1, 3]) else resolve_year(m[1, 3])
      return(tibble(incident_date_clean = make_date(yr, mo, 1L),
                    incident_date_raw   = txt,
                    incident_date_precision = "month"))
    }
  }

  m <- str_match(txt_lower, "^end\\s+(\\w+)\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    mo <- MONTH_LOOKUP[m[1, 2]]
    if (!is.na(mo)) {
      yr <- if (nchar(m[1, 3]) == 4) as.integer(m[1, 3]) else resolve_year(m[1, 3])
      last_day <- as.Date(ceiling_date(make_date(yr, mo, 1L), "month")) - 1L
      return(tibble(incident_date_clean = last_day,
                    incident_date_raw   = txt,
                    incident_date_precision = "month"))
    }
  }

  m <- str_match(txt_lower, "^in\\s+(\\w+)\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    mo <- MONTH_LOOKUP[m[1, 2]]
    if (!is.na(mo)) {
      yr <- if (nchar(m[1, 3]) == 4) as.integer(m[1, 3]) else resolve_year(m[1, 3])
      return(tibble(incident_date_clean = make_date(yr, mo, 15L),
                    incident_date_raw   = txt,
                    incident_date_precision = "month"))
    }
    season <- SEASON_MID[[m[1, 2]]]
    if (!is.null(season)) {
      yr <- if (nchar(m[1, 3]) == 4) as.integer(m[1, 3]) else resolve_year(m[1, 3])
      return(tibble(incident_date_clean = make_date(yr, season$month, season$day),
                    incident_date_raw   = txt,
                    incident_date_precision = "imprecise"))
    }
  }

  m <- str_match(txt_lower, "^(\\w+)/(\\w+)\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    mo1 <- MONTH_LOOKUP[m[1, 2]]
    if (!is.na(mo1)) {
      yr <- if (nchar(m[1, 4]) == 4) as.integer(m[1, 4]) else resolve_year(m[1, 4])
      return(tibble(incident_date_clean = make_date(yr, mo1, 15L),
                    incident_date_raw   = txt,
                    incident_date_precision = "imprecise"))
    }
  }

  m <- str_match(txt_lower, "^begin\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    yr <- if (nchar(m[1, 2]) == 4) as.integer(m[1, 2]) else resolve_year(m[1, 2])
    return(tibble(incident_date_clean = make_date(yr, 1L, 1L),
                  incident_date_raw   = txt,
                  incident_date_precision = "year_only"))
  }
  m <- str_match(txt_lower, "^end\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    yr <- if (nchar(m[1, 2]) == 4) as.integer(m[1, 2]) else resolve_year(m[1, 2])
    return(tibble(incident_date_clean = make_date(yr, 12L, 31L),
                  incident_date_raw   = txt,
                  incident_date_precision = "year_only"))
  }

  m <- str_match(txt_lower, "^in\\s+mid\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    yr <- if (nchar(m[1, 2]) == 4) as.integer(m[1, 2]) else resolve_year(m[1, 2])
    return(tibble(incident_date_clean = make_date(yr, 7L, 1L),
                  incident_date_raw   = txt,
                  incident_date_precision = "year_only"))
  }

  m <- str_match(txt_lower, "^in\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    yr <- if (nchar(m[1, 2]) == 4) as.integer(m[1, 2]) else resolve_year(m[1, 2])
    return(tibble(incident_date_clean = make_date(yr, 7L, 1L),
                  incident_date_raw   = txt,
                  incident_date_precision = "year_only"))
  }

  m <- str_match(txt, "^(\\d{4})$")
  if (!is.na(m[1, 1])) {
    yr <- as.integer(m[1, 2])
    return(tibble(incident_date_clean = make_date(yr, 7L, 1L),
                  incident_date_raw   = txt,
                  incident_date_precision = "year_only"))
  }

  if (!is.na(sort_date_posix)) {
    sd <- tryCatch(as.Date(as.numeric(sort_date_posix),
                           origin = "1899-12-30"),
                   error = function(e) NA)
    if (!is.na(sd)) {
      return(tibble(incident_date_clean = sd,
                    incident_date_raw   = txt,
                    incident_date_precision = "imprecise"))
    }
  }

  result_na
}

collapse_flags <- function(df, flag_cols) {
  map_chr(seq_len(nrow(df)), function(i) {
    vals <- as.numeric(df[i, flag_cols, drop = TRUE])
    active <- flag_cols[!is.na(vals) & vals > 0]
    if (length(active) == 0L) NA_character_
    else if (length(active) == 1L) active
    else "multiple"
  })
}

# ── 1. Read raw ─────────────────────────────────────────────────────────────
raw_path <- file.path(BASE_DIR, "data", "raw", "united",
                      "UNITED_ListOfRefugeeDeaths.xlsx")
df_raw <- read_excel(raw_path, sheet = "LIST", skip = 1, col_types = "text")

# ── 2. Assign clean column names ────────────────────────────────────────────
raw_names <- names(df_raw)
n_cols <- length(raw_names)

core_names <- c(
  "record_id", "sort_date", "spacer_1", "found_dead_text", "n_deaths",
  "name_gender_age", "region_of_origin", "cause_of_death_text",
  "source_text", "spacer_2", "country_of_death", "place_of_death",
  "latitude", "longitude", "crossing_border"
)
country_flag_names <- paste0("cflag_", raw_names[16:67])

manner_raw    <- raw_names[69:80]
manner_names  <- paste0("manner_",
                        str_to_lower(str_remove(str_replace_all(manner_raw, "[/ ]+", "_"),
                                                "\\.\\.\\..+$")))

suicide_raw   <- raw_names[82:91]
suicide_names <- paste0("suicide_",
                        str_to_lower(str_remove(str_replace_all(suicide_raw, "[/ ]+", "_"),
                                                "\\.\\.\\..+$")))

transport_org_raw   <- raw_names[93:96]
transport_org_names <- paste0("torg_",
                              str_to_lower(str_remove(str_replace_all(transport_org_raw, "[/ ]+", "_"),
                                                      "\\.\\.\\..+$")))

transport_means_raw   <- raw_names[98:103]
transport_means_names <- paste0("tmeans_",
                                str_to_lower(str_remove(str_replace_all(transport_means_raw, "[/ ]+", "_"),
                                                        "\\.\\.\\..+$")))

state_services_raw   <- raw_names[105:108]
state_services_names <- paste0("state_",
                               str_to_lower(str_replace_all(state_services_raw, "[/ ]+", "_")))

where_raw   <- raw_names[110:114]
where_names <- paste0("where_",
                      str_to_lower(str_replace_all(where_raw, "[/ ]+", "_")))

tail_names <- c("event_group", "spacer_7", "weblink", "spacer_8",
                "long_description")

spacers_between <- c("spacer_3", "spacer_4", "spacer_5", "spacer_6")

new_names <- c(
  core_names,
  country_flag_names,
  spacers_between[1],
  manner_names,
  spacers_between[2],
  suicide_names,
  spacers_between[3],
  transport_org_names,
  spacers_between[4],
  transport_means_names,
  "spacer_transport_end",
  state_services_names,
  "spacer_state_end",
  where_names,
  tail_names
)

if (length(new_names) < n_cols) {
  new_names <- c(new_names, paste0("extra_", seq_len(n_cols - length(new_names))))
} else if (length(new_names) > n_cols) {
  new_names <- new_names[seq_len(n_cols)]
}

names(df_raw) <- new_names
df <- df_raw |> select(-all_of(grep("^spacer_", names(df_raw), value = TRUE)))

# ── 3. Drop summary / junk rows ─────────────────────────────────────────────
df <- df |>
  filter(toupper(str_trim(found_dead_text)) != "TOTAL" | is.na(found_dead_text))
df$n_deaths_num <- as.numeric(df$n_deaths)
df <- df |> filter(is.na(n_deaths_num) | n_deaths_num <= 10000)
df <- df |> filter(!is.na(record_id) & str_trim(record_id) != "")

# ── 4. Type conversion ──────────────────────────────────────────────────────
df <- df |>
  mutate(
    record_id       = as.integer(record_id),
    n_deaths        = as.numeric(n_deaths),
    latitude        = as.numeric(latitude),
    longitude       = as.numeric(longitude),
    crossing_border = as.numeric(crossing_border)
  )

flag_cols_all <- grep("^cflag_|^manner_|^suicide_|^torg_|^tmeans_|^state_|^where_",
                       names(df), value = TRUE)
for (fc in flag_cols_all) {
  df[[fc]] <- suppressWarnings(as.numeric(df[[fc]]))
}

# ── 5. Parse dates ──────────────────────────────────────────────────────────
date_parsed <- map2_dfr(df$found_dead_text, df$sort_date, parse_united_date)

df <- df |>
  bind_cols(date_parsed) |>
  mutate(
    incident_year  = year(incident_date_clean),
    incident_month = month(incident_date_clean)
  )

# ── 6. Collapse flag columns ────────────────────────────────────────────────
df$manner_of_death <- str_remove(collapse_flags(df, grep("^manner_",  names(df), value = TRUE)), "^manner_")
df$suicide_manner  <- str_remove(collapse_flags(df, grep("^suicide_", names(df), value = TRUE)), "^suicide_")
df$transport_org   <- str_remove(collapse_flags(df, grep("^torg_",    names(df), value = TRUE)), "^torg_")
df$transport_means <- str_remove(collapse_flags(df, grep("^tmeans_",  names(df), value = TRUE)), "^tmeans_")
df$state_services  <- str_remove(collapse_flags(df, grep("^state_",   names(df), value = TRUE)), "^state_")
df$where_died      <- str_remove(collapse_flags(df, grep("^where_",   names(df), value = TRUE)), "^where_")

cflag_cols <- grep("^cflag_", names(df), value = TRUE)
df$crossing_countries <- map_chr(seq_len(nrow(df)), function(i) {
  vals <- as.numeric(df[i, cflag_cols, drop = TRUE])
  active <- cflag_cols[!is.na(vals) & vals > 0]
  if (length(active) == 0L) NA_character_
  else paste(str_remove(active, "^cflag_"), collapse = ";")
})

# ── 7. CMR indicator ────────────────────────────────────────────────────────
CMR_COUNTRIES_DEATH <- c("italy", "libya", "malta", "tunisia", "algeria",
                          "mediterranean")
CMR_PLACE_KEYWORDS  <- "mediterranean|lampedusa|sicil|malta|libya|tunis|channel of sicily|pantelleria"
CMR_COUNTRY_CODES   <- c("IT", "LY", "MT", "TN", "DZ")

df <- df |>
  mutate(
    is_cmr = (
      tolower(country_of_death) %in% CMR_COUNTRIES_DEATH |
      str_detect(tolower(place_of_death), CMR_PLACE_KEYWORDS) |
      map_lgl(crossing_countries, function(cc) {
        if (is.na(cc)) return(FALSE)
        any(str_split(cc, ";")[[1]] %in% CMR_COUNTRY_CODES)
      })
    )
  )

# ── 8. Select final columns ─────────────────────────────────────────────────
united_clean <- df |>
  transmute(
    record_id, incident_date_clean, incident_date_raw,
    incident_date_precision, incident_year, incident_month,
    n_deaths, name_gender_age, region_of_origin, cause_of_death_text,
    source_text, country_of_death, place_of_death,
    latitude, longitude, crossing_countries,
    manner_of_death, suicide_manner,
    transport_means, transport_org, state_services, where_died,
    event_group, is_cmr, weblink, long_description
  )

saveRDS(united_clean,
        file.path(BASE_DIR, "data", "processed", "united_incidents.RDS"))
