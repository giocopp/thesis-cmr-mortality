# 06_clean_united.R
# =================
# Clean UNITED List of Refugee Deaths dataset.
#
# The UNITED dataset records refugee/migrant deaths in and around Europe
# since 1993, compiled from media reports by the UNITED network. Unlike
# IOM MMP (which starts in 2014), UNITED covers 1993-2026, providing the
# only incident-level source for pre-2014 deaths.
#
# Fixes applied:
#   - Skip merged header row (skip=1), assign clean column names positionally
#   - Filter total/summary rows (death count > 10000, "TOTAL" text)
#   - Parse 340 unique fuzzy date patterns ("in mid Feb 26", "end Dec 25",
#     etc.) with precision tracking matching IOM convention
#   - Collapse 50+ country-flag columns into single crossing_countries field
#   - Collapse manner-of-death, transport, etc. flag columns into categoricals
#
# Input:
#   data/raw/united/UNITED_ListOfRefugeeDeaths.xlsx  (sheet "LIST")
#
# Output:
#   data/processed/united_incidents.RDS

library(dplyr)
library(readxl)
library(stringr)
library(lubridate)
library(purrr)
library(tidyr)

BASE_DIR <- here::here()

cat("============================================================\n")
cat("CLEAN UNITED DATA\n")
cat("============================================================\n\n")


# ====================================================================
# Month and season lookup tables
# ====================================================================
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


# ====================================================================
# Helper: resolve short year (1-3 digits)
# ====================================================================
resolve_year <- function(yy) {
  yy <- as.integer(yy)
  # Single digit (e.g., "6" = 2006) or <= 30: treat as 2000s
  ifelse(yy <= 30L, 2000L + yy, 1900L + yy)
}


# ====================================================================
# Helper: parse a single UNITED date text value
# ====================================================================
parse_united_date <- function(raw_text, sort_date_posix = NA) {
  # Returns tibble row: incident_date_clean, incident_date_raw, incident_date_precision

  result_na <- tibble(incident_date_clean = as.Date(NA),
                      incident_date_raw   = as.character(raw_text),
                      incident_date_precision = "unknown")

  if (is.null(raw_text) || is.na(raw_text) || str_trim(raw_text) == "") {
    result_na$incident_date_raw <- ""
    return(result_na)
  }

  txt <- str_trim(raw_text)

  # --- 1. Standard DD/MM/YY or DD-MM-YY (1-4 digit year) ---
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

  # --- 2. "in mid Mon YY" ---
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

  # --- 3. "mid Mon YY" (no "in" prefix) ---
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

  # --- 4. "begin Mon YY" / "beg Mon YY" ---
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

  # --- 5. "end Mon YY" ---
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

  # --- 6. "in Mon YY" / "in Month YY" ---
  m <- str_match(txt_lower, "^in\\s+(\\w+)\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    mo <- MONTH_LOOKUP[m[1, 2]]
    if (!is.na(mo)) {
      yr <- if (nchar(m[1, 3]) == 4) as.integer(m[1, 3]) else resolve_year(m[1, 3])
      return(tibble(incident_date_clean = make_date(yr, mo, 15L),
                    incident_date_raw   = txt,
                    incident_date_precision = "month"))
    }
    # Check for season in this same pattern (e.g., "in spring 25")
    season <- SEASON_MID[[m[1, 2]]]
    if (!is.null(season)) {
      yr <- if (nchar(m[1, 3]) == 4) as.integer(m[1, 3]) else resolve_year(m[1, 3])
      return(tibble(incident_date_clean = make_date(yr, season$month, season$day),
                    incident_date_raw   = txt,
                    incident_date_precision = "imprecise"))
    }
  }

  # --- 7. "Mon/Mon YY" (month range, e.g., "Jan/Feb 19") ---
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

  # --- 8. "begin YYYY" / "end YYYY" / "begin YY" / "end YY" ---
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

  # --- 9. "in mid YY" (year only, no month) ---
  m <- str_match(txt_lower, "^in\\s+mid\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    yr <- if (nchar(m[1, 2]) == 4) as.integer(m[1, 2]) else resolve_year(m[1, 2])
    return(tibble(incident_date_clean = make_date(yr, 7L, 1L),
                  incident_date_raw   = txt,
                  incident_date_precision = "year_only"))
  }

  # --- 10. "in YYYY" / "in YY" ---
  m <- str_match(txt_lower, "^in\\s+(\\d{2,4})$")
  if (!is.na(m[1, 1])) {
    yr <- if (nchar(m[1, 2]) == 4) as.integer(m[1, 2]) else resolve_year(m[1, 2])
    return(tibble(incident_date_clean = make_date(yr, 7L, 1L),
                  incident_date_raw   = txt,
                  incident_date_precision = "year_only"))
  }

  # --- 11. Bare 4-digit year ---
  m <- str_match(txt, "^(\\d{4})$")
  if (!is.na(m[1, 1])) {
    yr <- as.integer(m[1, 2])
    return(tibble(incident_date_clean = make_date(yr, 7L, 1L),
                  incident_date_raw   = txt,
                  incident_date_precision = "year_only"))
  }

  # --- 12. Fallback: use sort_date column (Excel serial number) ---
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

  # --- 13. Final fallback ---
  return(result_na)
}


# ====================================================================
# Helper: collapse flag columns into a single categorical
# ====================================================================
collapse_flags <- function(df, flag_cols) {
  # For each row, return the name of the single active flag column.
  # If multiple active: "multiple". If none: NA.
  map_chr(seq_len(nrow(df)), function(i) {
    vals <- as.numeric(df[i, flag_cols, drop = TRUE])
    active <- flag_cols[!is.na(vals) & vals > 0]
    if (length(active) == 0L) NA_character_
    else if (length(active) == 1L) active
    else "multiple"
  })
}


# ====================================================================
# 1. Read raw data
# ====================================================================
cat("--- 1. Reading raw data ---\n")

raw_path <- file.path(BASE_DIR, "data", "raw", "united",
                      "UNITED_ListOfRefugeeDeaths.xlsx")

df_raw <- read_excel(raw_path, sheet = "LIST", skip = 1, col_types = "text")
cat(sprintf("  Raw dimensions: %d rows x %d cols\n", nrow(df_raw), ncol(df_raw)))


# ====================================================================
# 2. Assign clean column names
# ====================================================================
cat("\n--- 2. Assigning clean column names ---\n")

raw_names <- names(df_raw)
n_cols <- length(raw_names)
cat(sprintf("  Columns read: %d\n", n_cols))

# Core fields (positions 1-15)
core_names <- c(
  "record_id",           # 1
  "sort_date",           # 2
  "spacer_1",            # 3
  "found_dead_text",     # 4
  "n_deaths",            # 5
  "name_gender_age",     # 6
  "region_of_origin",    # 7
  "cause_of_death_text", # 8
  "source_text",         # 9
  "spacer_2",            # 10
  "country_of_death",    # 11
  "place_of_death",      # 12
  "latitude",            # 13
  "longitude",           # 14
  "crossing_border"      # 15
)

# Country flag columns (positions 16-67)
country_flag_names <- paste0("cflag_", raw_names[16:67])

# Spacer at 68
# Manner of death (positions 69-80)
manner_raw <- raw_names[69:80]
manner_clean <- str_replace_all(manner_raw, "[/ ]+", "_") |>
  str_remove("\\.\\.\\..+$") |>
  str_to_lower()
manner_names <- paste0("manner_", manner_clean)

# Spacer at 81
# Suicide manner (positions 82-91)
suicide_raw <- raw_names[82:91]
suicide_clean <- str_replace_all(suicide_raw, "[/ ]+", "_") |>
  str_remove("\\.\\.\\..+$") |>
  str_to_lower()
suicide_names <- paste0("suicide_", suicide_clean)

# Spacer at 92
# Transport organisation (positions 93-96)
transport_org_raw <- raw_names[93:96]
transport_org_clean <- str_replace_all(transport_org_raw, "[/ ]+", "_") |>
  str_remove("\\.\\.\\..+$") |>
  str_to_lower()
transport_org_names <- paste0("torg_", transport_org_clean)

# Spacer at 97
# Transport means (positions 98-103)
transport_means_raw <- raw_names[98:103]
transport_means_clean <- str_replace_all(transport_means_raw, "[/ ]+", "_") |>
  str_remove("\\.\\.\\..+$") |>
  str_to_lower()
transport_means_names <- paste0("tmeans_", transport_means_clean)

# Spacer at 104
# State services (positions 105-108)
state_services_raw <- raw_names[105:108]
state_services_clean <- str_replace_all(state_services_raw, "[/ ]+", "_") |>
  str_to_lower()
state_services_names <- paste0("state_", state_services_clean)

# Spacer at 109
# Where (positions 110-114)
where_raw <- raw_names[110:114]
where_clean <- str_replace_all(where_raw, "[/ ]+", "_") |>
  str_to_lower()
where_names <- paste0("where_", where_clean)

# Spacer at 115
# Tail columns
tail_names <- c("event_group",    # 116
                "spacer_7",       # 117
                "weblink",        # 118
                "spacer_8",       # 119
                "long_description") # 120

# Build full name vector
spacers_between <- c("spacer_3", "spacer_4", "spacer_5", "spacer_6")
# Positions: 68, 81, 92, 97, 104, 109, 115

new_names <- c(
  core_names,                # 1-15
  country_flag_names,        # 16-67
  spacers_between[1],        # 68
  manner_names,              # 69-80
  spacers_between[2],        # 81
  suicide_names,             # 82-91
  spacers_between[3],        # 92
  transport_org_names,       # 93-96
  spacers_between[4],        # 97
  transport_means_names,     # 98-103
  "spacer_transport_end",    # 104
  state_services_names,      # 105-108
  "spacer_state_end",        # 109
  where_names,               # 110-114
  tail_names                 # 115(=spacer in tail), 116-120
)

# Handle case where Excel adds extra trailing columns
if (length(new_names) < n_cols) {
  new_names <- c(new_names, paste0("extra_", seq_len(n_cols - length(new_names))))
} else if (length(new_names) > n_cols) {
  new_names <- new_names[seq_len(n_cols)]
}

names(df_raw) <- new_names
cat(sprintf("  Assigned %d clean names\n", length(new_names)))

# Drop spacer columns
spacer_cols <- grep("^spacer_", names(df_raw), value = TRUE)
df <- df_raw |> select(-all_of(spacer_cols))
cat(sprintf("  After dropping %d spacer columns: %d cols remain\n",
            length(spacer_cols), ncol(df)))


# ====================================================================
# 3. Filter junk rows
# ====================================================================
cat("\n--- 3. Filtering junk rows ---\n")

n_before <- nrow(df)

# Drop TOTAL rows
total_idx <- which(toupper(str_trim(df$found_dead_text)) == "TOTAL")
if (length(total_idx) > 0) {
  df <- df[-total_idx, ]
  cat(sprintf("  Dropped %d TOTAL rows\n", length(total_idx)))
}

# Drop rows with extreme death counts (summary rows)
df$n_deaths_num <- as.numeric(df$n_deaths)
extreme_idx <- which(df$n_deaths_num > 10000)
if (length(extreme_idx) > 0) {
  df <- df[-extreme_idx, ]
  cat(sprintf("  Dropped %d rows with n_deaths > 10000\n", length(extreme_idx)))
}

# Drop rows without a record ID
na_id <- sum(is.na(df$record_id) | str_trim(df$record_id) == "")
if (na_id > 0) {
  df <- df |> filter(!is.na(record_id) & str_trim(record_id) != "")
  cat(sprintf("  Dropped %d rows with missing record_id\n", na_id))
}

cat(sprintf("  Rows: %d -> %d (dropped %d)\n", n_before, nrow(df), n_before - nrow(df)))


# ====================================================================
# 4. Type conversion
# ====================================================================
cat("\n--- 4. Converting column types ---\n")

df <- df |>
  mutate(
    record_id      = as.integer(record_id),
    n_deaths       = as.numeric(n_deaths),
    latitude       = as.numeric(latitude),
    longitude      = as.numeric(longitude),
    crossing_border = as.numeric(crossing_border)
  )

# Convert all flag columns to numeric
flag_col_patterns <- c("^cflag_", "^manner_", "^suicide_", "^torg_",
                        "^tmeans_", "^state_", "^where_")
flag_cols_all <- grep(paste(flag_col_patterns, collapse = "|"),
                       names(df), value = TRUE)
for (fc in flag_cols_all) {
  df[[fc]] <- suppressWarnings(as.numeric(df[[fc]]))
}
cat(sprintf("  Converted %d flag columns to numeric\n", length(flag_cols_all)))

cat(sprintf("  record_id NAs: %d\n", sum(is.na(df$record_id))))
cat(sprintf("  n_deaths NAs: %d\n", sum(is.na(df$n_deaths))))
cat(sprintf("  latitude NAs: %d / %d\n", sum(is.na(df$latitude)), nrow(df)))
cat(sprintf("  longitude NAs: %d / %d\n", sum(is.na(df$longitude)), nrow(df)))


# ====================================================================
# 5. Parse dates
# ====================================================================
cat("\n--- 5. Parsing dates ---\n")

date_parsed <- map2_dfr(df$found_dead_text, df$sort_date, parse_united_date)

df <- df |>
  bind_cols(date_parsed) |>
  mutate(
    incident_year  = year(incident_date_clean),
    incident_month = month(incident_date_clean)
  )

precision_dist <- table(df$incident_date_precision, useNA = "ifany")
cat("  Date precision distribution:\n")
for (p in names(precision_dist)) {
  cat(sprintf("    %-12s: %d\n", p, precision_dist[p]))
}

date_range <- range(df$incident_date_clean, na.rm = TRUE)
cat(sprintf("  Date range: %s to %s\n",
            as.character(date_range[1]), as.character(date_range[2])))
cat(sprintf("  Unparseable dates (NA): %d / %d\n",
            sum(is.na(df$incident_date_clean)), nrow(df)))


# ====================================================================
# 6. Collapse flag columns into tidy categoricals
# ====================================================================
cat("\n--- 6. Collapsing flag columns ---\n")

# --- Manner of death ---
manner_cols <- grep("^manner_", names(df), value = TRUE)
df$manner_of_death <- collapse_flags(df, manner_cols)
# Clean up names: remove prefix
df$manner_of_death <- str_remove(df$manner_of_death, "^manner_")
cat(sprintf("  manner_of_death: %d non-NA\n", sum(!is.na(df$manner_of_death))))

# --- Suicide manner ---
suicide_cols <- grep("^suicide_", names(df), value = TRUE)
df$suicide_manner <- collapse_flags(df, suicide_cols)
df$suicide_manner <- str_remove(df$suicide_manner, "^suicide_")
cat(sprintf("  suicide_manner: %d non-NA\n", sum(!is.na(df$suicide_manner))))

# --- Transport organisation ---
torg_cols <- grep("^torg_", names(df), value = TRUE)
df$transport_org <- collapse_flags(df, torg_cols)
df$transport_org <- str_remove(df$transport_org, "^torg_")
cat(sprintf("  transport_org: %d non-NA\n", sum(!is.na(df$transport_org))))

# --- Transport means ---
tmeans_cols <- grep("^tmeans_", names(df), value = TRUE)
df$transport_means <- collapse_flags(df, tmeans_cols)
df$transport_means <- str_remove(df$transport_means, "^tmeans_")
cat(sprintf("  transport_means: %d non-NA\n", sum(!is.na(df$transport_means))))

# --- State services ---
state_cols <- grep("^state_", names(df), value = TRUE)
df$state_services <- collapse_flags(df, state_cols)
df$state_services <- str_remove(df$state_services, "^state_")
cat(sprintf("  state_services: %d non-NA\n", sum(!is.na(df$state_services))))

# --- Where died ---
where_cols <- grep("^where_", names(df), value = TRUE)
df$where_died <- collapse_flags(df, where_cols)
df$where_died <- str_remove(df$where_died, "^where_")
cat(sprintf("  where_died: %d non-NA\n", sum(!is.na(df$where_died))))

# --- Crossing countries (semicolon-separated list) ---
cflag_cols <- grep("^cflag_", names(df), value = TRUE)
df$crossing_countries <- map_chr(seq_len(nrow(df)), function(i) {
  vals <- as.numeric(df[i, cflag_cols, drop = TRUE])
  active <- cflag_cols[!is.na(vals) & vals > 0]
  if (length(active) == 0L) NA_character_
  else paste(str_remove(active, "^cflag_"), collapse = ";")
})
cat(sprintf("  crossing_countries: %d non-NA\n",
            sum(!is.na(df$crossing_countries))))


# ====================================================================
# 7. Create CMR indicator
# ====================================================================
cat("\n--- 7. Creating CMR indicator ---\n")

CMR_COUNTRIES_DEATH <- c("italy", "libya", "malta", "tunisia", "algeria",
                          "mediterranean")
CMR_PLACE_KEYWORDS <- "mediterranean|lampedusa|sicil|malta|libya|tunis|channel of sicily|pantelleria"
CMR_COUNTRY_CODES <- c("IT", "LY", "MT", "TN", "DZ")

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

n_cmr <- sum(df$is_cmr, na.rm = TRUE)
cat(sprintf("  CMR records: %d / %d (%.1f%%)\n",
            n_cmr, nrow(df), 100 * n_cmr / nrow(df)))

# CMR by country of death
cmr_by_country <- df |>
  filter(is_cmr) |>
  count(country_of_death, sort = TRUE)
cat("  CMR breakdown by country_of_death:\n")
for (i in seq_len(min(10, nrow(cmr_by_country)))) {
  cat(sprintf("    %-20s: %d\n", cmr_by_country$country_of_death[i],
              cmr_by_country$n[i]))
}


# ====================================================================
# 8. Select final columns
# ====================================================================
cat("\n--- 8. Selecting final columns ---\n")

united_clean <- df |>
  transmute(
    record_id,
    incident_date_clean,
    incident_date_raw,
    incident_date_precision,
    incident_year,
    incident_month,
    n_deaths,
    name_gender_age,
    region_of_origin,
    cause_of_death_text,
    source_text,
    country_of_death,
    place_of_death,
    latitude,
    longitude,
    crossing_countries,
    manner_of_death,
    suicide_manner,
    transport_means,
    transport_org,
    state_services,
    where_died,
    event_group,
    is_cmr,
    weblink,
    long_description
  )

cat(sprintf("  Final dimensions: %d rows x %d cols\n",
            nrow(united_clean), ncol(united_clean)))


# ====================================================================
# 9. Validation
# ====================================================================
cat("\n--- 9. Validation ---\n")

# Null record IDs
null_ids <- sum(is.na(united_clean$record_id))
cat(sprintf("  Null record IDs: %d  %s\n", null_ids,
            ifelse(null_ids == 0, "PASS", "WARN")))

# Duplicate record IDs
dup_ids <- sum(duplicated(united_clean$record_id))
cat(sprintf("  Duplicate record IDs: %d  %s\n", dup_ids,
            ifelse(dup_ids == 0, "PASS", "WARN")))

# Date coverage
na_dates <- sum(is.na(united_clean$incident_date_clean))
cat(sprintf("  Unparseable dates: %d / %d (%.2f%%)  %s\n",
            na_dates, nrow(united_clean),
            100 * na_dates / nrow(united_clean),
            ifelse(na_dates < 10, "PASS", "WARN")))

# Coordinate coverage
na_coords <- sum(is.na(united_clean$latitude) | is.na(united_clean$longitude))
cat(sprintf("  Missing coordinates: %d / %d (%.2f%%)\n",
            na_coords, nrow(united_clean),
            100 * na_coords / nrow(united_clean)))

# Year range
cat(sprintf("  Year range: %d to %d\n",
            min(united_clean$incident_year, na.rm = TRUE),
            max(united_clean$incident_year, na.rm = TRUE)))

# Total deaths
cat(sprintf("  Total deaths: %s\n",
            format(sum(united_clean$n_deaths, na.rm = TRUE), big.mark = ",")))

# CMR summary
cat(sprintf("  CMR records: %d\n", sum(united_clean$is_cmr, na.rm = TRUE)))
cat(sprintf("  CMR deaths: %s\n",
            format(sum(united_clean$n_deaths[united_clean$is_cmr], na.rm = TRUE),
                   big.mark = ",")))

# Yearly death counts
cat("\n  Yearly death totals:\n")
yearly <- united_clean |>
  group_by(incident_year) |>
  summarise(incidents = n(),
            deaths = sum(n_deaths, na.rm = TRUE),
            .groups = "drop") |>
  arrange(incident_year)
for (i in seq_len(nrow(yearly))) {
  cat(sprintf("    %d: %4d incidents, %6s deaths\n",
              yearly$incident_year[i],
              yearly$incidents[i],
              format(yearly$deaths[i], big.mark = ",")))
}


# ====================================================================
# 10. Save
# ====================================================================
cat("\n--- 10. Saving ---\n")

out_path <- file.path(BASE_DIR, "data", "processed", "united_incidents.RDS")
saveRDS(united_clean, out_path)
cat(sprintf("  Saved to: %s\n", out_path))
cat(sprintf("  File size: %.1f KB\n", file.size(out_path) / 1024))

cat("\n============================================================\n")
cat("DONE — UNITED cleaning complete\n")
cat("============================================================\n")
