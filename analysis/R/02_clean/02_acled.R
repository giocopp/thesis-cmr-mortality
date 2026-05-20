# 02_clean_acled.R
# ================
# Process ACLED weekly aggregated data for Libya and Tunisia.
# Produces a daily panel where each day carries its week's conflict values.
#
# Following Rodriguez-Sanchez et al. (2023):
#   conflict = Battles + Explosions/Remote violence + Violence against civilians
#
# Input:  data/raw/acled/Africa_aggregated_data_up_to_week_of-2026-03-21.xlsx
# Output: data/processed/acled_daily.RDS

library(dplyr)
library(tidyr)
library(lubridate)
library(readxl)

BASE_DIR <- here::here()
CONFLICT_TYPES <- c("Battles", "Explosions/Remote violence",
                     "Violence against civilians")

cat("============================================================\n")
cat("BUILD ACLED WEEKLY CONFLICT DATA (LIBYA + TUNISIA)\n")
cat("============================================================\n\n")

# ── 1. Load and filter ────────────────────────────────────
cat("--- 1. Loading ACLED data ---\n")

raw <- read_excel(file.path(BASE_DIR, "data", "raw", "acled",
                             "Africa_aggregated_data_up_to_week_of-2026-03-21.xlsx"),
                  sheet = "Sheet1") |>
  mutate(week_date = as.Date(WEEK))

cat(sprintf("  Raw: %d rows, %s to %s\n", nrow(raw),
    min(raw$week_date), max(raw$week_date)))

# Filter to Libya and Tunisia
lt <- raw |>
  filter(COUNTRY %in% c("Libya", "Tunisia"))
cat(sprintf("  Libya + Tunisia: %d rows\n", nrow(lt)))

# ── 2. Aggregate to country-week level ────────────────────
cat("\n--- 2. Aggregating to country-week ---\n")

# Separate event types for each country
country_week <- lt |>
  group_by(week_date, COUNTRY, EVENT_TYPE) |>
  summarise(events = sum(EVENTS, na.rm = TRUE),
            fatalities = sum(FATALITIES, na.rm = TRUE),
            .groups = "drop")

# Pivot: one row per week with columns for each country × event type
make_weekly <- function(country_name, prefix) {
  cw <- country_week |>
    filter(COUNTRY == country_name) |>
    select(-COUNTRY)

  # Conflict composite (R-S definition)
  conflict <- cw |>
    filter(EVENT_TYPE %in% CONFLICT_TYPES) |>
    group_by(week_date) |>
    summarise(conflict = sum(events),
              conflict_fatalities = sum(fatalities), .groups = "drop")

  # Separated event types
  battles <- cw |> filter(EVENT_TYPE == "Battles") |>
    select(week_date, battles = events, battles_fat = fatalities)
  expvio <- cw |> filter(EVENT_TYPE == "Explosions/Remote violence") |>
    select(week_date, expvio = events, expvio_fat = fatalities)
  violciv <- cw |> filter(EVENT_TYPE == "Violence against civilians") |>
    select(week_date, violciv = events, violciv_fat = fatalities)
  protests <- cw |> filter(EVENT_TYPE == "Protests") |>
    select(week_date, protests = events)
  riots <- cw |> filter(EVENT_TYPE == "Riots") |>
    select(week_date, riots = events)

  # Merge all for this country
  out <- conflict |>
    full_join(battles, by = "week_date") |>
    full_join(expvio, by = "week_date") |>
    full_join(violciv, by = "week_date") |>
    full_join(protests, by = "week_date") |>
    full_join(riots, by = "week_date")

  # Replace NAs with 0 (weeks with no events of a given type)
  out <- out |> mutate(across(-week_date, ~ replace_na(.x, 0L)))

  # Prefix column names
  names(out)[-1] <- paste0(prefix, "_", names(out)[-1])
  out
}

libya_weekly <- make_weekly("Libya", "libya")
tunisia_weekly <- make_weekly("Tunisia", "tunisia")

# Merge Libya and Tunisia
weekly <- libya_weekly |>
  full_join(tunisia_weekly, by = "week_date") |>
  mutate(across(-week_date, ~ replace_na(.x, 0L))) |>
  arrange(week_date)

cat(sprintf("  Weekly panel: %d weeks (%s to %s)\n",
    nrow(weekly), min(weekly$week_date), max(weekly$week_date)))

# ── 3. Expand to daily (option 4: each day gets its week's value) ─
cat("\n--- 3. Expanding to daily ---\n")

# Create daily spine
daily_spine <- tibble(date = seq(as.Date("2014-01-01"), as.Date("2025-12-31"), by = "day")) |>
  mutate(week_date = floor_date(date, "week", week_start = 6))
# ACLED weeks appear to start on Saturday based on the data; let me verify

# Check what day of week the ACLED dates fall on
acled_dow <- unique(wday(weekly$week_date, label = TRUE))
cat(sprintf("  ACLED week dates fall on: %s\n", paste(acled_dow, collapse = ", ")))

# Use floor_date with the matching week start
# wday: 1=Sun, 2=Mon, ..., 7=Sat
acled_wday <- unique(wday(weekly$week_date))
cat(sprintf("  ACLED wday numbers: %s\n", paste(acled_wday, collapse = ", ")))

# Assign each day to its ACLED week by finding the most recent ACLED week date
# This is more robust than assuming a specific week start
daily_spine <- daily_spine |>
  select(date) |>
  mutate(week_date = as.Date(cut(date, breaks = sort(unique(weekly$week_date)),
                                  include.lowest = TRUE, right = FALSE)))

# Simpler approach: for each day, find the ACLED week it belongs to
# by matching to the closest preceding week_date
daily_spine <- tibble(date = seq(as.Date("2014-01-01"), as.Date("2025-12-31"), by = "day"))

# For each day, find the ACLED week: the largest week_date <= date
weekly_dates <- sort(unique(weekly$week_date))
daily_spine <- daily_spine |>
  rowwise() |>
  mutate(week_date = max(weekly_dates[weekly_dates <= date])) |>
  ungroup()

# Merge weekly data onto daily spine
acled_daily <- daily_spine |>
  left_join(weekly, by = "week_date")

# Check for NAs (days before first ACLED week or after last)
n_na <- sum(is.na(acled_daily$libya_conflict))
cat(sprintf("  Daily panel: %d days\n", nrow(acled_daily)))
cat(sprintf("  Days with NA conflict (outside ACLED range): %d\n", n_na))

# ── 4. Validation ─────────────────────────────────────────
cat("\n--- 4. Validation ---\n")

cat("\n  Libya conflict by year (weekly totals):\n")
acled_daily |>
  mutate(year = year(date)) |>
  filter(year >= 2014, year <= 2025) |>
  group_by(year) |>
  summarise(
    conflict_events = sum(libya_conflict, na.rm = TRUE) / 7,  # undo daily repeat
    fatalities = sum(libya_conflict_fatalities, na.rm = TRUE) / 7,
    .groups = "drop"
  ) |>
  mutate(conflict_events = round(conflict_events),
         fatalities = round(fatalities)) |>
  print(n = 15)

cat("\n  Tunisia conflict by year (weekly totals):\n")
acled_daily |>
  mutate(year = year(date)) |>
  filter(year >= 2014, year <= 2025) |>
  group_by(year) |>
  summarise(
    conflict_events = sum(tunisia_conflict, na.rm = TRUE) / 7,
    fatalities = sum(tunisia_conflict_fatalities, na.rm = TRUE) / 7,
    .groups = "drop"
  ) |>
  mutate(conflict_events = round(conflict_events),
         fatalities = round(fatalities)) |>
  print(n = 15)

# Verify daily repetition is correct: all days in the same week should match
check <- acled_daily |>
  filter(date >= "2017-01-01", date <= "2017-01-14") |>
  select(date, week_date, libya_conflict, tunisia_conflict)
cat("\n  Sample: first 2 weeks of 2017:\n")
print(check, n = 14)

# ── 5. Save ───────────────────────────────────────────────
cat("\n--- 5. Saving ---\n")

saveRDS(acled_daily, file.path(BASE_DIR, "data", "processed", "acled_daily.RDS"))
cat(sprintf("Saved: data/processed/acled_daily.RDS (%d rows x %d cols)\n",
    nrow(acled_daily), ncol(acled_daily)))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
