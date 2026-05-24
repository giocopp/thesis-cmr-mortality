# ACLED weekly conflict data for Libya and Tunisia, broadcast to daily.

library(dplyr)
library(tidyr)
library(lubridate)
library(readxl)

BASE_DIR <- here::here()
CONFLICT_TYPES <- c("Battles", "Explosions/Remote violence",
                     "Violence against civilians")

# ── 1. Load and filter ──────────────────────────────────────────────────────
raw <- read_excel(file.path(BASE_DIR, "data", "raw", "acled",
                             "Africa_aggregated_data_up_to_week_of-2026-03-21.xlsx"),
                  sheet = "Sheet1") |>
  mutate(week_date = as.Date(WEEK))

lt <- raw |>
  filter(COUNTRY %in% c("Libya", "Tunisia"))

# ── 2. Aggregate to country-week ────────────────────────────────────────────
country_week <- lt |>
  group_by(week_date, COUNTRY, EVENT_TYPE) |>
  summarise(events = sum(EVENTS, na.rm = TRUE),
            fatalities = sum(FATALITIES, na.rm = TRUE),
            .groups = "drop")

make_weekly <- function(country_name, prefix) {
  cw <- country_week |>
    filter(COUNTRY == country_name) |>
    select(-COUNTRY)

  conflict <- cw |>
    filter(EVENT_TYPE %in% CONFLICT_TYPES) |>
    group_by(week_date) |>
    summarise(conflict = sum(events),
              conflict_fatalities = sum(fatalities), .groups = "drop")

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

  out <- conflict |>
    full_join(battles, by = "week_date") |>
    full_join(expvio, by = "week_date") |>
    full_join(violciv, by = "week_date") |>
    full_join(protests, by = "week_date") |>
    full_join(riots, by = "week_date") |>
    mutate(across(-week_date, ~ replace_na(.x, 0L)))

  names(out)[-1] <- paste0(prefix, "_", names(out)[-1])
  out
}

libya_weekly   <- make_weekly("Libya", "libya")
tunisia_weekly <- make_weekly("Tunisia", "tunisia")

weekly <- libya_weekly |>
  full_join(tunisia_weekly, by = "week_date") |>
  mutate(across(-week_date, ~ replace_na(.x, 0L))) |>
  arrange(week_date)

# ── 3. Expand to daily (each day inherits its week's value) ─────────────────
daily_spine  <- tibble(date = seq(as.Date("2014-01-01"), as.Date("2025-12-31"), by = "day"))
weekly_dates <- sort(unique(weekly$week_date))

daily_spine <- daily_spine |>
  rowwise() |>
  mutate(week_date = max(weekly_dates[weekly_dates <= date])) |>
  ungroup()

acled_daily <- daily_spine |>
  left_join(weekly, by = "week_date")

# ── 4. Save ─────────────────────────────────────────────────────────────────
saveRDS(acled_daily, file.path(BASE_DIR, "data", "processed", "acled_daily.RDS"))
