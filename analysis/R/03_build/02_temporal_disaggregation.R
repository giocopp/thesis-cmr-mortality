# Proportional Denton disaggregation of monthly LCG/TCG interceptions to
# daily, using Frontex daily departures as the high-frequency indicator.

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()
CMR_DEPARTURES <- c("Libya", "Tunisia", "Algeria")

# ── 1. Monthly interceptions ────────────────────────────────────────────────
monthly_raw <- readRDS(file.path(BASE_DIR, "data", "processed",
                                  "iom_med_crossings_monthly.RDS"))

monthly <- monthly_raw |>
  transmute(
    ym          = as.Date(date),
    lcg_monthly = replace_na(as.numeric(interceptions_by_libyan_coast_guard), 0),
    tcg_monthly = replace_na(as.numeric(interceptions_by_tunisian_coast_guard), 0)
  )

# ── 2. Daily Frontex indicator ──────────────────────────────────────────────
frx_raw <- readRDS(file.path(BASE_DIR, "data", "processed",
                              "frontex_incidents.RDS"))

frx <- frx_raw |> filter(country_of_departure %in% CMR_DEPARTURES)

FRX_END <- max(frx$date)

frx_daily <- frx |>
  group_by(date) |>
  summarise(
    indicator_lcg = sum(num_persons[country_of_departure == "Libya"], na.rm = TRUE),
    indicator_tcg = sum(num_persons[country_of_departure == "Tunisia"], na.rm = TRUE),
    .groups = "drop"
  )

# ── 3. Date spine ───────────────────────────────────────────────────────────
last_full_month <- floor_date(FRX_END, "month")
days_in_last    <- as.integer(FRX_END - last_full_month) + 1
panel_end_month <- if (days_in_last < 15) last_full_month - months(1) else last_full_month
panel_end_date  <- min(FRX_END, ceiling_date(panel_end_month, "month") - days(1))

spine <- tibble(date = seq(as.Date("2014-01-01"), panel_end_date, by = "day")) |>
  mutate(ym = floor_date(date, "month")) |>
  left_join(frx_daily, by = "date") |>
  replace_na(list(indicator_lcg = 0, indicator_tcg = 0)) |>
  left_join(monthly, by = "ym") |>
  replace_na(list(lcg_monthly = 0, tcg_monthly = 0))

# ── 4. Proportional Denton disaggregation ───────────────────────────────────
disagg <- spine |>
  group_by(ym) |>
  mutate(
    S_lcg  = sum(indicator_lcg),
    S_tcg  = sum(indicator_tcg),
    n_days = n(),
    lcg_pushbacks = case_when(
      lcg_monthly == 0 ~ 0,
      S_lcg > 0        ~ lcg_monthly * (indicator_lcg / S_lcg),
      TRUE             ~ lcg_monthly / n_days
    ),
    disagg_method_lcg = case_when(
      lcg_monthly == 0 ~ "zero",
      S_lcg > 0        ~ "denton",
      TRUE             ~ "uniform"
    ),
    tcg_pushbacks = case_when(
      tcg_monthly == 0 ~ 0,
      S_tcg > 0        ~ tcg_monthly * (indicator_tcg / S_tcg),
      TRUE             ~ tcg_monthly / n_days
    ),
    disagg_method_tcg = case_when(
      tcg_monthly == 0 ~ "zero",
      S_tcg > 0        ~ "denton",
      TRUE             ~ "uniform"
    )
  ) |>
  ungroup()

# Largest-remainder rounding keeps monthly totals exact for integer counts.
lr_round_month <- function(x, target_total) {
  target_total <- as.integer(round(target_total))
  if (target_total == 0L || length(x) == 0L) return(rep(0L, length(x)))
  floor_x  <- floor(x)
  residual <- target_total - sum(floor_x)
  if (residual > 0) {
    idx <- order(x - floor_x, decreasing = TRUE)[seq_len(residual)]
    floor_x[idx] <- floor_x[idx] + 1
  } else if (residual < 0) {
    idx <- order(x - floor_x)[seq_len(-residual)]
    floor_x[idx] <- pmax(floor_x[idx] - 1, 0)
  }
  as.integer(floor_x)
}

disagg <- disagg |>
  group_by(ym) |>
  mutate(
    lcg_pushbacks = lr_round_month(lcg_pushbacks, first(lcg_monthly)),
    tcg_pushbacks = lr_round_month(tcg_pushbacks, first(tcg_monthly))
  ) |>
  ungroup() |>
  mutate(lcg_tcg_pushbacks = lcg_pushbacks + tcg_pushbacks)

# ── 5. Validation ───────────────────────────────────────────────────────────
monthly_check <- disagg |>
  group_by(ym) |>
  summarise(
    lcg_diff = abs(sum(lcg_pushbacks) - first(lcg_monthly)),
    tcg_diff = abs(sum(tcg_pushbacks) - first(tcg_monthly)),
    .groups = "drop"
  )

stopifnot(
  "LCG daily sums do not match monthly totals" = all(monthly_check$lcg_diff < 0.01),
  "TCG daily sums do not match monthly totals" = all(monthly_check$tcg_diff < 0.01),
  "Negative LCG daily values found"             = all(disagg$lcg_pushbacks >= 0),
  "Negative TCG daily values found"             = all(disagg$tcg_pushbacks >= 0),
  "NA values in lcg_pushbacks"                  = !any(is.na(disagg$lcg_pushbacks)),
  "NA values in tcg_pushbacks"                  = !any(is.na(disagg$tcg_pushbacks)),
  "NA values in lcg_tcg_pushbacks"              = !any(is.na(disagg$lcg_tcg_pushbacks)),
  "Duplicate dates found"                       = !any(duplicated(disagg$date))
)

# ── 6. Save ─────────────────────────────────────────────────────────────────
output <- disagg |>
  select(date, lcg_pushbacks, tcg_pushbacks, lcg_tcg_pushbacks,
         disagg_method_lcg, disagg_method_tcg)

out_dir <- file.path(BASE_DIR, "analysis", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(output, file.path(out_dir, "interceptions_daily_disagg.RDS"))
