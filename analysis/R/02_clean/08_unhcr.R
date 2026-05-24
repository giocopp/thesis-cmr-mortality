# UNHCR daily Italian arrivals; fill calendar gaps as zero-arrival days.

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()

raw <- read_csv(file.path(BASE_DIR, "data", "raw", "unhcr",
                            "unhcr_daily_arrivals_italy.csv"),
                 show_col_types = FALSE) |>
  transmute(date = as.Date(data_date), arrivals = individuals)

all_dates <- tibble(date = seq(min(raw$date), max(raw$date), by = "day"))

arrivals_clean <- all_dates |>
  left_join(raw, by = "date") |>
  replace_na(list(arrivals = 0L))

saveRDS(arrivals_clean, file.path(BASE_DIR, "data", "processed",
                                    "unhcr_daily_arrivals.RDS"))
