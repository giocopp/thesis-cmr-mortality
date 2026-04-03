# 05_clean_unhcr.R
# ================
# Clean UNHCR daily arrivals data: fill missing dates with zero arrivals.
#
# The raw UNHCR CSV drops days with zero arrivals — those days simply
# don't appear in the file. This script creates a complete daily time
# series with explicit zeros for all calendar days.
#
# Input:  data/raw/unhcr/unhcr_daily_arrivals_italy.csv
# Output: data/processed/unhcr_daily_arrivals.RDS

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()

cat("============================================================\n")
cat("CLEAN UNHCR DAILY ARRIVALS\n")
cat("============================================================\n\n")

# Load raw data
raw <- read_csv(file.path(BASE_DIR, "data", "raw", "unhcr",
                            "unhcr_daily_arrivals_italy.csv"),
                 show_col_types = FALSE) %>%
  transmute(date = as.Date(data_date), arrivals = individuals)

cat(sprintf("Raw file: %d rows, %s to %s\n",
    nrow(raw), min(raw$date), max(raw$date)))
cat(sprintf("Days with arrivals > 0: %d\n", sum(raw$arrivals > 0)))
cat(sprintf("Days with arrivals = 0: %d\n", sum(raw$arrivals == 0)))

# Create complete daily grid from first to last date
all_dates <- tibble(date = seq(min(raw$date), max(raw$date), by = "day"))

cat(sprintf("Calendar days in range: %d\n", nrow(all_dates)))
cat(sprintf("Missing from raw: %d (these are zero-arrival days)\n",
    nrow(all_dates) - nrow(raw)))

# Merge and fill missing with 0
arrivals_clean <- all_dates %>%
  left_join(raw, by = "date") %>%
  replace_na(list(arrivals = 0L))

# Verify
cat(sprintf("\nClean data: %d rows\n", nrow(arrivals_clean)))
cat(sprintf("  Zero-arrival days: %d (%.1f%%)\n",
    sum(arrivals_clean$arrivals == 0),
    100 * mean(arrivals_clean$arrivals == 0)))
cat(sprintf("  Total arrivals: %d\n", sum(arrivals_clean$arrivals)))
cat(sprintf("  Date range: %s to %s\n",
    min(arrivals_clean$date), max(arrivals_clean$date)))

# Save
saveRDS(arrivals_clean, file.path(BASE_DIR, "data", "processed",
                                    "unhcr_daily_arrivals.RDS"))
cat("\nSaved: data/processed/unhcr_daily_arrivals.RDS\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
