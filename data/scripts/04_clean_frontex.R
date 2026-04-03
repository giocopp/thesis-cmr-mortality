# 04_clean_frontex.R
# ==================
# Clean Frontex Themis incident-level data.
#
# Reads the raw Excel, parses dates, classifies boat types,
# standardizes column names to snake_case.
# Keeps ALL incidents (CMR filtering is done in analysis scripts).
#
# Input:  data/raw/frontex/pad-194_themis_2014_2023.xlsx
# Output: data/processed/frontex_incidents.RDS

library(dplyr)
library(readxl)
library(lubridate)
library(stringr)

BASE_DIR <- here::here()

cat("============================================================\n")
cat("CLEAN FRONTEX THEMIS DATA\n")
cat("============================================================\n\n")

# ── 1. Load raw data ─────────────────────────────────────
cat("--- 1. Loading Frontex Themis ---\n")

frx_raw <- read_excel(
  file.path(BASE_DIR, "data", "raw", "frontex",
            "pad-194_themis_2014_2023.xlsx"),
  sheet = "Sheet1"
)

cat(sprintf("  Raw: %d incidents, %d columns\n", nrow(frx_raw), ncol(frx_raw)))
cat(sprintf("  Columns: %s\n", paste(names(frx_raw), collapse = ", ")))

# ── 2. Clean and standardize ─────────────────────────────
cat("\n--- 2. Cleaning ---\n")

frx <- frx_raw %>%
  transmute(
    date                 = as.Date(DetectionDate),
    country_of_departure = CountryOfDeparture,
    transport_type       = TransportType,
    boat_category        = case_when(
      grepl("inflatable|rubber|zodiac|dinghy", TransportType, ignore.case = TRUE) ~ "Inflatable",
      grepl("wooden|wood", TransportType, ignore.case = TRUE) ~ "Wooden",
      grepl("metal|makeshift", TransportType, ignore.case = TRUE) ~ "Metal",
      TRUE ~ "Other"
    ),
    num_persons          = num_total_persons,
    num_deaths           = num_DeathCases,
    num_migrants         = num_total_irreg_migrants,
    sar_flag             = (SAR == "Yes"),
    in_op_area           = (ReferenceToOpArea == "in"),
    operation_name       = OperationName
  )

cat(sprintf("  Date range: %s to %s\n", min(frx$date), max(frx$date)))

# ── 3. Summary and validation ────────────────────────────
cat("\n--- 3. Validation ---\n")

cat(sprintf("  Total incidents: %d\n", nrow(frx)))
cat(sprintf("  Departure countries:\n"))
dep_counts <- sort(table(frx$country_of_departure), decreasing = TRUE)
for (nm in names(dep_counts)) {
  cat(sprintf("    %s: %d\n", nm, dep_counts[[nm]]))
}

cat(sprintf("\n  Boat type distribution:\n"))
boat_counts <- table(frx$boat_category)
for (nm in names(boat_counts)) {
  cat(sprintf("    %s: %d\n", nm, boat_counts[[nm]]))
}

cat(sprintf("\n  SAR incidents: %d (%.1f%%)\n",
    sum(frx$sar_flag, na.rm = TRUE),
    100 * mean(frx$sar_flag, na.rm = TRUE)))

na_dates <- sum(is.na(frx$date))
if (na_dates > 0) {
  cat(sprintf("  WARNING: %d rows with NA date\n", na_dates))
} else {
  cat("  No NA dates — PASSED\n")
}

# ── 4. Save ──────────────────────────────────────────────
saveRDS(frx, file.path(BASE_DIR, "data", "processed", "frontex_incidents.RDS"))
cat(sprintf("\nSaved: data/processed/frontex_incidents.RDS (%d rows)\n", nrow(frx)))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
