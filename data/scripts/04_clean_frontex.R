# 04_clean_frontex.R
# ==================
# Clean Frontex Themis incident-level data.
#
# Reads the raw Excel, parses dates, classifies boat types,
# standardizes column names to snake_case.
# Keeps ALL incidents (CMR filtering is done in analysis scripts).
#
# Known source-data issues:
#   - 1,127 incidents (8%) have num_persons = 0. These are mostly
#     "Not SAR: Coast Guard" (571) and "Not SAR: Self-arrived" (220).
#     Possible reasons: boat detected with no one aboard, failed crossings,
#     or data entry gaps. This is a quirk of the source data.
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
    sar_flag             = case_when(SAR == "Yes" ~ TRUE,
                                     SAR == "No" ~ FALSE,
                                     TRUE ~ NA),
    in_op_area           = (ReferenceToOpArea == "in"),
    operation_name       = OperationName,
    detected_by          = TypeOfDetectedBy,
    intercepted_by       = TypeOfInterceptedBy,
    ngo_involved         = grepl("NGO vessel", TypeOfDetectedBy, ignore.case = TRUE) |
                           grepl("NGO vessel", TypeOfInterceptedBy, ignore.case = TRUE)
  )

# ── 2b. Classify event type ─────────────────────────────
# Based on SAR flag + primary interceptor (first entry before ";")
frx <- frx %>%
  mutate(
    primary_intercept = sub(";.*", "", intercepted_by),
    event_type = case_when(
      # SAR events (SAR=Yes)
      sar_flag & grepl("NGO vessel", intercepted_by)        ~ "SAR: NGO",
      sar_flag & grepl("EUNAVFOR", intercepted_by)           ~ "SAR: EU operations (IRINI)",
      sar_flag & grepl("Commercial|fishing|Merchant",
                       intercepted_by, ignore.case = TRUE)   ~ "SAR: Commercial vessels",
      sar_flag & grepl("CPV|CPB|OPV|Marina Militare|Mare Sicuro|Mare Nostrum",
                       intercepted_by)                       ~ "SAR: Italian authorities",
      sar_flag                                               ~ "SAR: Other",
      # Non-SAR events (SAR=No)
      !sar_flag & grepl("Land patrol", intercepted_by)       ~ "Not SAR: Land patrol",
      !sar_flag & grepl("No interception", intercepted_by)   ~ "Not SAR: Self-arrived",
      !sar_flag & grepl("CPV|CPB|OPV|Marina Militare",
                        intercepted_by)                      ~ "Not SAR: Coast Guard",
      !sar_flag                                              ~ "Not SAR: Other",
      # SAR flag NA — classify by interceptor
      grepl("Land patrol", intercepted_by)                   ~ "Not SAR: Land patrol",
      grepl("No interception", intercepted_by)               ~ "Not SAR: Self-arrived",
      grepl("NGO vessel", intercepted_by)                    ~ "SAR: NGO",
      grepl("CPV|CPB|OPV|Marina Militare",
            intercepted_by)                                  ~ "Not SAR: Coast Guard",
      TRUE                                                   ~ "Not SAR: Other"
    ),
    event_type_agg = case_when(
      startsWith(event_type, "SAR:")     ~ "SAR",
      startsWith(event_type, "Not SAR:") ~ "Not SAR",
      TRUE                               ~ "Unknown"
    )
  ) %>%
  select(-primary_intercept)

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
cat(sprintf("  NGO-involved incidents: %d (%.1f%%)\n",
    sum(frx$ngo_involved, na.rm = TRUE),
    100 * mean(frx$ngo_involved, na.rm = TRUE)))

cat("\n  Event type distribution:\n")
et_counts <- sort(table(frx$event_type), decreasing = TRUE)
for (nm in names(et_counts)) {
  cat(sprintf("    %5d (%.1f%%)  %s\n", et_counts[[nm]],
      100 * et_counts[[nm]] / nrow(frx), nm))
}

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
