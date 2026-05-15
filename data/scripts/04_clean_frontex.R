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
    incident_id          = IncidentNumber,
    date                 = as.Date(DetectionDate),
    country_of_departure = CountryOfDeparture,
    transport_type       = TransportType,
    boat_category        = case_when(
      grepl("inflatable|rubber|zodiac|dinghy", TransportType, ignore.case = TRUE) ~ "Inflatable",
      grepl("wooden|wood", TransportType, ignore.case = TRUE) ~ "Wooden",
      grepl("metal|makeshift", TransportType, ignore.case = TRUE) ~ "Metal",
      grepl("fibre glass|fiber glass|fibreglass|fiberglass",
            TransportType, ignore.case = TRUE) ~ "Fibre glass",
      TRUE ~ "Other"
    ),
    num_persons          = num_total_persons,
    num_deaths           = num_DeathCases,
    num_migrants         = num_total_irreg_migrants,
    n_transport_means    = num_TransportMeansNumber,
    sar_ops             = case_when(SAR == "Yes" ~ TRUE,
                                     SAR == "No" ~ FALSE,
                                     TRUE ~ NA),
    in_op_area           = (ReferenceToOpArea == "in"),
    detection_by_frx_asset     = case_when(DetectionByFrxAsset == "Yes" ~ TRUE,
                                           DetectionByFrxAsset == "No"  ~ FALSE,
                                           TRUE                          ~ NA),
    interception_by_frx_asset  = case_when(InterceptionByFrxAsset == "Yes" ~ TRUE,
                                           InterceptionByFrxAsset == "No"  ~ FALSE,
                                           TRUE                              ~ NA),
    other_frx_asset_involvement = case_when(OtherFrontexAssetInvolvement == "Yes" ~ TRUE,
                                            OtherFrontexAssetInvolvement == "No"  ~ FALSE,
                                            TRUE                                   ~ NA),
    operation_name       = OperationName,
    detected_by          = TypeOfDetectedBy,
    intercepted_by       = TypeOfInterceptedBy,
    ngo_involved         = grepl("NGO vessel", TypeOfDetectedBy, ignore.case = TRUE) |
                           grepl("NGO vessel", TypeOfInterceptedBy, ignore.case = TRUE)
  )

# ── 2b. Classify interceptor and detector (independent of SAR flag) ─────
# Priority order handles compound "A;B" values by first match wins.
# Token-level classifiers used for multi-actor flag below.
classify_int_token <- function(tok) {
  tok <- trimws(tok)
  dplyr::case_when(
    tok == "" | is.na(tok)                                     ~ NA_character_,
    grepl("NGO vessel", tok)                                   ~ "NGO",
    grepl("EUNAVFOR", tok)                                     ~ "EU_ops",
    grepl("Marina|Mare", tok)                                  ~ "Ita_ops",
    grepl("Commercial|fishing|Merchant", tok,
          ignore.case = TRUE)                                  ~ "Commercial",
    grepl("CPV|CPB|OPV", tok)                                  ~ "EU_Coast_Guard",
    grepl("Land", tok)                                         ~ "Land_patrol",
    grepl("No interception", tok)                              ~ "No_intercept",
    TRUE                                                       ~ "Other"
  )
}
classify_det_token <- function(tok) {
  tok <- trimws(tok)
  dplyr::case_when(
    tok == "" | is.na(tok)                                     ~ NA_character_,
    grepl("FWA|HELO|RPAS|MAS", tok)                            ~ "Aerial",
    grepl("Call-", tok)                                        ~ "Call_Ext_report",
    grepl("NGO vessel", tok)                                   ~ "NGO",
    grepl("Commercial|fishing|Merchant", tok,
          ignore.case = TRUE)                                  ~ "Commercial",
    grepl("EUNAVFOR", tok)                                     ~ "EU_ops",
    grepl("Marina|Mare", tok)                                  ~ "Ita_ops",
    grepl("CPV|CPB|OPV", tok)                                  ~ "EU_Coast_Guard",
    grepl("Coastal Surveillance", tok)                         ~ "Coastal_surv",
    grepl("Land", tok)                                         ~ "Land_patrol",
    TRUE                                                       ~ "Other"
  )
}
n_distinct_cats <- function(s, classifier) {
  if (is.na(s) || s == "") return(0L)
  toks <- strsplit(s, ";")[[1]]
  cats <- classifier(toks)
  length(unique(cats[!is.na(cats)]))
}

frx <- frx %>%
  mutate(
    interceptor_type = case_when(
      is.na(intercepted_by)                                      ~ "NA",
      grepl("NGO vessel", intercepted_by)                        ~ "NGO",
      grepl("EUNAVFOR", intercepted_by)                          ~ "EU_ops",
      grepl("Marina|Mare", intercepted_by)                       ~ "Ita_ops",
      grepl("Commercial|fishing|Merchant", intercepted_by,
            ignore.case = TRUE)                                  ~ "Commercial",
      grepl("CPV|CPB|OPV", intercepted_by)                       ~ "EU_Coast_Guard",
      grepl("Land", intercepted_by)                              ~ "Land_patrol",
      grepl("No interception", intercepted_by)                   ~ "No_intercept",
      TRUE                                                       ~ "Other"
    ),
    detector_type = case_when(
      is.na(detected_by)                                         ~ "NA",
      grepl("FWA|HELO|RPAS|MAS", detected_by)                    ~ "Aerial",
      grepl("Call-", detected_by)                                ~ "Call_Ext_report",
      grepl("NGO vessel", detected_by)                           ~ "NGO",
      grepl("Commercial|fishing|Merchant", detected_by,
            ignore.case = TRUE)                                  ~ "Commercial",
      grepl("EUNAVFOR", detected_by)                             ~ "EU_ops",
      grepl("Marina|Mare", detected_by)                          ~ "Ita_ops",
      grepl("CPV|CPB|OPV", detected_by)                          ~ "EU_Coast_Guard",
      grepl("Coastal Surveillance", detected_by)                 ~ "Coastal_surv",
      grepl("Land", detected_by)                                 ~ "Land_patrol",
      TRUE                                                       ~ "Other"
    ),
    multi_actors_inv = vapply(intercepted_by,
                              function(s) n_distinct_cats(s, classify_int_token),
                              integer(1)) >= 2 |
                       vapply(detected_by,
                              function(s) n_distinct_cats(s, classify_det_token),
                              integer(1)) >= 2
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
    sum(frx$sar_ops, na.rm = TRUE),
    100 * mean(frx$sar_ops, na.rm = TRUE)))
cat(sprintf("  NGO-involved incidents: %d (%.1f%%)\n",
    sum(frx$ngo_involved, na.rm = TRUE),
    100 * mean(frx$ngo_involved, na.rm = TRUE)))

cat("\n  Interceptor type distribution:\n")
it_counts <- sort(table(frx$interceptor_type), decreasing = TRUE)
for (nm in names(it_counts)) {
  cat(sprintf("    %5d (%.1f%%)  %s\n", it_counts[[nm]],
      100 * it_counts[[nm]] / nrow(frx), nm))
}

cat("\n  Detector type distribution:\n")
dt_counts <- sort(table(frx$detector_type), decreasing = TRUE)
for (nm in names(dt_counts)) {
  cat(sprintf("    %5d (%.1f%%)  %s\n", dt_counts[[nm]],
      100 * dt_counts[[nm]] / nrow(frx), nm))
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
