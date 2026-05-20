# 07_crosscheck_united_iom.R
# ===========================
# Cross-check UNITED refugee deaths against IOM Missing Migrants Project.
#
# Purpose: Understand the overlap between the two data sources before
# deciding how to use UNITED data in the analysis. This script is
# DESCRIPTIVE — it reports what matches and what does not, without
# choosing a deduplication strategy.
#
# Matching strategy:
#   1. Filter both datasets to CMR + overlapping period (2014-2025)
#   2. Fuzzy join: date window (+-3 days) + geographic proximity
#      (haversine < 50 km) + death count similarity (ratio in [0.5, 2.0])
#   3. Classify: matched, UNITED-only, IOM-only
#   4. Aggregate comparison by year
#
# Input:
#   data/processed/united_incidents.RDS
#   data/processed/iom_mmp_incidents.RDS
#
# Output:
#   data/processed/crosscheck_united_iom.RDS
#   output/tables/crosscheck_united_iom_summary.txt

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(purrr)

BASE_DIR <- here::here()

cat("============================================================\n")
cat("CROSS-CHECK: UNITED vs IOM\n")
cat("============================================================\n\n")


# ====================================================================
# Matching parameters
# ====================================================================
DATE_WINDOW  <- 3L       # days
DIST_MAX_KM  <- 50       # km
COUNT_RATIO  <- c(0.5, 2.0)  # death count ratio bounds


# ====================================================================
# Haversine distance (reused from 04b_merge_triton_coords.R)
# ====================================================================
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; rd <- pi / 180
  dlat <- (lat2 - lat1) * rd
  dlon <- (lon2 - lon1) * rd
  a <- sin(dlat / 2)^2 + cos(lat1 * rd) * cos(lat2 * rd) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}


# ====================================================================
# 1. Load and filter both datasets
# ====================================================================
cat("--- 1. Loading datasets ---\n")

# UNITED: CMR records in overlap period
united_all <- readRDS(file.path(BASE_DIR, "data", "processed",
                                 "united_incidents.RDS"))
united <- united_all |>
  filter(is_cmr,
         !is.na(incident_date_clean),
         incident_year >= 2014L,
         incident_year <= 2025L)

cat(sprintf("  UNITED: %d total -> %d CMR in 2014-2025\n",
            nrow(united_all), nrow(united)))

# IOM: Central Mediterranean, matching _helpers.R conventions
iom_all <- readRDS(file.path(BASE_DIR, "data", "processed",
                              "iom_mmp_incidents.RDS"))
iom <- iom_all |>
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) %in% c("incident", "split incident"),
         `Country of Incident` %in% c("Algeria", "Italy", "Libya",
                                       "Malta", "Tunisia")) |>
  transmute(
    iom_id     = `Main ID`,
    date       = as.Date(incident_date_clean),
    n_dead     = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE),
    lat        = as.numeric(Latitude),
    lon        = as.numeric(Longitude),
    cause_cat  = `Cause of death (category)`,
    location   = `Location of death`
  ) |>
  filter(!is.na(date),
         year(date) >= 2014L,
         year(date) <= 2025L)

cat(sprintf("  IOM: %d total -> %d CMR in 2014-2025\n",
            nrow(iom_all), nrow(iom)))

cat(sprintf("\n  UNITED deaths in period: %s\n",
            format(sum(united$n_deaths, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  IOM deaths in period: %s\n",
            format(sum(iom$n_dead, na.rm = TRUE), big.mark = ",")))


# ====================================================================
# 2. Fuzzy matching
# ====================================================================
cat("\n--- 2. Fuzzy matching ---\n")
cat(sprintf("  Parameters: date_window=%d days, dist_max=%d km, count_ratio=[%.1f, %.1f]\n",
            DATE_WINDOW, DIST_MAX_KM, COUNT_RATIO[1], COUNT_RATIO[2]))

# Split UNITED into records with and without coordinates
united_geo <- united |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  arrange(incident_date_clean)

united_nogeo <- united |>
  filter(is.na(latitude) | is.na(longitude))

cat(sprintf("  UNITED with coords: %d, without: %d\n",
            nrow(united_geo), nrow(united_nogeo)))

# Pre-sort IOM by date for efficient windowing
iom_sorted <- iom |>
  filter(!is.na(lat), !is.na(lon)) |>
  arrange(date)

# --- Geographic + date + count matching ---
matches <- vector("list", nrow(united_geo))
match_count <- 0L

pb_step <- max(1L, nrow(united_geo) %/% 20L)

for (i in seq_len(nrow(united_geo))) {
  u <- united_geo[i, ]

  # Date window candidates
  date_diff_days <- as.integer(iom_sorted$date - u$incident_date_clean)
  candidates_idx <- which(abs(date_diff_days) <= DATE_WINDOW)

  if (length(candidates_idx) == 0L) next

  candidates <- iom_sorted[candidates_idx, ]
  candidates$date_diff <- abs(date_diff_days[candidates_idx])

  # Geographic filter
  candidates$dist_km <- haversine_km(u$latitude, u$longitude,
                                      candidates$lat, candidates$lon)
  candidates <- candidates |> filter(!is.na(dist_km), dist_km <= DIST_MAX_KM)

  if (nrow(candidates) == 0L) next

  # Death count similarity
  candidates$count_ratio <- u$n_deaths / pmax(candidates$n_dead, 1)
  candidates <- candidates |>
    filter(count_ratio >= COUNT_RATIO[1], count_ratio <= COUNT_RATIO[2])

  if (nrow(candidates) == 0L) next

  # Pick best match: minimize combined score
  best <- candidates |>
    mutate(score = date_diff + dist_km / DIST_MAX_KM) |>
    slice_min(score, n = 1L, with_ties = FALSE)

  match_count <- match_count + 1L
  matches[[i]] <- tibble(
    united_id   = u$record_id,
    iom_id      = best$iom_id,
    date_diff   = best$date_diff,
    dist_km     = best$dist_km,
    count_ratio = best$count_ratio,
    match_type  = ifelse(best$date_diff == 0L & best$dist_km < 5,
                         "strong", "fuzzy")
  )

  if (i %% pb_step == 0L) {
    cat(sprintf("    Progress: %d / %d (%d matches so far)\n",
                i, nrow(united_geo), match_count))
  }
}

match_df <- bind_rows(matches)

cat(sprintf("\n  Geographic matches: %d / %d UNITED records (%.1f%%)\n",
            nrow(match_df), nrow(united_geo),
            100 * nrow(match_df) / nrow(united_geo)))

# --- Date-only matching for records without coordinates ---
date_only_matches <- vector("list", nrow(united_nogeo))
date_only_count <- 0L

for (i in seq_len(nrow(united_nogeo))) {
  u <- united_nogeo[i, ]

  date_diff_days <- as.integer(iom$date - u$incident_date_clean)
  candidates_idx <- which(abs(date_diff_days) <= DATE_WINDOW)

  if (length(candidates_idx) == 0L) next

  candidates <- iom[candidates_idx, ]
  candidates$date_diff <- abs(date_diff_days[candidates_idx])

  # Count similarity only
  candidates$count_ratio <- u$n_deaths / pmax(candidates$n_dead, 1)
  candidates <- candidates |>
    filter(count_ratio >= COUNT_RATIO[1], count_ratio <= COUNT_RATIO[2])

  if (nrow(candidates) == 0L) next

  best <- candidates |>
    slice_min(date_diff, n = 1L, with_ties = FALSE)

  date_only_count <- date_only_count + 1L
  date_only_matches[[i]] <- tibble(
    united_id   = u$record_id,
    iom_id      = best$iom_id,
    date_diff   = best$date_diff,
    dist_km     = NA_real_,
    count_ratio = best$count_ratio,
    match_type  = "date_only"
  )
}

date_only_df <- bind_rows(date_only_matches)

if (nrow(date_only_df) > 0L) {
  match_df <- bind_rows(match_df, date_only_df)
  cat(sprintf("  Date-only matches (no coords): %d\n", nrow(date_only_df)))
}

# Check for IOM records matched more than once
iom_match_freq <- table(match_df$iom_id)
multi_matched <- sum(iom_match_freq > 1)
cat(sprintf("  IOM records matched >1 time: %d\n", multi_matched))


# ====================================================================
# 3. Match quality summary
# ====================================================================
cat("\n--- 3. Match quality ---\n")

if (nrow(match_df) > 0L) {
  type_dist <- table(match_df$match_type)
  for (t in names(type_dist)) {
    cat(sprintf("  %-12s: %d\n", t, type_dist[t]))
  }

  geo_matches <- match_df |> filter(match_type != "date_only")
  if (nrow(geo_matches) > 0L) {
    cat(sprintf("\n  Distance (geo matches): median=%.1f km, mean=%.1f km, max=%.1f km\n",
                median(geo_matches$dist_km, na.rm = TRUE),
                mean(geo_matches$dist_km, na.rm = TRUE),
                max(geo_matches$dist_km, na.rm = TRUE)))
    cat(sprintf("  Date diff (all matches): median=%d days, mean=%.1f days\n",
                median(match_df$date_diff),
                mean(match_df$date_diff)))
  }
}


# ====================================================================
# 4. Classify all records
# ====================================================================
cat("\n--- 4. Classifying records ---\n")

# UNITED classification
united_result <- united |>
  mutate(match_status = ifelse(record_id %in% match_df$united_id,
                               "matched", "united_only"))

# IOM classification
iom_matched_ids <- unique(match_df$iom_id)
iom_result <- iom |>
  mutate(match_status = ifelse(iom_id %in% iom_matched_ids,
                               "matched", "iom_only"))

n_united_matched <- sum(united_result$match_status == "matched")
n_united_only    <- sum(united_result$match_status == "united_only")
n_iom_matched    <- sum(iom_result$match_status == "matched")
n_iom_only       <- sum(iom_result$match_status == "iom_only")

cat(sprintf("  UNITED: %d matched, %d UNITED-only\n",
            n_united_matched, n_united_only))
cat(sprintf("  IOM:    %d matched, %d IOM-only\n",
            n_iom_matched, n_iom_only))
cat(sprintf("  UNITED match rate: %.1f%%\n",
            100 * n_united_matched / nrow(united_result)))
cat(sprintf("  IOM match rate: %.1f%%\n",
            100 * n_iom_matched / nrow(iom_result)))


# ====================================================================
# 5. Yearly comparison
# ====================================================================
cat("\n--- 5. Yearly comparison ---\n")

yearly_united <- united |>
  group_by(year = incident_year) |>
  summarise(united_incidents = n(),
            united_deaths = sum(n_deaths, na.rm = TRUE),
            .groups = "drop")

yearly_iom <- iom |>
  group_by(year = year(date)) |>
  summarise(iom_incidents = n(),
            iom_deaths = sum(n_dead, na.rm = TRUE),
            .groups = "drop")

yearly_united_matched <- united_result |>
  filter(match_status == "matched") |>
  group_by(year = incident_year) |>
  summarise(united_matched = n(),
            united_matched_deaths = sum(n_deaths, na.rm = TRUE),
            .groups = "drop")

yearly_comparison <- yearly_united |>
  full_join(yearly_iom, by = "year") |>
  left_join(yearly_united_matched, by = "year") |>
  replace_na(list(united_incidents = 0L, united_deaths = 0,
                  iom_incidents = 0L, iom_deaths = 0,
                  united_matched = 0L, united_matched_deaths = 0)) |>
  mutate(
    united_only_incidents = united_incidents - united_matched,
    match_rate_pct = ifelse(united_incidents > 0,
                            round(100 * united_matched / united_incidents, 1), 0),
    death_ratio = ifelse(iom_deaths > 0,
                         round(united_deaths / iom_deaths, 2), NA_real_)
  ) |>
  arrange(year)

cat("\n  Year  | UNITED inc/deaths | IOM inc/deaths | Matched | Match% | Death ratio\n")
cat("  ------|-------------------|----------------|---------|--------|------------\n")
for (i in seq_len(nrow(yearly_comparison))) {
  r <- yearly_comparison[i, ]
  cat(sprintf("  %d | %4d / %5s      | %4d / %5s   | %4d    | %5.1f%% | %s\n",
              r$year,
              r$united_incidents, format(r$united_deaths, big.mark = ","),
              r$iom_incidents, format(r$iom_deaths, big.mark = ","),
              r$united_matched, r$match_rate_pct,
              ifelse(is.na(r$death_ratio), "N/A",
                     sprintf("%.2f", r$death_ratio))))
}


# ====================================================================
# 6. Death count comparison for matched pairs
# ====================================================================
cat("\n--- 6. Matched pair death count comparison ---\n")

if (nrow(match_df) > 0L) {
  paired <- match_df |>
    left_join(united |> select(record_id, n_deaths, country_of_death),
              by = c("united_id" = "record_id")) |>
    left_join(iom |> select(iom_id, n_dead),
              by = "iom_id")

  cat(sprintf("  Matched pairs: %d\n", nrow(paired)))
  cat(sprintf("  Sum of UNITED deaths (matched): %s\n",
              format(sum(paired$n_deaths, na.rm = TRUE), big.mark = ",")))
  cat(sprintf("  Sum of IOM deaths (matched): %s\n",
              format(sum(paired$n_dead, na.rm = TRUE), big.mark = ",")))
  cat(sprintf("  Mean count ratio (UNITED/IOM): %.2f\n",
              mean(paired$count_ratio, na.rm = TRUE)))

  # Exact count agreement
  exact_match <- sum(paired$n_deaths == paired$n_dead, na.rm = TRUE)
  cat(sprintf("  Exact death count agreement: %d / %d (%.1f%%)\n",
              exact_match, nrow(paired),
              100 * exact_match / nrow(paired)))
}


# ====================================================================
# 7. UNITED-only records summary
# ====================================================================
cat("\n--- 7. UNITED-only records (not matched to IOM) ---\n")

united_only <- united_result |> filter(match_status == "united_only")
cat(sprintf("  Total UNITED-only: %d incidents, %s deaths\n",
            nrow(united_only),
            format(sum(united_only$n_deaths, na.rm = TRUE), big.mark = ",")))

if (nrow(united_only) > 0L) {
  cat("\n  By country of death:\n")
  uo_country <- united_only |>
    count(country_of_death, sort = TRUE)
  for (i in seq_len(min(10, nrow(uo_country)))) {
    cat(sprintf("    %-20s: %d\n", uo_country$country_of_death[i],
                uo_country$n[i]))
  }

  cat("\n  By manner of death:\n")
  uo_manner <- united_only |>
    count(manner_of_death, sort = TRUE)
  for (i in seq_len(min(8, nrow(uo_manner)))) {
    cat(sprintf("    %-30s: %d\n",
                ifelse(is.na(uo_manner$manner_of_death[i]), "(NA)",
                       uo_manner$manner_of_death[i]),
                uo_manner$n[i]))
  }

  cat("\n  By date precision:\n")
  uo_prec <- united_only |> count(incident_date_precision, sort = TRUE)
  for (i in seq_len(nrow(uo_prec))) {
    cat(sprintf("    %-12s: %d\n", uo_prec$incident_date_precision[i],
                uo_prec$n[i]))
  }
}


# ====================================================================
# 8. Save results
# ====================================================================
cat("\n--- 8. Saving results ---\n")

crosscheck_result <- list(
  matches            = match_df,
  united_classified  = united_result,
  iom_classified     = iom_result,
  yearly_comparison  = yearly_comparison,
  params             = list(date_window = DATE_WINDOW,
                            dist_max_km = DIST_MAX_KM,
                            count_ratio = COUNT_RATIO)
)

rds_path <- file.path(BASE_DIR, "data", "processed",
                       "crosscheck_united_iom.RDS")
saveRDS(crosscheck_result, rds_path)
cat(sprintf("  Saved RDS: %s\n", rds_path))

# Text summary
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))
txt_path <- tbl_path("02_clean", "10_crosscheck_united_iom_summary.txt")
sink(txt_path)
cat("UNITED vs IOM Cross-Check Summary\n")
cat("==================================\n\n")
cat(sprintf("Parameters: date_window=%d days, dist_max=%d km, count_ratio=[%.1f, %.1f]\n\n",
            DATE_WINDOW, DIST_MAX_KM, COUNT_RATIO[1], COUNT_RATIO[2]))

cat(sprintf("UNITED CMR records (2014-2025): %d incidents, %s deaths\n",
            nrow(united),
            format(sum(united$n_deaths, na.rm = TRUE), big.mark = ",")))
cat(sprintf("IOM CMR records (2014-2025):    %d incidents, %s deaths\n\n",
            nrow(iom),
            format(sum(iom$n_dead, na.rm = TRUE), big.mark = ",")))

cat(sprintf("Matched:      %d UNITED records -> %d unique IOM records\n",
            n_united_matched, n_iom_matched))
cat(sprintf("UNITED-only:  %d records (%s deaths)\n",
            n_united_only,
            format(sum(united_only$n_deaths, na.rm = TRUE), big.mark = ",")))
cat(sprintf("IOM-only:     %d records\n\n", n_iom_only))

if (nrow(match_df) > 0L) {
  cat("Match type breakdown:\n")
  print(table(match_df$match_type))
  cat("\n")
}

cat("\nYearly comparison:\n")
print(as.data.frame(yearly_comparison), row.names = FALSE)
sink()
cat(sprintf("  Saved text summary: %s\n", txt_path))

cat("\n============================================================\n")
cat("DONE — Cross-check complete\n")
cat("============================================================\n")
