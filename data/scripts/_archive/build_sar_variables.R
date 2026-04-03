# build_sar_variables.R
# ====================
# Build daily SAR (Search and Rescue) operation variables.
#
# Source: Rodriguez-Sanchez et al. (2023), Table S2 — dates of NGO-led
# and government SAR operations in the Central Mediterranean.
# Original code: replication/rodriguez-sanchez/original/code/
#   dates_search_and_rescue_EU_NGOS_dates.R
#
# Output variables:
#   n_ngo_vessels:    count of active NGO SAR vessels on each day
#   n_gov_operations: count of active government operations
#   mare_nostrum:     binary indicator for Mare Nostrum period
#   coastguard_libya: binary indicator for Libyan Coast Guard operations
#
# Output: data/processed/sar_ngo_ops_daily_RS.RDS

library(dplyr)
library(lubridate)

BASE_DIR <- here::here()

cat("============================================================\n")
cat("BUILD DAILY SAR VARIABLES\n")
cat("============================================================\n\n")

# ── 1. Create daily date spine ───────────────────────────────
dates <- tibble(date = seq(as.Date("2008-01-01"), as.Date("2025-12-31"), by = "day"))

# ── 2. NGO vessel indicators ────────────────────────────────
# Adapted from Rodriguez-Sanchez et al. (2023)
# Each vessel is 1 when actively operating, 0 otherwise
# Pause/break periods set to 0

in_range <- function(d, start, end) d >= as.Date(start) & d <= as.Date(end)

sar <- dates %>% mutate(
  # --- SeaWatch vessels ---
  SAR_SeaWatch123 = case_when(
    in_range(date, "2015-06-20", "2018-07-01") ~ 1L,
    in_range(date, "2018-10-21", "2019-01-30") ~ 1L,
    in_range(date, "2019-02-23", "2019-05-17") ~ 1L,
    in_range(date, "2019-06-02", "2019-06-28") ~ 1L,
    in_range(date, "2019-12-31", "2020-02-27") ~ 1L,
    in_range(date, "2020-06-07", "2020-06-16") ~ 1L,
    in_range(date, "2021-02-25", "2021-03-26") ~ 1L,
    TRUE ~ 0L),
  SAR_SeaWatch4 = case_when(
    in_range(date, "2020-08-15", "2020-09-19") ~ 1L,
    in_range(date, "2021-03-03", "2021-10-31") ~ 1L,
    TRUE ~ 0L),

  # --- Sea-Eye vessels ---
  SAR_SeaEye = as.integer(in_range(date, "2016-04-19", "2018-06-21")),
  SAR_Seefuchs = as.integer(in_range(date, "2017-05-18", "2018-06-06")),
  SAR_AlanKurdi = case_when(
    in_range(date, "2018-12-21", "2020-05-04") ~ 1L,
    in_range(date, "2020-09-12", "2020-09-25") ~ 1L,
    TRUE ~ 0L),

  # --- Open Arms vessels ---
  SAR_Astral = as.integer(in_range(date, "2016-06-01", "2021-10-31")),
  SAR_GolfoAzzuro = as.integer(in_range(date, "2016-12-01", "2017-12-01")),
  SAR_OpenArms = case_when(
    in_range(date, "2017-12-01", "2018-03-16") ~ 1L,
    in_range(date, "2018-04-17", "2019-01-13") ~ 1L,
    in_range(date, "2019-04-24", "2019-08-20") ~ 1L,
    in_range(date, "2019-08-30", "2021-04-17") ~ 1L,
    TRUE ~ 0L),

  # --- Other NGO vessels ---
  SAR_MareLiberum = as.integer(in_range(date, "2018-08-26", "2021-10-31")),
  SAR_Mediterranea = case_when(
    in_range(date, "2018-10-03", "2019-03-18") ~ 1L,
    in_range(date, "2019-03-28", "2019-05-09") ~ 1L,
    in_range(date, "2019-07-02", "2019-07-06") ~ 1L,
    in_range(date, "2019-08-23", "2019-09-01") ~ 1L,
    in_range(date, "2020-02-05", "2020-03-18") ~ 1L,
    in_range(date, "2020-06-10", "2020-09-25") ~ 1L,
    TRUE ~ 0L),
  SAR_SMH = case_when(
    in_range(date, "2018-10-01", "2019-01-17") ~ 1L,
    in_range(date, "2019-04-18", "2020-05-07") ~ 1L,
    in_range(date, "2020-12-09", "2021-10-31") ~ 1L,
    TRUE ~ 0L),
  SAR_LouiseMichel = as.integer(in_range(date, "2020-08-22", "2020-10-22")),
  SAR_RefugeeRescue = as.integer(in_range(date, "2016-01-15", "2020-08-14")),
  SAR_MOAS = as.integer(in_range(date, "2014-08-26", "2017-09-06")),
  SAR_JugendRettet = as.integer(in_range(date, "2016-07-24", "2017-08-01")),
  SAR_Lifeline = case_when(
    in_range(date, "2017-09-14", "2018-06-26") ~ 1L,
    in_range(date, "2019-08-26", "2019-09-02") ~ 1L,
    TRUE ~ 0L),
  SAR_MSFandSOS = case_when(
    in_range(date, "2015-05-09", "2018-11-19") ~ 1L,
    in_range(date, "2019-07-21", "2021-10-31") ~ 1L,
    TRUE ~ 0L),
  SAR_MSF_BourbonArgos = case_when(
    in_range(date, "2015-05-09", "2015-08-16") ~ 1L,
    in_range(date, "2015-10-03", "2016-01-14") ~ 1L,
    in_range(date, "2016-05-06", "2016-11-20") ~ 1L,
    TRUE ~ 0L),
  SAR_MSF_Dignity1 = case_when(
    in_range(date, "2015-06-13", "2015-11-04") ~ 1L,
    in_range(date, "2016-04-22", "2016-10-04") ~ 1L,
    TRUE ~ 0L),
  SAR_MSF_VosPrudence = as.integer(in_range(date, "2017-03-20", "2017-10-05")),
  SAR_Resqship = as.integer(in_range(date, "2019-04-01", "2019-10-28")),
  SAR_Lifeboat = as.integer(in_range(date, "2016-07-01", "2017-09-22")),
  SAR_SavetheChildren = as.integer(in_range(date, "2016-09-08", "2017-10-23"))
)

# ── 3. Government/EU operation indicators ────────────────────
sar <- sar %>% mutate(
  mare_nostrum = as.integer(in_range(date, "2013-10-18", "2014-10-27")),
  frontex_triton = as.integer(in_range(date, "2014-11-01", "2018-02-01")),
  frontex_themis = as.integer(in_range(date, "2018-02-01", "2025-12-31")),
  eunavfor_sophia = as.integer(in_range(date, "2015-06-22", "2020-03-31")),
  eunavfor_irini = as.integer(in_range(date, "2020-03-31", "2025-12-31")),
  coastguard_libya = as.integer(in_range(date, "2017-02-02", "2025-12-31"))
)

# ── 4. Aggregate counts ─────────────────────────────────────
ngo_cols <- grep("^SAR_", names(sar), value = TRUE)
gov_cols <- c("mare_nostrum", "frontex_triton", "frontex_themis",
              "eunavfor_sophia", "eunavfor_irini")

sar <- sar %>% mutate(
  n_ngo_vessels = rowSums(across(all_of(ngo_cols))),
  n_gov_operations = rowSums(across(all_of(gov_cols)))
)

cat(sprintf("Daily SAR panel: %d days (%s to %s)\n",
    nrow(sar), min(sar$date), max(sar$date)))
cat(sprintf("NGO vessels tracked: %d\n", length(ngo_cols)))
cat(sprintf("Peak active NGO vessels: %d\n", max(sar$n_ngo_vessels)))

# Summary by year
cat("\nNGO vessels active (annual mean):\n")
sar %>%
  mutate(year = year(date)) %>%
  filter(year >= 2014, year <= 2024) %>%
  group_by(year) %>%
  summarise(mean_ngo = round(mean(n_ngo_vessels), 1),
            max_ngo = max(n_ngo_vessels), .groups = "drop") %>%
  print()

# ── 5. Save ──────────────────────────────────────────────────
# Keep only the summary columns + key indicators
sar_out <- sar %>%
  select(date, n_ngo_vessels, n_gov_operations,
         mare_nostrum, coastguard_libya)

saveRDS(sar_out, file.path(BASE_DIR, "data", "processed", "sar_ngo_ops_daily_RS.RDS"))
cat("\nSaved: data/processed/sar_ngo_ops_daily_RS.RDS\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
