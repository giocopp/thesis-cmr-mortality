# clean_migrant_files.R
# ====================
# Clean Migrant Files (2000-2016) for pre-IOM period (2008-2013).
#
# The Migrant Files is a researcher-compiled database of migrant deaths
# in Europe, assembled from media and NGO reports. Less systematic than
# IOM MMP but the only source for pre-2014 incident-level data.
#
# We use Migrant Files for 2008-2013 only. IOM MMP covers 2014+.
# No overlap: IOM is used for all years it covers.
#
# Cause filter (mapped to IOM equivalents):
#   "drowning or exhaustion related death"       -> IOM "Drowning"
#   "unknown - supposedly exhaustion related death" -> IOM "Mixed or unknown"
#
# Input:  data/raw/migrant files/Migrant_Files_2000-2016.xlsx
# Output: data/processed/migrant_files_cmr_pre_iom.RDS

library(readxl)
library(dplyr)
library(lubridate)

BASE_DIR <- here::here()

cat("============================================================\n")
cat("CLEAN MIGRANT FILES (pre-IOM, up to 2014-02-16)\n")
cat("============================================================\n\n")

# ── 1. Load ──────────────────────────────────────────────────
raw <- read_excel(file.path(BASE_DIR, "data", "raw", "migrant files",
                             "Migrant_Files_2000-2016.xlsx"),
                  sheet = "Events")

cat(sprintf("Raw events: %d\n", nrow(raw)))

# ── 2. Filter and clean ─────────────────────────────────────
# CMR route
cmr <- raw %>%
  filter(grepl("central med", `route (Frontex)`, ignore.case = TRUE))
cat(sprintf("CMR events: %d\n", nrow(cmr)))

# Drowning/exhaustion causes (equivalent to IOM Drowning + Mixed/unknown)
sea_causes <- c("drowning or exhaustion related death",
                "unknown - supposedly exhaustion related death")

IOM_START <- as.Date("2014-02-17")  # first IOM MMP event on CMR

mf <- cmr %>%
  filter(CartoDB_Cause_of_death %in% sea_causes) %>%
  mutate(
    date = as.Date(date, format = "%Y-%m-%dT%H:%M:%SZ"),
    dead_missing = as.numeric(dead_and_missing),
    dead_missing = if_else(is.na(dead_missing), 0, dead_missing),
    lat = as.numeric(latitude),
    lon = as.numeric(longitude),
    cause_mf = CartoDB_Cause_of_death,
    source = "Migrant Files"
  ) %>%
  filter(!is.na(date), date < IOM_START) %>%  # MF only for pre-IOM period
  select(date, dead_missing, lat, lon, cause_mf, source, Year)

cat(sprintf("CMR drowning 2008-2013: %d events, %.0f dead+missing\n",
    nrow(mf), sum(mf$dead_missing)))

# ── 3. Summary ──────────────────────────────────────────────
cat("\nBy year:\n")
for (yr in 2008:2013) {
  sub <- mf %>% filter(Year == yr)
  cat(sprintf("  %d: %3d events, %5.0f dead+missing\n",
      yr, nrow(sub), sum(sub$dead_missing)))
}

cat(sprintf("\nCoordinates: lat [%.1f, %.1f], lon [%.1f, %.1f]\n",
    min(mf$lat, na.rm = TRUE), max(mf$lat, na.rm = TRUE),
    min(mf$lon, na.rm = TRUE), max(mf$lon, na.rm = TRUE)))

# ── 4. Save ──────────────────────────────────────────────────
saveRDS(mf, file.path(BASE_DIR, "data", "processed",
                       "migrant_files_cmr_pre_iom.RDS"))
cat("\nSaved: data/processed/migrant_files_cmr_pre_iom.RDS\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
