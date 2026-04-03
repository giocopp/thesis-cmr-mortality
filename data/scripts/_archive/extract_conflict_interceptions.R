# extract_conflict_interceptions.R
# ================================
# Extract Libya conflict indicators (ACLED) and coast guard interceptions
# from Rodriguez-Sanchez et al. (2023) extension dataset.
#
# Source: replication/rodriguez-sanchez/extension/data/df_extended.RDS
# Original ACLED data: https://acleddata.com
# Original interception data: IOM Missing Migrants Project
#
# Output: data/processed/conflicts_interceptions_monthly_RS.RDS
#   Monthly panel with:
#   - Libya conflict events (battles, explosions, violence against civilians)
#   - Libyan Coast Guard interceptions (pushbacks)
#   - Tunisian Coast Guard interceptions (pushbacks)

library(dplyr)
library(lubridate)

BASE_DIR <- here::here()

cat("============================================================\n")
cat("EXTRACT CONFLICT & INTERCEPTION DATA\n")
cat("============================================================\n\n")

df <- readRDS(file.path(BASE_DIR, "replication", "rodriguez-sanchez",
                         "extension", "data", "df_extended.RDS"))

# Extract relevant columns
monthly <- df %>%
  transmute(
    date = as.Date(date),
    # ACLED conflict events in Libya (contemporaneous)
    battles_libya     = num_battles_Libya,
    expvio_libya      = num_expvio_Libya,     # explosions & remote violence
    violciv_libya     = num_violciv_Libya,     # violence against civilians
    protests_libya    = num_protest_Libya,
    riots_libya       = num_riots_Libya,
    strdev_libya      = num_strdev_Libya,      # strategic developments
    # Total conflict intensity (battles + explosions + violence)
    conflict_libya    = num_battles_Libya + num_expvio_Libya + num_violciv_Libya,
    # ACLED conflict events in Tunisia
    battles_tunisia   = num_battles_Tunisia,
    expvio_tunisia    = num_expvio_Tunisia,
    violciv_tunisia   = num_violciv_Tunisia,
    protests_tunisia  = num_protest_Tunisia,
    riots_tunisia     = num_riots_Tunisia,
    conflict_tunisia  = num_battles_Tunisia + num_expvio_Tunisia + num_violciv_Tunisia,
    # Coast guard interceptions
    lcg_interceptions = LCG_pushbacks_count,   # Libyan Coast Guard
    tcg_interceptions = TCG_pushbacks_count    # Tunisian Coast Guard
  )

cat(sprintf("Monthly panel: %d months (%s to %s)\n",
    nrow(monthly), min(monthly$date), max(monthly$date)))

cat("\nConflict indicators (Libya, ACLED):\n")
cat(sprintf("  Battles:     mean=%.1f  max=%.0f\n",
    mean(monthly$battles_libya, na.rm=TRUE), max(monthly$battles_libya, na.rm=TRUE)))
cat(sprintf("  Explosions:  mean=%.1f  max=%.0f\n",
    mean(monthly$expvio_libya, na.rm=TRUE), max(monthly$expvio_libya, na.rm=TRUE)))
cat(sprintf("  Viol. civ.:  mean=%.1f  max=%.0f\n",
    mean(monthly$violciv_libya, na.rm=TRUE), max(monthly$violciv_libya, na.rm=TRUE)))
cat(sprintf("  Total conflict: mean=%.1f  max=%.0f\n",
    mean(monthly$conflict_libya, na.rm=TRUE), max(monthly$conflict_libya, na.rm=TRUE)))

cat("\nInterceptions:\n")
cat(sprintf("  LCG (Libya):  mean=%.0f  max=%.0f\n",
    mean(monthly$lcg_interceptions), max(monthly$lcg_interceptions)))
cat(sprintf("  TCG (Tunisia): mean=%.0f  max=%.0f\n",
    mean(monthly$tcg_interceptions), max(monthly$tcg_interceptions)))

# Coverage check
cat("\nNon-NA coverage:\n")
cat(sprintf("  Conflict: %d / %d months\n",
    sum(!is.na(monthly$conflict_libya)), nrow(monthly)))
cat(sprintf("  LCG interceptions: %d / %d months\n",
    sum(!is.na(monthly$lcg_interceptions)), nrow(monthly)))

saveRDS(monthly, file.path(BASE_DIR, "data", "processed",
                            "conflicts_interceptions_monthly_RS.RDS"))
cat("\nSaved: data/processed/conflicts_interceptions_monthly_RS.RDS\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
