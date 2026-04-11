# 057_crossings_check.R
# =====================
# MECHANISM TEST: did rough weather deter departures pre-MoU and stop
# deterring post-MoU?
#
# Hypothesis (user-formulated):
#   The flat post-MoU SWH-mortality gradient is downstream of a flat
#   post-MoU SWH-CROSSINGS gradient. Pre-MoU, smugglers waited for calm
#   weather to dispatch boats (and an NGO SAR fleet was available to
#   rescue boats that did depart). Post-MoU, with reduced SAR and an
#   active Libyan/Tunisian coast guard pushback regime, smugglers
#   dispatch regardless of weather. If true, then:
#     - Pre-MoU slope of crossings on SWH should be NEGATIVE
#     - Post-MoU slope should be near zero (or positive)
#
# This is consistent with Camarena (2024), who finds a negative wave-
# height coefficient on Italy arrivals, and with Deiana, Maheshri &
# Mastrobuoni (2024), who model weather sensitivity as a function of
# SAR intensity.
#
# Outcomes tested:
#   (1) crossing_attempts  : composite (frx_persons + lcg_tcg_pushbacks + n_dead_missing)
#   (2) frx_persons        : people in Frontex-recorded events (EU side)
#   (3) frx_incidents      : number of Frontex-recorded boats/incidents
#   (4) lcg_tcg_pushbacks  : people pushed back by Libyan/Tunisian coast guards
#                            (mostly post-MoU; pre-MoU is sparse)
#   (5) arrivals           : UNHCR Italy daily arrivals (NB: many days NA)
#
# Spec: same as 05_reduced_form_primary.R primary
#   outcome ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year
#   NegBin (fenegbin), NW(28) SEs, panel.id = ~unit + date
#
# Periods: 2014-2020 (3y pre, 3y post) and 2014-2023 (full)
#
# Output: output/tables/057_crossings_check.txt

library(tidyverse)
library(fixest)
library(lubridate)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
YEAR_START <- 2014
PERIODS <- c(2020, 2023)

cat("============================================================\n")
cat("057  WEATHER ELASTICITY OF CROSSINGS (mechanism test)\n")
cat("============================================================\n\n")

# ── 1. Load daily-agg panel ─────────────────────────────────
cat("--- 1. Loading panel ---\n")

da <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_complete.RDS")) %>%
  arrange(date) %>%
  mutate(unit = 1L)

cat(sprintf("  daily-agg: %d rows (%s to %s)\n",
    nrow(da), min(da$date), max(da$date)))

# ── 2. Helper to fit and report ─────────────────────────────
fit_and_report <- function(outcome, data, label) {
  d <- data %>% filter(!is.na(.data[[outcome]]), !is.na(swh_prevweek)) %>%
    mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)))

  if (nrow(d) < 200 || sum(d[[outcome]] > 0, na.rm = TRUE) < 30) {
    cat(sprintf("    %-30s SKIPPED (N=%d, n_nonzero=%d)\n",
        label, nrow(d), sum(d[[outcome]] > 0, na.rm = TRUE)))
    return(invisible(NULL))
  }

  f <- as.formula(paste0(outcome,
    " ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year"))

  m <- tryCatch(
    fenegbin(f, data = d, vcov = NW(28), panel.id = ~unit + date),
    error = function(e) { cat("    err:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(invisible(NULL))

  ct <- coeftable(m, vcov = NW(28))
  b1_row <- which(rownames(ct) == "swh_prevweek_z")
  b3_row <- grep(":post_mou$", rownames(ct))
  b1 <- ct[b1_row, 1]; se1 <- ct[b1_row, 2]
  b3 <- ct[b3_row, 1]; se3 <- ct[b3_row, 2]

  # Post-MoU slope = b1 + b3, with SE via delta method (assuming
  # cov(b1, b3) is small relative to vars; use full vcov for accuracy).
  V <- vcov(m, vcov = NW(28))
  c_b1 <- "swh_prevweek_z"
  c_b3 <- rownames(ct)[b3_row]
  if (all(c(c_b1, c_b3) %in% rownames(V))) {
    var_post <- V[c_b1, c_b1] + V[c_b3, c_b3] + 2 * V[c_b1, c_b3]
    se_post <- sqrt(var_post)
  } else {
    se_post <- NA_real_
  }
  b_post <- b1 + b3
  p1 <- 2 * pnorm(-abs(b1 / se1))
  p3 <- 2 * pnorm(-abs(b3 / se3))
  p_post <- 2 * pnorm(-abs(b_post / se_post))

  cat(sprintf("    %-25s  N = %5d  n_nonzero = %4d\n",
      label, nrow(d), sum(d[[outcome]] > 0)))
  cat(sprintf("      pre-MoU slope (b1)        = %+.3f (SE %.3f)  IRR = %.3f  p = %.4f\n",
      b1, se1, exp(b1), p1))
  cat(sprintf("      post-MoU shift (b3)       = %+.3f (SE %.3f)  IRR_mult = %.3f  p = %.4f\n",
      b3, se3, exp(b3), p3))
  cat(sprintf("      post-MoU slope (b1 + b3)  = %+.3f (SE %.3f)  IRR = %.3f  p = %.4f\n",
      b_post, se_post, exp(b_post), p_post))
}

# ── 3. Run by period and outcome ────────────────────────────
sink_file <- file.path(BASE_DIR, "output", "tables", "057_crossings_check.txt")
sink(sink_file)

cat("057  WEATHER ELASTICITY OF CROSSINGS (mechanism test)\n")
cat("=====================================================\n")
cat("Outcomes: composite + Frontex + pushbacks + UNHCR arrivals\n")
cat("Spec: outcome ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year\n")
cat("NegBin (fenegbin), NW(28) SEs.\n\n")
cat("KEY: pre-MoU slope NEGATIVE = rough seas deter departures\n")
cat("     post-MoU slope NEAR ZERO = departure decoupled from weather\n\n")

OUTCOMES <- list(
  list(col = "crossing_attempts",  label = "(1) crossing_attempts"),
  list(col = "frx_persons",        label = "(2) frx_persons"),
  list(col = "frx_incidents",      label = "(3) frx_incidents"),
  list(col = "lcg_tcg_pushbacks",  label = "(4) lcg_tcg_pushbacks"),
  list(col = "arrivals",           label = "(5) UNHCR arrivals")
)

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("\n=== %s ===\n", label))

  d_period <- da %>% filter(year(date) >= YEAR_START, year(date) <= ye)
  cat(sprintf("  N rows = %d\n\n", nrow(d_period)))

  for (oc in OUTCOMES) {
    fit_and_report(oc$col, d_period, oc$label)
    cat("\n")
  }
}

sink()
cat(sprintf("\nSaved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
