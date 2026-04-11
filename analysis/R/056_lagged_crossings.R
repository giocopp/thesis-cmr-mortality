# 056_lagged_crossings.R
# ======================
# Enhancement #7: LAGGED DEPARTURE PROXY as a pre-determined control.
#
# The primary model deliberately avoids conditioning on crossings because
# same-day crossings are post-treatment. A LAGGED version (7-day lag) is
# pre-determined relative to today's deaths, so it can serve as a legitimate
# control for "exposure" — the number of people at sea.
#
# Control: log1p(crossing_attempts_lag7)
#   crossing_attempts = frx_persons + lcg_tcg_pushbacks + n_dead_missing
#   (from analysis/R/02_build_daily_panel.R)
#
# Sample: 2014 through **2023-06-09** (Frontex data ends there).
#         So for both "period labels" we cap at 2023 regardless.
#
# Three flavors:
#   (A) daily-agg:  direct use of daily_panel_complete.RDS
#   (B) 2-bloc:     zone panel collapsed to AFR/EU blocs with national
#                   crossings control (same across blocs on the same date)
#   (C) 4-country:  zone panel at country level with national crossings control
#                   (zone-level crossings don't exist)
#
# Output: output/tables/056_lagged_crossings.txt

library(tidyverse)
library(fixest)
library(lubridate)

BASE_DIR <- here::here()
YEAR_START <- 2014
FRX_END_YEAR <- 2023
CROSS_LAG <- 7

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("056  LAGGED DEPARTURE PROXY (crossing_attempts_lag7)\n")
cat("============================================================\n\n")

# ── 1. Load data ────────────────────────────────────────────
cat("--- 1. Loading panels ---\n")

# Drop the panel's broad n_dead_missing and replace with the analytical
# series via the shared helper. Default = incident-only, core corridor,
# all causes. Change the call to test sensitivity variants.
da <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_complete.RDS")) %>%
  select(-n_dead_missing) %>%
  left_join(build_iom_daily(), by = "date") %>%
  replace_na(list(n_dead_missing = 0)) %>%
  arrange(date) %>%
  mutate(log1p_crossings_lag7 = log1p(dplyr::lag(crossing_attempts, CROSS_LAG)))

zp <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_zone.RDS"))

# Pull national crossings_lag7 onto the zone panel (common control per date)
cross_daily <- da %>% select(date, log1p_crossings_lag7)
zp <- zp %>% left_join(cross_daily, by = "date")

# Collapse to 2 blocs, then add national crossings control
bloc <- zp %>%
  group_by(date, sar_bloc) %>%
  summarise(
    n_dead_missing       = sum(n_dead_missing),
    swh_prevweek         = mean(swh_prevweek, na.rm = TRUE),
    log1p_crossings_lag7 = first(log1p_crossings_lag7),  # national, same on both blocs
    .groups = "drop"
  ) %>%
  mutate(
    post_mou   = as.integer(date >= as.Date("2017-07-01")),
    year       = year(date),
    month_year = factor(format(date, "%Y-%m")),
    iso_week   = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))
  )
dim(bloc$date) <- NULL

cat(sprintf("  daily-agg: %d rows (%s to %s)\n",
    nrow(da), min(da$date), max(da$date)))
cat(sprintf("  2-bloc:    %d rows\n", nrow(bloc)))
cat(sprintf("  4-country: %d rows (after crossing merge)\n", nrow(zp)))
cat(sprintf("  crossings column coverage: %d non-NA in zone panel\n",
    sum(!is.na(zp$log1p_crossings_lag7))))

# ── 2. Estimation ───────────────────────────────────────────
cat("\n--- 2. Estimation ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables", "056_lagged_crossings.txt")
sink(sink_file)

cat("056  LAGGED DEPARTURE PROXY as control\n")
cat("======================================\n")
cat(sprintf("Control: log1p(crossing_attempts_lag%d) — pre-determined (lag 7 days)\n",
    CROSS_LAG))
cat("Sample: 2014 to 2023-06-09 (Frontex data end date)\n\n")
cat("Each cell compares (i) baseline model without the control to\n")
cat("(ii) model with the lagged crossings control added.\n\n")

# Data frames for analysis — single sample (2014-2023)
d_da <- da %>%
  filter(year(date) >= YEAR_START, year(date) <= FRX_END_YEAR,
         !is.na(swh_prevweek), !is.na(log1p_crossings_lag7)) %>%
  mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)),
         unit = 1L)

d_bl <- bloc %>%
  filter(year >= YEAR_START, year <= FRX_END_YEAR,
         !is.na(swh_prevweek), !is.na(log1p_crossings_lag7)) %>%
  mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)))

d_zp <- zp %>%
  filter(year >= YEAR_START, year <= FRX_END_YEAR,
         !is.na(swh_prevweek), !is.na(log1p_crossings_lag7)) %>%
  mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)))

cat(sprintf("daily-agg  N = %d\n", nrow(d_da)))
cat(sprintf("2-bloc     N = %d\n", nrow(d_bl)))
cat(sprintf("4-country  N = %d\n\n", nrow(d_zp)))

# Helper to extract the interaction coef + SE
get_b3 <- function(m, vcov_type = NW(28)) {
  ct <- coeftable(m, vcov = vcov_type)
  r <- grep(":post_mou$", rownames(ct))
  tibble(coef = ct[r, 1], se = ct[r, 2],
         p = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

# --- (A) daily-agg ---
cat("=== (A) daily-agg ===\n")

m_da_base <- fenegbin(
  n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year,
  data = d_da, vcov = NW(28), panel.id = ~unit + date
)
m_da_ctrl <- fenegbin(
  n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou +
    log1p_crossings_lag7 | month_year,
  data = d_da, vcov = NW(28), panel.id = ~unit + date
)

print(etable(m_da_base, m_da_ctrl,
              headers = c("baseline", "+ crossings_lag7"),
              vcov = NW(28), se.below = TRUE))

b3_da_base <- get_b3(m_da_base)
b3_da_ctrl <- get_b3(m_da_ctrl)
cat(sprintf("\n  baseline      : b3 = %+.3f (SE = %.3f)  p = %.4f\n",
    b3_da_base$coef, b3_da_base$se, b3_da_base$p))
cat(sprintf("  + control     : b3 = %+.3f (SE = %.3f)  p = %.4f\n",
    b3_da_ctrl$coef, b3_da_ctrl$se, b3_da_ctrl$p))

# Coefficient on the control itself
ct_ctrl <- coeftable(m_da_ctrl, vcov = NW(28))
if ("log1p_crossings_lag7" %in% rownames(ct_ctrl)) {
  cat(sprintf("  control itself: log1p_crossings_lag7 = %+.3f (SE = %.3f)  p = %.4f\n",
      ct_ctrl["log1p_crossings_lag7", 1],
      ct_ctrl["log1p_crossings_lag7", 2],
      2 * pnorm(-abs(ct_ctrl["log1p_crossings_lag7", 1] /
                      ct_ctrl["log1p_crossings_lag7", 2]))))
}

# --- (B) 2-bloc ---
cat("\n=== (B) 2-bloc ===\n")

m_bl_base <- fenegbin(
  n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou |
    month_year + sar_bloc,
  data = d_bl, vcov = NW(28), panel.id = ~sar_bloc + date
)
m_bl_ctrl <- fenegbin(
  n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou +
    log1p_crossings_lag7 | month_year + sar_bloc,
  data = d_bl, vcov = NW(28), panel.id = ~sar_bloc + date
)

print(etable(m_bl_base, m_bl_ctrl,
              headers = c("baseline", "+ crossings_lag7"),
              vcov = NW(28), se.below = TRUE))

b3_bl_base <- get_b3(m_bl_base)
b3_bl_ctrl <- get_b3(m_bl_ctrl)
cat(sprintf("\n  baseline      : b3 = %+.3f (SE = %.3f)  p = %.4f\n",
    b3_bl_base$coef, b3_bl_base$se, b3_bl_base$p))
cat(sprintf("  + control     : b3 = %+.3f (SE = %.3f)  p = %.4f\n",
    b3_bl_ctrl$coef, b3_bl_ctrl$se, b3_bl_ctrl$p))

ct_ctrl_bl <- coeftable(m_bl_ctrl, vcov = NW(28))
if ("log1p_crossings_lag7" %in% rownames(ct_ctrl_bl)) {
  cat(sprintf("  control itself: log1p_crossings_lag7 = %+.3f (SE = %.3f)  p = %.4f\n",
      ct_ctrl_bl["log1p_crossings_lag7", 1],
      ct_ctrl_bl["log1p_crossings_lag7", 2],
      2 * pnorm(-abs(ct_ctrl_bl["log1p_crossings_lag7", 1] /
                      ct_ctrl_bl["log1p_crossings_lag7", 2]))))
}

# --- (C) 4-country ---
cat("\n=== (C) 4-country ===\n")

m_zp_base <- fenegbin(
  n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou |
    month_year + country,
  data = d_zp, vcov = NW(28), panel.id = ~country + date
)
m_zp_ctrl <- fenegbin(
  n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou +
    log1p_crossings_lag7 | month_year + country,
  data = d_zp, vcov = NW(28), panel.id = ~country + date
)

print(etable(m_zp_base, m_zp_ctrl,
              headers = c("baseline", "+ crossings_lag7"),
              vcov = NW(28), se.below = TRUE))

b3_zp_base <- get_b3(m_zp_base)
b3_zp_ctrl <- get_b3(m_zp_ctrl)
cat(sprintf("\n  baseline      : b3 = %+.3f (SE = %.3f)  p = %.4f\n",
    b3_zp_base$coef, b3_zp_base$se, b3_zp_base$p))
cat(sprintf("  + control     : b3 = %+.3f (SE = %.3f)  p = %.4f\n",
    b3_zp_ctrl$coef, b3_zp_ctrl$se, b3_zp_ctrl$p))

ct_ctrl_zp <- coeftable(m_zp_ctrl, vcov = NW(28))
if ("log1p_crossings_lag7" %in% rownames(ct_ctrl_zp)) {
  cat(sprintf("  control itself: log1p_crossings_lag7 = %+.3f (SE = %.3f)  p = %.4f\n",
      ct_ctrl_zp["log1p_crossings_lag7", 1],
      ct_ctrl_zp["log1p_crossings_lag7", 2],
      2 * pnorm(-abs(ct_ctrl_zp["log1p_crossings_lag7", 1] /
                      ct_ctrl_zp["log1p_crossings_lag7", 2]))))
}

cat("\nLimitation: both zone flavors use NATIONAL crossings_lag7 common\n")
cat("across zones on the same date (zone-level crossings don't exist).\n")

sink()
cat(sprintf("\nSaved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
