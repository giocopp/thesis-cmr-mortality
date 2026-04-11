# diag_swh_weighted_compare.R
# ===========================
# Run the primary reduced-form spec and the 051 zone-panel spec with
# BOTH the polygon-mean SWH (swh_prevweek) and the death-weighted SWH
# (swh_w_prevweek), and print the headline coefficients side-by-side.
#
# If swh_w_prevweek gives a materially different (tighter / larger) β₃,
# the weighted measure is worth adopting. Otherwise the 0.95 correlation
# between the two means they're telling the same story.
#
# Input:  analysis/data/daily_panel_complete.RDS  (with swh_w columns)
#         analysis/data/daily_panel_zone.RDS      (with swh_w columns)
#         data/processed/iom_mmp_incidents.RDS    (via build_iom_daily)
# Output: output/tables/diag_swh_weighted_compare.txt

library(tidyverse)
library(fixest)
library(lubridate)

BASE_DIR      <- here::here()
START_DATE    <- as.Date("2014-01-01")
SYMMETRIC_END <- as.Date("2020-12-31")

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("diag_swh_weighted_compare: polygon-mean vs death-weighted SWH\n")
cat("============================================================\n\n")

# ── 1. Daily-agg panel (primary sample) ───────────────────
cat("--- 1. Loading daily-agg panel ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                             "daily_panel_complete.RDS")) %>%
  select(-n_dead_missing) %>%
  arrange(date) %>%
  mutate(year = year(date), unit = 1L)

PANEL_END <- max(panel$date)
cat(sprintf("  panel: %d days, %s to %s\n",
            nrow(panel), min(panel$date), PANEL_END))

# Build primary death series via the shared helper (same as 05)
daily_primary <- build_iom_daily()
panel <- panel %>%
  left_join(daily_primary, by = "date") %>%
  replace_na(list(n_dead_missing = 0))

cat(sprintf("  total deaths in primary sample: %.0f\n",
            sum(panel$n_dead_missing)))

# ── 2. Zone panel ──────────────────────────────────────────
cat("\n--- 2. Loading zone panel ---\n")

zp <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_zone.RDS"))
dim(zp$date) <- NULL  # fixest bug: NW + panel.id chokes on dim attr
cat(sprintf("  zone panel: %d rows (%d days x 4 zones)\n",
            nrow(zp), length(unique(zp$date))))

# Pre-collapse to 2-bloc (matches 051_zone_panel.R pattern)
bloc <- zp %>%
  group_by(date, sar_bloc) %>%
  summarise(
    n_dead_missing  = sum(n_dead_missing),
    swh             = mean(swh, na.rm = TRUE),
    swh_prev3days   = mean(swh_prev3days, na.rm = TRUE),
    swh_prevweek    = mean(swh_prevweek, na.rm = TRUE),
    swh_w           = mean(swh_w, na.rm = TRUE),
    swh_w_prev3days = mean(swh_w_prev3days, na.rm = TRUE),
    swh_w_prevweek  = mean(swh_w_prevweek, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    post_mou   = as.integer(date >= as.Date("2017-07-01")),
    year       = year(date),
    month_year = factor(format(date, "%Y-%m"))
  )
dim(bloc$date) <- NULL
cat(sprintf("  2-bloc panel: %d rows\n", nrow(bloc)))

# ── 3. Run each spec with both SWH variables ──────────────
cat("\n--- 3. Running primary + 051 with both SWH variables ---\n")

fit_period_primary <- function(swh_var, end_date) {
  d <- panel %>%
    filter(between(date, START_DATE, end_date)) %>%
    mutate(swh_z       = as.numeric(scale(.data[[swh_var]])),
           month_year  = factor(format(date, "%Y-%m")))

  fenegbin(n_dead_missing ~ swh_z + swh_z:post_mou | month_year,
           data = d, vcov = NW(28), panel.id = ~unit + date)
}

fit_period_zone <- function(swh_var, end_year, fe_spec) {
  source_tbl <- if (fe_spec == "2bloc") bloc else zp
  d <- source_tbl %>%
    filter(year >= 2014, year <= end_year,
           !is.na(.data[[swh_var]])) %>%
    mutate(swh_z = as.numeric(scale(.data[[swh_var]]))) %>%
    as.data.frame()
  dim(d$date) <- NULL

  if (fe_spec == "2bloc") {
    f <- as.formula("n_dead_missing ~ swh_z + swh_z:post_mou | month_year + sar_bloc")
    p_id <- ~sar_bloc + date
  } else {
    f <- as.formula("n_dead_missing ~ swh_z + swh_z:post_mou | month_year + country")
    p_id <- ~country + date
  }
  fenegbin(f, data = d, vcov = NW(28), panel.id = p_id)
}

extract_row <- function(m, label, vcov_type = NW(28)) {
  ct <- coeftable(m, vcov = vcov_type)
  b1_row <- which(rownames(ct) == "swh_z")
  b3_row <- grep(":post_mou$", rownames(ct))
  b1 <- ct[b1_row, 1]; se1 <- ct[b1_row, 2]
  b3 <- ct[b3_row, 1]; se3 <- ct[b3_row, 2]
  V <- vcov(m, vcov = vcov_type)
  var_post <- V["swh_z", "swh_z"] +
              V[rownames(ct)[b3_row], rownames(ct)[b3_row]] +
              2 * V["swh_z", rownames(ct)[b3_row]]
  b_post <- b1 + b3
  se_post <- sqrt(var_post)

  tibble(
    spec    = label,
    b1      = b1,      se1      = se1,
    b3      = b3,      se3      = se3,
    b_post  = b_post,  se_post  = se_post,
    p_b3    = 2 * pnorm(-abs(b3 / se3))
  )
}

cat("  Fitting primary (daily-agg) ...\n")
m_primary_mean_full <- fit_period_primary("swh_prevweek",    PANEL_END)
m_primary_w_full    <- fit_period_primary("swh_w_prevweek",  PANEL_END)
m_primary_mean_sym  <- fit_period_primary("swh_prevweek",    SYMMETRIC_END)
m_primary_w_sym     <- fit_period_primary("swh_w_prevweek",  SYMMETRIC_END)

cat("  Fitting zone 4-country and 2-bloc ...\n")
m_zone4_mean_full <- fit_period_zone("swh_prevweek",    2023, "4country")
m_zone4_w_full    <- fit_period_zone("swh_w_prevweek",  2023, "4country")
m_zone2_mean_full <- fit_period_zone("swh_prevweek",    2023, "2bloc")
m_zone2_w_full    <- fit_period_zone("swh_w_prevweek",  2023, "2bloc")
m_zone4_mean_sym  <- fit_period_zone("swh_prevweek",    2020, "4country")
m_zone4_w_sym     <- fit_period_zone("swh_w_prevweek",  2020, "4country")
m_zone2_mean_sym  <- fit_period_zone("swh_prevweek",    2020, "2bloc")
m_zone2_w_sym     <- fit_period_zone("swh_w_prevweek",  2020, "2bloc")

results <- bind_rows(
  extract_row(m_primary_mean_full, "Primary 2014..full   | mean"),
  extract_row(m_primary_w_full,    "Primary 2014..full   | weighted"),
  extract_row(m_primary_mean_sym,  "Primary 2014-2020    | mean"),
  extract_row(m_primary_w_sym,     "Primary 2014-2020    | weighted"),
  extract_row(m_zone4_mean_full,   "Zone 4-country 2014-2023 | mean"),
  extract_row(m_zone4_w_full,      "Zone 4-country 2014-2023 | weighted"),
  extract_row(m_zone2_mean_full,   "Zone 2-bloc 2014-2023    | mean"),
  extract_row(m_zone2_w_full,      "Zone 2-bloc 2014-2023    | weighted"),
  extract_row(m_zone4_mean_sym,    "Zone 4-country 2014-2020 | mean"),
  extract_row(m_zone4_w_sym,       "Zone 4-country 2014-2020 | weighted"),
  extract_row(m_zone2_mean_sym,    "Zone 2-bloc 2014-2020    | mean"),
  extract_row(m_zone2_w_sym,       "Zone 2-bloc 2014-2020    | weighted")
)

# ── 4. Report ─────────────────────────────────────────────
sink_file <- file.path(BASE_DIR, "output", "tables",
                        "diag_swh_weighted_compare.txt")
sink(sink_file)

cat("POLYGON-MEAN vs DEATH-WEIGHTED SWH — head-to-head\n")
cat("==================================================\n\n")
cat("Compared variables:\n")
cat("  swh_prevweek     (polygon-mean, current primary)\n")
cat("  swh_w_prevweek   (death-weighted, static weights from historical incidents)\n\n")
cat(sprintf("Cor(swh, swh_w) in the panel: %.4f\n",
    cor(panel$swh, panel$swh_w, use = "pairwise.complete.obs")))
cat(sprintf("Cor(swh_prevweek, swh_w_prevweek): %.4f\n\n",
    cor(panel$swh_prevweek, panel$swh_w_prevweek,
        use = "pairwise.complete.obs")))

cat("All specs use NegBin (fenegbin), NW(28) SEs.\n")
cat("b1 = pre-MoU slope; b3 = interaction; b_post = b1 + b3 = post-MoU slope.\n\n")

print_row <- function(r) {
  cat(sprintf(
    "  %-40s  b1=%+.3f (%.3f)  b3=%+.3f (%.3f, p=%.4f)  post=%+.3f\n",
    r$spec, r$b1, r$se1, r$b3, r$se3, r$p_b3, r$b_post))
}
for (i in seq_len(nrow(results))) print_row(results[i, ])

cat("\n--- PAIRED DIFFS (weighted minus mean) ---\n")
paired <- results %>%
  mutate(base = sub(" \\| .*$", "", spec),
         kind = sub("^.* \\| ", "", spec)) %>%
  select(base, kind, b3, se3) %>%
  pivot_wider(names_from = kind, values_from = c(b3, se3))

for (i in seq_len(nrow(paired))) {
  r <- paired[i, ]
  d_b3 <- r$b3_weighted - r$b3_mean
  d_se <- r$se3_weighted - r$se3_mean
  cat(sprintf("  %-30s  Δb3 = %+.3f   Δse = %+.3f\n", r$base, d_b3, d_se))
}

sink()
cat(sprintf("\nSaved: %s\n", sink_file))
cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
