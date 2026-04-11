# 054_wild_bootstrap.R
# ====================
# Enhancement #5: WILD CLUSTER BOOTSTRAP at the ISO-week level.
#
# Robustness check for the SE of the primary interaction coefficient
# (swh_prevweek_z x post_mou), complementing NW(28) and cluster(iso_week).
#
# fwildclusterboot::boottest() ONLY supports OLS via feols(). It does not
# support fenegbin or fepois. Two parallel bootstrap approaches:
#
#   (1) Pairs cluster bootstrap on fepois (primary):
#       Sample iso_week clusters with replacement, refit fepois B times,
#       compute bootstrap SE and percentile p-value. Matches the count
#       structure of the NegBin headline; standard approach since
#       Cameron-Gelbach-Miller (2008).
#
#   (2) Wild cluster bootstrap on feols with log1p(deaths) (secondary):
#       Rigorous wild-cluster bootstrap via fwildclusterboot::boottest
#       on a linear approximation. Different functional form but the
#       "gold standard" inference method for clustered errors.
#
# Flavors:
#   (A) daily-agg:   FE = month_year,              panel.id = unit+date
#   (B) 2-bloc:      FE = month_year + sar_bloc,   panel.id = sar_bloc+date
#   (C) 4-country:   FE = month_year + country,    panel.id = country+date
#
# Output: output/tables/054_wild_bootstrap.txt

library(tidyverse)
library(fixest)
library(fwildclusterboot)
library(lubridate)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
YEAR_START <- 2014
PERIODS <- c(2020, 2023)
B_BOOT <- 999

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

set.seed(42)

cat("============================================================\n")
cat("054  WILD / PAIRS CLUSTER BOOTSTRAP (iso_week)\n")
cat("============================================================\n\n")

# ── 1. Load panels ──────────────────────────────────────────
cat("--- 1. Loading panels ---\n")

# Drop the panel's broad n_dead_missing and replace with the analytical
# series via the shared helper. Default = incident-only, core corridor,
# all causes. Change the call to test sensitivity variants.
da <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_complete.RDS")) %>%
  select(-n_dead_missing) %>%
  left_join(build_iom_daily(), by = "date") %>%
  replace_na(list(n_dead_missing = 0)) %>%
  arrange(date)

zp <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_zone.RDS"))

# Collapse to 2 blocs (AFR/EU)
bloc <- zp %>%
  group_by(date, sar_bloc) %>%
  summarise(
    n_dead_missing = sum(n_dead_missing),
    swh_prevweek   = mean(swh_prevweek, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    post_mou   = as.integer(date >= MOU_DATE),
    year       = year(date),
    month_year = factor(format(date, "%Y-%m")),
    iso_week   = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))
  )
dim(bloc$date) <- NULL

cat(sprintf("  daily-agg:  %d rows\n", nrow(da)))
cat(sprintf("  2-bloc:     %d rows\n", nrow(bloc)))
cat(sprintf("  4-country:  %d rows\n", nrow(zp)))

# ── 2. Bootstrap function: pairs cluster on fepois ──────────
# Sample iso_week clusters with replacement, refit fepois, store coefficient.
pairs_cluster_boot_fepois <- function(formula, data, panel_id_formula,
                                       cluster_col, focal_coef,
                                       B = B_BOOT) {
  data[[cluster_col]] <- as.character(data[[cluster_col]])
  clusters <- unique(data[[cluster_col]])
  n_clust <- length(clusters)

  # Original fit
  m0 <- fepois(formula, data = data, vcov = NW(28),
                panel.id = panel_id_formula)
  beta_hat <- coef(m0)[[focal_coef]]

  # Pre-split by cluster for speed
  idx_by_clust <- split(seq_len(nrow(data)), data[[cluster_col]])

  boot_betas <- numeric(B)
  ok <- 0L
  for (b in seq_len(B)) {
    samp_clusters <- sample(clusters, n_clust, replace = TRUE)
    rows <- unlist(idx_by_clust[samp_clusters], use.names = FALSE)
    d_b <- data[rows, , drop = FALSE]

    m_b <- tryCatch(
      suppressMessages(fepois(formula, data = d_b)),
      error = function(e) NULL
    )
    if (!is.null(m_b) && focal_coef %in% names(coef(m_b))) {
      boot_betas[b] <- coef(m_b)[[focal_coef]]
      ok <- ok + 1L
    } else {
      boot_betas[b] <- NA_real_
    }
  }
  valid <- !is.na(boot_betas)
  boot_betas <- boot_betas[valid]
  if (length(boot_betas) < 100)
    warning(sprintf("Only %d valid bootstrap reps", length(boot_betas)))

  beta_sd <- sd(boot_betas)
  ci_lo <- quantile(boot_betas, 0.025)
  ci_hi <- quantile(boot_betas, 0.975)

  # Studentized / Percentile-t p-value
  t_obs <- beta_hat / beta_sd
  t_boot <- (boot_betas - beta_hat) / beta_sd
  p_pairs <- (1 + sum(abs(t_boot) >= abs(t_obs))) / (length(boot_betas) + 1)

  list(
    beta_hat = beta_hat,
    boot_se  = beta_sd,
    ci_lo    = ci_lo,
    ci_hi    = ci_hi,
    p_pairs  = p_pairs,
    n_boot   = length(boot_betas),
    n_ok     = ok
  )
}

# ── 3. Run bootstraps ───────────────────────────────────────
cat("\n--- 3. Estimation and bootstrap ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables", "054_wild_bootstrap.txt")
sink(sink_file)

cat("054  WILD / PAIRS CLUSTER BOOTSTRAP AT iso_week\n")
cat("================================================\n")
cat("Target: swh_prevweek_z:post_mou\n\n")
cat("Approaches:\n")
cat("  (1) fepois pairs cluster bootstrap on iso_week (B =",
    B_BOOT, "reps)\n")
cat("  (2) feols log1p(deaths) wild cluster bootstrap via boottest\n")
cat("\nfwildclusterboot::boottest() only supports OLS. The fepois pairs\n")
cat("bootstrap is the primary count-model robustness; the feols wild\n")
cat("bootstrap is a linear-approximation robustness using the 'proper' wild\n")
cat("cluster procedure.\n\n")
cat("Flavors:\n")
cat("  [A] daily-agg   FE=month_year\n")
cat("  [B] 2-bloc      FE=month_year + sar_bloc\n")
cat("  [C] 4-country   FE=month_year + country\n\n")

results <- list()

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("\n=== %s ===\n", label))

  # --- prepare data ---
  d_da <- da %>%
    filter(year(date) >= YEAR_START, year(date) <= ye,
           !is.na(swh_prevweek)) %>%
    mutate(
      swh_prevweek_z = as.numeric(scale(swh_prevweek)),
      unit = 1L,
      log1p_deaths = log1p(n_dead_missing)
    ) %>%
    as.data.frame()

  d_bl <- bloc %>%
    filter(year >= YEAR_START, year <= ye, !is.na(swh_prevweek)) %>%
    mutate(
      swh_prevweek_z = as.numeric(scale(swh_prevweek)),
      log1p_deaths = log1p(n_dead_missing)
    ) %>%
    as.data.frame()

  d_zp <- zp %>%
    filter(year >= YEAR_START, year <= ye, !is.na(swh_prevweek)) %>%
    mutate(
      swh_prevweek_z = as.numeric(scale(swh_prevweek)),
      log1p_deaths = log1p(n_dead_missing)
    ) %>%
    as.data.frame()

  cat(sprintf("  daily-agg  N = %d | iso_weeks = %d\n",
      nrow(d_da), length(unique(d_da$iso_week))))
  cat(sprintf("  2-bloc     N = %d | iso_weeks = %d\n",
      nrow(d_bl), length(unique(d_bl$iso_week))))
  cat(sprintf("  4-country  N = %d | iso_weeks = %d\n",
      nrow(d_zp), length(unique(d_zp$iso_week))))

  # ---------------- (A) daily-agg ----------------
  cat("\n  [A] daily-agg\n")

  t0 <- Sys.time()
  res_da_pairs <- pairs_cluster_boot_fepois(
    formula = n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou |
      month_year,
    data = d_da,
    panel_id_formula = ~unit + date,
    cluster_col = "iso_week",
    focal_coef = "swh_prevweek_z:post_mou"
  )
  cat(sprintf("    pairs cluster boot: %d valid reps in %.1f s\n",
      res_da_pairs$n_boot,
      as.numeric(Sys.time() - t0, units = "secs")))
  cat(sprintf("      b3 = %+.3f, boot SE = %.3f, p = %.4f, 95%% CI [%.3f, %.3f]\n",
      res_da_pairs$beta_hat, res_da_pairs$boot_se, res_da_pairs$p_pairs,
      res_da_pairs$ci_lo, res_da_pairs$ci_hi))

  # feols wild cluster bootstrap (log1p)
  m_ols_da <- feols(log1p_deaths ~ swh_prevweek_z + swh_prevweek_z:post_mou |
                      month_year, data = d_da)
  bt_da <- tryCatch(
    boottest(m_ols_da,
              clustid = "iso_week",
              param = "swh_prevweek_z:post_mou",
              B = 1999, type = "rademacher"),
    error = function(e) { cat("    boottest err:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(bt_da)) {
    cat(sprintf("    feols log1p wild cluster boot: b = %+.3f, p = %.4f, 95%% CI [%.3f, %.3f]\n",
        bt_da$point_estimate, bt_da$p_val,
        bt_da$conf_int[1], bt_da$conf_int[2]))
  }

  # ---------------- (B) 2-bloc ----------------
  cat("\n  [B] 2-bloc\n")

  t0 <- Sys.time()
  res_bl_pairs <- pairs_cluster_boot_fepois(
    formula = n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou |
      month_year + sar_bloc,
    data = d_bl,
    panel_id_formula = ~sar_bloc + date,
    cluster_col = "iso_week",
    focal_coef = "swh_prevweek_z:post_mou"
  )
  cat(sprintf("    pairs cluster boot: %d valid reps in %.1f s\n",
      res_bl_pairs$n_boot,
      as.numeric(Sys.time() - t0, units = "secs")))
  cat(sprintf("      b3 = %+.3f, boot SE = %.3f, p = %.4f, 95%% CI [%.3f, %.3f]\n",
      res_bl_pairs$beta_hat, res_bl_pairs$boot_se, res_bl_pairs$p_pairs,
      res_bl_pairs$ci_lo, res_bl_pairs$ci_hi))

  m_ols_bl <- feols(log1p_deaths ~ swh_prevweek_z + swh_prevweek_z:post_mou |
                      month_year + sar_bloc, data = d_bl)
  bt_bl <- tryCatch(
    boottest(m_ols_bl,
              clustid = "iso_week",
              param = "swh_prevweek_z:post_mou",
              B = 1999, type = "rademacher"),
    error = function(e) { cat("    boottest err:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(bt_bl)) {
    cat(sprintf("    feols log1p wild cluster boot: b = %+.3f, p = %.4f, 95%% CI [%.3f, %.3f]\n",
        bt_bl$point_estimate, bt_bl$p_val,
        bt_bl$conf_int[1], bt_bl$conf_int[2]))
  }

  # ---------------- (C) 4-country ----------------
  cat("\n  [C] 4-country\n")

  t0 <- Sys.time()
  res_zp_pairs <- pairs_cluster_boot_fepois(
    formula = n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou |
      month_year + country,
    data = d_zp,
    panel_id_formula = ~country + date,
    cluster_col = "iso_week",
    focal_coef = "swh_prevweek_z:post_mou"
  )
  cat(sprintf("    pairs cluster boot: %d valid reps in %.1f s\n",
      res_zp_pairs$n_boot,
      as.numeric(Sys.time() - t0, units = "secs")))
  cat(sprintf("      b3 = %+.3f, boot SE = %.3f, p = %.4f, 95%% CI [%.3f, %.3f]\n",
      res_zp_pairs$beta_hat, res_zp_pairs$boot_se, res_zp_pairs$p_pairs,
      res_zp_pairs$ci_lo, res_zp_pairs$ci_hi))

  m_ols_zp <- feols(log1p_deaths ~ swh_prevweek_z + swh_prevweek_z:post_mou |
                      month_year + country, data = d_zp)
  bt_zp <- tryCatch(
    boottest(m_ols_zp,
              clustid = "iso_week",
              param = "swh_prevweek_z:post_mou",
              B = 1999, type = "rademacher"),
    error = function(e) { cat("    boottest err:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(bt_zp)) {
    cat(sprintf("    feols log1p wild cluster boot: b = %+.3f, p = %.4f, 95%% CI [%.3f, %.3f]\n",
        bt_zp$point_estimate, bt_zp$p_val,
        bt_zp$conf_int[1], bt_zp$conf_int[2]))
  }

  results[[label]] <- list(
    da = res_da_pairs,
    bl = res_bl_pairs,
    zp = res_zp_pairs,
    bt_da = bt_da,
    bt_bl = bt_bl,
    bt_zp = bt_zp
  )
}

sink()
cat(sprintf("\nSaved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
