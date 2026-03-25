# 03b_placebo_stability.R
# =======================
# Placebo treatment dates and pre-period gradient stability tests.
# These test the gradient stability assumption (Assumption 3 in estimands.qmd).
#
# Runs all tests across three weather windows: day0, 3-day mean, 7-day mean.
#
# Test 1: Placebo treatment dates — re-estimate Weather x Post interaction at
#         many candidate dates. If the MoU date is special, beta_3 should
#         stand out.
#
# Test 2: Pre-period gradient by year — estimate the SWH slope separately
#         for each pre-MoU year.  Stability = the gradient was constant
#         before the intervention.
#
# Input:  data/processed/cmr_events_with_weather.RDS
# Output: output/figures/placebo_beta3.pdf
#         output/figures/placebo_beta3.png
#         output/tables/placebo_results.csv
#         printed diagnostics

library(fixest)
library(data.table)
library(ggplot2)

BASE_DIR <- here::here()
d <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                      "cmr_events_with_weather.RDS")))
d[, grid_1deg := paste0(sprintf("%.0f", round(grid_lat)), "_",
                          sprintf("%.0f", round(grid_lon)))]

MOU_DATE <- as.Date("2017-02-01")

# Candidate dates: quarterly from 2015-Q1 to 2018-Q3
placebo_dates <- as.Date(c(
  "2015-01-01", "2015-04-01", "2015-07-01", "2015-10-01",
  "2016-01-01", "2016-04-01", "2016-07-01", "2016-10-01",
  "2017-02-01",  # actual MoU
  "2017-07-01", "2017-10-01",
  "2018-01-01", "2018-04-01", "2018-07-01"
))

# Weather windows to test
weather_windows <- list(
  day0 = list(
    label = "Day-0",
    swh = "swh_day0", wind = "wind_day0", gust = "i10fg_day0"
  ),
  d3 = list(
    label = "3-day mean",
    swh = "swh_mean_3d", wind = "wind_mean_3d", gust = "i10fg_mean_3d"
  ),
  d7 = list(
    label = "7-day mean",
    swh = "swh_mean_7d", wind = "wind_mean_7d", gust = "i10fg_mean_7d"
  )
)

# ============================================================
# Test 1: Placebo treatment dates (SWH x Post) — all windows
# ============================================================
cat("============================================================\n")
cat("TEST 1: PLACEBO TREATMENT DATES — SWH x Post\n")
cat("============================================================\n\n")

all_placebo <- list()

for (wname in names(weather_windows)) {
  w <- weather_windows[[wname]]
  cat(sprintf("--- Window: %s (var: %s) ---\n", w$label, w$swh))

  results_swh <- rbindlist(lapply(placebo_dates, function(pd) {
    d[, post_placebo := as.integer(date >= pd)]

    swh_var <- w$swh
    wind_var <- w$wind
    n_pre  <- sum(!is.na(d[[swh_var]]) & d$post_placebo == 0)
    n_post <- sum(!is.na(d[[swh_var]]) & d$post_placebo == 1)
    if (n_pre < 30 || n_post < 30) return(NULL)

    fml <- as.formula(sprintf(
      "dead_missing ~ %s * post_placebo + %s | grid_1deg + year_fac + month_fac",
      swh_var, wind_var))

    m <- tryCatch(fenegbin(fml, data = d, vcov = "hetero"), error = function(e) NULL)
    if (is.null(m)) return(NULL)

    ct <- summary(m, vcov = "hetero")$coeftable
    int_pat <- paste0(swh_var, ":post_placebo")
    int_row <- grep(int_pat, rownames(ct), fixed = TRUE)
    if (length(int_row) == 0) return(NULL)

    data.table(
      date   = pd, is_mou = pd == MOU_DATE,
      beta   = ct[int_row, 1], se = ct[int_row, 2], p = ct[int_row, 4],
      n_pre  = n_pre, n_post = n_post
    )
  }))

  if (nrow(results_swh) > 0) {
    results_swh[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]
    results_swh[, `:=`(variable = "SWH", window = w$label)]
    cat(sprintf("  MoU estimate: beta=%+.4f, SE=%.4f, p=%.4f\n",
        results_swh[is_mou == TRUE, beta],
        results_swh[is_mou == TRUE, se],
        results_swh[is_mou == TRUE, p]))
    all_placebo[[length(all_placebo) + 1]] <- results_swh
  }
  cat("\n")
}

# ============================================================
# Test 2: Placebo treatment dates (Gust x Post) — all windows
# ============================================================
cat("============================================================\n")
cat("TEST 2: PLACEBO TREATMENT DATES — Gust x Post\n")
cat("============================================================\n\n")

for (wname in names(weather_windows)) {
  w <- weather_windows[[wname]]
  cat(sprintf("--- Window: %s (var: %s) ---\n", w$label, w$gust))

  results_gust <- rbindlist(lapply(placebo_dates, function(pd) {
    d[, post_placebo := as.integer(date >= pd)]

    gust_var <- w$gust
    swh_var  <- w$swh
    wind_var <- w$wind
    n_pre  <- sum(!is.na(d[[gust_var]]) & d$post_placebo == 0)
    n_post <- sum(!is.na(d[[gust_var]]) & d$post_placebo == 1)
    if (n_pre < 30 || n_post < 30) return(NULL)

    fml <- as.formula(sprintf(
      "dead_missing ~ %s * post_placebo + %s + %s | grid_1deg + year_fac + month_fac",
      gust_var, swh_var, wind_var))

    m <- tryCatch(fenegbin(fml, data = d, vcov = "hetero"), error = function(e) NULL)
    if (is.null(m)) return(NULL)

    ct <- summary(m, vcov = "hetero")$coeftable
    int_pat <- paste0(gust_var, ":post_placebo")
    int_row <- grep(int_pat, rownames(ct), fixed = TRUE)
    if (length(int_row) == 0) return(NULL)

    data.table(
      date   = pd, is_mou = pd == MOU_DATE,
      beta   = ct[int_row, 1], se = ct[int_row, 2], p = ct[int_row, 4],
      n_pre  = n_pre, n_post = n_post
    )
  }))

  if (nrow(results_gust) > 0) {
    results_gust[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]
    results_gust[, `:=`(variable = "Gust", window = w$label)]
    cat(sprintf("  MoU estimate: beta=%+.4f, SE=%.4f, p=%.4f\n",
        results_gust[is_mou == TRUE, beta],
        results_gust[is_mou == TRUE, se],
        results_gust[is_mou == TRUE, p]))
    all_placebo[[length(all_placebo) + 1]] <- results_gust
  }
  cat("\n")
}

# ============================================================
# Test 3: Pre-period gradient stability — all windows
# ============================================================
cat("============================================================\n")
cat("TEST 3: PRE-PERIOD GRADIENT STABILITY\n")
cat("============================================================\n\n")

for (wname in names(weather_windows)) {
  w <- weather_windows[[wname]]
  swh_var  <- w$swh
  wind_var <- w$wind

  d_pre <- d[post_mou == 0 & !is.na(d[[swh_var]])]
  cat(sprintf("--- Window: %s --- (n_pre = %d)\n", w$label, nrow(d_pre)))

  if (nrow(d_pre) < 30) { cat("  Too few obs\n\n"); next }

  # Wald test: SWH × year interactions jointly zero?
  fml_pre <- as.formula(sprintf(
    "dead_missing ~ %s * factor(year) + %s | month_fac",
    swh_var, wind_var))

  m_pre <- tryCatch(
    fenegbin(fml_pre, data = d_pre, vcov = "hetero"),
    error = function(e) { cat("  fenegbin failed:", e$message, "\n"); NULL }
  )

  if (!is.null(m_pre)) {
    ct_pre <- summary(m_pre, vcov = "hetero")$coeftable
    swh_rows <- grep(swh_var, rownames(ct_pre), fixed = TRUE)
    cat("  SWH gradient by year:\n")
    print(round(ct_pre[swh_rows, , drop = FALSE], 4))

    int_pat <- paste0(swh_var, ":factor")
    int_rows <- grep(int_pat, rownames(ct_pre))
    if (length(int_rows) >= 2) {
      w_test <- tryCatch(
        wald(m_pre, keep = int_pat, vcov = "hetero"),
        error = function(e) NULL
      )
      if (!is.null(w_test)) {
        cat("\n  Wald test (H0: gradient constant across pre-MoU years):\n")
        print(w_test)
      }
    }
  }

  # Year-by-year slopes
  cat(sprintf("\n  Year-by-year %s slopes (separate Poisson, no grid FE):\n", swh_var))
  for (yr in sort(unique(d_pre$year))) {
    dsub <- d_pre[year == yr]
    if (nrow(dsub) < 10) {
      cat(sprintf("    %d: n=%d (too few)\n", yr, nrow(dsub)))
      next
    }
    fml_yr <- as.formula(sprintf("dead_missing ~ %s + %s", swh_var, wind_var))
    m_yr <- tryCatch(glm(fml_yr, data = dsub, family = poisson), error = function(e) NULL)
    if (!is.null(m_yr)) {
      b  <- coef(m_yr)[swh_var]
      se <- summary(m_yr)$coefficients[swh_var, 2]
      cat(sprintf("    %d: n=%d, beta=%+.4f (SE=%.4f), exp=%.4f\n",
          yr, nrow(dsub), b, se, exp(b)))
    }
  }
  cat("\n")
}

# ============================================================
# Plot: placebo coefficient paths — faceted by window × variable
# ============================================================
cat("============================================================\n")
cat("GENERATING PLOTS\n")
cat("============================================================\n\n")

results_all <- rbindlist(all_placebo, fill = TRUE)

# Order windows for plotting
results_all[, window := factor(window, levels = c("Day-0", "3-day mean", "7-day mean"))]

p <- ggplot(results_all, aes(x = date, y = beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red",
             linewidth = 0.6) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "steelblue") +
  geom_point(aes(shape = is_mou), size = 2.5) +
  geom_line(linewidth = 0.4) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  facet_grid(variable ~ window, scales = "free_y") +
  labs(
    title = "Placebo treatment dates: Weather x Post interaction",
    subtitle = "Red dashed line = actual MoU date (Feb 2017). Diamond = true estimate.",
    x = "Placebo treatment date",
    y = expression(hat(beta)[3] ~ "(interaction coefficient)")
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(BASE_DIR, "output", "figures", "placebo_beta3.pdf"),
       p, width = 14, height = 8)
ggsave(file.path(BASE_DIR, "output", "figures", "placebo_beta3.png"),
       p, width = 14, height = 8, dpi = 200)
cat("Saved: output/figures/placebo_beta3.pdf + .png\n")

# Save table
fwrite(results_all, file.path(BASE_DIR, "output", "tables", "placebo_results.csv"))
cat("Saved: output/tables/placebo_results.csv\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
