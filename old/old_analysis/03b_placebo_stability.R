# 03b_placebo_stability.R
# =======================
# Placebo treatment dates and gradient evolution over time.
#
# These tests validate the primary specification from 03c:
#   daily panel, core SWH day-0, NegBin, week-year FE.
#
# Test 1: Placebo treatment dates
#   Re-estimate SWH_core x Post at many candidate break dates.
#   If the MoU date is special, beta_3 at Jul 2017 should stand out
#   relative to other dates.
#
# Input:  data/processed/cmr_daily_weather_panel.RDS
# Output: output/figures/placebo_beta3.pdf
#         output/tables/placebo_results.csv
#         printed diagnostics

library(fixest)
library(data.table)
library(ggplot2)

BASE_DIR <- here::here()
d <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                      "cmr_daily_weather_panel.RDS")))

MOU_DATE <- as.Date("2017-07-01")

cat("============================================================\n")
cat("PLACEBO AND STABILITY TESTS (daily panel, primary spec)\n")
cat("============================================================\n\n")

cat(sprintf("Panel: %d days, outcome mean=%.2f\n\n", nrow(d), mean(d$n_dead_missing)))


# ============================================================
# Test 1: Placebo treatment dates — SWH_core x Post
# ============================================================
cat("============================================================\n")
cat("TEST 1: PLACEBO TREATMENT DATES\n")
cat("============================================================\n\n")

# Candidate dates: quarterly from 2015-Q1 to 2019-Q1
placebo_dates <- as.Date(c(
  "2015-01-01", "2015-04-01", "2015-07-01", "2015-10-01",
  "2016-01-01", "2016-04-01", "2016-07-01", "2016-10-01",
  "2017-07-01",  # actual MoU
  "2017-07-01", "2017-10-01",
  "2018-01-01", "2018-04-01", "2018-07-01",
  "2019-01-01"
))

placebo_results <- rbindlist(lapply(placebo_dates, function(pd) {
  d[, post_placebo := as.integer(date >= pd)]

  # Need sufficient pre/post observations for week-year FE estimation
  n_pre  <- sum(d$post_placebo == 0)
  n_post <- sum(d$post_placebo == 1)
  if (n_pre < 60 || n_post < 60) return(NULL)

  # Week-year FE: some cells will be split by the placebo cut.
  # This is expected — it's the within-week variation that identifies.
  # Capture warnings to detect convergence failures.
  warn_msg <- NULL
  m <- tryCatch(
    withCallingHandlers(
      fenegbin(n_dead_missing ~ swh_core + swh_core:post_placebo | week_year_fac,
               data = d, vcov = "hetero"),
      warning = function(w) { warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
    ),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)

  # Check for convergence failure (singular information matrix → NA SEs)
  converged <- is.null(warn_msg) || !grepl("did not converge|singular", warn_msg)

  ct <- summary(m, vcov = "hetero")$coeftable
  int_row <- grep("swh_core:post_placebo|post_placebo:swh_core",
                  rownames(ct), fixed = FALSE)
  if (length(int_row) == 0) return(NULL)

  data.table(
    date      = pd,
    is_mou    = pd == MOU_DATE,
    beta      = ct[int_row, 1],
    se        = if (converged) ct[int_row, 2] else NA_real_,
    p         = if (converged) ct[int_row, 4] else NA_real_,
    irr       = exp(ct[int_row, 1]),
    n_obs     = nobs(m),
    converged = converged
  )
}))

if (nrow(placebo_results) > 0) {
  placebo_results[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]

  n_converged <- sum(placebo_results$converged)
  n_failed    <- sum(!placebo_results$converged)
  cat(sprintf("Placebo estimates (SWH_core x Post, NegBin, week-year FE):\n"))
  cat(sprintf("  %d converged, %d failed (marked [NC])\n\n", n_converged, n_failed))
  cat(sprintf("  %-12s %+8s %7s %7s %8s %3s\n",
      "Date", "Beta", "SE", "p", "IRR", ""))
  cat(paste0("  ", paste(rep("-", 60), collapse = "")), "\n")
  for (i in seq_len(nrow(placebo_results))) {
    r <- placebo_results[i]
    marker <- ifelse(r$is_mou, " <-- MoU", "")
    if (!r$converged) {
      cat(sprintf("  %-12s %+8.4f %7s %7s %8.4f [NC]%s\n",
          as.character(r$date), r$beta, "  NA", "  NA", r$irr, marker))
    } else {
      stars <- ifelse(r$p < 0.01, "***",
               ifelse(r$p < 0.05, "**",
               ifelse(r$p < 0.1, "*", "")))
      cat(sprintf("  %-12s %+8.4f %7.4f %7.4f %8.4f %3s%s\n",
          as.character(r$date), r$beta, r$se, r$p, r$irr, stars, marker))
    }
  }
}


# ============================================================
# Plot: placebo coefficient path
# ============================================================
cat("\n============================================================\n")
cat("GENERATING PLOTS\n")
cat("============================================================\n\n")

if (nrow(placebo_results) > 0) {

  # Plot only converged estimates (those with CIs)
  plot_dt <- placebo_results[converged == TRUE]

  p <- ggplot(placebo_results, aes(x = date, y = beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red",
               linewidth = 0.6) +
    geom_ribbon(data = plot_dt, aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15,
                fill = "steelblue") +
    geom_line(linewidth = 0.4) +
    geom_point(aes(shape = is_mou, size = is_mou,
                   colour = converged)) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4), guide = "none") +
    scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey60"),
                        guide = "none") +
    labs(
      title = expression("Placebo treatment dates: " * hat(beta)[3] *
                          " (SWH"[core] * " × Post)"),
      subtitle = "Daily panel, NegBin, week-year FE. Red line = actual MoU (Jul 2017).\nDiamond = true estimate. Grey dots = NegBin did not converge (no CI).",
      x = "Placebo treatment date",
      y = expression(hat(beta)[3])
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(BASE_DIR, "output", "figures", "placebo_beta3.pdf"),
         p, width = 10, height = 6)
  ggsave(file.path(BASE_DIR, "output", "figures", "placebo_beta3.png"),
         p, width = 10, height = 6, dpi = 200)
  cat("Saved: output/figures/placebo_beta3.pdf + .png\n")
}

# Save table
fwrite(placebo_results, file.path(BASE_DIR, "output", "tables",
                                   "placebo_results.csv"))
cat("Saved: output/tables/placebo_results.csv\n")


# ============================================================
# Test 3: Gradient evolution over time
# ============================================================
cat("\n============================================================\n")
cat("TEST 3: GRADIENT EVOLUTION (rolling window + year-by-year)\n")
cat("============================================================\n\n")

# a) Rolling 2-year window: estimate SWH coefficient from NegBin
#    with week-year FE, stepping in 6-month increments.
all_years <- sort(unique(d$year))
window_half <- 365  # ±1 year = 2-year window

# Window centers: every 6 months from mid-2015 to mid-2024
window_centers <- seq.Date(as.Date("2015-01-01"), as.Date("2024-07-01"),
                           by = "6 months")

cat(sprintf("Rolling 2-year window: %d centers\n", length(window_centers)))

rolling_results <- rbindlist(lapply(window_centers, function(center) {
  w_start <- center - window_half
  w_end   <- center + window_half
  dsub <- d[date >= w_start & date <= w_end]

  # Need enough event-days for estimation
  n_events <- sum(dsub$n_dead_missing > 0)
  if (n_events < 20) return(NULL)

  # Recode FE within window
  dsub[, wyfac := factor(week_year)]

  roll_warn <- NULL
  m <- tryCatch(
    withCallingHandlers(
      fenegbin(n_dead_missing ~ swh_core | wyfac,
               data = dsub, vcov = "hetero"),
      warning = function(w) { roll_warn <<- conditionMessage(w); invokeRestart("muffleWarning") }
    ),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)

  converged <- is.null(roll_warn) || !grepl("did not converge|singular", roll_warn)

  ct <- summary(m, vcov = "hetero")$coeftable
  swh_row <- which(rownames(ct) == "swh_core")
  if (length(swh_row) == 0) return(NULL)

  data.table(
    center    = center,
    beta      = ct[swh_row, 1],
    se        = if (converged) ct[swh_row, 2] else NA_real_,
    p         = if (converged) ct[swh_row, 4] else NA_real_,
    n_days    = nrow(dsub),
    n_events  = n_events,
    converged = converged
  )
}))

if (nrow(rolling_results) > 0) {
  rolling_results[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]

  cat("\nRolling SWH gradient:\n")
  for (i in seq_len(nrow(rolling_results))) {
    r <- rolling_results[i]
    flag <- if (!r$converged) " [NC]" else ""
    cat(sprintf("  %s: beta=%+.4f, n_events=%d%s\n",
        as.character(r$center), r$beta, r$n_events, flag))
  }
}

# b) Year-by-year SWH coefficient (full sample, each year separately)
cat("\nYear-by-year SWH coefficient (NegBin with week-year FE):\n")

yearly_results <- rbindlist(lapply(all_years, function(yr) {
  dsub <- d[year == yr]
  n_events <- sum(dsub$n_dead_missing > 0)
  if (n_events < 10) {
    cat(sprintf("  %d: n_events=%d (too few)\n", yr, n_events))
    return(NULL)
  }

  dsub[, wyfac := factor(week_year)]

  yr_warn <- NULL
  m <- tryCatch(
    withCallingHandlers(
      fenegbin(n_dead_missing ~ swh_core | wyfac,
               data = dsub, vcov = "hetero"),
      warning = function(w) { yr_warn <<- conditionMessage(w); invokeRestart("muffleWarning") }
    ),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)

  converged <- is.null(yr_warn) || !grepl("did not converge|singular", yr_warn)

  ct <- summary(m, vcov = "hetero")$coeftable
  swh_row <- which(rownames(ct) == "swh_core")
  if (length(swh_row) == 0) return(NULL)

  res <- data.table(
    year      = yr,
    beta      = ct[swh_row, 1],
    se        = if (converged) ct[swh_row, 2] else NA_real_,
    p         = if (converged) ct[swh_row, 4] else NA_real_,
    n_events  = n_events,
    converged = converged
  )

  flag <- if (!converged) " [NC]" else ""
  cat(sprintf("  %d: beta=%+.4f, SE=%.4f, n_events=%d%s\n",
      yr, res$beta, ifelse(converged, res$se, NA), n_events, flag))
  res
}))

if (nrow(yearly_results) > 0) {
  yearly_results[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]
}

# c) Gradient evolution plot: rolling window (blue) + year-by-year (red)
cat("\n--- Generating gradient evolution plots ---\n")

if (nrow(rolling_results) > 0) {
  roll_conv <- rolling_results[converged == TRUE]

  p_roll <- ggplot(roll_conv, aes(x = center, y = beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red",
               linewidth = 0.6) +
    annotate("text", x = MOU_DATE, y = max(roll_conv$ci_hi, na.rm = TRUE),
             label = "MoU", colour = "red", hjust = -0.1, size = 3) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15,
                fill = "steelblue") +
    geom_line(colour = "steelblue", linewidth = 0.5) +
    geom_point(size = 2) +
    labs(
      title = expression("Rolling 2-year window: SWH-mortality gradient over time"),
      subtitle = "NegBin with week-year FE within each window. Core geography.",
      x = "Window center",
      y = expression("SWH coefficient " * hat(beta)[SWH])
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(BASE_DIR, "output", "figures", "gradient_rolling.pdf"),
         p_roll, width = 10, height = 6)
  ggsave(file.path(BASE_DIR, "output", "figures", "gradient_rolling.png"),
         p_roll, width = 10, height = 6, dpi = 200)
  cat("Saved: output/figures/gradient_rolling.pdf + .png\n")
}

if (nrow(yearly_results) > 0) {
  yr_conv <- yearly_results[converged == TRUE]

  p_yr <- ggplot(yr_conv, aes(x = year, y = beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 2017, linetype = "dashed", colour = "red",
               linewidth = 0.6) +
    annotate("text", x = 2017, y = max(yr_conv$ci_hi, na.rm = TRUE),
             label = "MoU", colour = "red", hjust = -0.1, size = 3) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                  colour = "firebrick") +
    geom_point(size = 3, colour = "firebrick") +
    scale_x_continuous(breaks = all_years) +
    labs(
      title = "SWH-mortality gradient by year",
      subtitle = "Year-specific SWH coefficient from NegBin with week-year FE. Core geography.",
      x = "Year",
      y = expression("SWH coefficient " * hat(beta)[SWH])
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(BASE_DIR, "output", "figures", "gradient_by_year.pdf"),
         p_yr, width = 10, height = 6)
  ggsave(file.path(BASE_DIR, "output", "figures", "gradient_by_year.png"),
         p_yr, width = 10, height = 6, dpi = 200)
  cat("Saved: output/figures/gradient_by_year.pdf + .png\n")
}

# Combined evolution plot
if (nrow(rolling_results) > 0 && nrow(yearly_results) > 0) {
  roll_conv <- rolling_results[converged == TRUE]
  yr_conv   <- yearly_results[converged == TRUE]
  yr_conv[, center := as.Date(paste0(year, "-07-01"))]

  p_evo <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red",
               linewidth = 0.6) +
    annotate("text", x = MOU_DATE,
             y = max(c(roll_conv$ci_hi, yr_conv$ci_hi), na.rm = TRUE),
             label = "MoU (Jul 2017)", colour = "red", hjust = -0.1, size = 3) +
    # Rolling window
    geom_ribbon(data = roll_conv, aes(x = center, ymin = ci_lo, ymax = ci_hi),
                alpha = 0.12, fill = "steelblue") +
    geom_line(data = roll_conv, aes(x = center, y = beta),
              colour = "steelblue", linewidth = 0.5) +
    geom_point(data = roll_conv, aes(x = center, y = beta),
               colour = "steelblue", size = 2) +
    # Year-by-year
    geom_errorbar(data = yr_conv, aes(x = center, ymin = ci_lo, ymax = ci_hi),
                  width = 40, colour = "firebrick", linewidth = 0.4) +
    geom_point(data = yr_conv, aes(x = center, y = beta),
               colour = "firebrick", size = 3) +
    labs(
      title = "Evolution of the SWH-mortality gradient over time",
      subtitle = "Blue: rolling 2-year window. Red: year-by-year NegBin. Week-year FE, core geography.",
      x = NULL,
      y = expression("SWH coefficient " * hat(beta)[SWH])
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(BASE_DIR, "output", "figures", "gradient_evolution.pdf"),
         p_evo, width = 10, height = 6)
  ggsave(file.path(BASE_DIR, "output", "figures", "gradient_evolution.png"),
         p_evo, width = 10, height = 6, dpi = 200)
  cat("Saved: output/figures/gradient_evolution.pdf + .png\n")
}


cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
