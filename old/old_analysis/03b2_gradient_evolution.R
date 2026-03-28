# 03b2_gradient_evolution.R
# =========================
# Two diagnostics for the gradient stability assumption:
#
# Part 1 — Gradient LEVEL over time:
#   Rolling 2-year window estimates of the SWH lag-1 coefficient
#   (deaths ~ SWH | FE). Shows the trajectory of the gradient.
#
# Part 2 — Gradient CHANGE (beta_3) as a function of sample width:
#   Expanding windows centered on the MoU date, estimating the
#   SWH x Post interaction (deaths ~ SWH + SWH:Post | FE).
#   Shows whether the estimated change is stable as more data
#   is included, or is driven by specific periods.
#
# Input:  data/processed/cmr_daily_weather_panel.RDS
# Output: output/figures/gradient_evolution_lag1.pdf / .png
#         output/figures/expanding_window_beta3.pdf / .png
#         output/tables/expanding_window_beta3.csv

library(fixest)
library(data.table)
library(ggplot2)

BASE_DIR <- here::here()
d <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                      "cmr_daily_weather_panel.RDS")))

MOU_DATE <- as.Date("2017-07-01")

cat("============================================================\n")
cat("GRADIENT EVOLUTION: SWH lag-1, rolling 2-year window\n")
cat("============================================================\n\n")

# ============================================================
# Rolling 2-year window, stepping in 6-month increments
# ============================================================
window_half <- 365  # ±1 year = 2-year window
window_centers <- seq.Date(as.Date("2015-01-01"), as.Date("2024-07-01"),
                           by = "6 months")

cat(sprintf("Rolling 2-year window: %d centers\n\n", length(window_centers)))

estimate_gradient <- function(dsub, fe_type = "weekly") {
  n_events <- sum(dsub$n_dead_missing > 0)
  if (n_events < 15) return(NULL)

  if (fe_type == "weekly") {
    dsub[, fe_fac := factor(week_year)]
  } else {
    dsub[, fe_fac := factor(month_year)]
  }

  warn_msg <- NULL
  m <- tryCatch(
    withCallingHandlers(
      fenegbin(n_dead_missing ~ swh_core_lag1 | fe_fac,
               data = dsub[!is.na(swh_core_lag1)], vcov = "hetero"),
      warning = function(w) { warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
    ),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)

  converged <- is.null(warn_msg) || !grepl("did not converge|singular", warn_msg)

  ct <- summary(m, vcov = "hetero")$coeftable
  swh_row <- which(rownames(ct) == "swh_core_lag1")
  if (length(swh_row) == 0) return(NULL)

  data.table(
    beta      = ct[swh_row, 1],
    se        = if (converged) ct[swh_row, 2] else NA_real_,
    p         = if (converged) ct[swh_row, 4] else NA_real_,
    n_events  = n_events,
    converged = converged
  )
}

# --- Weekly FE ---
cat("--- Weekly FE ---\n")
roll_weekly <- rbindlist(lapply(window_centers, function(center) {
  dsub <- d[date >= (center - window_half) & date <= (center + window_half)]
  res <- estimate_gradient(dsub, "weekly")
  if (!is.null(res)) res[, center := center]
  res
}))

if (nrow(roll_weekly) > 0) {
  roll_weekly[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se,
                      fe = "Week-year FE")]
  for (i in seq_len(nrow(roll_weekly))) {
    r <- roll_weekly[i]
    cat(sprintf("  %s: beta=%+.3f (SE=%.3f), n_events=%d%s\n",
        as.character(r$center), r$beta, ifelse(r$converged, r$se, NA),
        r$n_events, ifelse(r$converged, "", " [NC]")))
  }
}

# --- Monthly FE ---
cat("\n--- Monthly FE ---\n")
roll_monthly <- rbindlist(lapply(window_centers, function(center) {
  dsub <- d[date >= (center - window_half) & date <= (center + window_half)]
  res <- estimate_gradient(dsub, "monthly")
  if (!is.null(res)) res[, center := center]
  res
}))

if (nrow(roll_monthly) > 0) {
  roll_monthly[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se,
                        fe = "Month-year FE")]
  for (i in seq_len(nrow(roll_monthly))) {
    r <- roll_monthly[i]
    cat(sprintf("  %s: beta=%+.3f (SE=%.3f), n_events=%d%s\n",
        as.character(r$center), r$beta, ifelse(r$converged, r$se, NA),
        r$n_events, ifelse(r$converged, "", " [NC]")))
  }
}

# ============================================================
# Year-by-year estimates (week-year FE)
# ============================================================
cat("\n--- Year-by-year (week-year FE) ---\n")
all_years <- sort(unique(d$year))

yearly_wk <- rbindlist(lapply(all_years, function(yr) {
  dsub <- d[year == yr]
  res <- estimate_gradient(dsub, "weekly")
  if (!is.null(res)) {
    res[, center := as.Date(paste0(yr, "-07-01"))]
    res[, fe := "Year-by-year (week FE)"]
    cat(sprintf("  %d: beta=%+.3f, n_events=%d%s\n",
        yr, res$beta, res$n_events, ifelse(res$converged, "", " [NC]")))
  }
  res
}))

if (nrow(yearly_wk) > 0) {
  yearly_wk[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]
}

cat("\n--- Year-by-year (month-year FE) ---\n")

yearly_mo <- rbindlist(lapply(all_years, function(yr) {
  dsub <- d[year == yr]
  res <- estimate_gradient(dsub, "monthly")
  if (!is.null(res)) {
    res[, center := as.Date(paste0(yr, "-07-01"))]
    res[, fe := "Year-by-year (month FE)"]
    cat(sprintf("  %d: beta=%+.3f, n_events=%d%s\n",
        yr, res$beta, res$n_events, ifelse(res$converged, "", " [NC]")))
  }
  res
}))

if (nrow(yearly_mo) > 0) {
  yearly_mo[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]
}

# ============================================================
# Plot: rolling + year-by-year, both FE structures
# ============================================================
cat("\n--- Generating plots ---\n")

# --- Plot 1: Week-year FE (rolling + year-by-year) ---
roll_wk_conv <- roll_weekly[converged == TRUE]
yr_wk_conv   <- yearly_wk[converged == TRUE]

p_wk <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red",
             linewidth = 0.6) +
  annotate("text", x = MOU_DATE,
           y = max(c(roll_wk_conv$ci_hi, yr_wk_conv$ci_hi), na.rm = TRUE) * 0.9,
           label = "MoU (Jul 2017)", colour = "red", hjust = -0.1, size = 3) +
  # Rolling window
  geom_ribbon(data = roll_wk_conv, aes(x = center, ymin = ci_lo, ymax = ci_hi),
              alpha = 0.12, fill = "steelblue") +
  geom_line(data = roll_wk_conv, aes(x = center, y = beta),
            colour = "steelblue", linewidth = 0.5) +
  geom_point(data = roll_wk_conv, aes(x = center, y = beta),
             colour = "steelblue", size = 2) +
  # Year-by-year
  geom_errorbar(data = yr_wk_conv, aes(x = center, ymin = ci_lo, ymax = ci_hi),
                width = 40, colour = "firebrick", linewidth = 0.4) +
  geom_point(data = yr_wk_conv, aes(x = center, y = beta),
             colour = "firebrick", size = 3) +
  labs(
    title = "SWH lag-1 mortality gradient over time (week-year FE)",
    subtitle = "Blue: rolling 2-year window. Red: year-by-year. NegBin, drowning/suspected, core corridor.",
    x = NULL,
    y = expression("SWH lag-1 coefficient " * hat(beta)[SWH])
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(BASE_DIR, "output", "figures", "gradient_evolution_lag1_weekly.pdf"),
       p_wk, width = 10, height = 6)
ggsave(file.path(BASE_DIR, "output", "figures", "gradient_evolution_lag1_weekly.png"),
       p_wk, width = 10, height = 6, dpi = 200)
cat("Saved: gradient_evolution_lag1_weekly.pdf + .png\n")

# --- Plot 2: Month-year FE (rolling + year-by-year) ---
roll_mo_conv <- roll_monthly[converged == TRUE]
yr_mo_conv   <- yearly_mo[converged == TRUE]

p_mo <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red",
             linewidth = 0.6) +
  annotate("text", x = MOU_DATE,
           y = max(c(roll_mo_conv$ci_hi, yr_mo_conv$ci_hi), na.rm = TRUE) * 0.9,
           label = "MoU (Jul 2017)", colour = "red", hjust = -0.1, size = 3) +
  # Rolling window
  geom_ribbon(data = roll_mo_conv, aes(x = center, ymin = ci_lo, ymax = ci_hi),
              alpha = 0.12, fill = "darkorange") +
  geom_line(data = roll_mo_conv, aes(x = center, y = beta),
            colour = "darkorange", linewidth = 0.5) +
  geom_point(data = roll_mo_conv, aes(x = center, y = beta),
             colour = "darkorange", size = 2) +
  # Year-by-year
  geom_errorbar(data = yr_mo_conv, aes(x = center, ymin = ci_lo, ymax = ci_hi),
                width = 40, colour = "firebrick", linewidth = 0.4) +
  geom_point(data = yr_mo_conv, aes(x = center, y = beta),
             colour = "firebrick", size = 3) +
  labs(
    title = "SWH lag-1 mortality gradient over time (month-year FE)",
    subtitle = "Orange: rolling 2-year window. Red: year-by-year. NegBin, drowning/suspected, core corridor.",
    x = NULL,
    y = expression("SWH lag-1 coefficient " * hat(beta)[SWH])
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(BASE_DIR, "output", "figures", "gradient_evolution_lag1_monthly.pdf"),
       p_mo, width = 10, height = 6)
ggsave(file.path(BASE_DIR, "output", "figures", "gradient_evolution_lag1_monthly.png"),
       p_mo, width = 10, height = 6, dpi = 200)
cat("Saved: gradient_evolution_lag1_monthly.pdf + .png\n")

# --- Plot 3: Combined (both FE rolling, no year-by-year for clarity) ---
plot_dt <- rbind(roll_wk_conv, roll_mo_conv)

p_both <- ggplot(plot_dt, aes(x = center, y = beta, colour = fe, fill = fe)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red",
             linewidth = 0.6) +
  annotate("text", x = MOU_DATE, y = max(plot_dt$ci_hi, na.rm = TRUE) * 0.9,
           label = "MoU (Jul 2017)", colour = "red", hjust = -0.1, size = 3) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("Week-year FE" = "steelblue",
                                  "Month-year FE" = "darkorange")) +
  scale_fill_manual(values = c("Week-year FE" = "steelblue",
                                "Month-year FE" = "darkorange")) +
  labs(
    title = "SWH lag-1 mortality gradient over time",
    subtitle = "Rolling 2-year window, NegBin. Drowning/suspected drowning, core corridor.",
    x = NULL,
    y = expression("SWH lag-1 coefficient " * hat(beta)[SWH]),
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(BASE_DIR, "output", "figures", "gradient_evolution_lag1.pdf"),
       p_both, width = 10, height = 6)
ggsave(file.path(BASE_DIR, "output", "figures", "gradient_evolution_lag1.png"),
       p_both, width = 10, height = 6, dpi = 200)
cat("Saved: gradient_evolution_lag1.pdf + .png\n")

# ============================================================
# Part 2: Expanding-window beta_3 around MoU
# ============================================================
# The rolling window above estimates the LEVEL of the gradient.
# This section estimates the CHANGE (beta_3 = SWH x Post interaction)
# from windows that expand symmetrically around the MoU date.
#
# A stable beta_3 across window widths means the result does not
# depend on which years are included. A beta_3 that shrinks as
# distant data is added means the effect is local to the MoU period.

cat("\n============================================================\n")
cat("EXPANDING-WINDOW BETA_3\n")
cat("============================================================\n\n")

expand_halfdays <- seq(180, 1800, by = 90)

estimate_interaction <- function(dsub, fe_type = "weekly") {
  n_events_pre  <- sum(dsub$n_dead_missing[dsub$post_mou == 0] > 0)
  n_events_post <- sum(dsub$n_dead_missing[dsub$post_mou == 1] > 0)
  if (n_events_pre < 10 || n_events_post < 10) return(NULL)

  if (fe_type == "weekly") {
    dsub[, fe_fac := factor(week_year)]
  } else {
    dsub[, fe_fac := factor(month_year)]
  }

  warn_msg <- NULL
  m <- tryCatch(
    withCallingHandlers(
      fenegbin(n_dead_missing ~ swh_core_lag1 + swh_core_lag1:post_mou | fe_fac,
               data = dsub[!is.na(swh_core_lag1)], vcov = "hetero"),
      warning = function(w) { warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
    ),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)

  converged <- is.null(warn_msg) || !grepl("did not converge|singular", warn_msg)
  ct <- summary(m, vcov = "hetero")$coeftable
  int_row <- grep("post_mou", rownames(ct))
  if (length(int_row) == 0) return(NULL)

  data.table(
    beta3     = ct[int_row, 1],
    se        = if (converged) ct[int_row, 2] else NA_real_,
    p         = if (converged) ct[int_row, 4] else NA_real_,
    n_pre     = sum(dsub$post_mou == 0),
    n_post    = sum(dsub$post_mou == 1),
    n_events_pre  = n_events_pre,
    n_events_post = n_events_post,
    converged = converged
  )
}

# --- Week-year FE ---
cat("--- Expanding window, week-year FE ---\n")
expand_wk <- rbindlist(lapply(expand_halfdays, function(hd) {
  w_start <- max(MOU_DATE - hd, min(d$date))
  w_end   <- min(MOU_DATE + hd, max(d$date))
  dsub <- d[date >= w_start & date <= w_end]
  res <- estimate_interaction(dsub, "weekly")
  if (!is.null(res)) {
    res[, `:=`(halfdays = hd, halfyears = round(hd / 365.25, 2),
               w_start = w_start, w_end = w_end, fe = "Week-year FE")]
  }
  res
}))

if (nrow(expand_wk) > 0) {
  expand_wk[, `:=`(ci_lo = beta3 - 1.96 * se, ci_hi = beta3 + 1.96 * se)]
  for (i in seq_len(nrow(expand_wk))) {
    r <- expand_wk[i]
    cat(sprintf("  +/-%4.1f yr: beta3=%+.3f (SE=%.3f, p=%.3f), events pre=%d post=%d%s\n",
        r$halfyears, r$beta3, ifelse(r$converged, r$se, NA),
        ifelse(r$converged, r$p, NA),
        r$n_events_pre, r$n_events_post,
        ifelse(r$converged, "", " [NC]")))
  }
}

# --- Month-year FE ---
cat("\n--- Expanding window, month-year FE ---\n")
expand_mo <- rbindlist(lapply(expand_halfdays, function(hd) {
  w_start <- max(MOU_DATE - hd, min(d$date))
  w_end   <- min(MOU_DATE + hd, max(d$date))
  dsub <- d[date >= w_start & date <= w_end]
  res <- estimate_interaction(dsub, "monthly")
  if (!is.null(res)) {
    res[, `:=`(halfdays = hd, halfyears = round(hd / 365.25, 2),
               w_start = w_start, w_end = w_end, fe = "Month-year FE")]
  }
  res
}))

if (nrow(expand_mo) > 0) {
  expand_mo[, `:=`(ci_lo = beta3 - 1.96 * se, ci_hi = beta3 + 1.96 * se)]
  for (i in seq_len(nrow(expand_mo))) {
    r <- expand_mo[i]
    cat(sprintf("  +/-%4.1f yr: beta3=%+.3f (SE=%.3f, p=%.3f), events pre=%d post=%d%s\n",
        r$halfyears, r$beta3, ifelse(r$converged, r$se, NA),
        ifelse(r$converged, r$p, NA),
        r$n_events_pre, r$n_events_post,
        ifelse(r$converged, "", " [NC]")))
  }
}

# --- Save expanding-window results ---
expand_all <- rbind(expand_wk, expand_mo, fill = TRUE)
fwrite(expand_all, file.path(BASE_DIR, "output", "tables",
                              "expanding_window_beta3.csv"))
cat("\nSaved: output/tables/expanding_window_beta3.csv\n")

# --- Expanding-window plot ---
cat("\n--- Generating expanding-window plot ---\n")

expand_conv <- expand_all[converged == TRUE]

if (nrow(expand_conv) > 0) {
  p_expand <- ggplot(expand_conv, aes(x = halfyears, y = beta3,
                                       colour = fe, fill = fe)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    scale_colour_manual(values = c("Week-year FE" = "steelblue",
                                    "Month-year FE" = "darkorange")) +
    scale_fill_manual(values = c("Week-year FE" = "steelblue",
                                  "Month-year FE" = "darkorange")) +
    labs(
      title = expression("Expanding-window " * hat(beta)[3] *
                          " (SWH lag-1 × Post)"),
      subtitle = "Windows expand symmetrically from MoU (Jul 2017). NegBin, hetero-robust SE.",
      x = "Window halfwidth (years from MoU)",
      y = expression(hat(beta)[3] * " (SWH × Post interaction)"),
      colour = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")

  ggsave(file.path(BASE_DIR, "output", "figures", "expanding_window_beta3.pdf"),
         p_expand, width = 10, height = 6)
  ggsave(file.path(BASE_DIR, "output", "figures", "expanding_window_beta3.png"),
         p_expand, width = 10, height = 6, dpi = 200)
  cat("Saved: expanding_window_beta3.pdf + .png\n")
}

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
