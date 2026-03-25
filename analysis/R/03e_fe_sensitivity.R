# 03e_fe_sensitivity.R
# ====================
# FE sensitivity analysis: how does beta_3 (SWH x Post) change with FE granularity?
#
# The key empirical finding is that coarse FE (year+month) give negative beta_3,
# while fine FE (week-year) give positive beta_3.  This reflects two channels:
#   - Across-week variation (coarse FE): deterrence — rough weeks → fewer crossings
#   - Within-week variation (fine FE):   danger    — rough days → more deaths per crossing
#
# This script systematically traces beta_3 across FE granularity:
#   Part 1: Daily panel (core SWH)
#   Part 2: Event level (incident-location SWH)
#   Part 3: Comparison figure
#
# Input:  data/processed/cmr_daily_weather_panel.RDS
#         data/processed/cmr_events_with_weather.RDS
# Output: output/tables/fe_sensitivity.csv
#         output/figures/fe_sensitivity.pdf

library(fixest)
library(data.table)
library(ggplot2)

BASE_DIR <- here::here()

# ============================================================
# Helpers
# ============================================================
extract_interaction <- function(model, spec_label, model_type, fe_label, n_fe_cells,
                                vcov_type = "hetero") {
  ct <- summary(model, vcov = vcov_type)$coeftable
  int_rows <- grep(":post_mou|post_mou:", rownames(ct))
  if (length(int_rows) == 0) return(NULL)
  # Take first interaction only (the SWH one)
  r <- int_rows[1]
  data.table(
    model_type = model_type,
    fe = fe_label,
    spec = spec_label,
    n_fe_cells = n_fe_cells,
    beta = ct[r, 1], se = ct[r, 2], p = ct[r, 4],
    irr = exp(ct[r, 1]),
    ci_lo = ct[r, 1] - 1.96 * ct[r, 2],
    ci_hi = ct[r, 1] + 1.96 * ct[r, 2],
    n_obs = nobs(model)
  )
}


# ============================================================
# PART 1: DAILY PANEL — Core SWH x Post across FE levels
# ============================================================
cat("============================================================\n")
cat("PART 1: DAILY PANEL FE SENSITIVITY (core SWH day-0)\n")
cat("============================================================\n\n")

d <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                      "cmr_daily_weather_panel.RDS")))

# Ensure all FE variables exist
d[, month_fac := factor(month)]
d[, year_fac := factor(year)]
d[, half_year := paste0(year, ifelse(month <= 6, "H1", "H2"))]
d[, half_year_fac := factor(half_year)]
d[, quarter := paste0(year, "Q", ceiling(month / 3))]
d[, quarter_fac := factor(quarter)]
d[, year_month := paste0(year, "_", sprintf("%02d", month))]
d[, year_month_fac := factor(year_month)]
# week_year_fac already exists

cat(sprintf("Panel: %d days, outcome mean=%.2f, var/mean=%.1f\n\n",
    nrow(d), mean(d$n_dead_missing), var(d$n_dead_missing) / mean(d$n_dead_missing)))

# FE ladder (coarsest to finest)
fe_specs_daily <- list(
  list(label = "Month-of-year", fml = n_dead_missing ~ swh_core * post_mou + wind_core | month_fac,
       cells_var = "month_fac"),
  list(label = "Year", fml = n_dead_missing ~ swh_core * post_mou + wind_core | year_fac,
       cells_var = "year_fac"),
  list(label = "Year + Month", fml = n_dead_missing ~ swh_core * post_mou + wind_core | year_fac + month_fac,
       cells_var = c("year_fac", "month_fac")),
  list(label = "Half-year", fml = n_dead_missing ~ swh_core * post_mou + wind_core | half_year_fac,
       cells_var = "half_year_fac"),
  list(label = "Quarter-year", fml = n_dead_missing ~ swh_core * post_mou + wind_core | quarter_fac,
       cells_var = "quarter_fac"),
  list(label = "Year x Month", fml = n_dead_missing ~ swh_core * post_mou + wind_core | year_month_fac,
       cells_var = "year_month_fac"),
  list(label = "Week-year", fml = n_dead_missing ~ swh_core * post_mou + wind_core | week_year_fac,
       cells_var = "week_year_fac")
)

daily_results <- list()

for (sp in fe_specs_daily) {
  n_cells <- if (length(sp$cells_var) == 1) {
    uniqueN(d[[sp$cells_var]])
  } else {
    # Additive FE: report sum of unique levels
    sum(sapply(sp$cells_var, function(v) uniqueN(d[[v]])))
  }

  m <- tryCatch(
    fenegbin(sp$fml, data = d, vcov = "hetero"),
    error = function(e) { cat(sprintf("  %s: FAILED — %s\n", sp$label, e$message)); NULL }
  )

  if (!is.null(m)) {
    res <- extract_interaction(m, sp$label, "Daily panel", sp$label, n_cells)
    if (!is.null(res)) {
      daily_results[[length(daily_results) + 1]] <- res
      stars <- ifelse(res$p < 0.01, "***", ifelse(res$p < 0.05, "**",
               ifelse(res$p < 0.1, "*", "")))
      cat(sprintf("  %-20s (%3d cells): beta=%+7.4f  SE=%6.4f  p=%6.4f %3s  N=%d\n",
          sp$label, n_cells, res$beta, res$se, res$p, stars, res$n_obs))
    }
  }
}

daily_dt <- rbindlist(daily_results)
cat(sprintf("\nSign flip: coarse FE → negative, fine FE → positive\n"))
cat(sprintf("  Year+Month: beta=%+.4f\n", daily_dt[fe == "Year + Month", beta]))
cat(sprintf("  Week-year:  beta=%+.4f\n", daily_dt[fe == "Week-year", beta]))


# ============================================================
# PART 2: EVENT LEVEL — Incident-location SWH x Post
# ============================================================
cat("\n============================================================\n")
cat("PART 2: EVENT-LEVEL FE SENSITIVITY (incident-location SWH day-0)\n")
cat("============================================================\n\n")

df <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                       "cmr_events_with_weather.RDS")))
df[, grid_1deg := paste0(sprintf("%.0f", round(grid_lat)), "_",
                          sprintf("%.0f", round(grid_lon)))]
df[, month_fac := factor(month(date))]
df[, year_fac := factor(year)]
df[, half_year := paste0(year, ifelse(month(date) <= 6, "H1", "H2"))]
df[, half_year_fac := factor(half_year)]
df[, quarter := paste0(year, "Q", ceiling(month(date) / 3))]
df[, quarter_fac := factor(quarter)]
df[, year_month := paste0(year, "_", sprintf("%02d", month(date)))]
df[, year_month_fac := factor(year_month)]

cat(sprintf("Events: %d (pre=%d, post=%d), grid cells (1deg): %d\n\n",
    nrow(df), sum(df$post_mou == 0), sum(df$post_mou == 1), uniqueN(df$grid_1deg)))

# FE ladder with grid (coarsest to finest)
fe_specs_event <- list(
  list(label = "Grid + Month", fml = dead_missing ~ swh_day0 * post_mou + wind_day0 | grid_1deg + month_fac,
       cells_expr = quote(uniqueN(df$grid_1deg) + uniqueN(df$month_fac))),
  list(label = "Grid + Year", fml = dead_missing ~ swh_day0 * post_mou + wind_day0 | grid_1deg + year_fac,
       cells_expr = quote(uniqueN(df$grid_1deg) + uniqueN(df$year_fac))),
  list(label = "Grid + Year + Month", fml = dead_missing ~ swh_day0 * post_mou + wind_day0 | grid_1deg + year_fac + month_fac,
       cells_expr = quote(uniqueN(df$grid_1deg) + uniqueN(df$year_fac) + uniqueN(df$month_fac))),
  list(label = "Grid + Half-year", fml = dead_missing ~ swh_day0 * post_mou + wind_day0 | grid_1deg + half_year_fac,
       cells_expr = quote(uniqueN(df$grid_1deg) + uniqueN(df$half_year_fac))),
  list(label = "Grid + Quarter", fml = dead_missing ~ swh_day0 * post_mou + wind_day0 | grid_1deg + quarter_fac,
       cells_expr = quote(uniqueN(df$grid_1deg) + uniqueN(df$quarter_fac))),
  list(label = "Grid + Year x Month", fml = dead_missing ~ swh_day0 * post_mou + wind_day0 | grid_1deg + year_month_fac,
       cells_expr = quote(uniqueN(df$grid_1deg) + uniqueN(df$year_month_fac)))
)

event_results <- list()

for (sp in fe_specs_event) {
  n_cells <- eval(sp$cells_expr)

  m <- tryCatch(
    fenegbin(sp$fml, data = df, vcov = "hetero"),
    error = function(e) { cat(sprintf("  %s: FAILED — %s\n", sp$label, e$message)); NULL }
  )

  if (!is.null(m)) {
    res <- extract_interaction(m, sp$label, "Event level", sp$label, n_cells)
    if (!is.null(res)) {
      event_results[[length(event_results) + 1]] <- res
      stars <- ifelse(res$p < 0.01, "***", ifelse(res$p < 0.05, "**",
               ifelse(res$p < 0.1, "*", "")))
      cat(sprintf("  %-25s (%3d FE): beta=%+7.4f  SE=%6.4f  p=%6.4f %3s  N=%d\n",
          sp$label, n_cells, res$beta, res$se, res$p, stars, res$n_obs))
    }
  }
}

event_dt <- rbindlist(event_results)


# ============================================================
# PART 3: Combined results and figure
# ============================================================
cat("\n============================================================\n")
cat("PART 3: COMBINED RESULTS\n")
cat("============================================================\n\n")

all_results <- rbind(daily_dt, event_dt, fill = TRUE)

# Order FE from coarsest to finest
fe_order <- c("Month-of-year", "Year", "Grid + Month", "Grid + Year",
              "Year + Month", "Grid + Year + Month",
              "Half-year", "Grid + Half-year",
              "Quarter-year", "Grid + Quarter",
              "Year x Month", "Grid + Year x Month",
              "Week-year")
all_results[, fe := factor(fe, levels = fe_order)]

# Print combined table
cat(sprintf("%-12s %-25s %3s %+8s %7s %7s %8s %5s\n",
    "Model", "FE", "Cells", "Beta", "SE", "p", "IRR", "N"))
cat(paste(rep("-", 85), collapse = ""), "\n")
for (i in seq_len(nrow(all_results))) {
  r <- all_results[i]
  stars <- ifelse(r$p < 0.01, "***", ifelse(r$p < 0.05, "**",
           ifelse(r$p < 0.1, "*", "")))
  cat(sprintf("%-12s %-25s %3d %+8.4f %7.4f %7.4f %8.4f %5d %s\n",
      r$model_type, as.character(r$fe), r$n_fe_cells,
      r$beta, r$se, r$p, r$irr, r$n_obs, stars))
}

# Save
fwrite(all_results, file.path(BASE_DIR, "output", "tables", "fe_sensitivity.csv"))
cat("\nSaved: output/tables/fe_sensitivity.csv\n")


# ============================================================
# Figure: beta_3 across FE granularity
# ============================================================
cat("\n--- Generating figure ---\n")

# For plotting, use n_fe_cells as x (continuous measure of FE granularity)
plot_dt <- copy(all_results)
plot_dt[, fe_label := paste0(fe, "\n(", n_fe_cells, " cells)")]

# Separate panels for daily vs event
p <- ggplot(plot_dt, aes(x = n_fe_cells, y = beta, colour = model_type,
                          shape = model_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi),
                  size = 0.6, position = position_dodge(width = 15)) +
  geom_line(aes(group = model_type), linewidth = 0.4, alpha = 0.5) +
  scale_colour_manual(values = c("Daily panel" = "#2166AC", "Event level" = "#B2182B")) +
  scale_shape_manual(values = c("Daily panel" = 16, "Event level" = 17)) +
  geom_text(aes(label = as.character(fe)), size = 2.3, vjust = -1.5, hjust = 0.5,
            show.legend = FALSE) +
  labs(
    title = expression("FE sensitivity: " * hat(beta)[3] * " (SWH × Post) across fixed-effect granularity"),
    subtitle = "Daily panel: core SWH [11-15,32-36]. Event level: incident-location SWH.\n95% CI, hetero-robust SEs. NegBin.",
    x = "Number of FE cells (coarser → finer)",
    y = expression(hat(beta)[3] * " (SWH × Post interaction)"),
    colour = "Model", shape = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(BASE_DIR, "output", "figures", "fe_sensitivity.pdf"),
       p, width = 12, height = 7)
ggsave(file.path(BASE_DIR, "output", "figures", "fe_sensitivity.png"),
       p, width = 12, height = 7, dpi = 200)
cat("Saved: output/figures/fe_sensitivity.pdf + .png\n")

# Cleaner version: faceted by model type, FE labels on x-axis
plot_dt2 <- copy(all_results)
plot_dt2[, fe_short := fe]

p2 <- ggplot(plot_dt2, aes(x = fe, y = beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.5,
                  colour = "#2166AC") +
  facet_wrap(~ model_type, scales = "free_x") +
  labs(
    title = expression("FE sensitivity: " * hat(beta)[3] * " (SWH × Post interaction)"),
    subtitle = "Daily panel: core SWH. Event level: incident-location SWH. NegBin, hetero-robust 95% CI.",
    x = "Fixed effects (coarser → finer)",
    y = expression(hat(beta)[3])
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(BASE_DIR, "output", "figures", "fe_sensitivity_faceted.pdf"),
       p2, width = 14, height = 6)
ggsave(file.path(BASE_DIR, "output", "figures", "fe_sensitivity_faceted.png"),
       p2, width = 14, height = 6, dpi = 200)
cat("Saved: output/figures/fe_sensitivity_faceted.pdf + .png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
