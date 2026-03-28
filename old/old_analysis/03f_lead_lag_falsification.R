# 03f_lead_lag_falsification.R
# ============================
# Lead-lag falsification + 6 main lag specifications.
#
# Part 1: Lead-lag coefficient path (Camarena et al. style)
#   Estimate n_dead_missing ~ FE + sum_{k=-3}^{+3} beta_k * SWH_{t-k}
#   Leads (future weather) should be null; lags 1-2 significant, 3+ decaying.
#   Two versions: baseline (no interaction) and with SWH x Post interaction.
#
# Part 2: Six main lag specifications (each a separate model)
#   lag0, lag1, lag2, lag3, avg(0-3), avg(1-3)
#   Each: n_dead_missing ~ SWH_var + SWH_var:post_mou | FE
#   Both week-year and month-year FE; Newey-West SEs.
#
# Input:  data/processed/cmr_daily_weather_panel.RDS
# Output: output/figures/lead_lag_coefpath.pdf / .png
#         output/figures/lag_spec_comparison.pdf / .png
#         output/tables/lead_lag_results.csv
#         output/tables/lag_spec_results.csv

library(fixest)
library(data.table)
library(ggplot2)

BASE_DIR <- here::here()
d <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                      "cmr_daily_weather_panel.RDS")))

MOU_DATE <- as.Date("2017-07-01")

# Recompute post_mou with correct date
d[, post_mou := as.integer(date >= MOU_DATE)]

cat("============================================================\n")
cat("LEAD-LAG FALSIFICATION + LAG SPECIFICATIONS\n")
cat("============================================================\n\n")
cat(sprintf("Panel: %d days, MoU date: %s\n", nrow(d), MOU_DATE))
cat(sprintf("Pre: %d days, Post: %d days\n", sum(d$post_mou == 0), sum(d$post_mou == 1)))
cat(sprintf("Outcome mean: %.3f, non-zero days: %d (%.1f%%)\n\n",
    mean(d$n_dead_missing), sum(d$n_dead_missing > 0),
    100 * mean(d$n_dead_missing > 0)))

# ============================================================
# 0. Create lead and average variables
# ============================================================
d <- d[order(date)]

# Leads (future weather)
d[, swh_core_lead1 := shift(swh_core, n = 1, type = "lead")]
d[, swh_core_lead2 := shift(swh_core, n = 2, type = "lead")]
d[, swh_core_lead3 := shift(swh_core, n = 3, type = "lead")]

# Average lag 0-3
d[, swh_core_avg03 := rowMeans(.SD, na.rm = FALSE),
  .SDcols = c("swh_core", "swh_core_lag1", "swh_core_lag2", "swh_core_lag3")]

cat("Variables created:\n")
cat(sprintf("  Leads: %d/%d/%d non-NA (lead1/2/3)\n",
    sum(!is.na(d$swh_core_lead1)), sum(!is.na(d$swh_core_lead2)),
    sum(!is.na(d$swh_core_lead3))))
cat(sprintf("  Avg 0-3: %d non-NA\n", sum(!is.na(d$swh_core_avg03))))
cat(sprintf("  Prev3d (avg 1-3): %d non-NA\n\n", sum(!is.na(d$swh_core_prev3d))))

# ============================================================
# Part 1: LEAD-LAG FALSIFICATION
# ============================================================
cat("============================================================\n")
cat("PART 1: LEAD-LAG COEFFICIENT PATH\n")
cat("============================================================\n\n")

# Subset to complete cases for lead-lag model
d_ll <- d[!is.na(swh_core_lead3) & !is.na(swh_core_lag3)]
cat(sprintf("Lead-lag sample: %d days (dropped %d for leads/lags)\n\n",
    nrow(d_ll), nrow(d) - nrow(d_ll)))

# --- Model A: Baseline (no interaction) ---
cat("--- Model A: Baseline lead-lag path ---\n")
m_base <- fenegbin(
  n_dead_missing ~ swh_core_lead3 + swh_core_lead2 + swh_core_lead1 +
    swh_core + swh_core_lag1 + swh_core_lag2 + swh_core_lag3 | week_year_fac,
  data = d_ll, vcov = NW ~ date
)
cat("Converged. Coefficients:\n")
ct_base <- summary(m_base)$coeftable
print(round(ct_base[1:7, ], 4))

# --- Model B: With Post interaction ---
cat("\n--- Model B: Lead-lag path x Post ---\n")
m_inter <- fenegbin(
  n_dead_missing ~ swh_core_lead3 + swh_core_lead2 + swh_core_lead1 +
    swh_core + swh_core_lag1 + swh_core_lag2 + swh_core_lag3 +
    swh_core_lead3:post_mou + swh_core_lead2:post_mou + swh_core_lead1:post_mou +
    swh_core:post_mou + swh_core_lag1:post_mou + swh_core_lag2:post_mou +
    swh_core_lag3:post_mou | week_year_fac,
  data = d_ll, vcov = NW ~ date
)
cat("Converged. Baseline coefficients:\n")
ct_inter <- summary(m_inter)$coeftable
print(round(ct_inter[1:7, ], 4))
cat("\nInteraction coefficients:\n")
print(round(ct_inter[8:14, ], 4))

# --- Extract for plotting ---
extract_path <- function(ct, rows, label) {
  k <- c(3, 2, 1, 0, -1, -2, -3)  # lead3..lag3 -> k=+3..−3
  data.table(
    k = k[seq_along(rows)],
    beta = ct[rows, 1],
    se = ct[rows, 2],
    p = ct[rows, 4],
    model = label
  )
}

path_base <- extract_path(ct_base, 1:7, "Baseline")
path_base[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]

path_inter_main <- extract_path(ct_inter, 1:7, "Baseline (with interaction model)")
path_inter_main[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]

path_inter_int <- extract_path(ct_inter, 8:14, "Interaction (SWH x Post)")
path_inter_int[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]

# --- Plot ---
plot_base <- rbind(path_base)
plot_inter <- rbind(path_inter_main, path_inter_int)

p_base <- ggplot(plot_base, aes(x = k, y = beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey70") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "steelblue") +
  geom_line(colour = "steelblue", linewidth = 0.6) +
  geom_point(colour = "steelblue", size = 3) +
  scale_x_continuous(breaks = -3:3,
    labels = c("Lead 3", "Lead 2", "Lead 1", "Day 0",
               "Lag 1", "Lag 2", "Lag 3")) +
  annotate("text", x = 1.5, y = max(plot_base$ci_hi) * 0.9,
           label = "Leads (future)\nshould be null", size = 3, colour = "grey40") +
  annotate("text", x = -1.5, y = max(plot_base$ci_hi) * 0.9,
           label = "Lags (past)\nshould matter", size = 3, colour = "grey40") +
  labs(
    title = expression("Lead-lag coefficient path: " * hat(beta)[k] *
                        " for SWH at t-k"),
    subtitle = "NegBin, week-year FE, Newey-West SEs. Baseline model (no Post interaction).",
    x = NULL, y = expression(hat(beta)[k])
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_inter <- ggplot(plot_inter, aes(x = k, y = beta, colour = model, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey70") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = -3:3,
    labels = c("Lead 3", "Lead 2", "Lead 1", "Day 0",
               "Lag 1", "Lag 2", "Lag 3")) +
  scale_colour_manual(values = c("Baseline (with interaction model)" = "steelblue",
                                  "Interaction (SWH x Post)" = "firebrick")) +
  scale_fill_manual(values = c("Baseline (with interaction model)" = "steelblue",
                                "Interaction (SWH x Post)" = "firebrick")) +
  labs(
    title = expression("Lead-lag coefficient path with " * hat(gamma)[k] *
                        " (SWH x Post interaction)"),
    subtitle = "NegBin, week-year FE, Newey-West SEs.",
    x = NULL, y = "Coefficient",
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

library(patchwork)
p_ll <- p_base / p_inter

ggsave(file.path(BASE_DIR, "output", "figures", "lead_lag_coefpath.pdf"),
       p_ll, width = 10, height = 10)
ggsave(file.path(BASE_DIR, "output", "figures", "lead_lag_coefpath.png"),
       p_ll, width = 10, height = 10, dpi = 200)
cat("\nSaved: output/figures/lead_lag_coefpath.pdf + .png\n")

# Save lead-lag table
ll_all <- rbind(path_base, path_inter_main, path_inter_int)
fwrite(ll_all, file.path(BASE_DIR, "output", "tables", "lead_lag_results.csv"))
cat("Saved: output/tables/lead_lag_results.csv\n")

# ============================================================
# Part 2: SIX LAG SPECIFICATIONS
# ============================================================
cat("\n============================================================\n")
cat("PART 2: LAG SPECIFICATION COMPARISON\n")
cat("============================================================\n\n")

spec_list <- list(
  "Lag 0"     = "swh_core",
  "Lag 1"     = "swh_core_lag1",
  "Lag 2"     = "swh_core_lag2",
  "Lag 3"     = "swh_core_lag3",
  "Avg 0-3"   = "swh_core_avg03",
  "Avg 1-3"   = "swh_core_prev3d"
)

fe_list <- list(
  "Week-year FE"  = "week_year_fac",
  "Month-year FE" = "month_year_fac"
)

results <- rbindlist(lapply(names(spec_list), function(spec_name) {
  swh_var <- spec_list[[spec_name]]

  rbindlist(lapply(names(fe_list), function(fe_name) {
    fe_var <- fe_list[[fe_name]]
    fml <- as.formula(paste0(
      "n_dead_missing ~ ", swh_var, " + ", swh_var, ":post_mou | ", fe_var))

    m <- tryCatch(
      fenegbin(fml, data = d[!is.na(get(swh_var))], vcov = NW ~ date),
      error = function(e) NULL
    )
    if (is.null(m)) return(NULL)

    ct <- summary(m)$coeftable
    # Main effect row
    main_row <- which(rownames(ct) == swh_var)
    # Interaction row
    int_row <- grep(paste0(swh_var, ":post_mou|post_mou:", swh_var), rownames(ct))

    if (length(main_row) == 0 || length(int_row) == 0) return(NULL)

    data.table(
      spec      = spec_name,
      fe        = fe_name,
      beta1     = ct[main_row, 1],
      beta1_se  = ct[main_row, 2],
      beta1_p   = ct[main_row, 4],
      beta3     = ct[int_row, 1],
      beta3_se  = ct[int_row, 2],
      beta3_p   = ct[int_row, 4],
      beta3_irr = exp(ct[int_row, 1]),
      n_obs     = nobs(m)
    )
  }))
}))

results[, `:=`(
  beta3_ci_lo = beta3 - 1.96 * beta3_se,
  beta3_ci_hi = beta3 + 1.96 * beta3_se
)]

# --- Print ---
for (fe_label in names(fe_list)) {
  cat(sprintf("\n%s:\n", fe_label))
  cat(sprintf("  %-10s %+8s %7s %7s  |  %+8s %7s %7s %8s\n",
      "Spec", "beta1", "SE", "p", "beta3", "SE", "p", "IRR"))
  cat(paste0("  ", paste(rep("-", 80), collapse = "")), "\n")
  sub <- results[fe == fe_label]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i]
    stars <- ifelse(r$beta3_p < 0.01, "***",
             ifelse(r$beta3_p < 0.05, "**",
             ifelse(r$beta3_p < 0.1, "*", "")))
    cat(sprintf("  %-10s %+8.4f %7.4f %7.4f  |  %+8.4f %7.4f %7.4f %8.3f %s\n",
        r$spec, r$beta1, r$beta1_se, r$beta1_p,
        r$beta3, r$beta3_se, r$beta3_p, r$beta3_irr, stars))
  }
}

# --- Save ---
fwrite(results, file.path(BASE_DIR, "output", "tables", "lag_spec_results.csv"))
cat("\nSaved: output/tables/lag_spec_results.csv\n")

# --- Plot: beta3 comparison ---
results[, spec_f := factor(spec, levels = rev(c("Lag 0", "Lag 1", "Lag 2",
                                                 "Lag 3", "Avg 0-3", "Avg 1-3")))]

p_specs <- ggplot(results, aes(y = spec_f, x = beta3, colour = fe)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = beta3_ci_lo, xmax = beta3_ci_hi),
                  position = position_dodge(width = 0.4), size = 0.5) +
  scale_colour_manual(values = c("Week-year FE" = "steelblue",
                                  "Month-year FE" = "darkorange")) +
  labs(
    title = expression(hat(beta)[3] * " (SWH x Post) across lag specifications"),
    subtitle = "NegBin, Newey-West SEs. Daily panel, full sea zone.",
    x = expression(hat(beta)[3]),
    y = NULL, colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(BASE_DIR, "output", "figures", "lag_spec_comparison.pdf"),
       p_specs, width = 9, height = 5)
ggsave(file.path(BASE_DIR, "output", "figures", "lag_spec_comparison.png"),
       p_specs, width = 9, height = 5, dpi = 200)
cat("Saved: output/figures/lag_spec_comparison.pdf + .png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
