# Pure descriptive (model-free) statistics relating prev-week SWH to the four
# outcomes of interest, stratified by YEAR (not pre/post MoU), so the
# evolution across the institutional timeline (Mare Nostrum end 2014-11,
# NGO peak 2015-17, MoU + Minniti 2017-07, NGO crackdown 2018-19, LCG/TCG
# era 2020+) is visible directly.
#
# Outcomes:
#   (1) total deaths          (sum of n_dead_missing)
#   (2) deadly-event days     (count and Pr of n_dead_missing > 0)
#   (3) people crossing       (sum of frx_persons — Frontex-engaged persons)
#   (4) Frontex incidents     (sum of frx_incidents — Frontex events)
#
# SWH deciles are computed on the FULL 2014-2023 sample so cells are
# directly comparable across years (same physical SWH range per decile).
#
# Three views:
#   (a) year x SWH decile heatmap of share of each outcome (one panel per
#       outcome). Each year row sums to 100% within the panel.
#   (b) annual outcome-weighted SWH mean: where does each outcome's centre
#       of mass sit each year? 4 lines on one panel.
#   (c) annual probability of a deadly day, plus mean deaths per death-day
#       (a measure of how concentrated mortality is).
#
# This script reads the same panel and primary death series as 05/058. It
# does NOT fit any model. NB: 2023 is partial (data ends 2023-05-31), so
# its annual aggregates cover Jan-May only — flagged in the output.
#
# In:  analysis/data/daily_panel_complete.RDS
# Out: output/tables/11_descriptives_swh.txt
#      output/figures/11_descriptives_swh_heatmap.png      (year x decile share, 4 panels)
#      output/figures/11_descriptives_swh_centre_of_mass.png  (annual centre of mass)
#      output/figures/11_descriptives_swh_pr_deadly_day.png    (annual extensive-margin)

library(tidyverse)
library(lubridate)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")

# Institutional reference dates (used for vlines on annual plots).
INST <- tibble(
  date  = as.Date(c("2014-11-01", "2017-07-01", "2018-06-01")),
  year  = c(2014.83, 2017.5, 2018.42),
  label = c("Mare Nostrum end", "MoU + Minniti", "Salvini crackdown"),
  col   = c("#1B7837", "#D6604D", "#762A83")
)

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("059  DESCRIPTIVE: SWH vs DEATHS / EVENTS / CROSSINGS BY YEAR\n")
cat("============================================================\n\n")

# ── 1. Load panel, rebuild primary death series ────────────
panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS")) |>
  select(-n_dead_missing) |>
  left_join(build_iom_daily(), by = "date") |>
  replace_na(list(n_dead_missing = 0))

d <- panel |>
  filter(year(date) >= 2014, year(date) <= 2023, !is.na(swh_prev5days)) |>
  mutate(
    yr          = year(date),
    deadly      = as.integer(n_dead_missing > 0),
    swh_decile  = ntile(swh_prev5days, 10)   # full-sample deciles
  )

# Decile boundaries (for the report)
decile_brk <- d |>
  group_by(swh_decile) |>
  summarise(swh_min = min(swh_prev5days),
            swh_max = max(swh_prev5days),
            swh_mid = round(median(swh_prev5days), 3),
            .groups = "drop")
cat("Decile boundaries (m):\n")
print(decile_brk)

# ── 2. Year-level summary table ────────────────────────────
cat("\n--- 1. Year-level summary ---\n")
year_tab <- d |>
  group_by(yr) |>
  summarise(
    n_days          = n(),
    n_death_days    = sum(deadly),
    pr_death_day    = round(n_death_days / n_days, 3),
    total_deaths    = sum(n_dead_missing),
    deaths_per_day  = round(total_deaths / n_days, 2),
    deaths_per_ddday = round(total_deaths / pmax(n_death_days, 1), 1),
    total_frx_pers  = sum(frx_persons,   na.rm = TRUE),
    total_frx_inc   = sum(frx_incidents, na.rm = TRUE),
    .groups = "drop"
  )
print(year_tab, n = Inf, width = Inf)

# ── 3. Year × decile cross-tab ─────────────────────────────
cat("\n--- 2. Year x SWH decile cross-tab (long format) ---\n")
yd_tab <- d |>
  group_by(yr, swh_decile) |>
  summarise(
    n_days          = n(),
    n_death_days    = sum(deadly),
    pr_death_day    = round(n_death_days / n_days, 3),
    total_deaths    = sum(n_dead_missing),
    total_frx_pers  = sum(frx_persons,   na.rm = TRUE),
    total_frx_inc   = sum(frx_incidents, na.rm = TRUE),
    .groups = "drop"
  )
# Add within-year shares for each outcome
yd_share <- yd_tab |>
  group_by(yr) |>
  mutate(
    share_days   = n_days         / sum(n_days),
    share_deaths = total_deaths   / pmax(sum(total_deaths),   1),
    share_ddays  = n_death_days   / pmax(sum(n_death_days),   1),
    share_pers   = total_frx_pers / pmax(sum(total_frx_pers), 1),
    share_inc    = total_frx_inc  / pmax(sum(total_frx_inc),  1)
  ) |>
  ungroup()

# Print one year as a sanity check
cat("\nExample (2017):\n")
print(yd_share |> filter(yr == 2017) |>
        select(swh_decile, n_days, total_deaths,
               share_deaths, share_pers, share_inc),
      n = Inf)

# ── 4. Outcome-weighted centre-of-mass by year ─────────────
wmean <- function(x, w) {
  w <- pmax(w, 0)
  if (sum(w) == 0) NA_real_ else sum(x * w) / sum(w)
}

centre_yr <- d |>
  group_by(yr) |>
  summarise(
    n_days                  = n(),
    `unweighted mean`       = mean(swh_prev5days),
    `deaths-weighted`       = wmean(swh_prev5days, n_dead_missing),
    `death-day SWH mean`    = if (sum(deadly) > 0) mean(swh_prev5days[deadly == 1]) else NA_real_,
    `frx-pers-weighted`     = wmean(swh_prev5days, frx_persons),
    `frx-inc-weighted`      = wmean(swh_prev5days, frx_incidents),
    .groups = "drop"
  )
cat("\n--- 3. SWH centre-of-mass by year ---\n")
print(centre_yr |>
        mutate(across(where(is.numeric), \(x) round(x, 3))),
      width = Inf)

# Long format for the centre-of-mass plot
centre_long <- centre_yr |>
  pivot_longer(c(`unweighted mean`, `deaths-weighted`, `death-day SWH mean`,
                 `frx-pers-weighted`, `frx-inc-weighted`),
               names_to = "weight", values_to = "swh_m") |>
  mutate(weight = factor(weight, levels = c(
    "unweighted mean", "death-day SWH mean", "deaths-weighted",
    "frx-pers-weighted", "frx-inc-weighted")))

# ── 5. Plot 1: heatmap of decile shares by year, faceted by outcome ─
cat("\n--- 4. Heatmap figure ---\n")

heat_long <- yd_share |>
  select(yr, swh_decile,
         `(1) deaths`            = share_deaths,
         `(2) death-days`        = share_ddays,
         `(3) Frontex persons`   = share_pers,
         `(4) Frontex incidents` = share_inc) |>
  pivot_longer(c(-yr, -swh_decile), names_to = "outcome", values_to = "share") |>
  mutate(outcome = factor(outcome, levels = c(
    "(1) deaths", "(2) death-days",
    "(3) Frontex persons", "(4) Frontex incidents")))

p_heat <- ggplot(heat_long, aes(x = swh_decile, y = factor(yr), fill = share)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = scales::percent(share, accuracy = 1)),
            size = 2.6, colour = "grey15") +
  scale_x_continuous(breaks = 1:10, labels = function(x) paste0("D", x),
                     expand = c(0, 0)) +
  scale_y_discrete(limits = rev) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       labels = scales::percent_format(accuracy = 1),
                       name = "Share within year") +
  facet_wrap(~ outcome, ncol = 2) +
  labs(
    title    = "Share of each outcome by SWH decile, by year",
    subtitle = paste("SWH deciles computed on the full 2014-2023 sample.",
                     "Each row sums to 100% within each panel.",
                     "2023 is partial (Jan-May)."),
    x = "Prev-week SWH decile (D1 = calmest, D10 = roughest)",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid       = element_blank(),
        legend.position  = "right",
        strip.text       = element_text(face = "bold"))

ggsave(fig_path("04_descriptive", "02_swh_heatmap.png"),
       p_heat, width = 11, height = 9, dpi = 200)

# ── 6. Plot 2: annual centre-of-mass time series ───────────
weight_cols <- c("unweighted mean"     = "#999999",
                 "death-day SWH mean"  = "#000000",
                 "deaths-weighted"     = "#B2182B",
                 "frx-pers-weighted"   = "#2166AC",
                 "frx-inc-weighted"    = "#4393C3")

p_com <- centre_long |>
  ggplot(aes(yr, swh_m, colour = weight, group = weight)) +
  geom_vline(data = INST, aes(xintercept = year),
             linetype = "dotted", colour = "grey50", linewidth = 0.5) +
  geom_text(data = INST,
            aes(x = year, y = 1.25, label = label),
            inherit.aes = FALSE,
            angle = 90, hjust = 1, vjust = -0.3, size = 2.8, colour = "grey50") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = weight_cols) +
  scale_x_continuous(breaks = 2014:2023) +
  labs(
    title    = "Annual SWH centre of mass, by outcome",
    subtitle = paste("Each point is the SWH value at which the outcome's",
                     "mass is centred in that year. 2023 is partial."),
    x = NULL, y = "Prev-week SWH (m)", colour = "Weighting"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "top")

ggsave(fig_path("04_descriptive", "02_swh_centre_of_mass.png"),
       p_com, width = 10, height = 6, dpi = 200)

# ── 7. Plot 3: annual probability of a deadly day + intensity ─
# Annual Pr(deadly day) and annual deaths-per-deadly-day (mass-casualty
# intensity). Two y-axes on one plot via geom_line with separate scales,
# rendered as two stacked panels for clarity.

intensity_df <- year_tab |>
  select(yr, pr_death_day, deaths_per_ddday)

p_pr <- intensity_df |>
  ggplot(aes(yr, pr_death_day)) +
  geom_vline(data = INST, aes(xintercept = year),
             linetype = "dotted", colour = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 0.7, colour = "#2166AC") +
  geom_point(size = 2.5, colour = "#2166AC") +
  scale_x_continuous(breaks = 2014:2023) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Annual Pr(deadly day): share of days with at least one death",
       x = NULL, y = "Pr(deadly day)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_int <- intensity_df |>
  ggplot(aes(yr, deaths_per_ddday)) +
  geom_vline(data = INST, aes(xintercept = year),
             linetype = "dotted", colour = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 0.7, colour = "#B2182B") +
  geom_point(size = 2.5, colour = "#B2182B") +
  scale_x_continuous(breaks = 2014:2023) +
  labs(title = "Annual mean deaths per deadly day (mass-casualty intensity)",
       x = NULL, y = "deaths / death-day") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

library(patchwork)
p_combined <- p_pr / p_int
ggsave(fig_path("04_descriptive", "02_swh_pr_deadly_day.png"),
       p_combined, width = 10, height = 7, dpi = 200)

# ── 8. Save text output ────────────────────────────────────
sink_file <- tbl_path("04_descriptive", "02_swh.txt")
sink(sink_file)
old_opts <- options(tibble.width = Inf, tibble.print_max = Inf)
on.exit(options(old_opts), add = TRUE)

cat("059  SWH vs OUTCOMES — annual descriptive (no model)\n")
cat("=====================================================\n")
cat(sprintf("Sample: 2014-01-01 to 2023-05-31  N=%d days\n", nrow(d)))
cat("NB: 2023 is partial (Jan-May only).\n\n")

cat("=== Decile boundaries (m) ===\n")
print(decile_brk)
cat("\n")

cat("=== Annual summary ===\n")
print(year_tab, n = Inf)
cat("\n")

cat("=== SWH centre-of-mass by year ===\n")
print(centre_yr |>
        mutate(across(where(is.numeric), \(x) round(x, 3))))
cat("\nRead: 'deaths-weighted' is sum(SWH * deaths) / sum(deaths) for that year.\n")
cat("If it rises across years, the centre of the deaths distribution\n")
cat("is shifting toward rougher weather.\n\n")

cat("=== Year x decile shares (long format) ===\n")
print(yd_share |>
        select(yr, swh_decile, n_days, total_deaths, share_deaths,
               share_ddays, share_pers, share_inc) |>
        mutate(across(starts_with("share"), \(x) round(x, 3))),
      n = Inf)

sink()
cat(sprintf("\nSaved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
