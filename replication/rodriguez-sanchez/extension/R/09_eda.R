# 09_eda.R
# =======
# Exploratory data analysis for the mortality extension.
#
# Part 1: Time trends — mortality rate, death counts, crossing attempts
#         for CMR and other Mediterranean routes (EMR, WMR).
# Part 2: Correlations — three CMR outcomes vs. 19 exogenous predictors.
#
# Outputs: output/figures/eda/
# Requires: targets cache (run tar_make() first).

library(targets)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(lubridate)
library(sf)
library(ggspatial)

OUT_DIR <- "output/figures/eda"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

theme_eda <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40")
  )


# =========================================================================
# 0. Load data
# =========================================================================

df_model <- tar_read(df_model)
df_full  <- df_model$df_full

# CMR outcomes (raw + log)
df_cmr <- df_full %>%
  filter(date >= as.Date("2011-02-01") & date <= as.Date("2021-09-01")) %>%
  mutate(
    crossings     = crossings_CMR,
    deaths        = dead_and_missing_Central_Mediterranean,
    rate_100      = mortality_rate_100,
    log_crossings = log(crossings_CMR),
    log_deaths    = log(deaths + 1),
    log_rate      = log(rate_100 + 0.01),
    period = factor(case_when(
      date < as.Date("2013-10-01") ~ "Pre-MN",
      date < as.Date("2014-11-01") ~ "Mare Nostrum",
      date < as.Date("2017-02-01") ~ "NGO SAR",
      TRUE                         ~ "Post-MoU"
    ), levels = c("Pre-MN", "Mare Nostrum", "NGO SAR", "Post-MoU"))
  )

# Other routes from original dataset
df_orig <- readRDS("data/raw/df.RDS") %>%
  filter(date >= as.Date("2011-02-01") & date <= as.Date("2021-09-01"))

df_routes <- df_orig %>%
  transmute(
    date,
    deaths_CMR   = ifelse(is.na(dead_and_missing_Central_Mediterranean), 0,
                          dead_and_missing_Central_Mediterranean),
    arrivals_CMR = arrivals_CMR,
    deaths_EMR   = ifelse(is.na(dead_and_missing_Eastern_Mediterranean), 0,
                          dead_and_missing_Eastern_Mediterranean),
    arrivals_EMR = arrivals_EMR,
    deaths_WMR   = ifelse(is.na(dead_and_missing_Western_Mediterranean), 0,
                          dead_and_missing_Western_Mediterranean),
    arrivals_WMR = arrivals_WMR
  ) %>%
  mutate(
    rate_CMR = ifelse(arrivals_CMR > 0, deaths_CMR / arrivals_CMR * 100, NA),
    rate_EMR = ifelse(arrivals_EMR > 0, deaths_EMR / arrivals_EMR * 100, NA),
    rate_WMR = ifelse(arrivals_WMR > 0, deaths_WMR / arrivals_WMR * 100, NA)
  )

# Intervention dates
interventions <- data.frame(
  date  = as.Date(c("2013-10-01", "2014-11-01", "2017-02-01")),
  label = c("Mare Nostrum", "NGO SAR", "EU-Libya MoU")
)

# 19 base exogenous predictors (lag 0)
base_vars <- c(
  "wave_height_central_med", "wave_period_central_med",
  "wave_direction_central_med", "wind_speed_central_med",
  "sst_central_med", "sst_anomaly_central_med",
  "cloud_cover_central_med", "low_cloud_central_med",
  "dewpoint_depression_central_med", "temperature_central_med",
  "wind_speed_departure_coast", "cloud_cover_departure_coast",
  "temperature_departure_coast", "precipitation_departure_coast",
  "wave_max_central_med", "wave_sd_central_med",
  "wave_days_above_2m",
  "current_speed_central_med", "current_against_route"
)

var_labels <- c(
  wave_height_central_med       = "Wave height",
  wave_period_central_med       = "Wave period",
  wave_direction_central_med    = "Wave direction",
  wind_speed_central_med        = "Wind speed",
  sst_central_med               = "SST",
  sst_anomaly_central_med       = "SST anomaly",
  cloud_cover_central_med       = "Cloud cover",
  low_cloud_central_med         = "Low cloud",
  dewpoint_depression_central_med = "Dewpoint dep.",
  temperature_central_med       = "Air temp.",
  wind_speed_departure_coast    = "Wind (coast)",
  cloud_cover_departure_coast   = "Cloud (coast)",
  temperature_departure_coast   = "Temp. (coast)",
  precipitation_departure_coast = "Precip. (coast)",
  wave_max_central_med          = "Wave max",
  wave_sd_central_med           = "Wave SD",
  wave_days_above_2m            = "Days >2m waves",
  current_speed_central_med     = "Current speed",
  current_against_route         = "Current (opposing)"
)

period_colors <- c(
  "Pre-MN"       = "#E69F00",
  "Mare Nostrum"  = "#56B4E9",
  "NGO SAR"      = "#009E73",
  "Post-MoU"     = "#D55E00"
)

route_colors <- c(
  "Central Mediterranean" = "firebrick",
  "Eastern Mediterranean" = "steelblue",
  "Western Mediterranean" = "forestgreen"
)


# =========================================================================
# PART 1: TIME TRENDS
# =========================================================================

cat("=== PART 1: Time trends ===\n")

# --- 1a. CMR: raw scale (3 panels) ---

make_cmr_panel <- function(df, yvar, ylab, title, color) {
  ggplot(df, aes(date, .data[[yvar]])) +
    geom_line(color = color, linewidth = 0.6) +
    geom_vline(data = interventions, aes(xintercept = date),
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_text(data = interventions, aes(x = date, y = Inf, label = label),
              vjust = 1.5, hjust = 0.5, size = 2.5, color = "grey40") +
    labs(x = NULL, y = ylab, title = title) +
    theme_eda
}

p1a <- (
  make_cmr_panel(df_cmr, "crossings", "Crossing attempts",
                 "Crossing attempts", "steelblue") +
    scale_y_continuous(labels = scales::comma)
) /
  make_cmr_panel(df_cmr, "deaths", "Deaths & missing",
                 "Deaths & missing", "firebrick") /
  make_cmr_panel(df_cmr, "rate_100", "Rate (%)",
                 "Mortality rate (deaths / crossings)", "darkorange") +
  plot_annotation(
    title    = "Central Mediterranean Route",
    subtitle = "Feb 2011 - Sep 2021. Dashed lines = intervention dates.",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(OUT_DIR, "01_cmr_time_trends.png"), p1a,
       width = 10, height = 10, dpi = 200, bg = "white")
cat("  Saved 01_cmr_time_trends.png\n")


# --- 1b. CMR: log scale (as modeled) ---

p1b <- (
  make_cmr_panel(df_cmr, "log_crossings", "log(crossings)",
                 "log(crossings)", "steelblue") /
  make_cmr_panel(df_cmr, "log_deaths", "log(deaths + 1)",
                 "log(deaths + 1)", "firebrick") /
  make_cmr_panel(df_cmr, "log_rate", "log(rate + 0.01)",
                 "log(mortality rate + 0.01)", "darkorange")
) +
  plot_annotation(
    title    = "Central Mediterranean Route (log scale, as modeled)",
    subtitle = "These are the actual outcome variables entering CausalImpact.",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(OUT_DIR, "02_cmr_time_trends_log.png"), p1b,
       width = 10, height = 10, dpi = 200, bg = "white")
cat("  Saved 02_cmr_time_trends_log.png\n")


# --- 1c. Cross-route: deaths ---

df_deaths_long <- df_routes %>%
  select(date, deaths_CMR, deaths_EMR, deaths_WMR) %>%
  pivot_longer(-date, names_to = "route", values_to = "deaths") %>%
  mutate(route = recode(route,
    deaths_CMR = "Central Mediterranean",
    deaths_EMR = "Eastern Mediterranean",
    deaths_WMR = "Western Mediterranean"
  ))

p1c <- ggplot(df_deaths_long, aes(date, deaths, color = route)) +
  geom_line(linewidth = 0.6) +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = route_colors) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Deaths & missing", color = "Route",
       title = "Monthly Deaths by Mediterranean Route",
       subtitle = "Feb 2011 - Sep 2021") +
  theme_eda

ggsave(file.path(OUT_DIR, "03_deaths_by_route.png"), p1c,
       width = 10, height = 5, dpi = 200, bg = "white")
cat("  Saved 03_deaths_by_route.png\n")


# --- 1d. Cross-route: arrivals ---

df_arrivals_long <- df_routes %>%
  select(date, arrivals_CMR, arrivals_EMR, arrivals_WMR) %>%
  pivot_longer(-date, names_to = "route", values_to = "arrivals") %>%
  mutate(route = recode(route,
    arrivals_CMR = "Central Mediterranean",
    arrivals_EMR = "Eastern Mediterranean",
    arrivals_WMR = "Western Mediterranean"
  ))

p1d <- ggplot(df_arrivals_long, aes(date, arrivals, color = route)) +
  geom_line(linewidth = 0.6) +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = route_colors) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Arrivals", color = "Route",
       title = "Monthly Arrivals by Mediterranean Route",
       subtitle = "Arrivals only (no pushback data for EMR/WMR).") +
  theme_eda

ggsave(file.path(OUT_DIR, "04_arrivals_by_route.png"), p1d,
       width = 10, height = 5, dpi = 200, bg = "white")
cat("  Saved 04_arrivals_by_route.png\n")


# --- 1e. Cross-route: mortality rate ---

df_rate_long <- df_routes %>%
  select(date, rate_CMR, rate_EMR, rate_WMR) %>%
  pivot_longer(-date, names_to = "route", values_to = "rate") %>%
  mutate(route = recode(route,
    rate_CMR = "Central Mediterranean",
    rate_EMR = "Eastern Mediterranean",
    rate_WMR = "Western Mediterranean"
  )) %>%
  filter(!is.na(rate))

p1e <- ggplot(df_rate_long, aes(date, rate, color = route)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = route_colors) +
  labs(x = NULL, y = "Mortality rate (%)", color = "Route",
       title = "Monthly Mortality Rate by Mediterranean Route",
       subtitle = "Deaths / arrivals (approx. for EMR/WMR, no pushback data).") +
  theme_eda

ggsave(file.path(OUT_DIR, "05_mortality_rate_by_route.png"), p1e,
       width = 10, height = 5, dpi = 200, bg = "white")
cat("  Saved 05_mortality_rate_by_route.png\n")


# --- 1f. Deaths vs crossings scatter (colored by period) ---

p1f <- ggplot(df_cmr, aes(crossings, deaths, color = period)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed",
              linewidth = 0.5, aes(group = 1), color = "grey40") +
  scale_color_manual(values = period_colors) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Crossing attempts", y = "Deaths & missing", color = "Period",
       title = "Deaths vs. Crossings (CMR)",
       subtitle = sprintf("r = %.2f (overall)",
                          cor(df_cmr$crossings, df_cmr$deaths, use = "complete.obs"))) +
  theme_eda

ggsave(file.path(OUT_DIR, "06_deaths_vs_crossings.png"), p1f,
       width = 8, height = 6, dpi = 200, bg = "white")
cat("  Saved 06_deaths_vs_crossings.png\n")


# =========================================================================
# PART 2: CORRELATIONS
# =========================================================================

cat("\n=== PART 2: Correlations with exogenous predictors ===\n")

available_vars <- base_vars[base_vars %in% names(df_cmr)]
if (length(available_vars) < length(base_vars)) {
  cat("  Missing:", paste(setdiff(base_vars, available_vars), collapse = ", "), "\n")
}

# --- 2a. Compute correlations ---

cor_results <- data.frame(variable = available_vars, label = var_labels[available_vars],
                          stringsAsFactors = FALSE)

for (v in available_vars) {
  x <- df_cmr[[v]]
  ct_c <- cor.test(x, df_cmr$log_crossings, use = "complete.obs")
  ct_d <- cor.test(x, df_cmr$log_deaths,    use = "complete.obs")
  ct_r <- cor.test(x, df_cmr$log_rate,      use = "complete.obs")

  cor_results[cor_results$variable == v, "r_crossings"] <- ct_c$estimate
  cor_results[cor_results$variable == v, "p_crossings"] <- ct_c$p.value
  cor_results[cor_results$variable == v, "r_deaths"]    <- ct_d$estimate
  cor_results[cor_results$variable == v, "p_deaths"]    <- ct_d$p.value
  cor_results[cor_results$variable == v, "r_rate"]      <- ct_r$estimate
  cor_results[cor_results$variable == v, "p_rate"]      <- ct_r$p.value
}

# Print summary
cat(sprintf("\n  N = %d months\n\n", nrow(df_cmr)))
star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "   ")))

cat(sprintf("  %-20s  %8s  %8s  %8s\n", "Variable", "r(cross)", "r(death)", "r(rate)"))
cat(paste0("  ", paste(rep("-", 60), collapse = ""), "\n"))
for (i in seq_len(nrow(cor_results))) {
  r <- cor_results[i, ]
  cat(sprintf("  %-20s  %+.3f%s  %+.3f%s  %+.3f%s\n",
              r$label, r$r_crossings, star(r$p_crossings),
              r$r_deaths, star(r$p_deaths), r$r_rate, star(r$p_rate)))
}
cat(sprintf("\n  Sig. with crossings: %d / %d\n",
            sum(cor_results$p_crossings < 0.05), nrow(cor_results)))
cat(sprintf("  Sig. with deaths:    %d / %d\n",
            sum(cor_results$p_deaths < 0.05), nrow(cor_results)))
cat(sprintf("  Sig. with rate:      %d / %d\n",
            sum(cor_results$p_rate < 0.05), nrow(cor_results)))


# --- 2b. Dot plot ---

label_order <- cor_results %>% arrange(r_crossings) %>% pull(label)

cor_long <- cor_results %>%
  select(label, r_crossings, r_deaths, r_rate) %>%
  pivot_longer(-label, names_to = "outcome", values_to = "r") %>%
  mutate(
    outcome = factor(recode(outcome,
      r_crossings = "log(crossings)",
      r_deaths    = "log(deaths + 1)",
      r_rate      = "log(rate + 0.01)"),
      levels = c("log(crossings)", "log(deaths + 1)", "log(rate + 0.01)")),
    label = factor(label, levels = label_order)
  )

p2b <- ggplot(cor_long, aes(r, label, color = outcome)) +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c(
    "log(crossings)"   = "steelblue",
    "log(deaths + 1)"  = "firebrick",
    "log(rate + 0.01)" = "darkorange"
  )) +
  labs(x = "Pearson r", y = NULL, color = "Outcome",
       title = "Exogenous Predictors: Correlations with CMR Outcomes",
       subtitle = "Lag 0 only. Blue = crossings, red = deaths, orange = rate.") +
  theme_eda +
  theme(axis.text.y = element_text(size = 8))

ggsave(file.path(OUT_DIR, "07_correlation_dotplot.png"), p2b,
       width = 9, height = 7, dpi = 200, bg = "white")
cat("  Saved 07_correlation_dotplot.png\n")


# --- 2c. Heatmap ---

cor_heat <- cor_results %>%
  select(label, r_crossings, r_deaths, r_rate) %>%
  pivot_longer(-label, names_to = "outcome", values_to = "r") %>%
  mutate(
    outcome = factor(recode(outcome,
      r_crossings = "Crossings", r_deaths = "Deaths", r_rate = "Rate"),
      levels = c("Crossings", "Deaths", "Rate")),
    label = factor(label, levels = label_order)
  )

p2c <- ggplot(cor_heat, aes(outcome, label, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.8) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-0.6, 0.6)) +
  labs(x = NULL, y = NULL, fill = "r",
       title = "Correlation Heatmap: Predictors vs. CMR Outcomes",
       subtitle = "Deaths track crossings (volume channel). Rate is near zero.") +
  theme_eda +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10, face = "bold"),
        panel.grid  = element_blank())

ggsave(file.path(OUT_DIR, "08_correlation_heatmap.png"), p2c,
       width = 7, height = 8, dpi = 200, bg = "white")
cat("  Saved 08_correlation_heatmap.png\n")


# --- 2d. Scatter: top 6 predictors vs crossings and vs rate ---

top6 <- cor_results %>% arrange(desc(abs(r_crossings))) %>% slice_head(n = 6)

scatter_long <- df_cmr %>%
  select(all_of(c("log_crossings", "log_rate", "period", top6$variable))) %>%
  pivot_longer(cols = all_of(top6$variable),
               names_to = "variable", values_to = "x") %>%
  mutate(label = factor(var_labels[variable], levels = top6$label))

p2d_cross <- ggplot(scatter_long, aes(x, log_crossings, color = period)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30",
              linewidth = 0.5, aes(group = 1)) +
  facet_wrap(~ label, scales = "free_x", ncol = 3) +
  scale_color_manual(values = period_colors) +
  labs(x = "Predictor value", y = "log(crossings)", color = "Period",
       title = "Top 6 Predictors vs. log(crossings)",
       subtitle = "Strong correlations: these predict WHEN people cross.") +
  theme_eda + theme(strip.text = element_text(size = 9))

ggsave(file.path(OUT_DIR, "09_scatter_crossings.png"), p2d_cross,
       width = 10, height = 7, dpi = 200, bg = "white")
cat("  Saved 09_scatter_crossings.png\n")

p2d_rate <- ggplot(scatter_long, aes(x, log_rate, color = period)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30",
              linewidth = 0.5, aes(group = 1)) +
  facet_wrap(~ label, scales = "free_x", ncol = 3) +
  scale_color_manual(values = period_colors) +
  labs(x = "Predictor value", y = "log(rate + 0.01)", color = "Period",
       title = "Same Top 6 Predictors vs. log(mortality rate)",
       subtitle = "Near-zero correlations: these do NOT predict per-crossing danger.") +
  theme_eda + theme(strip.text = element_text(size = 9))

ggsave(file.path(OUT_DIR, "10_scatter_rate.png"), p2d_rate,
       width = 10, height = 7, dpi = 200, bg = "white")
cat("  Saved 10_scatter_rate.png\n")


# =========================================================================
# PART 3: SUMMARY TABLES
# =========================================================================

cat("\n=== PART 3: Period-level summary ===\n\n")

period_summary <- df_cmr %>%
  group_by(period) %>%
  summarise(
    months      = n(),
    mean_cross  = mean(crossings, na.rm = TRUE),
    mean_deaths = mean(deaths, na.rm = TRUE),
    mean_rate   = mean(rate_100, na.rm = TRUE),
    sd_rate     = sd(rate_100, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("  %-15s  %5s  %10s  %8s  %8s  %8s\n",
            "Period", "N", "Cross/mo", "Death/mo", "Rate(%)", "SD(rate)"))
cat(paste0("  ", paste(rep("-", 58), collapse = ""), "\n"))
for (i in seq_len(nrow(period_summary))) {
  r <- period_summary[i, ]
  cat(sprintf("  %-15s  %5d  %10.0f  %8.0f  %8.2f  %8.2f\n",
              r$period, r$months, r$mean_cross, r$mean_deaths,
              r$mean_rate, r$sd_rate))
}

cat("\n=== Cross-route totals ===\n\n")

cat(sprintf("  %-30s  %8s  %10s  %10s\n",
            "Route", "Deaths", "Arrivals", "Mean rate(%)"))
cat(paste0("  ", paste(rep("-", 62), collapse = ""), "\n"))
for (rt in c("CMR", "EMR", "WMR")) {
  cat(sprintf("  %-30s  %8.0f  %10.0f  %10.2f\n",
              switch(rt, CMR = "Central Mediterranean",
                     EMR = "Eastern Mediterranean",
                     WMR = "Western Mediterranean"),
              sum(df_routes[[paste0("deaths_", rt)]], na.rm = TRUE),
              sum(df_routes[[paste0("arrivals_", rt)]], na.rm = TRUE),
              mean(df_routes[[paste0("rate_", rt)]], na.rm = TRUE)))
}


# =========================================================================
# PART 4: SEASONALITY
# =========================================================================

cat("\n=== PART 4: Seasonality ===\n")

df_cmr <- df_cmr %>%
  mutate(month_name = factor(month.abb[month(date)], levels = month.abb))

# --- 4a. Monthly boxplots for all 3 outcomes ---

make_month_box <- function(df, yvar, ylab, title, fill_col) {
  ggplot(df, aes(month_name, .data[[yvar]])) +
    geom_boxplot(fill = fill_col, alpha = 0.5, outlier.size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2,
                 color = "grey20") +
    labs(x = NULL, y = ylab, title = title) +
    theme_eda
}

p4a <- (
  make_month_box(df_cmr, "crossings", "Crossing attempts",
                 "Crossing attempts by month", "steelblue") +
    scale_y_continuous(labels = scales::comma)
) /
  make_month_box(df_cmr, "deaths", "Deaths & missing",
                 "Deaths by month", "firebrick") /
  make_month_box(df_cmr, "rate_100", "Rate (%)",
                 "Mortality rate by month", "darkorange") +
  plot_annotation(
    title    = "Seasonal Patterns (Central Mediterranean Route)",
    subtitle = "Diamond = mean. Strong summer peak in crossings/deaths, weaker for rate.",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(OUT_DIR, "11_seasonality_boxplots.png"), p4a,
       width = 10, height = 10, dpi = 200, bg = "white")
cat("  Saved 11_seasonality_boxplots.png\n")


# --- 4b. Monthly means with SE ribbons (seasonal profile) ---

monthly_profile <- df_cmr %>%
  group_by(month_name) %>%
  summarise(
    mean_cross  = mean(crossings, na.rm = TRUE),
    se_cross    = sd(crossings, na.rm = TRUE) / sqrt(n()),
    mean_deaths = mean(deaths, na.rm = TRUE),
    se_deaths   = sd(deaths, na.rm = TRUE) / sqrt(n()),
    mean_rate   = mean(rate_100, na.rm = TRUE),
    se_rate     = sd(rate_100, na.rm = TRUE) / sqrt(n()),
    .groups     = "drop"
  )

make_seasonal_line <- function(df, yvar, se_var, ylab, title, color) {
  ggplot(df, aes(month_name, .data[[yvar]], group = 1)) +
    geom_ribbon(aes(ymin = .data[[yvar]] - 1.96 * .data[[se_var]],
                    ymax = .data[[yvar]] + 1.96 * .data[[se_var]]),
                fill = color, alpha = 0.2) +
    geom_line(color = color, linewidth = 0.8) +
    geom_point(color = color, size = 2) +
    labs(x = NULL, y = ylab, title = title) +
    theme_eda
}

p4b <- (
  make_seasonal_line(monthly_profile, "mean_cross", "se_cross",
                     "Crossings/mo", "Crossing attempts", "steelblue") +
    scale_y_continuous(labels = scales::comma)
) /
  make_seasonal_line(monthly_profile, "mean_deaths", "se_deaths",
                     "Deaths/mo", "Deaths & missing", "firebrick") /
  make_seasonal_line(monthly_profile, "mean_rate", "se_rate",
                     "Rate (%)", "Mortality rate", "darkorange") +
  plot_annotation(
    title    = "Average Seasonal Profile (CMR, 2011-2021)",
    subtitle = "Mean +/- 95% CI. Note: rate is nearly flat — seasonality is volume-driven.",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(OUT_DIR, "12_seasonal_profile.png"), p4b,
       width = 10, height = 10, dpi = 200, bg = "white")
cat("  Saved 12_seasonal_profile.png\n")


# =========================================================================
# PART 5: PERIOD-LEVEL DISTRIBUTIONS
# =========================================================================

cat("\n=== PART 5: Period distributions ===\n")

# --- 5a. Violin + jitter plots for each outcome ---

make_period_violin <- function(df, yvar, ylab, title) {
  ggplot(df, aes(period, .data[[yvar]], fill = period)) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_jitter(aes(color = period), width = 0.15, size = 1.5, alpha = 0.6) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4,
                 fatten = 2, color = "grey20") +
    scale_fill_manual(values = period_colors) +
    scale_color_manual(values = period_colors) +
    labs(x = NULL, y = ylab, title = title) +
    theme_eda + theme(legend.position = "none")
}

p5a <- (
  make_period_violin(df_cmr, "crossings", "Crossing attempts",
                     "Crossing attempts") +
    scale_y_continuous(labels = scales::comma)
) /
  make_period_violin(df_cmr, "deaths", "Deaths & missing",
                     "Deaths & missing") /
  make_period_violin(df_cmr, "rate_100", "Rate (%)",
                     "Mortality rate") +
  plot_annotation(
    title    = "Outcome Distributions by Policy Period (CMR)",
    subtitle = "Crossbar = mean. Volume surges in MN/SAR; rate stays remarkably stable.",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(OUT_DIR, "13_period_violins.png"), p5a,
       width = 10, height = 10, dpi = 200, bg = "white")
cat("  Saved 13_period_violins.png\n")


# =========================================================================
# PART 6: CROSS-ROUTE MORTALITY (ZOOMED)
# =========================================================================

cat("\n=== PART 6: Cross-route mortality (zoomed) ===\n")

# The 78% CMR spike in 2019 crushes the y-axis. Cap at 30%.
p6 <- ggplot(df_rate_long %>% filter(!is.na(rate)),
             aes(date, pmin(rate, 30), color = route)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_text(data = interventions, aes(x = date, y = 30, label = label),
            vjust = 1.5, hjust = 0.5, size = 2.5, color = "grey40",
            inherit.aes = FALSE) +
  scale_color_manual(values = route_colors) +
  coord_cartesian(ylim = c(0, 30)) +
  labs(x = NULL, y = "Mortality rate (%)", color = "Route",
       title = "Monthly Mortality Rate by Route (y-axis capped at 30%)",
       subtitle = "CMR consistently most dangerous. EMR spikes in 2012; WMR in 2011-2013.") +
  theme_eda

ggsave(file.path(OUT_DIR, "14_mortality_rate_by_route_zoomed.png"), p6,
       width = 10, height = 5, dpi = 200, bg = "white")
cat("  Saved 14_mortality_rate_by_route_zoomed.png\n")


# =========================================================================
# PART 7: KEY PREDICTOR TIME SERIES OVERLAID WITH OUTCOMES
# =========================================================================

cat("\n=== PART 7: Predictor time series ===\n")

# Pick 4 key predictors: SST anomaly, wave height, wind speed, wave max
key_preds <- c("sst_anomaly_central_med", "wave_height_central_med",
               "wind_speed_central_med", "wave_max_central_med")
key_labels <- c("SST anomaly", "Wave height", "Wind speed", "Wave max")

# Scale predictors to [0,1] for overlay
rescale01 <- function(x) (x - min(x, na.rm = TRUE)) /
  diff(range(x, na.rm = TRUE))

df_pred_ts <- df_cmr %>%
  select(date, all_of(key_preds)) %>%
  pivot_longer(-date, names_to = "variable", values_to = "value") %>%
  mutate(label = factor(var_labels[variable],
                        levels = key_labels))

p7a <- ggplot(df_pred_ts, aes(date, value)) +
  geom_line(color = "grey30", linewidth = 0.5) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE,
              color = "steelblue", linewidth = 0.8) +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  facet_wrap(~ label, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Value",
       title = "Key Environmental Predictors Over Time",
       subtitle = "Blue = LOESS trend. These variables are strongly seasonal.") +
  theme_eda

ggsave(file.path(OUT_DIR, "15_predictor_time_series.png"), p7a,
       width = 10, height = 6, dpi = 200, bg = "white")
cat("  Saved 15_predictor_time_series.png\n")

# --- 7b. Dual-axis-style: standardized predictor vs. standardized rate ---

df_overlay <- df_cmr %>%
  mutate(
    z_rate     = scale(log_rate)[, 1],
    z_wave_ht  = scale(wave_height_central_med)[, 1],
    z_sst_anom = scale(sst_anomaly_central_med)[, 1]
  ) %>%
  select(date, z_rate, z_wave_ht, z_sst_anom) %>%
  pivot_longer(-date, names_to = "series", values_to = "z") %>%
  mutate(series = factor(recode(series,
    z_rate     = "Mortality rate (z)",
    z_wave_ht  = "Wave height (z)",
    z_sst_anom = "SST anomaly (z)"),
    levels = c("Mortality rate (z)", "Wave height (z)", "SST anomaly (z)")))

p7b <- ggplot(df_overlay, aes(date, z, color = series)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = c("Mortality rate (z)" = "darkorange",
                                "Wave height (z)"    = "steelblue",
                                "SST anomaly (z)"    = "firebrick")) +
  labs(x = NULL, y = "Standardized value (z-score)", color = NULL,
       title = "Mortality Rate vs. Key Predictors (Standardized)",
       subtitle = "Rate does not co-move with environmental covariates — different signal.") +
  theme_eda

ggsave(file.path(OUT_DIR, "16_rate_vs_predictors_overlay.png"), p7b,
       width = 10, height = 5, dpi = 200, bg = "white")
cat("  Saved 16_rate_vs_predictors_overlay.png\n")


# =========================================================================
# PART 8: CORRELATION STABILITY (PRE vs. POST-MoU)
# =========================================================================

cat("\n=== PART 8: Correlation stability ===\n")

mou_date <- as.Date("2017-02-01")

cor_by_period <- function(df, period_label) {
  out <- data.frame(variable = available_vars, label = var_labels[available_vars],
                    period = period_label, stringsAsFactors = FALSE)
  for (v in available_vars) {
    x <- df[[v]]
    if (sum(!is.na(x) & !is.na(df$log_crossings)) < 5) {
      out[out$variable == v, c("r_crossings", "r_deaths", "r_rate")] <- NA
      next
    }
    out[out$variable == v, "r_crossings"] <- cor(x, df$log_crossings, use = "complete.obs")
    out[out$variable == v, "r_deaths"]    <- cor(x, df$log_deaths, use = "complete.obs")
    out[out$variable == v, "r_rate"]      <- cor(x, df$log_rate, use = "complete.obs")
  }
  out
}

cor_pre  <- cor_by_period(df_cmr %>% filter(date < mou_date), "Pre-MoU")
cor_post <- cor_by_period(df_cmr %>% filter(date >= mou_date), "Post-MoU")

cor_stability <- bind_rows(cor_pre, cor_post)

# Dot plot comparing pre vs. post correlations for rate
cor_rate_stab <- cor_stability %>%
  select(label, period, r_rate) %>%
  pivot_wider(names_from = period, values_from = r_rate) %>%
  mutate(label = factor(label, levels = label_order))

p8 <- ggplot(cor_rate_stab, aes(y = label)) +
  geom_segment(aes(x = `Pre-MoU`, xend = `Post-MoU`, yend = label),
               color = "grey70", linewidth = 0.4) +
  geom_point(aes(x = `Pre-MoU`), color = "steelblue", size = 2.5) +
  geom_point(aes(x = `Post-MoU`), color = "firebrick", size = 2.5) +
  geom_vline(xintercept = 0, color = "grey50") +
  labs(x = "Pearson r with log(mortality rate)",
       y = NULL,
       title = "Correlation Stability: Pre-MoU vs. Post-MoU",
       subtitle = "Blue = pre-MoU, red = post-MoU. Rate correlations are weak and unstable.") +
  theme_eda +
  theme(axis.text.y = element_text(size = 8))

ggsave(file.path(OUT_DIR, "17_correlation_stability.png"), p8,
       width = 9, height = 7, dpi = 200, bg = "white")
cat("  Saved 17_correlation_stability.png\n")

# Print stability summary
cat("\n  Pre vs. Post-MoU correlations with log(rate):\n")
cat(sprintf("  %-20s  %8s  %8s  %8s\n", "Variable", "Pre", "Post", "Diff"))
cat(paste0("  ", paste(rep("-", 50), collapse = ""), "\n"))
for (i in seq_len(nrow(cor_rate_stab))) {
  r <- cor_rate_stab[i, ]
  cat(sprintf("  %-20s  %+.3f  %+.3f  %+.3f\n",
              r$label, r$`Pre-MoU`, r$`Post-MoU`,
              r$`Post-MoU` - r$`Pre-MoU`))
}


# =========================================================================
# PART 9: EXPORT SUMMARY TABLES
# =========================================================================

cat("\n=== PART 9: Exporting tables ===\n")

TABLE_DIR <- "output/tables/eda"
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 9a. Period summary ---
write.csv(period_summary, file.path(TABLE_DIR, "period_summary.csv"),
          row.names = FALSE)
cat("  Saved period_summary.csv\n")

# --- 9b. Correlation table ---
write.csv(cor_results, file.path(TABLE_DIR, "correlations_lag0.csv"),
          row.names = FALSE)
cat("  Saved correlations_lag0.csv\n")

# --- 9c. Full descriptive statistics ---
desc_stats <- df_cmr %>%
  summarise(across(
    c(crossings, deaths, rate_100, log_crossings, log_deaths, log_rate,
      all_of(available_vars)),
    list(
      n    = ~ sum(!is.na(.x)),
      mean = ~ mean(.x, na.rm = TRUE),
      sd   = ~ sd(.x, na.rm = TRUE),
      min  = ~ min(.x, na.rm = TRUE),
      q25  = ~ quantile(.x, 0.25, na.rm = TRUE),
      med  = ~ median(.x, na.rm = TRUE),
      q75  = ~ quantile(.x, 0.75, na.rm = TRUE),
      max  = ~ max(.x, na.rm = TRUE)
    ),
    .names = "{.col}__{.fn}"
  )) %>%
  pivot_longer(everything(), names_to = "stat", values_to = "value") %>%
  separate(stat, into = c("variable", "statistic"), sep = "__") %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  mutate(label = ifelse(variable %in% names(var_labels),
                        var_labels[variable], variable))

write.csv(desc_stats, file.path(TABLE_DIR, "descriptive_statistics.csv"),
          row.names = FALSE)
cat("  Saved descriptive_statistics.csv\n")

# --- 9d. Cross-route summary ---
route_summary <- data.frame(
  route = c("Central Mediterranean", "Eastern Mediterranean",
            "Western Mediterranean"),
  total_deaths   = sapply(c("CMR", "EMR", "WMR"),
    function(r) sum(df_routes[[paste0("deaths_", r)]], na.rm = TRUE)),
  total_arrivals = sapply(c("CMR", "EMR", "WMR"),
    function(r) sum(df_routes[[paste0("arrivals_", r)]], na.rm = TRUE)),
  mean_rate      = sapply(c("CMR", "EMR", "WMR"),
    function(r) mean(df_routes[[paste0("rate_", r)]], na.rm = TRUE)),
  median_rate    = sapply(c("CMR", "EMR", "WMR"),
    function(r) median(df_routes[[paste0("rate_", r)]], na.rm = TRUE)),
  row.names = NULL
)

write.csv(route_summary, file.path(TABLE_DIR, "cross_route_summary.csv"),
          row.names = FALSE)
cat("  Saved cross_route_summary.csv\n")

# --- 9e. Correlation stability table ---
write.csv(cor_stability, file.path(TABLE_DIR, "correlations_pre_post_mou.csv"),
          row.names = FALSE)
cat("  Saved correlations_pre_post_mou.csv\n")


# =========================================================================
# PART 10: MORTALITY RATE SPIKES BY PERIOD
# =========================================================================

cat("\n=== PART 10: Spike analysis ===\n")

# Define spike thresholds
spike_thresholds <- c(5, 10)

# --- 10a. Spike frequency table ---
spike_table <- df_cmr %>%
  group_by(period) %>%
  summarise(
    months       = n(),
    mean_rate    = mean(rate_100, na.rm = TRUE),
    median_rate  = median(rate_100, na.rm = TRUE),
    sd_rate      = sd(rate_100, na.rm = TRUE),
    max_rate     = max(rate_100, na.rm = TRUE),
    pct_above_5  = mean(rate_100 > 5, na.rm = TRUE) * 100,
    pct_above_10 = mean(rate_100 > 10, na.rm = TRUE) * 100,
    n_above_5    = sum(rate_100 > 5, na.rm = TRUE),
    n_above_10   = sum(rate_100 > 10, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n  Spike frequency by period:\n\n")
cat(sprintf("  %-15s  %5s  %8s  %8s  %8s  %8s  %8s\n",
            "Period", "N", "Mean(%)", "SD", "Max(%)", ">5%", ">10%"))
cat(paste0("  ", paste(rep("-", 65), collapse = ""), "\n"))
for (i in seq_len(nrow(spike_table))) {
  r <- spike_table[i, ]
  cat(sprintf("  %-15s  %5d  %8.2f  %8.2f  %8.1f  %3d (%2.0f%%)  %3d (%2.0f%%)\n",
              r$period, r$months, r$mean_rate, r$sd_rate, r$max_rate,
              r$n_above_5, r$pct_above_5, r$n_above_10, r$pct_above_10))
}

write.csv(spike_table, file.path(TABLE_DIR, "spike_frequency_by_period.csv"),
          row.names = FALSE)
cat("  Saved spike_frequency_by_period.csv\n")

# --- 10b. Spike frequency bar chart ---
spike_long <- spike_table %>%
  select(period, pct_above_5, pct_above_10) %>%
  pivot_longer(-period, names_to = "threshold", values_to = "pct") %>%
  mutate(threshold = factor(recode(threshold,
    pct_above_5  = "> 5%",
    pct_above_10 = "> 10%"),
    levels = c("> 5%", "> 10%")))

p10a <- ggplot(spike_long, aes(period, pct, fill = threshold)) +
  geom_col(position = "dodge", width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("> 5%" = "darkorange", "> 10%" = "firebrick")) +
  labs(x = NULL, y = "% of months", fill = "Rate exceeds",
       title = "Frequency of High-Mortality Months by Policy Period",
       subtitle = "MN/SAR periods have fewer extreme months; pre-MN and post-MoU have more.") +
  theme_eda

ggsave(file.path(OUT_DIR, "18_spike_frequency.png"), p10a,
       width = 8, height = 5, dpi = 200, bg = "white")
cat("  Saved 18_spike_frequency.png\n")

# --- 10c. Time series with spike highlighting ---
p10b <- ggplot(df_cmr, aes(date, rate_100)) +
  # Shaded period bands
  annotate("rect",
           xmin = as.Date(c("2011-02-01", "2013-10-01", "2014-11-01", "2017-02-01")),
           xmax = as.Date(c("2013-09-30", "2014-10-31", "2017-01-31", "2021-09-01")),
           ymin = -Inf, ymax = Inf,
           fill = unname(period_colors), alpha = 0.10) +
  # Threshold lines
  geom_hline(yintercept = 5, linetype = "dotted", color = "darkorange", linewidth = 0.4) +
  geom_hline(yintercept = 10, linetype = "dotted", color = "firebrick", linewidth = 0.4) +
  annotate("text", x = as.Date("2021-06-01"), y = 5.5, label = "5%",
           color = "darkorange", size = 2.5, hjust = 1) +
  annotate("text", x = as.Date("2021-06-01"), y = 10.5, label = "10%",
           color = "firebrick", size = 2.5, hjust = 1) +
  # Points colored by whether above threshold
  geom_line(color = "grey50", linewidth = 0.4) +
  geom_point(aes(color = case_when(
    rate_100 > 10 ~ "> 10%",
    rate_100 > 5  ~ "5-10%",
    TRUE          ~ "< 5%"
  )), size = 1.8, alpha = 0.8) +
  scale_color_manual(values = c("< 5%" = "grey60", "5-10%" = "darkorange",
                                "> 10%" = "firebrick"),
                     name = "Rate") +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_text(data = interventions, aes(x = date, y = max(df_cmr$rate_100) * 0.95,
                                       label = label),
            vjust = 1, hjust = 0.5, size = 2.5, color = "grey30") +
  coord_cartesian(ylim = c(0, max(df_cmr$rate_100) * 1.05)) +
  labs(x = NULL, y = "Mortality rate (%)",
       title = "Mortality Rate with Spike Highlighting",
       subtitle = "Background shading = policy periods. Spikes cluster in pre-MN and post-MoU.") +
  theme_eda

ggsave(file.path(OUT_DIR, "19_rate_spikes_highlighted.png"), p10b,
       width = 11, height = 5, dpi = 200, bg = "white")
cat("  Saved 19_rate_spikes_highlighted.png\n")


# =========================================================================
# PART 11: WITHIN-POST-MoU TREND
# =========================================================================

cat("\n=== PART 11: Within-post-MoU trend ===\n")

df_post <- df_cmr %>% filter(date >= as.Date("2017-02-01"))
df_post <- df_post %>%
  mutate(months_since_mou = as.numeric(difftime(date, as.Date("2017-02-01"),
                                                 units = "days")) / 30.44)

# --- 11a. Linear trend test ---
trend_rate   <- lm(rate_100   ~ months_since_mou, data = df_post)
trend_log    <- lm(log_rate   ~ months_since_mou, data = df_post)
trend_deaths <- lm(log_deaths ~ months_since_mou, data = df_post)
trend_cross  <- lm(log_crossings ~ months_since_mou, data = df_post)

cat("\n  Linear trend within post-MoU period (Feb 2017 - Sep 2021):\n\n")
print_trend <- function(label, model) {
  s <- summary(model)
  coef <- s$coefficients["months_since_mou", ]
  cat(sprintf("  %-25s  slope = %+.4f  (SE = %.4f, p = %.3f)\n",
              label, coef[1], coef[2], coef[4]))
}
print_trend("Mortality rate (%)", trend_rate)
print_trend("log(rate + 0.01)", trend_log)
print_trend("log(deaths + 1)", trend_deaths)
print_trend("log(crossings)", trend_cross)

# --- 11b. Rolling 6-month average within post-MoU ---
df_post <- df_post %>%
  arrange(date) %>%
  mutate(
    roll_rate_6m  = zoo::rollmean(rate_100, k = 6, fill = NA, align = "right"),
    roll_deaths_6m = zoo::rollmean(deaths, k = 6, fill = NA, align = "right"),
    roll_cross_6m  = zoo::rollmean(crossings, k = 6, fill = NA, align = "right")
  )

# --- 11c. Post-MoU trend plot ---
p11 <- ggplot(df_post, aes(date, rate_100)) +
  geom_point(color = "darkorange", size = 2, alpha = 0.6) +
  geom_line(aes(y = roll_rate_6m), color = "firebrick", linewidth = 0.9,
            na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30",
              linewidth = 0.6, fill = "grey80", alpha = 0.3) +
  # Mark key sub-events
  geom_vline(xintercept = as.Date("2018-06-01"), linetype = "dotted",
             color = "grey50", linewidth = 0.3) +
  annotate("text", x = as.Date("2018-06-01"), y = max(df_post$rate_100) * 0.9,
           label = "Salvini\nclosed ports", size = 2.3, color = "grey40",
           hjust = -0.1) +
  geom_vline(xintercept = as.Date("2020-03-01"), linetype = "dotted",
             color = "grey50", linewidth = 0.3) +
  annotate("text", x = as.Date("2020-03-01"), y = max(df_post$rate_100) * 0.9,
           label = "COVID-19", size = 2.3, color = "grey40",
           hjust = -0.1) +
  labs(x = NULL, y = "Mortality rate (%)",
       title = "Mortality Rate Trend Within Post-MoU Period",
       subtitle = sprintf(
         "Red = 6-month rolling mean. Linear trend: slope = %+.3f%%/mo (p = %.3f).",
         coef(trend_rate)["months_since_mou"],
         summary(trend_rate)$coefficients["months_since_mou", 4])) +
  theme_eda

ggsave(file.path(OUT_DIR, "20_post_mou_trend.png"), p11,
       width = 10, height = 5, dpi = 200, bg = "white")
cat("  Saved 20_post_mou_trend.png\n")

# --- 11d. Post-MoU: rate, deaths, crossings (3 panels) ---
p11d_rate <- ggplot(df_post, aes(date, rate_100)) +
  geom_point(color = "darkorange", size = 1.5, alpha = 0.5) +
  geom_line(aes(y = roll_rate_6m), color = "firebrick", linewidth = 0.8,
            na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30",
              linewidth = 0.5, linetype = "dashed") +
  labs(x = NULL, y = "Rate (%)", title = "Mortality rate") +
  theme_eda

p11d_deaths <- ggplot(df_post, aes(date, deaths)) +
  geom_point(color = "firebrick", size = 1.5, alpha = 0.5) +
  geom_line(aes(y = roll_deaths_6m), color = "firebrick4", linewidth = 0.8,
            na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30",
              linewidth = 0.5, linetype = "dashed") +
  labs(x = NULL, y = "Deaths & missing", title = "Death counts") +
  theme_eda

p11d_cross <- ggplot(df_post, aes(date, crossings)) +
  geom_point(color = "steelblue", size = 1.5, alpha = 0.5) +
  geom_line(aes(y = roll_cross_6m), color = "steelblue4", linewidth = 0.8,
            na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30",
              linewidth = 0.5, linetype = "dashed") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Crossing attempts", title = "Crossings") +
  theme_eda

p11d <- p11d_rate / p11d_deaths / p11d_cross +
  plot_annotation(
    title = "Post-MoU Trends: Rate, Deaths, and Crossings (Feb 2017 - Sep 2021)",
    subtitle = "Red/blue line = 6-month rolling mean. Dashed = linear trend.",
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(OUT_DIR, "21_post_mou_three_panels.png"), p11d,
       width = 10, height = 9, dpi = 200, bg = "white")
cat("  Saved 21_post_mou_three_panels.png\n")

# --- 11e. First vs. second half of post-MoU ---
df_post <- df_post %>%
  mutate(post_half = ifelse(date < as.Date("2019-05-01"),
                            "Early MoU (Feb 17 - Apr 19)",
                            "Late MoU (May 19 - Sep 21)"))

post_half_summary <- df_post %>%
  group_by(post_half) %>%
  summarise(
    months      = n(),
    mean_rate   = mean(rate_100, na.rm = TRUE),
    median_rate = median(rate_100, na.rm = TRUE),
    sd_rate     = sd(rate_100, na.rm = TRUE),
    pct_above_5 = mean(rate_100 > 5, na.rm = TRUE) * 100,
    mean_deaths = mean(deaths, na.rm = TRUE),
    mean_cross  = mean(crossings, na.rm = TRUE),
    .groups     = "drop"
  )

cat("\n  Post-MoU: early vs. late:\n\n")
cat(sprintf("  %-30s  %5s  %8s  %8s  %8s  %8s\n",
            "Sub-period", "N", "Rate(%)", "SD", ">5%", "Deaths/mo"))
cat(paste0("  ", paste(rep("-", 65), collapse = ""), "\n"))
for (i in seq_len(nrow(post_half_summary))) {
  r <- post_half_summary[i, ]
  cat(sprintf("  %-30s  %5d  %8.2f  %8.2f  %5.0f%%  %8.0f\n",
              r$post_half, r$months, r$mean_rate, r$sd_rate,
              r$pct_above_5, r$mean_deaths))
}

write.csv(post_half_summary, file.path(TABLE_DIR, "post_mou_early_vs_late.csv"),
          row.names = FALSE)
cat("  Saved post_mou_early_vs_late.csv\n")


# =========================================================================
# PART 12: EVENT-LEVEL INCIDENT ANALYSIS (IOM Missing Migrants Project)
# =========================================================================

cat("\n=== PART 12: Event-level incident analysis ===\n")

iom_path <- "../../Data/IOM Data/Clean/iom_mmp_2014_2025_all_types.csv"
attempts_path <- "../../Data/IOM Data/Clean/med_crossings_monthlyTS.csv"
if (!file.exists(iom_path) || !file.exists(attempts_path)) {
  if (!file.exists(iom_path)) {
    cat("  IOM data not found at", iom_path, ".\n")
  }
  if (!file.exists(attempts_path)) {
    cat("  Crossings TS not found at", attempts_path, ".\n")
  }
  cat("  Skipping Part 12.\n")
} else {

iom_raw <- readr::read_csv(iom_path, show_col_types = FALSE)
iom_cmr <- iom_raw %>% filter(Route == "Central Mediterranean")
attempts_raw <- readr::read_csv(attempts_path, show_col_types = FALSE)

attempts_cmr <- attempts_raw %>%
  mutate(
    date = as.Date(date),
    sea_arrivals_in_italy = suppressWarnings(as.numeric(sea_arrivals_in_italy)),
    sea_arrivals_in_malta = suppressWarnings(as.numeric(sea_arrivals_in_malta)),
    interceptions_by_libyan_coast_guard =
      suppressWarnings(as.numeric(interceptions_by_libyan_coast_guard)),
    interceptions_by_tunisian_coast_guard =
      suppressWarnings(as.numeric(interceptions_by_tunisian_coast_guard)),
    cmr = suppressWarnings(as.numeric(cmr)),
    arrivals_cmr = rowSums(across(c(sea_arrivals_in_italy, sea_arrivals_in_malta)),
                           na.rm = TRUE),
    interceptions_cmr = dplyr::coalesce(interceptions_by_libyan_coast_guard, 0) +
      dplyr::coalesce(interceptions_by_tunisian_coast_guard, 0),
    deaths_cmr = dplyr::coalesce(cmr, 0),
    attempts_cmr = arrivals_cmr + interceptions_cmr + deaths_cmr,
    year_month = floor_date(date, "month")
  ) %>%
  select(year_month, attempts_cmr) %>%
  filter(!is.na(year_month))

# Parse dates (multiple formats in the data)
iom_cmr <- iom_cmr %>%
  mutate(
    date_clean = as.Date(`Incident date`, format = "%Y-%m-%d"),
    date_clean = if_else(is.na(date_clean),
                         as.Date(`Incident date`, format = "%d.%m.%Y"),
                         date_clean),
    # Ranges like "06.07.2015-09.07.2015": take first date
    date_clean = if_else(is.na(date_clean),
                         as.Date(sub("[-/].*", "", `Incident date`), format = "%d.%m.%Y"),
                         date_clean)
  ) %>%
  filter(!is.na(date_clean), date_clean >= "2014-01-01")

max_iom_date <- max(iom_cmr$date_clean, na.rm = TRUE)
max_attempts_month <- max(attempts_cmr$year_month[!is.na(attempts_cmr$attempts_cmr)],
                          na.rm = TRUE)
analysis_end_month <- min(floor_date(max_iom_date, "month"), max_attempts_month)
analysis_end_date <- ceiling_date(analysis_end_month, "month") - days(1)
analysis_end_label <- format(analysis_end_month, "%b %Y")

iom_cmr <- iom_cmr %>%
  filter(date_clean <= analysis_end_date) %>%
  mutate(
    dead_missing = suppressWarnings(as.numeric(`No. dead/missing`)),
    survivors    = ifelse(is.na(`No. survivors`), NA, `No. survivors`),
    year_month   = floor_date(date_clean, "month"),
    period = factor(case_when(
      date_clean < as.Date("2014-11-01") ~ "Mare Nostrum",
      date_clean < as.Date("2017-02-01") ~ "NGO SAR",
      TRUE ~ "Post-MoU"
    ), levels = c("Mare Nostrum", "NGO SAR", "Post-MoU"))
  )

cat(sprintf("  Parsed %d CMR incidents (2014-%s, through %s)\n",
            nrow(iom_cmr), format(analysis_end_month, "%Y"), analysis_end_label))

# --- 12a. Period-level incident summary ---
period_incidents <- iom_cmr %>%
  group_by(period) %>%
  summarise(
    n_incidents       = n(),
    total_dead        = sum(dead_missing, na.rm = TRUE),
    mean_per_incident = mean(dead_missing, na.rm = TRUE),
    median_per_inc    = median(dead_missing, na.rm = TRUE),
    sd_per_inc        = sd(dead_missing, na.rm = TRUE),
    pct_large_50      = mean(dead_missing >= 50, na.rm = TRUE) * 100,
    pct_large_100     = mean(dead_missing >= 100, na.rm = TRUE) * 100,
    max_incident      = max(dead_missing, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n  Incident-level summary by period:\n\n")
cat(sprintf("  %-15s  %5s  %6s  %6s  %6s  %6s  %5s  %5s  %4s\n",
            "Period", "N", "Total", "Mean", "Med", "SD", ">50", ">100", "Max"))
cat(paste0("  ", paste(rep("-", 70), collapse = ""), "\n"))
for (i in seq_len(nrow(period_incidents))) {
  r <- period_incidents[i, ]
  cat(sprintf("  %-15s  %5d  %6.0f  %6.1f  %6.0f  %6.1f  %4.1f%%  %4.1f%%  %4.0f\n",
              r$period, r$n_incidents, r$total_dead,
              r$mean_per_incident, r$median_per_inc, r$sd_per_inc,
              r$pct_large_50, r$pct_large_100, r$max_incident))
}

write.csv(period_incidents, file.path(TABLE_DIR, "incident_summary_by_period.csv"),
          row.names = FALSE)
cat("  Saved incident_summary_by_period.csv\n")

# --- 12b. Deaths-per-incident distribution by period ---
p12a <- ggplot(iom_cmr %>% filter(dead_missing > 0),
               aes(dead_missing, fill = period)) +
  geom_histogram(binwidth = 5, alpha = 0.7, color = "white", linewidth = 0.2) +
  facet_wrap(~ period, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Mare Nostrum" = "#56B4E9",
                                "NGO SAR"     = "#009E73",
                                "Post-MoU"    = "#D55E00")) +
  scale_x_continuous(breaks = seq(0, 600, 50)) +
  labs(x = "Deaths & missing per incident", y = "Number of incidents",
       title = "Distribution of Deaths Per Incident by Period",
       subtitle = "MN/SAR have long right tails (mass casualty events). Post-MoU: many small incidents.") +
  theme_eda + theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "22_deaths_per_incident_hist.png"), p12a,
       width = 10, height = 8, dpi = 200, bg = "white")
cat("  Saved 22_deaths_per_incident_hist.png\n")

# --- 12c. Boxplot: deaths per incident by period ---
p12b <- ggplot(iom_cmr %>% filter(dead_missing > 0),
               aes(period, dead_missing, fill = period)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = period), width = 0.2, size = 0.8, alpha = 0.3) +
  scale_fill_manual(values = c("Mare Nostrum" = "#56B4E9",
                                "NGO SAR"     = "#009E73",
                                "Post-MoU"    = "#D55E00")) +
  scale_color_manual(values = c("Mare Nostrum" = "#56B4E9",
                                 "NGO SAR"     = "#009E73",
                                 "Post-MoU"    = "#D55E00")) +
  coord_cartesian(ylim = c(0, 200)) +
  labs(x = NULL, y = "Deaths & missing per incident (capped at 200)",
       title = "Deaths Per Incident by Period",
       subtitle = "MN: fewer but deadlier incidents. Post-MoU: many small incidents.") +
  theme_eda + theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "23_deaths_per_incident_box.png"), p12b,
       width = 8, height = 5, dpi = 200, bg = "white")
cat("  Saved 23_deaths_per_incident_box.png\n")

# --- 12d. Monthly: number of incidents over time ---
monthly_inc <- iom_cmr %>%
  group_by(year_month, period) %>%
  summarise(
    n_incidents    = n(),
    n_known_dead   = sum(!is.na(dead_missing)),
    total_dead     = ifelse(n_known_dead > 0, sum(dead_missing, na.rm = TRUE), NA_real_),
    max_single     = ifelse(n_known_dead > 0, max(dead_missing, na.rm = TRUE), NA_real_),
    mean_per_inc   = ifelse(n_known_dead > 0, mean(dead_missing, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    top1_share = ifelse(!is.na(total_dead) & total_dead > 0,
                        max_single / total_dead * 100, NA_real_)
  ) %>%
  select(-n_known_dead) %>%
  left_join(attempts_cmr, by = "year_month") %>%
  mutate(
    inc_per_1k_attempts = ifelse(attempts_cmr > 0,
                                 n_incidents / attempts_cmr * 1000, NA_real_)
  )

p12c_n <- ggplot(monthly_inc, aes(year_month, n_incidents, fill = period)) +
  geom_col(alpha = 0.7, width = 25) +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  scale_fill_manual(values = c("Mare Nostrum" = "#56B4E9",
                                "NGO SAR"     = "#009E73",
                                "Post-MoU"    = "#D55E00")) +
  labs(x = NULL, y = "Incidents / month",
       title = "Number of recorded incidents per month") +
  theme_eda + theme(legend.position = "none")

p12c_rate <- ggplot(monthly_inc %>% filter(!is.na(inc_per_1k_attempts)),
                    aes(year_month, inc_per_1k_attempts, fill = period)) +
  geom_col(alpha = 0.7, width = 25) +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  scale_fill_manual(values = c("Mare Nostrum" = "#56B4E9",
                                "NGO SAR"     = "#009E73",
                                "Post-MoU"    = "#D55E00")) +
  labs(x = NULL, y = "Incidents / 1,000 attempts",
       title = "Incident frequency per crossing attempts (monthly)") +
  theme_eda + theme(legend.position = "none")

p12c_mean <- ggplot(monthly_inc, aes(year_month, mean_per_inc, color = period)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(data = monthly_inc,
              aes(year_month, mean_per_inc, group = 1),
              inherit.aes = FALSE,
              method = "loess", span = 0.4, se = FALSE,
              color = "grey20", linewidth = 0.9) +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  scale_color_manual(values = c("Mare Nostrum" = "#56B4E9",
                                 "NGO SAR"     = "#009E73",
                                 "Post-MoU"    = "#D55E00")) +
  labs(x = NULL, y = "Mean deaths / incident",
       title = "Average deaths per incident (monthly)") +
  theme_eda + theme(legend.position = "none")

p12c_top1 <- ggplot(monthly_inc %>% filter(!is.na(top1_share)),
                     aes(year_month, top1_share, color = period)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(data = monthly_inc %>% filter(!is.na(top1_share)),
              aes(year_month, top1_share, group = 1),
              inherit.aes = FALSE,
              method = "loess", span = 0.4, se = FALSE,
              color = "grey20", linewidth = 0.9) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey50") +
  geom_vline(data = interventions, aes(xintercept = date),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  scale_color_manual(values = c("Mare Nostrum" = "#56B4E9",
                                 "NGO SAR"     = "#009E73",
                                 "Post-MoU"    = "#D55E00")) +
  labs(x = NULL, y = "% of monthly deaths",
       title = "Share of monthly deaths from single largest incident") +
  theme_eda + theme(legend.position = "none")

x_range <- range(monthly_inc$year_month)

p12c <- (p12c_n / p12c_mean / p12c_top1) +
  plot_annotation(
    title = paste0("Incident Structure Over Time (IOM MMP, CMR 2014-",
                   format(analysis_end_month, "%Y"), ")"),
    subtitle = paste0("Through ", analysis_end_label,
                      ": post-MoU has smaller incidents and less concentration."),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  ) &
  scale_x_date(limits = x_range)

ggsave(file.path(OUT_DIR, "24_incident_structure_over_time.png"), p12c,
       width = 10, height = 9.5, dpi = 200, bg = "white")
cat("  Saved 24_incident_structure_over_time.png\n")

# --- 12e. Cause of death by period ---
cause_by_period <- iom_cmr %>%
  filter(dead_missing > 0, !is.na(`Cause of death (category)`)) %>%
  group_by(period, `Cause of death (category)`) %>%
  summarise(
    n_incidents = n(),
    total_dead  = sum(dead_missing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(period) %>%
  mutate(pct_dead = total_dead / sum(total_dead) * 100) %>%
  ungroup()

# Top causes overall
top_causes <- cause_by_period %>%
  group_by(`Cause of death (category)`) %>%
  summarise(total = sum(total_dead), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 6) %>%
  pull(`Cause of death (category)`)

p12d <- cause_by_period %>%
  filter(`Cause of death (category)` %in% top_causes) %>%
  mutate(`Cause of death (category)` = factor(`Cause of death (category)`,
                                               levels = rev(top_causes))) %>%
  ggplot(aes(pct_dead, `Cause of death (category)`, fill = period)) +
  geom_col(position = "dodge", alpha = 0.7, width = 0.6) +
  scale_fill_manual(values = c("Mare Nostrum" = "#56B4E9",
                                "NGO SAR"     = "#009E73",
                                "Post-MoU"    = "#D55E00")) +
  labs(x = "% of period's deaths", y = NULL, fill = "Period",
       title = "Cause of Death by Period (Top 6 Categories)",
       subtitle = "Drowning dominates all periods. Check for shifts in other causes.") +
  theme_eda

ggsave(file.path(OUT_DIR, "25_cause_of_death_by_period.png"), p12d,
       width = 9, height = 5, dpi = 200, bg = "white")
cat("  Saved 25_cause_of_death_by_period.png\n")

# Print cause summary
cat("\n  Cause of death breakdown:\n\n")
for (p in levels(iom_cmr$period)) {
  cat(sprintf("  --- %s ---\n", p))
  sub <- cause_by_period %>%
    filter(period == p) %>%
    arrange(desc(total_dead)) %>%
    head(5)
  for (j in seq_len(nrow(sub))) {
    cat(sprintf("    %-30s  %5.0f deaths  (%4.1f%%)\n",
                sub$`Cause of death (category)`[j],
                sub$total_dead[j], sub$pct_dead[j]))
  }
}

write.csv(cause_by_period, file.path(TABLE_DIR, "cause_of_death_by_period.csv"),
          row.names = FALSE)
cat("  Saved cause_of_death_by_period.csv\n")

# --- 12f. Geographic shift: where do deaths occur? ---
cat("\n  --- 12f. Geographic shift in death locations ---\n")

# Compute distance from departure coast (Libya + Tunisia) for each incident
iom_cmr <- iom_cmr %>%
  mutate(
    lat = as.numeric(Latitude),
    lon = as.numeric(Longitude)
  ) %>%
  filter(!is.na(lat))

# Extract Libya + Tunisia coastline for distance calculation
departure_countries <- rnaturalearth::ne_countries(
  scale = "medium", country = c("Libya", "Tunisia"), returnclass = "sf"
)
departure_coast <- st_union(departure_countries)

inc_pts <- st_as_sf(iom_cmr %>% filter(!is.na(lat), !is.na(lon)),
                    coords = c("lon", "lat"), crs = 4326, remove = FALSE)
dist_m <- as.numeric(st_distance(inc_pts, departure_coast))
iom_cmr$dist_coast_km <- NA_real_
iom_cmr$dist_coast_km[!is.na(iom_cmr$lat) & !is.na(iom_cmr$lon)] <- dist_m / 1000

# Classify into distance-based zones
iom_cmr <- iom_cmr %>%
  mutate(
    geo_zone = case_when(
      dist_coast_km < 100  ~ "Near departure coast (<100 km)",
      dist_coast_km < 250  ~ "Mid-sea (100\u2013250 km)",
      dist_coast_km >= 250 ~ "Near Italy/Malta (\u2265250 km)",
      TRUE                 ~ NA_character_
    ),
    geo_zone = factor(geo_zone, levels = c(
      "Near departure coast (<100 km)",
      "Mid-sea (100\u2013250 km)",
      "Near Italy/Malta (\u2265250 km)"
    ))
  )

# Summary by period
geo_summary <- iom_cmr %>%
  filter(!is.na(geo_zone)) %>%
  group_by(period, geo_zone) %>%
  summarise(
    n_incidents = n(),
    total_dead  = sum(dead_missing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(period) %>%
  mutate(
    pct_incidents = n_incidents / sum(n_incidents) * 100,
    pct_dead      = total_dead / sum(total_dead) * 100
  ) %>%
  ungroup()

cat("\n  Death locations by period (% of incidents):\n\n")
cat(sprintf("  %-15s  %-35s  %5s  %6s  %6s\n",
            "Period", "Zone", "N", "%inc", "%dead"))
cat(paste0("  ", paste(rep("-", 75), collapse = ""), "\n"))
for (i in seq_len(nrow(geo_summary))) {
  r <- geo_summary[i, ]
  cat(sprintf("  %-15s  %-35s  %5d  %5.1f%%  %5.1f%%\n",
              r$period, r$geo_zone, r$n_incidents,
              r$pct_incidents, r$pct_dead))
}

write.csv(geo_summary, file.path(TABLE_DIR, "death_location_by_period.csv"),
          row.names = FALSE)
cat("  Saved death_location_by_period.csv\n")

# Distance-from-coast summary by period
dist_by_period <- iom_cmr %>%
  filter(!is.na(dist_coast_km)) %>%
  group_by(period) %>%
  summarise(
    n              = n(),
    median_dist_km = median(dist_coast_km, na.rm = TRUE),
    mean_dist_km   = mean(dist_coast_km, na.rm = TRUE),
    pct_under_100  = mean(dist_coast_km < 100, na.rm = TRUE) * 100,
    pct_over_250   = mean(dist_coast_km >= 250, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n  Distance-from-coast summary by period:\n")
for (i in seq_len(nrow(dist_by_period))) {
  r <- dist_by_period[i, ]
  cat(sprintf("  %-15s  n=%d  median=%.0f km  <100km=%.1f%%  >=250km=%.1f%%\n",
              r$period, r$n, r$median_dist_km, r$pct_under_100, r$pct_over_250))
}

write.csv(dist_by_period, file.path(TABLE_DIR, "death_distance_by_period.csv"),
          row.names = FALSE)
cat("  Saved death_distance_by_period.csv\n")

# --- Figure: geographic zone distribution by period ---
p12f_zone <- geo_summary %>%
  ggplot(aes(period, pct_incidents, fill = geo_zone)) +
  geom_col(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c(
    "Near departure coast (<100 km)" = "#D55E00",
    "Mid-sea (100\u2013250 km)"      = "#E69F00",
    "Near Italy/Malta (\u2265250 km)" = "#56B4E9"
  )) +
  labs(x = NULL, y = "% of incidents",
       fill = "Zone (distance from Libya/Tunisia coast)",
       title = "Where Deaths Occur: Geographic Shift by Period",
       subtitle = "Distance from Libya/Tunisia coastline. Post-MoU: deaths shifted toward departure coast.") +
  theme_eda +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 9))

ggsave(file.path(OUT_DIR, "26_death_location_by_period.png"), p12f_zone,
       width = 8, height = 5, dpi = 200, bg = "white")
cat("  Saved 26_death_location_by_period.png\n")

# --- Figure: distance-from-coast distribution by period (density) ---
p12f_lat <- iom_cmr %>%
  filter(!is.na(dist_coast_km), dead_missing > 0) %>%
  ggplot(aes(dist_coast_km, fill = period, color = period)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  geom_vline(xintercept = c(100, 250), linetype = "dotted", color = "grey40") +
  scale_fill_manual(values = c("Mare Nostrum" = "#56B4E9",
                                "NGO SAR"     = "#009E73",
                                "Post-MoU"    = "#D55E00")) +
  scale_color_manual(values = c("Mare Nostrum" = "#56B4E9",
                                 "NGO SAR"     = "#009E73",
                                 "Post-MoU"    = "#D55E00")) +
  annotate("text", x = 50, y = Inf, label = "Near\ncoast",
           vjust = 1.5, hjust = 0.5, size = 3, color = "grey40") +
  annotate("text", x = 175, y = Inf, label = "Mid-\nsea",
           vjust = 1.5, hjust = 0.5, size = 3, color = "grey40") +
  annotate("text", x = 350, y = Inf, label = "Near\nItaly/Malta",
           vjust = 1.5, hjust = 0.5, size = 3, color = "grey40") +
  labs(x = "Distance from Libya/Tunisia coast (km)", y = "Density",
       fill = "Period", color = "Period",
       title = "Distance from Departure Coast by Period",
       subtitle = "Post-MoU deaths cluster near the departure coast; MN/SAR deaths spread further out") +
  theme_eda +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "27_death_latitude_density.png"), p12f_lat,
       width = 9, height = 5, dpi = 200, bg = "white")
cat("  Saved 27_death_latitude_density.png\n")

# --- Figure: incident map with distance-from-coast contours ---
iom_map_data <- iom_cmr %>%
  filter(!is.na(lat), !is.na(lon), dead_missing > 0,
         lon >= 5, lon <= 22, lat >= 30, lat <= 42)

inc_sf <- st_as_sf(iom_map_data, coords = c("lon", "lat"),
                   crs = 4326, remove = FALSE)

coast <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Distance buffer rings from Libya + Tunisia coast (100 km, 250 km)
departure_coast_union <- st_transform(departure_coast, 3857)
buf_100  <- st_transform(st_buffer(departure_coast_union, 100000), 4326)
buf_250  <- st_transform(st_buffer(departure_coast_union, 250000), 4326)

# Clip buffers to sea (remove land overlap) for cleaner display
bbox_clip <- st_bbox(c(xmin = 5, xmax = 22, ymin = 30, ymax = 42),
                     crs = st_crs(4326)) %>% st_as_sfc()
buf_100_ring  <- st_intersection(st_cast(buf_100, "MULTILINESTRING"), bbox_clip)
buf_250_ring  <- st_intersection(st_cast(buf_250, "MULTILINESTRING"), bbox_clip)

p12f_map <- ggplot() +
  geom_sf(data = coast, fill = "grey90", color = "grey50", linewidth = 0.3) +
  geom_sf(data = buf_100_ring, linetype = "dashed",
          color = "grey40", linewidth = 0.4, fill = NA) +
  geom_sf(data = buf_250_ring, linetype = "dashed",
          color = "grey40", linewidth = 0.4, fill = NA) +
  annotate("text", x = 11.5, y = 34.2, label = "100 km",
           size = 2.8, color = "grey30", fontface = "italic") +
  annotate("text", x = 10, y = 36.2, label = "250 km",
           size = 2.8, color = "grey30", fontface = "italic") +
  geom_sf(data = inc_sf %>% arrange(desc(period)),
          aes(color = period, size = dead_missing),
          alpha = 0.45) +
  scale_color_manual(values = c("Mare Nostrum" = "#0072B2",
                                 "NGO SAR"     = "#009E73",
                                 "Post-MoU"    = "#D55E00")) +
  scale_size_continuous(range = c(0.5, 7), breaks = c(1, 10, 50, 200, 500),
                        name = "Deaths") +
  coord_sf(xlim = c(5, 22), ylim = c(30, 42), expand = FALSE) +
  labs(color = "Period",
       title = "Death Incident Locations (Central Mediterranean, 2014\u20132025)",
       subtitle = "Dashed lines: 100 km and 250 km from Libya/Tunisia coast") +
  theme_eda +
  theme(legend.position = "right")

ggsave(file.path(OUT_DIR, "28_death_incident_map.png"), p12f_map,
       width = 10, height = 7, dpi = 200, bg = "white")
cat("  Saved 28_death_incident_map.png\n")

} # end if IOM data exists


# =========================================================================
# DONE
# =========================================================================

n_figs <- length(list.files(OUT_DIR, pattern = "\\.png$"))
n_tabs <- length(list.files(TABLE_DIR, pattern = "\\.csv$"))
cat(sprintf("\n=== EDA complete. %d figures saved to %s/ ===\n", n_figs, OUT_DIR))
cat(sprintf("=== %d tables saved to %s/ ===\n", n_tabs, TABLE_DIR))
