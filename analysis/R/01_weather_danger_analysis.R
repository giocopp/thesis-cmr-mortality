# 01_weather_danger_analysis.R
# ============================
# Descriptive analysis: weather conditions and crossing danger.
#
# Questions:
#   Q1) How dangerous is being at sea during different weather conditions?
#   Q2) Is being at sea during rough conditions more dangerous post-MoU?
#   Q3) Have people at sea during rough conditions decreased or increased post-MoU?
#
# SWH measure: spatial mean over core corridor, avg lag 1-3 (swh_prev3days)
# Sea state cutoffs: Calm < 0.5m | Medium 0.5-1.25m | Rough > 1.25m
# Fatality rate: deaths / crossings (with interceptions)
# Temporal aggregation: daily and weekly (compared)
#
# Input:
#   analysis/data/daily_panel.RDS
#
# Output:
#   output/figures/q1_q2_q3_main.png
#   output/figures/q1_q2_q3_yearly.png
#   output/tables/fatality_rate_by_sea_state.csv

library(tidyverse)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")

cat("============================================================\n")
cat("WEATHER & CROSSING DANGER ANALYSIS\n")
cat("============================================================\n\n")

# ============================================================
# 1. Load daily panel
# ============================================================
cat("--- 1. Loading daily panel ---\n")
daily_raw <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS"))

daily <- daily_raw %>%
  filter(!is.na(arrivals), !is.na(swh_prev3days),
         year(date) >= 2015, year(date) <= 2021) %>%
  mutate(
    sea      = cut(swh_prev3days,
                   breaks = c(-Inf, 0.5, 1.25, Inf),
                   labels = c("Calm", "Medium", "Rough"),
                   right  = FALSE),
    post_mou = if_else(date >= MOU_DATE, "Post-MoU", "Pre-MoU") %>%
                 factor(levels = c("Pre-MoU", "Post-MoU")),
    year     = year(date)
  )

cat(sprintf("  Days: %d | Sea state: %s\n", nrow(daily),
    paste(table(daily$sea), collapse = " / ")))

# ============================================================
# 2. Weekly panel
# ============================================================
cat("--- 2. Weekly panel ---\n")

weekly <- daily %>%
  group_by(iso_week) %>%
  summarise(
    deaths        = sum(deaths),
    arrivals      = sum(arrivals),
    interceptions = sum(intercept_per_day),
    crossings     = sum(crossings),
    n_incidents   = sum(n_incidents),
    swh_prev3days = mean(swh_prev3days),
    date_start    = min(date),
    n_days        = n(),
    .groups       = "drop"
  ) %>%
  filter(n_days >= 5) %>%
  mutate(
    sea      = cut(swh_prev3days,
                   breaks = c(-Inf, 0.5, 1.25, Inf),
                   labels = c("Calm", "Medium", "Rough"),
                   right  = FALSE),
    post_mou = if_else(date_start >= MOU_DATE, "Post-MoU", "Pre-MoU") %>%
                 factor(levels = c("Pre-MoU", "Post-MoU")),
    year     = isoyear(date_start)
  )

cat(sprintf("  Weeks: %d | Sea state: %s\n", nrow(weekly),
    paste(table(weekly$sea), collapse = " / ")))

# ============================================================
# 3. Compute tables
# ============================================================

fr_summary <- function(df, ...) {
  df %>%
    filter(crossings > 0) %>%
    group_by(...) %>%
    summarise(
      obs         = n(),
      deaths      = sum(deaths),
      n_incidents = sum(n_incidents),
      crossings   = sum(crossings),
      fat_rate    = round(100 * sum(deaths) / sum(crossings), 2),
      .groups     = "drop"
    )
}

# Q1: overall
q1 <- bind_rows(
  fr_summary(daily, sea) %>% mutate(level = "Daily"),
  fr_summary(weekly, sea) %>% mutate(level = "Weekly")
)

# Q2: by period
q2 <- bind_rows(
  fr_summary(daily, post_mou, sea) %>% mutate(level = "Daily"),
  fr_summary(weekly, post_mou, sea) %>% mutate(level = "Weekly")
)

# Q3: crossing share
share_fn <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(total_crossings = sum(crossings), obs = n(), .groups = "drop") %>%
    group_by(post_mou) %>%
    mutate(share_pct = round(100 * total_crossings / sum(total_crossings), 1)) %>%
    ungroup()
}
q3 <- bind_rows(
  share_fn(daily, post_mou, sea) %>% mutate(level = "Daily"),
  share_fn(weekly, post_mou, sea) %>% mutate(level = "Weekly")
)

# Year-by-year
yr_fr <- bind_rows(
  fr_summary(daily, year, sea) %>% mutate(level = "Daily"),
  fr_summary(weekly, year, sea) %>% mutate(level = "Weekly")
)

yr_share <- bind_rows(
  daily %>% group_by(year, sea) %>%
    summarise(total = sum(crossings), .groups = "drop") %>%
    group_by(year) %>% mutate(share = round(100 * total / sum(total), 1), level = "Daily") %>% ungroup(),
  weekly %>% group_by(year, sea) %>%
    summarise(total = sum(crossings), .groups = "drop") %>%
    group_by(year) %>% mutate(share = round(100 * total / sum(total), 1), level = "Weekly") %>% ungroup()
)

# ============================================================
# 4. Print results
# ============================================================
cat("\n============================================================\n")
cat("Q1: FATALITY RATE BY SEA STATE (overall)\n")
cat("============================================================\n")
print(q1 %>% select(level, sea, obs, deaths, n_incidents, crossings, fat_rate))

cat("\n============================================================\n")
cat("Q2: FATALITY RATE BY SEA STATE x PERIOD\n")
cat("============================================================\n")
print(q2 %>% select(level, post_mou, sea, obs, deaths, n_incidents, crossings, fat_rate))

cat("\n============================================================\n")
cat("Q3: CROSSING SHARE BY SEA STATE x PERIOD\n")
cat("============================================================\n")
print(q3 %>% select(level, post_mou, sea, total_crossings, share_pct, obs))

write_csv(q2, file.path(BASE_DIR, "output", "tables", "fatality_rate_by_sea_state.csv"))
cat("\nSaved: output/tables/fatality_rate_by_sea_state.csv\n")

# ============================================================
# 5. Plots — main (Q1, Q2, Q3)
# ============================================================
cat("\n--- Generating plots ---\n")

fill_period <- scale_fill_manual(values = c("Pre-MoU" = "#D4820E", "Post-MoU" = "#D32F2F"))
fill_sea    <- scale_fill_manual(values = c("Calm" = "#2166AC", "Medium" = "#D4820E", "Rough" = "#B2182B"))

p_q1 <- ggplot(q1, aes(x = sea, y = fat_rate, fill = sea)) +
  geom_col(alpha = 0.8, width = 0.6) +
  facet_wrap(~level) +
  fill_sea + guides(fill = "none") +
  labs(title = "Q1: Fatality rate by sea state",
       subtitle = "SWH = avg prior 3 days | Calm < 0.5m | Medium 0.5-1.25m | Rough > 1.25m",
       y = "Fatality rate (%)", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_q2 <- ggplot(q2, aes(x = sea, y = fat_rate, fill = post_mou)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  facet_wrap(~level) +
  fill_period +
  labs(title = "Q2: Fatality rate by sea state and period",
       y = "Fatality rate (%)", x = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_q3 <- ggplot(q3, aes(x = sea, y = share_pct, fill = post_mou)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  facet_wrap(~level) +
  fill_period +
  labs(title = "Q3: Share of crossings by sea state and period",
       y = "Share of crossings (%)", x = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_main <- p_q1 / p_q2 / p_q3
ggsave(file.path(BASE_DIR, "output", "figures", "q1_q2_q3_main.png"),
       p_main, width = 12, height = 12, dpi = 200)

# ============================================================
# 6. Plots — year by year
# ============================================================

p_yr_fr <- ggplot(yr_fr, aes(x = factor(year), y = fat_rate, fill = sea)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  facet_wrap(~level) +
  fill_sea +
  labs(title = "Fatality rate by sea state, year by year",
       subtitle = "SWH = avg prior 3 days | Calm < 0.5m | Medium 0.5-1.25m | Rough > 1.25m",
       y = "Fatality rate (%)", x = NULL, fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_yr_share <- ggplot(yr_share, aes(x = factor(year), y = share, fill = sea)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  facet_wrap(~level) +
  fill_sea +
  labs(title = "Share of crossings by sea state, year by year",
       y = "Share (%)", x = NULL, fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_yearly <- p_yr_fr / p_yr_share
ggsave(file.path(BASE_DIR, "output", "figures", "q1_q2_q3_yearly.png"),
       p_yearly, width = 12, height = 9, dpi = 200)

cat("Saved: output/figures/q1_q2_q3_main.png\n")
cat("Saved: output/figures/q1_q2_q3_yearly.png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
