# 33_iom_vs_united_distribution.R
# ================================
# Descriptive comparison of IOM and UNITED death series, using the SAME
# filters applied in 32_lag_iom_vs_united.R. Purely visual / summary; no model.
#
# Four views:
#   (1) Monthly deaths by source, 2014-2023
#   (2) Cumulative deaths over time
#   (3) Distribution of daily death counts (log1p density)
#   (4) Scatter: monthly IOM vs UNITED, with y=x reference
#
# Filters:
#   IOM   : build_iom_daily() defaults (incident only, split EXCLUDED; drown+mixed, central)
#   UNITED: country in CMR + Mediterranean; manner_of_death in
#           {drowned, other_unknown}; spatial join to core_corridor polygon.
#
# Out: output/figures/33_iom_vs_united_distribution.png
#      output/tables/33_iom_vs_united_distribution.txt

library(tidyverse)
library(lubridate)
library(sf)
library(patchwork)

BASE_DIR <- here::here()

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("33  IOM vs UNITED death-series comparison\n")
cat("============================================================\n\n")

# ── 1. Sample window ─────────────────────────────────────────
panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))
date_min <- min(panel$date)
date_max <- max(panel$date)

cat(sprintf("Panel window: %s to %s\n", date_min, date_max))

# ── 2. IOM daily ─────────────────────────────────────────────
iom_daily <- build_iom_daily() %>%
  rename(n_dead_iom = n_dead_missing) %>%
  filter(between(date, date_min, date_max))

# ── 3. UNITED daily (matched filter) ─────────────────────────
# UNITED daily via the shared builder. Defaults (corridor spatial join;
# country in CMR+Med; manner drowned/other_unknown) replicate the previous
# inline filter exactly — single source of truth, see _helpers.R. Date
# clamp mirrors the iom_daily pattern above.
united_daily <- build_united_daily() %>%
  filter(between(date, date_min, date_max))

# ── 4. Join on full date spine so zero-days show up ──────────
all_dates <- tibble(date = seq(date_min, date_max, by = "day"))
daily <- all_dates %>%
  left_join(iom_daily,    by = "date") %>%
  left_join(united_daily, by = "date") %>%
  replace_na(list(n_dead_iom = 0, n_dead_united = 0))

cat(sprintf("\nDaily panel: %d days\n", nrow(daily)))

# ── 5. Summary table ─────────────────────────────────────────
summarise_src <- function(x, name) {
  tibble(
    source       = name,
    total_deaths = sum(x),
    death_days   = sum(x > 0),
    zero_days    = sum(x == 0),
    mean_daily   = mean(x),
    median_daily = median(x),
    sd_daily     = sd(x),
    p75          = quantile(x, 0.75),
    p95          = quantile(x, 0.95),
    p99          = quantile(x, 0.99),
    max_daily    = max(x)
  )
}

summary_tbl <- bind_rows(
  summarise_src(daily$n_dead_iom,    "IOM"),
  summarise_src(daily$n_dead_united, "UNITED")
)
cat("\n--- Summary stats (daily counts) ---\n")
print(summary_tbl, width = Inf)

# Correlation between daily counts
r_daily   <- cor(daily$n_dead_iom, daily$n_dead_united)
monthly <- daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(iom_m    = sum(n_dead_iom),
            united_m = sum(n_dead_united), .groups = "drop")
r_monthly <- cor(monthly$iom_m, monthly$united_m)
r_monthly_log <- cor(log1p(monthly$iom_m), log1p(monthly$united_m))

cat(sprintf("\nCorrelation (daily counts):             r = %.3f\n", r_daily))
cat(sprintf("Correlation (monthly counts):           r = %.3f\n", r_monthly))
cat(sprintf("Correlation (log1p monthly counts):     r = %.3f\n", r_monthly_log))

# ── 6. Panel (1): monthly time series ────────────────────────
monthly_long <- monthly %>%
  pivot_longer(c(iom_m, united_m), names_to = "source", values_to = "deaths") %>%
  mutate(source = if_else(source == "iom_m", "IOM", "UNITED"))

p1 <- ggplot(monthly_long, aes(month, deaths, colour = source)) +
  geom_vline(xintercept = as.Date("2017-07-01"), linetype = "dotted",
             colour = "#D32F2F", linewidth = 0.6) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "(1) Monthly deaths", x = NULL, y = "deaths / month",
       colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

# ── 7. Panel (2): cumulative deaths ──────────────────────────
cum_df <- daily %>%
  arrange(date) %>%
  transmute(date,
            IOM    = cumsum(n_dead_iom),
            UNITED = cumsum(n_dead_united)) %>%
  pivot_longer(c(IOM, UNITED), names_to = "source", values_to = "cum_deaths")

p2 <- ggplot(cum_df, aes(date, cum_deaths, colour = source)) +
  geom_vline(xintercept = as.Date("2017-07-01"), linetype = "dotted",
             colour = "#D32F2F", linewidth = 0.6) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "(2) Cumulative deaths", x = NULL, y = NULL,
       colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

# ── 8. Panel (3): daily-count density (log1p) ─────────────────
density_df <- daily %>%
  select(n_dead_iom, n_dead_united) %>%
  pivot_longer(everything(), names_to = "source", values_to = "deaths") %>%
  mutate(source = if_else(source == "n_dead_iom", "IOM", "UNITED"),
         log1p_deaths = log1p(deaths))

p3 <- ggplot(density_df, aes(log1p_deaths, fill = source, colour = source)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  scale_colour_manual(values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  scale_fill_manual  (values = c("IOM" = "#2166AC", "UNITED" = "#B2182B")) +
  labs(title = "(3) Daily-count density (log1p)",
       x = "log(1 + daily deaths)", y = "density", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

# ── 9. Panel (4): IOM vs UNITED monthly scatter ──────────────
max_m <- max(c(monthly$iom_m, monthly$united_m))

p4 <- ggplot(monthly, aes(iom_m, united_m)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.6, colour = "#2166AC") +
  scale_x_continuous(limits = c(0, max_m), labels = scales::comma) +
  scale_y_continuous(limits = c(0, max_m), labels = scales::comma) +
  labs(title = sprintf("(4) IOM vs UNITED monthly (r = %.3f)", r_monthly),
       subtitle = "Dashed = y = x", x = "IOM deaths / month",
       y = "UNITED deaths / month") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

combined <- (p1 | p2) / (p3 | p4)

fig_path <- file.path(BASE_DIR, "output", "figures",
                       "33_iom_vs_united_distribution.png")
ggsave(fig_path, combined, width = 12, height = 8, dpi = 200)
cat(sprintf("\nSaved: %s\n", fig_path))

# ── 10. Save text output ─────────────────────────────────────
sink_file <- file.path(BASE_DIR, "output", "tables",
                        "33_iom_vs_united_distribution.txt")
sink(sink_file)
cat("33  IOM vs UNITED distribution comparison\n")
cat("=========================================\n")
cat(sprintf("Window: %s to %s (%d days)\n",
            date_min, date_max, nrow(daily)))
cat("IOM:    build_iom_daily() defaults (incident only, split EXCLUDED; drown+mixed, central)\n")
cat("UNITED: country in {Algeria,Italy,Libya,Malta,Tunisia,Mediterranean},\n")
cat("        manner_of_death in {drowned, other_unknown},\n")
cat("        spatial join to core corridor polygon.\n\n")

cat("Summary stats (daily counts):\n")
print(summary_tbl, width = Inf)

cat(sprintf("\nCorrelation (daily counts):         r = %.3f\n", r_daily))
cat(sprintf("Correlation (monthly counts):       r = %.3f\n", r_monthly))
cat(sprintf("Correlation (log1p monthly counts): r = %.3f\n", r_monthly_log))

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
