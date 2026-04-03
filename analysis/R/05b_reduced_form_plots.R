# 05b_reduced_form_plots.R
# ========================
# Diagnostic plots for the reduced-form model, for BOTH sample periods:
#   1. Gradient plot: NegBin model-implied SWH-deaths relationship
#   2. Cutoff sensitivity: beta_3 across placebo treatment dates
#
# Input:  analysis/data/daily_panel.RDS (weather)
#         data/processed/iom_mmp_incidents.RDS (incidents)
# Output: output/figures/reduced_form_gradient.png
#         output/figures/reduced_form_placebo.png

library(tidyverse)
library(lubridate)
library(fixest)
library(patchwork)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
YEAR_START <- 2014
PERIODS <- c(2021, 2024)
SEA_CAUSES <- c("Drowning", "Mixed or unknown")
CORE <- list(lon_min = 10.0, lon_max = 15.1, lat_min = 32.4, lat_max = 37.8)

cat("============================================================\n")
cat("REDUCED-FORM DIAGNOSTIC PLOTS\n")
cat("============================================================\n\n")

# ── 1. Data (build once) ────────────────────────────────────
cat("--- 1. Data preparation ---\n")

iom <- readRDS(file.path(BASE_DIR, "data", "processed",
                           "iom_mmp_incidents.RDS")) %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) == "incident",
         `Cause of death (category)` %in% SEA_CAUSES,
         Longitude >= CORE$lon_min, Longitude <= CORE$lon_max,
         Latitude  >= CORE$lat_min, Latitude  <= CORE$lat_max) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date))
daily_deaths <- iom %>%
  group_by(date) %>%
  summarise(n_dead_missing = sum(dead_missing), .groups = "drop")

daily_full <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS")) %>%
  select(date, swh, swh_prevweek, iso_week) %>%
  left_join(daily_deaths, by = "date") %>%
  replace_na(list(n_dead_missing = 0)) %>%
  arrange(date) %>%
  mutate(
    post_mou   = as.integer(date >= MOU_DATE),
    period     = if_else(date >= MOU_DATE, "Post-MoU", "Pre-MoU") %>%
                   factor(levels = c("Pre-MoU", "Post-MoU")),
    year       = year(date),
    month_year = factor(format(date, "%Y-%m"))
  ) %>%
  filter(!is.na(swh_prevweek), year >= YEAR_START)

cat(sprintf("  Full panel: %d days\n", nrow(daily_full)))

# ── 2. Gradient plots (both periods) ────────────────────────
cat("\n--- 2. Gradient plots ---\n")

col_period <- c("Pre-MoU" = "#D4820E", "Post-MoU" = "#D32F2F")
gradient_plots <- list()

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("  %s\n", label))

  d <- daily_full %>%
    filter(year <= ye) %>%
    mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)))

  swh_sd   <- sd(d$swh_prevweek)
  swh_mean <- mean(d$swh_prevweek)

  m_z <- fenegbin(n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year,
                  data = d, vcov = "hetero")
  b <- coef(m_z)
  V <- vcov(m_z)
  cat(sprintf("    b1 = %.3f, b3 = %.3f\n", b[1], b[2]))

  se_pre  <- function(z) abs(z) * sqrt(V[1, 1])
  se_post <- function(z) abs(z) * sqrt(V[1, 1] + V[2, 2] + 2 * V[1, 2])

  z_grid <- seq(-2, 2, length.out = 200)
  m_grid <- z_grid * swh_sd + swh_mean

  pred_df <- bind_rows(
    tibble(swh_m = m_grid, rel_deaths = exp(b[1] * z_grid),
           ci_lo = exp(b[1] * z_grid - 1.96 * se_pre(z_grid)),
           ci_hi = exp(b[1] * z_grid + 1.96 * se_pre(z_grid)),
           period = "Pre-MoU"),
    tibble(swh_m = m_grid, rel_deaths = exp((b[1] + b[2]) * z_grid),
           ci_lo = exp((b[1] + b[2]) * z_grid - 1.96 * se_post(z_grid)),
           ci_hi = exp((b[1] + b[2]) * z_grid + 1.96 * se_post(z_grid)),
           period = "Post-MoU")
  ) %>% mutate(period = factor(period, levels = c("Pre-MoU", "Post-MoU")))

  gradient_plots[[label]] <- ggplot(pred_df,
      aes(x = swh_m, y = rel_deaths, colour = period, fill = period)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = swh_mean, linetype = "dotted", colour = "grey50", linewidth = 0.3) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual(values = col_period) +
    scale_fill_manual(values = col_period, guide = "none") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    coord_cartesian(ylim = c(0, max(pred_df$rel_deaths) * 1.05)) +
    annotate("text", x = swh_mean + 0.03, y = 0.15,
             label = sprintf("mean = %.2fm", swh_mean),
             size = 3, hjust = 0, colour = "grey40") +
    labs(title = sprintf("SWH-mortality gradient (%s)", label),
         subtitle = "Expected deaths relative to mean-weather baseline (= 1.0)",
         x = "Previous-week average SWH (m)",
         y = "Relative expected deaths",
         colour = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
}

p_gradient <- gradient_plots[[1]] + gradient_plots[[2]]
ggsave(file.path(BASE_DIR, "output", "figures", "reduced_form_gradient.png"),
       p_gradient, width = 14, height = 6, dpi = 200)
cat("Saved: output/figures/reduced_form_gradient.png\n")

# ── 3. Placebo plots (both periods) ─────────────────────────
cat("\n--- 3. Cutoff sensitivity ---\n")

run_placebo <- function(d, cutoff) {
  d <- d %>% mutate(post_placebo = as.integer(date >= cutoff))
  tryCatch({
    m <- fenegbin(n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_placebo | month_year,
                  data = d, vcov = "hetero")
    b <- coef(m); s <- sqrt(diag(vcov(m)))
    tibble(cutoff = cutoff, beta3 = b[2], se = s[2],
           p = 2 * pnorm(-abs(b[2] / s[2])))
  }, error = function(e) {
    tibble(cutoff = cutoff, beta3 = NA_real_, se = NA_real_, p = NA_real_)
  })
}

placebo_plots <- list()

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("  %s\n", label))

  d <- daily_full %>%
    filter(year <= ye) %>%
    mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)))

  placebo_dates <- seq(as.Date(paste0(YEAR_START + 1, "-01-01")),
                       as.Date(paste0(ye - 1, "-12-01")), by = "month")
  cat(sprintf("    Testing %d placebo dates...\n", length(placebo_dates)))

  pr <- map_dfr(placebo_dates, ~run_placebo(d, .x)) %>%
    filter(!is.na(beta3)) %>%
    mutate(ci_lo = beta3 - 1.96 * se, ci_hi = beta3 + 1.96 * se,
           is_mou = cutoff == as.Date("2017-07-01"), sig = p < 0.05)

  mou_b3 <- pr$beta3[pr$is_mou]
  cat(sprintf("    MoU beta3 = %+.3f | |b3| >= MoU: %d / %d\n",
      mou_b3, sum(abs(pr$beta3) >= abs(mou_b3)), nrow(pr)))

  placebo_plots[[label]] <- ggplot(pr, aes(x = cutoff, y = beta3)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "grey40") +
    geom_line(colour = "grey30", linewidth = 0.4) +
    geom_point(aes(colour = sig), size = 1.5) +
    geom_point(data = pr %>% filter(is_mou),
               colour = "#D32F2F", size = 4, shape = 18) +
    annotate("text", x = as.Date("2017-07-01"), y = max(pr$ci_hi) * 0.95,
             label = "MoU", colour = "#D32F2F", size = 3.5, fontface = "bold") +
    scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey60"),
                        labels = c("TRUE" = "p < 0.05", "FALSE" = "p >= 0.05"),
                        name = NULL) +
    scale_x_date(breaks = as.Date(paste0(YEAR_START:ye, "-01-01")),
                 date_labels = "%Y") +
    labs(title = sprintf("Cutoff sensitivity (%s)", label),
         subtitle = "Prev-week SWH (std.) x post_cutoff | Red diamond = MoU",
         x = "Placebo treatment date", y = expression(beta[3])) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
}

p_placebo <- placebo_plots[[1]] / placebo_plots[[2]]
ggsave(file.path(BASE_DIR, "output", "figures", "reduced_form_placebo.png"),
       p_placebo, width = 11, height = 10, dpi = 200)
cat("Saved: output/figures/reduced_form_placebo.png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
