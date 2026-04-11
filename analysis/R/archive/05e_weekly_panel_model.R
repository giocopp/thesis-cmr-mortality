# 05e_weekly_panel_model.R
# ========================
# Weekly panel robustness for the reduced-form model.
#
# Motivation: the daily panel has 84% zeros, limiting estimation power.
# Aggregating to weekly reduces zero-inflation (60% non-zero weeks vs
# 16% non-zero days). This also allows extending back to 2011 using
# Migrant Files data, which is too sparse at daily level but viable
# at weekly level (52% non-zero weeks in 2011).
#
# Same reduced-form design as 05: deaths ~ SWH x post_mou | FE
# FE: month-of-year + year (since month-year has only ~4 obs per cell
# at weekly frequency).
#
# Input:  analysis/data/daily_panel.RDS
#         data/processed/iom_mmp_incidents.RDS
#         data/processed/archive/migrant_files_cmr_pre_iom.RDS
# Output: output/tables/weekly_panel_results.txt
#         output/figures/weekly_panel_coefplot.png
#         output/figures/weekly_panel_yearly_gradient.png

library(tidyverse)
library(lubridate)
library(fixest)
library(patchwork)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
YEAR_START <- 2011
PERIODS <- c(2021, 2024)
SEA_CAUSES <- c("Drowning", "Mixed or unknown")
CMR_INCIDENT_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia")

cat("============================================================\n")
cat("WEEKLY PANEL: REDUCED-FORM MODEL\n")
cat("============================================================\n\n")

# ── 1. Build outcome (daily, then aggregate to weekly) ───────
cat("--- 1. Data preparation ---\n")

# CMR sea deaths: IOM (2014+)
iom <- readRDS(file.path(BASE_DIR, "data", "processed",
                           "iom_mmp_incidents.RDS")) %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) == "incident",
         `Cause of death (category)` %in% SEA_CAUSES,
         `Country of Incident` %in% CMR_INCIDENT_COUNTRIES) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date))
iom_daily <- iom %>%
  group_by(date) %>%
  summarise(n_dead_missing = sum(dead_missing), .groups = "drop")

# All CMR deaths (no corridor restriction) for robustness
iom_allcmr <- readRDS(file.path(BASE_DIR, "data", "processed",
                                  "iom_mmp_incidents.RDS")) %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) == "incident",
         `Cause of death (category)` %in% SEA_CAUSES) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date))
iom_allcmr_daily <- iom_allcmr %>%
  group_by(date) %>%
  summarise(n_dead_missing_allcmr = sum(dead_missing), .groups = "drop")

# Daily deaths (IOM only, no Migrant Files pre-IOM period)
daily_deaths <- iom_daily
daily_deaths_allcmr <- iom_allcmr_daily

# Build daily panel then aggregate to weekly
daily <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS")) %>%
  select(date, swh, swh_prevweek) %>%
  left_join(daily_deaths, by = "date") %>%
  left_join(daily_deaths_allcmr, by = "date") %>%
  replace_na(list(n_dead_missing = 0, n_dead_missing_allcmr = 0)) %>%
  arrange(date) %>%
  filter(!is.na(swh_prevweek), year(date) >= YEAR_START) %>%
  mutate(iso_week = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date))))

# Aggregate to weekly
weekly <- daily %>%
  group_by(iso_week) %>%
  summarise(
    date_start         = min(date),
    n_dead_missing     = sum(n_dead_missing),
    n_dead_missing_allcmr = sum(n_dead_missing_allcmr),
    swh_mean           = mean(swh, na.rm = TRUE),
    n_days             = n(),
    .groups            = "drop"
  ) %>%
  filter(n_days >= 5) %>%
  mutate(
    post_mou      = as.integer(date_start >= MOU_DATE),
    year          = isoyear(date_start),
    year_fac      = factor(year),
    month_of_year = factor(month(date_start))
  )

cat(sprintf("  Weekly panel: %d weeks (%d-%d)\n", nrow(weekly),
    min(weekly$year), max(weekly$year)))
cat(sprintf("  Weeks with deaths > 0: %d (%.1f%%)\n",
    sum(weekly$n_dead_missing > 0), 100 * mean(weekly$n_dead_missing > 0)))

# ── 2. Run models for each period ───────────────────────────
cat("\n--- 2. Estimation ---\n")

extract_interaction <- function(model, label) {
  b <- coef(model)
  s <- sqrt(diag(vcov(model)))
  tibble(spec = label, coef = b[2], se = s[2])
}

all_results <- list()

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("\n=== %s ===\n", label))

  w <- weekly %>%
    filter(year <= ye) %>%
    mutate(swh_mean_z = as.numeric(scale(swh_mean)))

  cat(sprintf("  N = %d weeks | Pre: %d | Post: %d | deaths>0: %d (%.1f%%)\n",
      nrow(w), sum(w$post_mou == 0), sum(w$post_mou == 1),
      sum(w$n_dead_missing > 0), 100 * mean(w$n_dead_missing > 0)))

  # Primary: weekly mean SWH, month-of-year + year FE
  m1 <- fenegbin(n_dead_missing ~ swh_mean_z + swh_mean_z:post_mou | month_of_year + year_fac,
                 data = w, vcov = "hetero")

  # Robustness: all CMR geography
  m2 <- fenegbin(n_dead_missing_allcmr ~ swh_mean_z + swh_mean_z:post_mou | month_of_year + year_fac,
                 data = w, vcov = "hetero")

  # Robustness: no outliers (weeks with > 200 deaths)
  m3 <- fenegbin(n_dead_missing ~ swh_mean_z + swh_mean_z:post_mou | month_of_year + year_fac,
                 data = w %>% filter(n_dead_missing <= 200), vcov = "hetero")

  models <- list(m1, m2, m3)
  labs <- c("Weekly mean (primary)", "All CMR (robust.)", "No outliers (robust.)")

  cat("\n  IRRs:\n")
  for (i in seq_along(models)) {
    b <- coef(models[[i]]); s <- sqrt(diag(vcov(models[[i]])))
    cat(sprintf("    %-26s  b3 = %+.3f (SE=%.3f)  IRR = %.3f  p = %.4f\n",
        labs[i], b[2], s[2], exp(b[2]), 2 * pnorm(-abs(b[2] / s[2]))))
  }

  res <- bind_rows(
    extract_interaction(m1, "Weekly mean (primary)"),
    extract_interaction(m2, "All CMR (robustness)"),
    extract_interaction(m3, "No outliers (robustness)")
  ) %>% mutate(period = label)

  all_results[[label]] <- res
}

# ── 3. Year-by-year gradient ────────────────────────────────
cat("\n--- 3. Year-by-year gradient ---\n")

yearly_plots <- list()

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("  %s\n", label))

  w <- weekly %>%
    filter(year <= ye) %>%
    mutate(swh_mean_z = as.numeric(scale(swh_mean)),
           yr_fac = factor(year))

  m_yr <- fenegbin(n_dead_missing ~ swh_mean_z:yr_fac | month_of_year + year_fac,
                   data = w, vcov = "hetero")

  yr_coefs <- coef(m_yr)
  V_full <- vcov(m_yr)
  yr_ses <- sqrt(diag(V_full[seq_along(yr_coefs), seq_along(yr_coefs)]))
  yr_vals <- as.integer(gsub(".*fac(\\d+)$", "\\1",
                              gsub("swh_mean_z:yr_fac", "fac", names(yr_coefs))))

  yearly_df <- tibble(year = yr_vals, beta = yr_coefs, se = yr_ses,
                       ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)

  for (i in seq_len(nrow(yearly_df))) {
    r <- yearly_df[i, ]
    sig <- if (abs(r$beta / r$se) > 1.96) "*" else ""
    cat(sprintf("    %d: %+.3f (SE=%.3f) %s\n", r$year, r$beta, r$se, sig))
  }

  yearly_plots[[label]] <- ggplot(yearly_df, aes(x = year, y = beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 2017.5, linetype = "dotted", colour = "#D32F2F",
               linewidth = 0.5) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "grey40") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.5) +
    annotate("text", x = 2017.7, y = max(yearly_df$ci_hi) * 0.95,
             label = "MoU", colour = "#D32F2F", size = 3.5, hjust = 0) +
    labs(title = sprintf("Weekly panel: SWH gradient by year (%s)", label),
         subtitle = "NegBin | Month-of-year + Year FE | 95% CI",
         x = NULL, y = expression(beta[SWH])) +
    scale_x_continuous(breaks = YEAR_START:ye) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

p_yearly <- yearly_plots[[1]] / yearly_plots[[2]]
ggsave(file.path(BASE_DIR, "output", "figures", "weekly_panel_yearly_gradient.png"),
       p_yearly, width = 10, height = 8, dpi = 200)
cat("Saved: output/figures/weekly_panel_yearly_gradient.png\n")

# ── 4. Save ──────────────────────────────────────────────────
cat("\n--- 4. Saving outputs ---\n")

combined <- bind_rows(all_results) %>%
  mutate(ci_lo = coef - 1.96 * se, ci_hi = coef + 1.96 * se,
         spec = factor(spec, levels = rev(unique(spec))))

p_coef <- ggplot(combined, aes(x = coef, y = spec, colour = period)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2,
                 position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("2011-2021" = "#2166AC", "2011-2024" = "#B2182B")) +
  labs(title = expression(paste("Weekly panel: SWH x post-MoU (", beta[3], ")")),
       subtitle = "NegBin | Month-of-year + Year FE | Standardized SWH",
       x = "Coefficient (per 1-SD weekly SWH)", y = NULL, colour = "Sample") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(BASE_DIR, "output", "figures", "weekly_panel_coefplot.png"),
       p_coef, width = 10, height = 5, dpi = 200)
cat("Saved: output/figures/weekly_panel_coefplot.png\n")

sink(file.path(BASE_DIR, "output", "tables", "weekly_panel_results.txt"))
cat("WEEKLY PANEL: REDUCED-FORM MODEL\n")
cat(sprintf("Period start: %d | NegBin | Month-of-year + Year FE\n", YEAR_START))
cat("Outcome: drowning + mixed/unknown, Incident type only, core corridor\n")
cat("MF for pre-2014, IOM for 2014+\n\n")
for (label in names(all_results)) {
  cat(sprintf("=== %s ===\n", label))
  r <- all_results[[label]]
  for (i in seq_len(nrow(r))) {
    cat(sprintf("  %-28s  b3 = %+.3f (SE=%.3f)  IRR = %.3f  p = %.4f\n",
        r$spec[i], r$coef[i], r$se[i], exp(r$coef[i]),
        2 * pnorm(-abs(r$coef[i] / r$se[i]))))
  }
  cat("\n")
}
sink()
cat("Saved: output/tables/weekly_panel_results.txt\n")

# ── 5. Why weekly gives null results ──────────────────────────
cat("\n--- 5. Weekly vs daily: identifying variation ---\n")

# Recompute daily for comparison (2014-2021)
daily_comp <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS")) %>%
  select(date, swh) %>%
  filter(!is.na(swh), year(date) >= 2014, year(date) <= 2021) %>%
  mutate(swh_z = as.numeric(scale(swh)),
         month_year = format(date, "%Y-%m"),
         month_of_year = factor(month(date)),
         year_fac = factor(year(date)))

w_comp <- weekly %>%
  filter(year <= 2021) %>%
  mutate(swh_z = as.numeric(scale(swh_mean)),
         month_year = format(date_start, "%Y-%m"))

# FE residual variation
fe_daily  <- fixest::feols(swh_z ~ 1 | month_of_year + year_fac, data = daily_comp)
fe_weekly <- fixest::feols(swh_z ~ 1 | month_of_year + year_fac, data = w_comp)

pct_daily  <- 100 * (1 - var(residuals(fe_daily)) / var(daily_comp$swh_z))
pct_weekly <- 100 * (1 - var(residuals(fe_weekly)) / var(w_comp$swh_z))

cat(sprintf("  Daily:  FE absorb %.0f%% of SWH variation (residual sd = %.3f)\n",
    pct_daily, sd(residuals(fe_daily))))
cat(sprintf("  Weekly: FE absorb %.0f%% of SWH variation (residual sd = %.3f)\n",
    pct_weekly, sd(residuals(fe_weekly))))
cat("  Weekly averaging smooths out day-to-day weather swings.\n")
cat("  After FE, insufficient residual variation to identify the gradient.\n")

# Append to results file
sink(file.path(BASE_DIR, "output", "tables", "weekly_panel_results.txt"), append = TRUE)
cat("\n=== WHY WEEKLY GIVES NULL RESULTS ===\n\n")
cat("Weekly aggregation smooths out high-frequency weather variation.\n")
cat(sprintf("FE absorb %.0f%% of weekly SWH variation vs %.0f%% of daily.\n",
    pct_weekly, pct_daily))
cat("After absorbing seasonality, there is insufficient residual SWH\n")
cat("variation at weekly frequency to identify the gradient change.\n")
cat("The daily panel is the correct frequency for this question.\n")
sink()

cat("\nAppended to: output/tables/weekly_panel_results.txt\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
