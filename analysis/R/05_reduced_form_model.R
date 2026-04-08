# 05_reduced_form_model.R
# =======================
# Reduced-form: does the MoU change the SWH-mortality relationship?
#
# Estimand: beta_3 in  deaths_t ~ b1*SWH_t + b3*SWH_t*post_mou | FE
#
# Runs all specifications for TWO sample periods (2014-2021 and 2014-2024).
#
# Specifications:
#   Primary:       prev-week SWH, core corridor deaths
#   Robustness:    prev-3d SWH, core corridor deaths
#   Robustness:    prev-week SWH, all CMR deaths
#   Robustness:    prev-week SWH, core corridor, no outliers (>100 deaths)
#   Falsification: next-7d SWH (future weather cannot cause past deaths)
#
# Input:  analysis/data/daily_panel.RDS (weather)
#         data/processed/iom_mmp_incidents.RDS (incidents)
# Output: output/tables/reduced_form_results.txt
#         output/figures/reduced_form_coefplot.png
#         output/figures/reduced_form_yearly_gradient.png

library(tidyverse)
library(lubridate)
library(fixest)
library(patchwork)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
YEAR_START <- 2014
PERIODS <- c(2021, 2024)
SEA_CAUSES <- c("Drowning", "Mixed or unknown")
CMR_INCIDENT_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia")

cat("============================================================\n")
cat("REDUCED-FORM: SWH x POST-MOU -> DEATHS\n")
cat("No conditioning on crossings (post-treatment variable)\n")
cat("============================================================\n\n")

# ── 1. Data (build once, filter per period) ──────────────────
cat("--- 1. Data preparation ---\n")

# CMR sea deaths (drowning + mixed/unknown, filtered by country of incident)
iom <- readRDS(file.path(BASE_DIR, "data", "processed",
                           "iom_mmp_incidents.RDS")) %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) == "incident",
         `Cause of death (category)` %in% SEA_CAUSES,
         `Country of Incident` %in% CMR_INCIDENT_COUNTRIES) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date))
daily_deaths <- iom %>%
  group_by(date) %>%
  summarise(n_dead_missing = sum(dead_missing), .groups = "drop")

# All CMR deaths (no corridor restriction)
iom_allcmr <- readRDS(file.path(BASE_DIR, "data", "processed",
                                  "iom_mmp_incidents.RDS")) %>%
  filter(Route == "Central Mediterranean",
         tolower(`Incident Type`) == "incident",
         `Cause of death (category)` %in% SEA_CAUSES) %>%
  mutate(date = as.Date(incident_date_clean),
         dead_missing = pmax(as.numeric(`No. dead/missing`), 0, na.rm = TRUE)) %>%
  filter(!is.na(date))
daily_deaths_allcmr <- iom_allcmr %>%
  group_by(date) %>%
  summarise(n_dead_missing_allcmr = sum(dead_missing), .groups = "drop")

cat(sprintf("  CMR deaths: %d days with deaths\n", nrow(daily_deaths)))
cat(sprintf("  All CMR:       %d days with deaths\n", nrow(daily_deaths_allcmr)))

# 1d. Merge with daily weather panel
daily_full <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS")) %>%
  select(date, swh, swh_prev3days, swh_prevweek, iso_week) %>%
  left_join(daily_deaths, by = "date") %>%
  left_join(daily_deaths_allcmr, by = "date") %>%
  replace_na(list(n_dead_missing = 0, n_dead_missing_allcmr = 0)) %>%
  arrange(date) %>%
  mutate(
    swh_next7avg = rowMeans(sapply(1:7, function(k) dplyr::lead(swh, k))),
    post_mou     = as.integer(date >= MOU_DATE),
    year         = year(date),
    month_year   = factor(format(date, "%Y-%m"))
  ) %>%
  filter(!is.na(swh_prev3days), year >= YEAR_START)

# ── 2. Run models for each period ───────────────────────────
cat("\n--- 2. Estimation ---\n")

extract_interaction <- function(model, label, vcov_type = NW(28)) {
  ct <- coeftable(model, vcov = vcov_type)
  tibble(spec = label, coef = ct[2, 1], se = ct[2, 2])
}

all_results <- list()

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("\n=== %s ===\n", label))

  d <- daily_full %>%
    filter(year <= ye) %>%
    mutate(swh_prev3days_z = as.numeric(scale(swh_prev3days)),
           swh_prevweek_z  = as.numeric(scale(swh_prevweek)),
           swh_next7avg_z  = as.numeric(scale(swh_next7avg)))

  cat(sprintf("  N = %d days | Pre: %d | Post: %d\n",
      nrow(d), sum(d$post_mou == 0), sum(d$post_mou == 1)))

  d1 <- d %>% filter(!is.na(swh_prevweek_z)) %>% mutate(unit = 1L)
  d2 <- d %>% mutate(unit = 1L)
  d4 <- d %>% filter(!is.na(swh_prevweek_z), n_dead_missing <= 100) %>% mutate(unit = 1L)
  d5 <- d %>% filter(!is.na(swh_next7avg_z)) %>% mutate(unit = 1L)

  m1 <- fenegbin(n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year,
                 data = d1, vcov = NW(28), panel.id = ~unit + date)
  m2 <- fenegbin(n_dead_missing ~ swh_prev3days_z + swh_prev3days_z:post_mou | month_year,
                 data = d2, vcov = NW(28), panel.id = ~unit + date)
  m3 <- fenegbin(n_dead_missing_allcmr ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year,
                 data = d1, vcov = NW(28), panel.id = ~unit + date)
  m4 <- fenegbin(n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year,
                 data = d4, vcov = NW(28), panel.id = ~unit + date)
  m5 <- fenegbin(n_dead_missing ~ swh_next7avg_z + swh_next7avg_z:post_mou | month_year,
                 data = d5, vcov = NW(28), panel.id = ~unit + date)

  models <- list(m1, m2, m3, m4, m5)
  labs <- c("Prev-week (primary)", "Prev-3d (robust.)",
            "All CMR (robust.)", "No outliers (robust.)",
            "Next-7d (falsif.)")

  cat("\n  --- NW(28) SEs ---\n")
  etable(m1, m2, m3, m4, m5, headers = labs, se.below = TRUE, vcov = NW(28))

  cat("\n  --- NW(14) SEs ---\n")
  etable(m1, m2, m3, m4, m5, headers = labs, se.below = TRUE, vcov = NW(14))

  cat("\n  IRRs (NW(28) SEs):\n")
  for (i in seq_along(models)) {
    ct <- coeftable(models[[i]], vcov = NW(28))
    cat(sprintf("    %-24s  b3 = %+.3f (SE=%.3f)  IRR = %.3f  p = %.4f\n",
        labs[i], ct[2,1], ct[2,2], exp(ct[2,1]), 2 * pnorm(-abs(ct[2,1] / ct[2,2]))))
  }

  res <- bind_rows(
    extract_interaction(m1, "Prev-week (primary)"),
    extract_interaction(m2, "Prev-3d (robustness)"),
    extract_interaction(m3, "All CMR (robustness)"),
    extract_interaction(m4, "No outliers (robustness)"),
    extract_interaction(m5, "Next-7d (falsification)")
  ) %>% mutate(period = label)

  all_results[[label]] <- res
}

# ── 3. Combined coefficient plot ────────────────────────────
cat("\n--- 3. Combined coefficient plot ---\n")

combined <- bind_rows(all_results) %>%
  mutate(
    ci_lo = coef - 1.96 * se,
    ci_hi = coef + 1.96 * se,
    spec  = factor(spec, levels = rev(c(
      "Prev-week (primary)", "Prev-3d (robustness)",
      "All CMR (robustness)", "No outliers (robustness)",
      "Next-7d (falsification)")))
  )

p_coef <- ggplot(combined, aes(x = coef, y = spec, colour = period)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2,
                 position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("2014-2021" = "#2166AC", "2014-2024" = "#B2182B")) +
  labs(title = expression(paste("Reduced-form: SWH x post-MoU interaction (", beta[3], ")")),
       subtitle = "NegBin, Newey-West(28) SEs, 95% CI | No conditioning on crossings",
       x = "Coefficient (per 1-SD SWH)", y = NULL, colour = "Sample") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(BASE_DIR, "output", "figures", "reduced_form_coefplot.png"),
       p_coef, width = 10, height = 6, dpi = 200)
cat("Saved: output/figures/reduced_form_coefplot.png\n")

# ── 4. Year-by-year gradient (both periods) ─────────────────
cat("\n--- 4. Year-by-year gradient ---\n")

yearly_plots <- list()

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("  %s\n", label))

  d <- daily_full %>%
    filter(year <= ye, !is.na(swh_prevweek)) %>%
    mutate(swh_prevweek_z = as.numeric(scale(swh_prevweek)),
           year_fac = factor(year),
           unit = 1L)

  m_yr <- fenegbin(n_dead_missing ~ swh_prevweek_z:year_fac | month_year,
                   data = d, vcov = NW(28), panel.id = ~unit + date)

  yr_coefs <- coef(m_yr)
  V_full <- vcov(m_yr, vcov = NW(28))
  yr_ses <- sqrt(diag(V_full[seq_along(yr_coefs), seq_along(yr_coefs)]))
  yr_vals <- as.integer(gsub(".*fac(\\d+)$", "\\1",
                              gsub("swh_prevweek_z:year_fac", "fac", names(yr_coefs))))

  yearly_df <- tibble(year = yr_vals, beta = yr_coefs, se = yr_ses,
                       ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se,
                       period = label)

  for (i in seq_len(nrow(yearly_df))) {
    r <- yearly_df[i, ]
    sig <- if (abs(r$beta / r$se) > 1.96) "*" else ""
    cat(sprintf("    %d: %+.3f (SE=%.3f) %s\n", r$year, r$beta, r$se, sig))
  }

  p <- ggplot(yearly_df, aes(x = year, y = beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 2017.5, linetype = "dotted", colour = "#D32F2F",
               linewidth = 0.5) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "grey40") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.5) +
    annotate("text", x = 2017.7, y = max(yearly_df$ci_hi) * 0.95,
             label = "MoU", colour = "#D32F2F", size = 3.5, hjust = 0) +
    labs(title = sprintf("SWH-mortality gradient by year (%s)", label),
         subtitle = "NegBin coefficient on standardized prev-week SWH | Month-year FE | NW(28) SEs | 95% CI",
         x = NULL, y = expression(beta[SWH])) +
    scale_x_continuous(breaks = YEAR_START:ye) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  yearly_plots[[label]] <- p
}

p_yearly_combined <- yearly_plots[[1]] / yearly_plots[[2]]
ggsave(file.path(BASE_DIR, "output", "figures", "reduced_form_yearly_gradient.png"),
       p_yearly_combined, width = 10, height = 8, dpi = 200)
cat("Saved: output/figures/reduced_form_yearly_gradient.png\n")

# ── 5. Save results table ───────────────────────────────────
cat("\n--- 5. Saving results table ---\n")

sink(file.path(BASE_DIR, "output", "tables", "reduced_form_results.txt"))
cat("REDUCED-FORM MODEL: n_dead_missing ~ SWH x post_mou | month-year FE\n")
cat("NegBin | Newey-West(28) SEs | Standardized SWH (per 1-SD)\n")
cat("Outcome: drowning + mixed/unknown, Incident type only\n\n")
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
cat("Saved: output/tables/reduced_form_results.txt\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
