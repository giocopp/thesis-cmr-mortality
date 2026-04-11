# 06_sar_moderation.R
# ===================
# Does NGO SAR presence moderate the weather-mortality relationship?
#
# Logic: if rescue vessels provide a safety buffer during rough weather,
# the SWH-mortality gradient should be weaker (more negative / less positive)
# when more NGOs are active. The interaction SWH x ngo is identified from
# cross-month variation in NGO presence x within-month variation in SWH.
# Month-year FE absorb intercepts, not slopes.
#
# Specifications:
#   m0: deaths ~ swh_z * post_mou                       (baseline)
#   m1: deaths ~ swh_z * ngo_z                           (SAR replaces post_mou)
#   m2: deaths ~ swh_z * ngo_z + swh_z * ngo_z_sq        (quadratic/convex)
#   m3: deaths ~ swh_z * post_mou + swh_z * ngo_z        (horse race)
#
# All with month-year FE, NW(28) SEs.
# Sample: 2014-2021 (Rodriguez-Sanchez et al. 2023 NGO data ends 2021).
#
# Input:  analysis/data/daily_panel.RDS
#         data/processed/archive/sar_ngo_ops_daily_RS.RDS
#         data/processed/iom_mmp_incidents.RDS
# Output: output/tables/sar_moderation_results.txt
#         output/figures/sar_moderation_gradient.png

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
SEA_CAUSES <- c("Drowning", "Mixed or unknown")
CMR_INCIDENT_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia")

cat("============================================================\n")
cat("SAR MODERATION: do NGO vessels weaken the SWH-mortality link?\n")
cat("============================================================\n\n")

# ── 1. Data ────────────────────────────────────────────────
cat("--- 1. Data preparation ---\n")

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

weather <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel.RDS")) %>%
  select(date, swh, swh_prevweek)

sar <- readRDS(file.path(BASE_DIR, "data", "processed", "archive", "sar_ngo_ops_daily_RS.RDS"))

d <- weather %>%
  left_join(daily_deaths, by = "date") %>%
  left_join(sar, by = "date") %>%
  replace_na(list(n_dead_missing = 0)) %>%
  arrange(date) %>%
  mutate(
    post_mou   = as.integer(date >= MOU_DATE),
    year       = year(date),
    month_year = factor(format(date, "%Y-%m")),
    unit       = 1L
  ) %>%
  filter(!is.na(swh_prevweek), year >= 2014, year <= 2021) %>%
  mutate(
    swh_prevweek_z = as.numeric(scale(swh_prevweek)),
    ngo_z          = as.numeric(scale(n_ngo_vessels)),
    ngo_z_sq       = ngo_z^2
  )

cat(sprintf("  N = %d days (2014-2021)\n", nrow(d)))
cat(sprintf("  NGO vessels: mean=%.1f, SD=%.1f, range=%d-%d\n",
    mean(d$n_ngo_vessels), sd(d$n_ngo_vessels),
    min(d$n_ngo_vessels), max(d$n_ngo_vessels)))
cat(sprintf("  Cor(ngo, post_mou): %.3f\n", cor(d$n_ngo_vessels, d$post_mou)))
cat(sprintf("  Cor(ngo_z, swh_prevweek_z): %.3f\n", cor(d$ngo_z, d$swh_prevweek_z)))

# ── 2. Estimation ─────────────────────────────────────────
cat("\n--- 2. Estimation ---\n")

# m0: Baseline — SWH x post_mou (same as 05_reduced_form)
m0 <- fenegbin(n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou | month_year,
               data = d, vcov = NW(28), panel.id = ~unit + date)

# m1: SAR moderation — SWH x ngo (replaces post_mou with continuous SAR)
m1 <- fenegbin(n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:ngo_z | month_year,
               data = d, vcov = NW(28), panel.id = ~unit + date)

# m2: Quadratic — SWH x ngo + SWH x ngo^2 (convex/cumulative effect)
m2 <- fenegbin(n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:ngo_z +
                 swh_prevweek_z:ngo_z_sq | month_year,
               data = d, vcov = NW(28), panel.id = ~unit + date)

# m3: Horse race — SWH x post_mou + SWH x ngo
m3 <- fenegbin(n_dead_missing ~ swh_prevweek_z + swh_prevweek_z:post_mou +
                 swh_prevweek_z:ngo_z | month_year,
               data = d, vcov = NW(28), panel.id = ~unit + date)

cat("\n")
etable(m0, m1, m2, m3,
       headers = c("Baseline", "SAR linear", "SAR quadratic", "Horse race"),
       se.below = TRUE, vcov = NW(28))

# Print key coefficients
cat("\n  Key coefficients (NW(28) SEs):\n")
models <- list("m0: Baseline" = m0, "m1: SAR linear" = m1,
               "m2: SAR quadratic" = m2, "m3: Horse race" = m3)
for (nm in names(models)) {
  ct <- coeftable(models[[nm]], vcov = NW(28))
  cat(sprintf("\n  %s:\n", nm))
  for (i in seq_len(nrow(ct))) {
    cat(sprintf("    %-35s  b=%+.3f (SE=%.3f) p=%.4f\n",
        rownames(ct)[i], ct[i,1], ct[i,2],
        2 * pnorm(-abs(ct[i,1] / ct[i,2]))))
  }
}

# ── 3. Binned gradient plot ───────────────────────────────
cat("\n--- 3. Gradient by NGO vessel bins ---\n")

d <- d %>%
  mutate(ngo_bin = case_when(
    n_ngo_vessels == 0                          ~ "0 vessels",
    n_ngo_vessels >= 1 & n_ngo_vessels <= 4     ~ "1-4 vessels",
    n_ngo_vessels >= 5 & n_ngo_vessels <= 8     ~ "5-8 vessels",
    n_ngo_vessels >= 9                          ~ "9+ vessels"
  ),
  ngo_bin = factor(ngo_bin, levels = c("0 vessels", "1-4 vessels",
                                        "5-8 vessels", "9+ vessels")))

# Estimate SWH gradient separately in each bin (no FE interaction needed)
bin_results <- list()
for (b in levels(d$ngo_bin)) {
  dsub <- d %>% filter(ngo_bin == b)
  if (nrow(dsub) < 60) next
  m <- tryCatch(
    fenegbin(n_dead_missing ~ swh_prevweek_z | month_year,
             data = dsub, vcov = "hetero"),
    error = function(e) NULL
  )
  if (is.null(m)) next
  ct <- coeftable(m)
  bin_results[[b]] <- tibble(
    bin = b, n_days = nrow(dsub),
    beta = ct[1, 1], se = ct[1, 2]
  )
  cat(sprintf("  %-14s  N=%4d  beta=%+.3f (SE=%.3f)\n",
      b, nrow(dsub), ct[1,1], ct[1,2]))
}

bin_df <- bind_rows(bin_results) %>%
  mutate(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se,
         bin = factor(bin, levels = c("0 vessels", "1-4 vessels",
                                       "5-8 vessels", "9+ vessels")))

p <- ggplot(bin_df, aes(x = bin, y = beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_text(aes(label = sprintf("N=%d", n_days)), vjust = -1.5, size = 3) +
  labs(title = "SWH-mortality gradient by NGO rescue vessel count",
       subtitle = "NegBin | Month-year FE | HC SEs | 2014-2021",
       x = "Active NGO rescue vessels", y = expression(beta[SWH])) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(BASE_DIR, "output", "figures", "sar_moderation_gradient.png"),
       p, width = 8, height = 6, dpi = 200)
cat("\nSaved: output/figures/sar_moderation_gradient.png\n")

# ── 4. Save results ───────────────────────────────────────
cat("\n--- 4. Saving results ---\n")

sink(file.path(BASE_DIR, "output", "tables", "sar_moderation_results.txt"))
cat("SAR MODERATION: do NGO vessels weaken the SWH-mortality link?\n")
cat("NegBin | NW(28) SEs | Standardized SWH and NGO vessels | Month-year FE\n")
cat(sprintf("Sample: 2014-2021, N = %d\n", nrow(d)))
cat(sprintf("Cor(ngo, post_mou): %.3f\n\n", cor(d$n_ngo_vessels, d$post_mou)))

cat("m0: Baseline (SWH x post_mou)\n")
cat("m1: SAR linear (SWH x ngo_z)\n")
cat("m2: SAR quadratic (SWH x ngo_z + SWH x ngo_z^2)\n")
cat("m3: Horse race (SWH x post_mou + SWH x ngo_z)\n\n")

etable(m0, m1, m2, m3,
       headers = c("Baseline", "SAR linear", "SAR quadratic", "Horse race"),
       se.below = TRUE, vcov = NW(28))

cat("\nBinned SWH-mortality gradient by NGO count (HC SEs):\n")
for (i in seq_len(nrow(bin_df))) {
  r <- bin_df[i, ]
  cat(sprintf("  %-14s  N=%4d  beta=%+.3f (SE=%.3f)\n",
      r$bin, r$n_days, r$beta, r$se))
}
sink()
cat("Saved: output/tables/sar_moderation_results.txt\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
