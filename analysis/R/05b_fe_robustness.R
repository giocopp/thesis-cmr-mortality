# 05b_fe_robustness.R
# ===================
# Robustness of the primary reduced form (05_reduced_form_primary.R) to the
# choice of fixed effects. Same data, same outcome, same SEs, same SWH measure
# (prev-week, raw metres). Only the FE block changes.
#
# Three things this script answers:
#   1. Is month-year FE correctly specified?  -> inspect whether post_mou is
#      absorbed, how much within-FE SWH variation remains, how many levels
#      the FE has, and whether b3 is stable under alternative FE.
#   2. How does it compare to coarser FE?     -> b1, b3, SE, pseudo-R2, AIC,
#      variance-absorbed numbers across 7 FE specs.
#   3. Can variance be reduced on month-year FE without breaking identification?
#      -> cluster SE at month_year, NW bandwidths, drop zero-death days, drop
#      months with no SWH variation.
#
# Full panel only (2014-01-01 to PANEL_END), primary drown+mixed outcome.
#
# In:  analysis/data/daily_panel_complete.RDS
# Out: output/tables/05b_fe_robustness.txt
#      output/figures/05b_fe_robustness_coefplot.png

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR   <- here::here()
MOU_DATE   <- as.Date("2017-07-01")
START_DATE <- as.Date("2014-01-01")

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("05b  FE ROBUSTNESS FOR PRIMARY REDUCED FORM\n")
cat("============================================================\n\n")

# ── 1. Load panel and rebuild primary death series ─────────
cat("--- 1. Loading panel and primary death series ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS")) %>%
  select(-n_dead_missing) %>%
  arrange(date) %>%
  mutate(unit = 1L, year = year(date))

PANEL_END <- max(panel$date)

daily_primary <- build_iom_daily()
panel <- panel %>%
  left_join(daily_primary, by = "date") %>%
  replace_na(list(n_dead_missing = 0))

d <- panel %>%
  filter(between(date, START_DATE, PANEL_END)) %>%
  mutate(
    month_year       = factor(format(date, "%Y-%m")),
    month_of_year    = factor(month(date)),
    quarter_year     = factor(paste0(year(date), "Q", quarter(date))),
    week_of_year     = factor(isoweek(date)),
    year_fac         = factor(year),
    dow              = factor(wday(date, week_start = 1))
  )

cat(sprintf("  N = %d days (%s to %s)\n", nrow(d), min(d$date), max(d$date)))
cat(sprintf("  pre-MoU: %d | post-MoU: %d\n",
            sum(d$post_mou == 0), sum(d$post_mou == 1)))
cat(sprintf("  total deaths: %.0f\n", sum(d$n_dead_missing)))

# ── 2. Variance in SWH absorbed by each candidate FE ───────
# Auxiliary OLS: swh_prevweek ~ 1 | FE. R2 = share of SWH explained by FE,
# i.e., share removed from the within-FE identifying variation.
cat("\n--- 2. SWH variance absorbed by each FE ---\n")

fe_specs_aux <- list(
  "year"                 = "year_fac",
  "month-of-year"        = "month_of_year",
  "year + month"         = "year_fac + month_of_year",
  "quarter-year"         = "quarter_year",
  "month-year (current)" = "month_year",
  "month-year + dow"     = "month_year + dow"
)

absorbed <- map_dfr(names(fe_specs_aux), function(nm) {
  fe_rhs <- fe_specs_aux[[nm]]
  fml <- as.formula(paste("swh_prevweek ~ 1 |", fe_rhs))
  m <- feols(fml, data = d, notes = FALSE)
  tibble(fe_spec = nm,
         abs_r2 = unname(r2(m, "r2")),
         n_fe   = sum(vapply(fixef(m), length, integer(1))))
})
print(absorbed, n = Inf)

# ── 3. Main NegBin fits across FE specs ────────────────────
# Formula matches 05_reduced_form_primary.R exactly (no explicit post_mou main
# effect). For fine FE (month-year, month-year+dow) post_mou is fully absorbed
# within each FE cell, so it is identified through the interaction only. For
# coarse FE (year, month-of-year, year+month, quarter-year) post_mou is not
# fully absorbed — the coarse-FE b3 therefore also picks up any unexplained
# mean shift at 2017-07-01 that these FEs cannot absorb. That is a known
# feature of the comparison, not a bug; see write-up in section 7.
cat("\n--- 3. NegBin fits ---\n")

fe_specs_nb <- list(
  "year"                 = "year_fac",
  "month-of-year"        = "month_of_year",
  "year + month"         = "year_fac + month_of_year",
  "quarter-year"         = "quarter_year",
  "month-year (current)" = "month_year",
  "month-year + dow"     = "month_year + dow"
)

fit_nb <- function(fe_rhs, data = d, vcv = NW(14)) {
  fml <- as.formula(paste(
    "n_dead_missing ~ swh_prevweek + swh_prevweek:post_mou |", fe_rhs
  ))
  fenegbin(fml, data = data, vcov = vcv, panel.id = ~ unit + date)
}

models <- map(fe_specs_nb, fit_nb)

extract_row <- function(m, nm) {
  ct <- coeftable(m)
  rn <- rownames(ct)
  b1_row <- which(rn == "swh_prevweek")
  b3_row <- grep("swh_prevweek:post_mou|post_mou:swh_prevweek", rn)

  get <- function(r) if (length(r) == 1) unname(ct[r, 1:2]) else c(NA_real_, NA_real_)

  v1 <- get(b1_row)
  v3 <- get(b3_row)

  tibble(
    fe_spec = nm,
    b1      = v1[1],  b1_se = v1[2],
    b3      = v3[1],  b3_se = v3[2],
    ll      = as.numeric(logLik(m)),
    aic     = AIC(m),
    bic     = BIC(m),
    pr2     = tryCatch(unname(r2(m, "pr2")), error = \(e) NA_real_),
    nobs    = nobs(m),
    n_fe    = sum(vapply(fixef(m), length, integer(1)))
  )
}

results <- imap_dfr(models, extract_row) %>%
  left_join(absorbed %>% select(fe_spec, abs_r2), by = "fe_spec") %>%
  mutate(
    b1_p = 2 * pnorm(-abs(b1 / b1_se)),
    b3_p = 2 * pnorm(-abs(b3 / b3_se))
  )

cat("\nResults table (NegBin, NW(14)):\n")
print(results %>%
  transmute(
    fe_spec,
    `abs R2 SWH` = sprintf("%.1f%%", abs_r2 * 100),
    n_fe,
    `b1 (SWH)`   = sprintf("%+.3f (%.3f)", b1, b1_se),
    `b1 p`       = sprintf("%.3f", b1_p),
    `b3 (SWHxMoU)` = sprintf("%+.3f (%.3f)", b3, b3_se),
    `b3 p`       = sprintf("%.3f", b3_p),
    ll           = round(ll, 1),
    AIC          = round(aic, 1),
    `pseudo-R2`  = sprintf("%.3f", pr2)
  ), n = Inf)

# ── 4. Variance reduction on month-year FE ─────────────────
cat("\n--- 4. Variance reduction on month-year FE ---\n")

base_fml <- n_dead_missing ~ swh_prevweek + swh_prevweek:post_mou | month_year

m_base <- fenegbin(base_fml, data = d, vcov = NW(14),
                   panel.id = ~ unit + date)

# 4a. NW bandwidths around the new primary choice (NW(14)). NW(7) and NW(21)
# bracket the rule-of-thumb daily-data bandwidth (~N^(1/4) ≈ 8 for N≈3300, or
# 1-3 weeks). cluster(month_year) is the most conservative one-way clustering
# at the same granularity as the FE; cluster(year) is a coarser cluster; iid
# is the anti-conservative reference (residual autocorrelation ignored).
vcv_variants <- list(
  "NW(7)"               = NW(7),
  "NW(14) [current]"    = NW(14),
  "NW(21)"              = NW(21),
  "cluster: month_year" = ~ month_year,
  "cluster: year"       = ~ year_fac,
  "iid"                 = "iid"
)

vcov_rows <- map_dfr(names(vcv_variants), function(nm) {
  ct <- coeftable(m_base, vcov = vcv_variants[[nm]])
  rn <- rownames(ct)
  b3_r <- grep("swh_prevweek:post_mou|post_mou:swh_prevweek", rn)
  b1_r <- which(rn == "swh_prevweek")
  tibble(
    variant = nm,
    b1      = ct[b1_r, 1],  b1_se = ct[b1_r, 2],
    b3      = ct[b3_r, 1],  b3_se = ct[b3_r, 2],
    b1_p    = 2 * pnorm(-abs(ct[b1_r, 1] / ct[b1_r, 2])),
    b3_p    = 2 * pnorm(-abs(ct[b3_r, 1] / ct[b3_r, 2]))
  )
})

cat("\n  SE variants on month-year FE (same point estimates, different vcov):\n")
print(vcov_rows %>%
  transmute(variant,
            b1 = sprintf("%+.3f", b1),
            b1_se = sprintf("%.3f", b1_se),
            b1_p  = sprintf("%.3f", b1_p),
            b3 = sprintf("%+.3f", b3),
            b3_se = sprintf("%.3f", b3_se),
            b3_p  = sprintf("%.3f", b3_p)),
  n = Inf)

# 4b. Sample restriction: drop zero-death days; drop months with <2 unique SWH
cat("\n  Sample-restriction variants on month-year FE:\n")

sample_variants <- list(
  "full"                    = d,
  "drop zero-death days"    = d %>% filter(n_dead_missing > 0),
  "cap deaths at 100"       = d %>% filter(n_dead_missing <= 100),
  "drop FE singletons"      = d %>% group_by(month_year) %>%
                                     filter(n() >= 7) %>% ungroup()
)

samp_rows <- imap_dfr(sample_variants, function(dd, nm) {
  dd <- dd %>%
    mutate(month_year = droplevels(factor(month_year)))
  m <- tryCatch(
    fenegbin(base_fml, data = dd, vcov = NW(14),
             panel.id = ~ unit + date),
    error = function(e) NULL)
  if (is.null(m)) return(tibble(variant = nm, N = nrow(dd),
                                b1 = NA, b1_se = NA, b3 = NA, b3_se = NA,
                                b1_p = NA, b3_p = NA))
  ct <- coeftable(m)
  rn <- rownames(ct)
  b3_r <- grep("swh_prevweek:post_mou|post_mou:swh_prevweek", rn)
  b1_r <- which(rn == "swh_prevweek")
  tibble(
    variant = nm, N = nobs(m),
    b1 = ct[b1_r, 1], b1_se = ct[b1_r, 2],
    b3 = ct[b3_r, 1], b3_se = ct[b3_r, 2],
    b1_p = 2 * pnorm(-abs(ct[b1_r, 1] / ct[b1_r, 2])),
    b3_p = 2 * pnorm(-abs(ct[b3_r, 1] / ct[b3_r, 2]))
  )
})
print(samp_rows %>%
  transmute(variant, N,
            b1 = sprintf("%+.3f (%.3f, p=%.3f)", b1, b1_se, b1_p),
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p)),
  n = Inf)

# ── 5. Coefficient plot ────────────────────────────────────
cat("\n--- 5. Coefficient plot ---\n")

plot_df <- bind_rows(
  results %>% transmute(fe_spec, coef = "b1 (SWH, pre-MoU)",
                         est = b1, se = b1_se),
  results %>% transmute(fe_spec, coef = "b3 (SWH x post-MoU)",
                         est = b3, se = b3_se)
) %>%
  mutate(
    fe_spec = factor(fe_spec, levels = rev(names(fe_specs_nb))),
    ci_lo   = est - 1.96 * se,
    ci_hi   = est + 1.96 * se
  )

p <- ggplot(plot_df, aes(est, fe_spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.5, colour = "#2166AC") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, colour = "#2166AC") +
  facet_wrap(~ coef, scales = "free_x") +
  labs(title    = "Reduced-form coefficients across FE specifications",
       subtitle = sprintf("NegBin, NW(14) SEs, full panel 2014-%s",
                          format(PANEL_END, "%Y-%m")),
       x = "Coefficient (per 1 metre SWH)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(BASE_DIR, "output", "figures",
                  "05b_fe_robustness_coefplot.png"),
       p, width = 11, height = 5, dpi = 200)
cat("Saved: output/figures/05b_fe_robustness_coefplot.png\n")

# ── 6. Write text output ───────────────────────────────────
cat("\n--- 6. Saving text output ---\n")

sink(file.path(BASE_DIR, "output", "tables", "05b_fe_robustness.txt"))
old_opts <- options(tibble.width = Inf, tibble.print_max = Inf,
                    pillar.max_dec_width = 6)
on.exit(options(old_opts), add = TRUE)
cat("05b  FE ROBUSTNESS OF PRIMARY REDUCED FORM\n")
cat("==========================================\n")
cat(sprintf("Sample: %s to %s (N = %d days, %.0f deaths)\n",
            min(d$date), max(d$date), nrow(d), sum(d$n_dead_missing)))
cat("Model: fenegbin(n_dead_missing ~ swh_prevweek + swh_prevweek:post_mou | <FE>)\n")
cat("       (matches 05_reduced_form_primary.R; post_mou main effect omitted)\n")
cat("SEs:   Newey-West(14) unless stated otherwise.\n")
cat("SWH:   raw prev-week mean in metres (not standardised).\n\n")

cat("=== 1. SWH variance absorbed by FE ===\n")
print(absorbed %>%
  mutate(`abs R2` = sprintf("%.1f%%", abs_r2 * 100)) %>%
  select(fe_spec, `abs R2`, n_fe),
  n = Inf)
cat("\n")

cat("=== 2. NegBin estimates across FE specs ===\n")
print(results %>%
  transmute(fe_spec,
            `abs R2` = sprintf("%.1f%%", abs_r2 * 100),
            n_fe,
            b1 = sprintf("%+.3f (%.3f, p=%.3f)", b1, b1_se, b1_p),
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p),
            ll = round(ll, 1), AIC = round(aic, 1),
            pr2 = sprintf("%.3f", pr2)),
  n = Inf)
cat("\n")

cat("=== 3. SE variants on month-year FE (same point estimates) ===\n")
print(vcov_rows %>%
  transmute(variant,
            b1 = sprintf("%+.3f (%.3f, p=%.3f)", b1, b1_se, b1_p),
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p)),
  n = Inf)
cat("\n")

cat("=== 4. Sample restriction on month-year FE (re-fit) ===\n")
print(samp_rows %>%
  transmute(variant, N,
            b1 = sprintf("%+.3f (%.3f, p=%.3f)", b1, b1_se, b1_p),
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p)),
  n = Inf)
cat("\n")
sink()
cat("Saved: output/tables/05b_fe_robustness.txt\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
