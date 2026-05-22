# ── FE-spec robustness of the primary reduced form ─────────────────────────
# Compares b1/b3 across year, month-of-year, year+month, quarter-year,
# month-year (primary), and month-year + dow. NegBin + Poisson, NW(14).

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("24  FE ROBUSTNESS (aligned with 20_primary_model.R)\n")
cat("============================================================\n\n")

# ── 1. 05d data prep ─────────────────────────────────────────
cat("--- 1. Loading panel + 05d data prep ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

iom_daily <- build_iom_daily()

panel <- panel |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  replace_na(list(n_dead_iom = 0)) |>
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag7  = dplyr::lag(zoo::rollmeanr(living_crossings, k = 7,
                                           fill = NA, align = "right"), 1),
    lc_lag14 = dplyr::lag(zoo::rollmeanr(living_crossings, k = 7,
                                           fill = NA, align = "right"), 8),
    unit           = 1L,
    year_fac       = factor(year),
    month_year_fac = factor(month_year),
    month_of_year  = factor(month(date)),
    quarter_year   = factor(paste0(year(date), "Q", quarter(date))),
    week_of_year   = factor(isoweek(date)),
    dow            = factor(wday(date, week_start = 1))
  )

d <- panel |> filter(!is.na(lc_lag14), !is.na(swh_prev5days))

cat(sprintf("  N = %d days (%s to %s)\n",
            nrow(d), min(d$date), max(d$date)))
cat(sprintf("  pre-MoU: %d | post-MoU: %d\n",
            sum(d$post_mou == 0), sum(d$post_mou == 1)))
cat(sprintf("  total deaths (IOM primary): %.0f\n", sum(d$n_dead_iom)))

# ── 2. SWH variance absorbed by each candidate FE ────────────
# Auxiliary OLS: swh_prev5days ~ 1 | FE. R^2 = share of SWH explained by FE,
# i.e., share removed from the within-FE identifying variation.
cat("\n--- 2. SWH variance absorbed by each FE ---\n")

fe_specs <- list(
  "year"                       = "year_fac",
  "month-of-year"              = "month_of_year",
  "year + month"               = "year_fac + month_of_year",
  "quarter-year"               = "quarter_year",
  "month-year (primary/20)"    = "month_year_fac",
  "month-year + dow"           = "month_year_fac + dow"
)

absorbed <- map_dfr(names(fe_specs), function(nm) {
  fml <- as.formula(paste("swh_prev5days ~ 1 |", fe_specs[[nm]]))
  m <- feols(fml, data = d, notes = FALSE)
  tibble(fe_spec = nm,
         abs_r2 = unname(r2(m, "r2")),
         n_fe   = sum(vapply(fixef(m), length, integer(1))))
})
print(absorbed, n = Inf)

# ── 3. Main fits across FE specs: NegBin + Poisson ───────────
# For fine FE (month-year, month-year+dow), post_mou is fully absorbed within
# each FE cell, so it is identified through the interaction only. For coarse
# FE, post_mou is not fully absorbed — the coarse-FE b3 also picks up any
# unexplained mean shift at 2017-02-02 those FEs cannot absorb. Known
# feature of the comparison, not a bug.
cat("\n--- 3. Fits across FE specs (NegBin + Poisson, NW(14)) ---\n")

fit_fam <- function(fn, fe_rhs) {
  fml <- as.formula(paste(
    "n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou |", fe_rhs
  ))
  fn(fml, data = d, vcov = NW(14), panel.id = ~unit + date)
}

models_nb   <- map(fe_specs, \(fe) fit_fam(fenegbin, fe))
models_pois <- map(fe_specs, \(fe) fit_fam(fepois,   fe))

extract_row <- function(m, nm, family) {
  ct <- coeftable(m, vcov = NW(14))
  rn <- rownames(ct)
  b1_row <- which(rn == "swh_prev5days")
  b3_row <- grep(":post_mou$", rn)
  get <- function(r) if (length(r) == 1) unname(ct[r, 1:2]) else c(NA_real_, NA_real_)
  v1 <- get(b1_row); v3 <- get(b3_row)
  tibble(
    fe_spec = nm, family = family,
    b1 = v1[1], b1_se = v1[2],
    b3 = v3[1], b3_se = v3[2],
    ll   = as.numeric(logLik(m)),
    aic  = AIC(m), bic = BIC(m),
    pr2  = tryCatch(unname(r2(m, "pr2")), error = \(e) NA_real_),
    nobs = nobs(m),
    n_fe = sum(vapply(fixef(m), length, integer(1)))
  )
}

results <- bind_rows(
  imap_dfr(models_nb,   \(m, nm) extract_row(m, nm, "NegBin")),
  imap_dfr(models_pois, \(m, nm) extract_row(m, nm, "Poisson"))
) |>
  left_join(absorbed |> select(fe_spec, abs_r2), by = "fe_spec") |>
  mutate(
    b1_p = 2 * pnorm(-abs(b1 / b1_se)),
    b3_p = 2 * pnorm(-abs(b3 / b3_se))
  )

cat("\n  NegBin results:\n")
print(results |> filter(family == "NegBin") |>
  transmute(fe_spec,
            `abs R2 SWH` = sprintf("%.1f%%", abs_r2 * 100),
            n_fe,
            b1 = sprintf("%+.3f (%.3f, p=%.3f)", b1, b1_se, b1_p),
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p),
            AIC = round(aic, 1)),
  n = Inf)

cat("\n  Poisson results:\n")
print(results |> filter(family == "Poisson") |>
  transmute(fe_spec,
            `abs R2 SWH` = sprintf("%.1f%%", abs_r2 * 100),
            n_fe,
            b1 = sprintf("%+.3f (%.3f, p=%.3f)", b1, b1_se, b1_p),
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p)),
  n = Inf)

# ── 4. SE variants on primary month_year_fac FE ──────────────
cat("\n--- 4. SE variants on month_year_fac (NegBin + Poisson) ---\n")

m_nb_base   <- fit_fam(fenegbin, "month_year_fac")
m_pois_base <- fit_fam(fepois,   "month_year_fac")

vcv_variants <- list(
  "NW(7)"                  = NW(7),
  "NW(14) [primary]"       = NW(14),
  "NW(21)"                 = NW(21),
  "cluster: month_year"    = ~ month_year_fac,
  "cluster: year"          = ~ year_fac,
  "iid"                    = "iid"
)

vcov_rows <- bind_rows(
  map_dfr(names(vcv_variants), function(nm) {
    ct <- coeftable(m_nb_base, vcov = vcv_variants[[nm]])
    rn <- rownames(ct)
    b3_r <- grep(":post_mou$", rn)
    b1_r <- which(rn == "swh_prev5days")
    tibble(family = "NegBin", variant = nm,
           b1 = ct[b1_r, 1], b1_se = ct[b1_r, 2],
           b3 = ct[b3_r, 1], b3_se = ct[b3_r, 2],
           b1_p = 2 * pnorm(-abs(ct[b1_r, 1] / ct[b1_r, 2])),
           b3_p = 2 * pnorm(-abs(ct[b3_r, 1] / ct[b3_r, 2])))
  }),
  map_dfr(names(vcv_variants), function(nm) {
    ct <- coeftable(m_pois_base, vcov = vcv_variants[[nm]])
    rn <- rownames(ct)
    b3_r <- grep(":post_mou$", rn)
    b1_r <- which(rn == "swh_prev5days")
    tibble(family = "Poisson", variant = nm,
           b1 = ct[b1_r, 1], b1_se = ct[b1_r, 2],
           b3 = ct[b3_r, 1], b3_se = ct[b3_r, 2],
           b1_p = 2 * pnorm(-abs(ct[b1_r, 1] / ct[b1_r, 2])),
           b3_p = 2 * pnorm(-abs(ct[b3_r, 1] / ct[b3_r, 2])))
  })
)

cat("\n  NegBin:\n")
print(vcov_rows |> filter(family == "NegBin") |>
  transmute(variant,
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p)),
  n = Inf)

cat("\n  Poisson:\n")
print(vcov_rows |> filter(family == "Poisson") |>
  transmute(variant,
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p)),
  n = Inf)

# ── 5. Sample-restriction variants on month_year_fac FE ──────
cat("\n--- 5. Sample-restriction variants (NegBin, NW(14)) ---\n")

sample_variants <- list(
  "full"                    = d,
  "drop zero-death days"    = d |> filter(n_dead_iom > 0),
  "cap deaths at 100"       = d |> filter(n_dead_iom <= 100),
  "drop FE singletons"      = d |> group_by(month_year_fac) |>
                                     filter(n() >= 7) |> ungroup()
)

samp_rows <- imap_dfr(sample_variants, function(dd, nm) {
  dd <- dd |> mutate(month_year_fac = droplevels(factor(month_year_fac)))
  m <- tryCatch(
    fenegbin(n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
             data = dd, vcov = NW(14), panel.id = ~unit + date),
    error = function(e) NULL)
  if (is.null(m)) return(tibble(variant = nm, N = nrow(dd),
                                b1 = NA, b1_se = NA, b3 = NA, b3_se = NA,
                                b1_p = NA, b3_p = NA))
  ct <- coeftable(m, vcov = NW(14))
  rn <- rownames(ct)
  b3_r <- grep(":post_mou$", rn)
  b1_r <- which(rn == "swh_prev5days")
  tibble(
    variant = nm, N = nobs(m),
    b1 = ct[b1_r, 1], b1_se = ct[b1_r, 2],
    b3 = ct[b3_r, 1], b3_se = ct[b3_r, 2],
    b1_p = 2 * pnorm(-abs(ct[b1_r, 1] / ct[b1_r, 2])),
    b3_p = 2 * pnorm(-abs(ct[b3_r, 1] / ct[b3_r, 2]))
  )
})
print(samp_rows |>
  transmute(variant, N,
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p)),
  n = Inf)

# ── 6. Coefficient plot ──────────────────────────────────────
cat("\n--- 6. Coefficient plot ---\n")

plot_df <- results |>
  transmute(fe_spec, family, coef = "b3 (SWH x post-MoU)",
            est = b3, se = b3_se) |>
  bind_rows(results |>
              transmute(fe_spec, family, coef = "b1 (SWH, pre-MoU)",
                        est = b1, se = b1_se)) |>
  mutate(
    fe_spec = factor(fe_spec, levels = rev(names(fe_specs))),
    ci_lo   = est - 1.96 * se,
    ci_hi   = est + 1.96 * se
  )

p <- ggplot(plot_df, aes(est, fe_spec, colour = family)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("NegBin" = "#2166AC", "Poisson" = "#B2182B")) +
  facet_wrap(~ coef, scales = "free_x") +
  labs(title    = "FE robustness: coefficients across FE specifications",
       subtitle = sprintf("Sample: %s to %s, N = %d days. NW(14) SEs.",
                          min(d$date), max(d$date), nrow(d)),
       x = "Coefficient (per 1 metre SWH)", y = NULL, colour = "Family") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

fig_out <- fig_path("06_robustness", "02_fe_robustness_coefplot.png")
ggsave(fig_out, p, width = 11, height = 5, dpi = 200)
cat(sprintf("Saved: %s\n", fig_out))

# ── 7. Save text output ──────────────────────────────────────
cat("\n--- 7. Saving text output ---\n")

sink_file <- tbl_path("06_robustness", "02_fe_robustness.txt")
sink(sink_file)
old_opts <- options(tibble.width = Inf, tibble.print_max = Inf)
on.exit(options(old_opts), add = TRUE)

cat("24  FE ROBUSTNESS OF PRIMARY REDUCED FORM\n")
cat("Aligned with 20_primary_model.R\n")
cat("=========================================\n")
cat(sprintf("Sample: %s to %s (N = %d days, %.0f deaths)\n",
            min(d$date), max(d$date), nrow(d), sum(d$n_dead_iom)))
cat("Model: n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou | <FE>\n")
cat("       (matches 20_primary_model.R; post_mou main effect omitted)\n")
cat("Families: NegBin (fenegbin) + Poisson QMLE (fepois)\n")
cat("SEs:   Newey-West(14) unless stated otherwise.\n\n")

cat("=== 1. SWH variance absorbed by FE ===\n")
print(absorbed |>
  mutate(`abs R2` = sprintf("%.1f%%", abs_r2 * 100)) |>
  select(fe_spec, `abs R2`, n_fe),
  n = Inf)

cat("\n=== 2. NegBin estimates across FE specs ===\n")
print(results |> filter(family == "NegBin") |>
  transmute(fe_spec,
            `abs R2` = sprintf("%.1f%%", abs_r2 * 100),
            n_fe,
            b1 = sprintf("%+.3f (%.3f, p=%.3f)", b1, b1_se, b1_p),
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p),
            AIC = round(aic, 1), pr2 = sprintf("%.3f", pr2)),
  n = Inf)

cat("\n=== 3. Poisson estimates across FE specs ===\n")
print(results |> filter(family == "Poisson") |>
  transmute(fe_spec,
            `abs R2` = sprintf("%.1f%%", abs_r2 * 100),
            n_fe,
            b1 = sprintf("%+.3f (%.3f, p=%.3f)", b1, b1_se, b1_p),
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p),
            pr2 = sprintf("%.3f", pr2)),
  n = Inf)

cat("\n=== 4. SE variants on month_year_fac (same point estimate) ===\n")
cat("NegBin:\n")
print(vcov_rows |> filter(family == "NegBin") |>
  transmute(variant,
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p)),
  n = Inf)
cat("\nPoisson:\n")
print(vcov_rows |> filter(family == "Poisson") |>
  transmute(variant,
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p)),
  n = Inf)

cat("\n=== 5. Sample restriction on month_year_fac (re-fit NegBin) ===\n")
print(samp_rows |>
  transmute(variant, N,
            b1 = sprintf("%+.3f (%.3f, p=%.3f)", b1, b1_se, b1_p),
            b3 = sprintf("%+.3f (%.3f, p=%.3f)", b3, b3_se, b3_p)),
  n = Inf)

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
