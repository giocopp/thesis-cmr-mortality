# Primary count + volume-controlled models: SWH × post-MoU on daily CMR deaths.
# Produces tab-primary, tab-rate, tab-appx-exposure (UNITED + IOM, NegBin + Poisson).

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

# ── 1. Load and prepare panel ───────────────────────────────────────────────
panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

iom_daily    <- build_iom_daily()
united_daily <- build_united_daily()

panel <- panel |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  left_join(united_daily, by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) |>
  add_crossing_exposure() |>
  mutate(
    swh_next1day  = dplyr::lead(swh, 1),
    swh_next3days = zoo::rollmean(dplyr::lead(swh, 1), k = 3,
                                  fill = NA, align = "left"),
    swh_next5days = zoo::rollmean(dplyr::lead(swh, 1), k = 5,
                                  fill = NA, align = "left"),
    swh_next7days = zoo::rollmean(dplyr::lead(swh, 1), k = 7,
                                  fill = NA, align = "left"),
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    unit           = 1L,
    month_year_fac = factor(month_year)
  )

d <- panel |> filter(!is.na(lc_lag14), !is.na(swh_prev5days))

# ── 2. Exposure-window sensitivity grid ─────────────────────────────────────
exposures <- c("swh_lag1", "swh_prev3days", "swh_prev5days", "swh_prevweek",
               "swh_next1day", "swh_next3days", "swh_next5days", "swh_next7days")
window_lookup <- c(
  swh_lag1      = "lag 1d",
  swh_prev3days = "lag 1-3d",
  swh_prev5days = "lag 1-5d",
  swh_prevweek  = "lag 1-7d",
  swh_next1day  = "lead 1d",
  swh_next3days = "lead 1-3d",
  swh_next5days = "lead 1-5d",
  swh_next7days = "lead 1-7d"
)
timing_lookup <- c(
  swh_lag1      = "past",   swh_prev3days = "past",
  swh_prev5days = "past",   swh_prevweek  = "past",
  swh_next1day  = "future_placebo", swh_next3days = "future_placebo",
  swh_next5days = "future_placebo", swh_next7days = "future_placebo"
)
source_lookup <- c(UNITED = "n_dead_united", IOM = "n_dead_iom")

d_sens <- d |> filter(if_all(all_of(exposures), ~ !is.na(.x)))

fit_sens_one <- function(x, outcome, source_label) {
  f <- as.formula(sprintf("%s ~ %s + %s:post_mou | month_year_fac",
                          outcome, x, x))
  models <- list(
    NegBin  = fenegbin(f, data = d_sens, vcov = NW(14), panel.id = ~unit + date),
    Poisson = fepois (f, data = d_sens, vcov = NW(14), panel.id = ~unit + date)
  )
  imap_dfr(models, function(m, fam) {
    ct <- coeftable(m, vcov = NW(14))
    r_main <- which(rownames(ct) == x)
    r_int  <- grep(":post_mou$", rownames(ct))
    r_int  <- r_int[grepl(x, rownames(ct)[r_int], fixed = TRUE)]
    tibble(
      source   = source_label,
      timing   = timing_lookup[[x]],
      window   = window_lookup[[x]],
      family   = fam,
      exposure = x,
      n_obs    = nobs(m),
      b1       = ct[r_main, 1],
      se1      = ct[r_main, 2],
      p1       = 2 * pnorm(-abs(ct[r_main, 1] / ct[r_main, 2])),
      b3       = ct[r_int,  1],
      se3      = ct[r_int,  2],
      p3       = 2 * pnorm(-abs(ct[r_int, 1] / ct[r_int, 2]))
    )
  })
}

sens_tbl <- map_dfr(names(source_lookup), function(src) {
  map_dfr(exposures, function(x) fit_sens_one(x, source_lookup[[src]], src))
}) |>
  mutate(
    source = factor(source, levels = c("UNITED", "IOM")),
    timing = factor(timing, levels = c("past", "future_placebo")),
    exposure_order = match(exposure, exposures)
  ) |>
  arrange(source, timing, exposure_order, family)

# ── 3. Primary count models (NegBin + Poisson, UNITED + IOM) ────────────────
m_nb_u   <- fenegbin(n_dead_united ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                     data = d, vcov = NW(14), panel.id = ~unit + date)
m_pois_u <- fepois  (n_dead_united ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                     data = d, vcov = NW(14), panel.id = ~unit + date)
m_nb     <- fenegbin(n_dead_iom    ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                     data = d, vcov = NW(14), panel.id = ~unit + date)
m_pois   <- fepois  (n_dead_iom    ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
                     data = d, vcov = NW(14), panel.id = ~unit + date)

# ── 4. Volume-controlled Poisson with log(crossing_attempts) ────────────────
d_rate <- d |>
  filter(crossing_attempts > 0) |>
  mutate(log_crossing_attempts = log(crossing_attempts))

m_rate_u <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  log_crossing_attempts | month_year_fac,
  data = d_rate, vcov = NW(14), panel.id = ~unit + date)

m_rate <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               log_crossing_attempts | month_year_fac,
  data = d_rate, vcov = NW(14), panel.id = ~unit + date)

# ── 5. Slope decomposition + elasticity helpers ─────────────────────────────
slope_summary <- function(m, label, family, estimand) {
  ct <- coeftable(m, vcov = NW(14))
  co <- coef(m)
  V  <- vcov(m, vcov = NW(14))
  b_pre <- unname(co["swh_prev5days"])
  b_shift <- unname(co["swh_prev5days:post_mou"])
  se_pre <- ct["swh_prev5days", 2]
  se_shift <- ct["swh_prev5days:post_mou", 2]
  b_post <- b_pre + b_shift
  v_post <- V["swh_prev5days", "swh_prev5days"] +
    V["swh_prev5days:post_mou", "swh_prev5days:post_mou"] +
    2 * V["swh_prev5days", "swh_prev5days:post_mou"]
  se_post <- sqrt(v_post)
  tibble(
    source = label, family = family, estimand = estimand, n_obs = nobs(m),
    b_pre = b_pre, se_pre = se_pre,
    p_pre = 2 * pnorm(-abs(b_pre / se_pre)),
    b_shift = b_shift, se_shift = se_shift,
    p_shift = 2 * pnorm(-abs(b_shift / se_shift)),
    b_post = b_post, se_post = se_post,
    p_post = 2 * pnorm(-abs(b_post / se_post))
  )
}

count_slopes <- bind_rows(
  slope_summary(m_nb_u,   "UNITED", "NegBin",  "count"),
  slope_summary(m_pois_u, "UNITED", "Poisson", "count"),
  slope_summary(m_nb,     "IOM",    "NegBin",  "count"),
  slope_summary(m_pois,   "IOM",    "Poisson", "count")
)

rate_slopes <- bind_rows(
  slope_summary(m_rate_u, "UNITED", "Poisson", "rate-free-exposure"),
  slope_summary(m_rate,   "IOM",    "Poisson", "rate-free-exposure")
)

elast_test <- function(m, xname, label) {
  ct <- coeftable(m, vcov = NW(14))
  b  <- ct[xname, 1]
  se <- ct[xname, 2]
  z  <- (b - 1) / se
  tibble(source = label, b_exposure = b, se = se,
         z_vs_1 = z, p_vs_1 = 2 * pnorm(-abs(z)))
}

elast_tbl <- bind_rows(
  elast_test(m_rate_u, "log_crossing_attempts", "UNITED"),
  elast_test(m_rate,   "log_crossing_attempts", "IOM")
)

# ── 6. LaTeX helpers ────────────────────────────────────────────────────────
sig_stars <- function(p) {
  ifelse(p < 0.001, "^{***}",
  ifelse(p < 0.01,  "^{**}",
  ifelse(p < 0.05,  "^{*}", "")))
}
fcoef <- function(b, p) sprintf("$%+.3f%s$", b, sig_stars(p))
fse   <- function(se)   sprintf("(%.3f)", se)
fint  <- function(x)    formatC(round(x), format = "d", big.mark = ",")
fr3   <- function(x)    formatC(x, format = "f", digits = 3)
fp    <- function(p)    ifelse(p < 0.001, "$< 0.001$", sprintf("$%.3f$", p))

pr2_nb_u   <- fixest::r2(m_nb_u,   "pr2")
pr2_pois_u <- fixest::r2(m_pois_u, "pr2")
pr2_nb     <- fixest::r2(m_nb,     "pr2")
pr2_pois   <- fixest::r2(m_pois,   "pr2")

# ── 7. tab-primary ──────────────────────────────────────────────────────────
cs <- count_slopes
L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[h!]")
add("\\centering")
add("\\small")
add("\\caption{Primary count estimates of the SWH--mortality slope shift around the 2 February 2017 MoU.}")
add("\\label{tab:primary}")
add("\\begin{tabular}{lcccc}")
add("\\hline")
add("                              & \\multicolumn{2}{c}{UNITED} & \\multicolumn{2}{c}{IOM (comparison)} \\\\")
add("                              & NegBin & Poisson & NegBin & Poisson \\\\")
add("\\hline")
add("$\\beta_1$: SWH$_{t-1:t-5}$            & ",
    fcoef(cs$b_pre[1],   cs$p_pre[1]),   " & ", fcoef(cs$b_pre[2],   cs$p_pre[2]),   " & ",
    fcoef(cs$b_pre[3],   cs$p_pre[3]),   " & ", fcoef(cs$b_pre[4],   cs$p_pre[4]),   " \\\\")
add("                                      & ",
    fse(cs$se_pre[1]),   " & ", fse(cs$se_pre[2]),   " & ",
    fse(cs$se_pre[3]),   " & ", fse(cs$se_pre[4]),   " \\\\")
add("$\\beta_3$: SWH$_{t-1:t-5}$ $\\times$ Post-MoU & ",
    fcoef(cs$b_shift[1], cs$p_shift[1]), " & ", fcoef(cs$b_shift[2], cs$p_shift[2]), " & ",
    fcoef(cs$b_shift[3], cs$p_shift[3]), " & ", fcoef(cs$b_shift[4], cs$p_shift[4]), " \\\\")
add("                                      & ",
    fse(cs$se_shift[1]), " & ", fse(cs$se_shift[2]), " & ",
    fse(cs$se_shift[3]), " & ", fse(cs$se_shift[4]), " \\\\")
add("$\\beta_1+\\beta_3$: implied post-MoU slope    & ",
    fcoef(cs$b_post[1],  cs$p_post[1]),  " & ", fcoef(cs$b_post[2],  cs$p_post[2]),  " & ",
    fcoef(cs$b_post[3],  cs$p_post[3]),  " & ", fcoef(cs$b_post[4],  cs$p_post[4]),  " \\\\")
add("                                      & ",
    fse(cs$se_post[1]),  " & ", fse(cs$se_post[2]),  " & ",
    fse(cs$se_post[3]),  " & ", fse(cs$se_post[4]),  " \\\\")
add("\\hline")
add("Month--year fixed effects             & Yes            & Yes           & Yes            & Yes          \\\\")
add("Newey--West SEs (lag 14)              & Yes            & Yes           & Yes            & Yes          \\\\")
add("Observations (days)                   & ",
    fint(cs$n_obs[1]), " & ", fint(cs$n_obs[2]), " & ",
    fint(cs$n_obs[3]), " & ", fint(cs$n_obs[4]), " \\\\")
add("Pseudo $R^{2}$                        & ",
    fr3(pr2_nb_u),   " & ", fr3(pr2_pois_u), " & ",
    fr3(pr2_nb),     " & ", fr3(pr2_pois),   " \\\\")
add("\\hline")
add("\\multicolumn{5}{l}{\\footnotesize Stars denote two-sided $p$-values: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\multicolumn{5}{l}{\\footnotesize Newey--West standard errors in parentheses; delta-method SEs for $\\beta_1+\\beta_3$.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
writeLines(L, tbl_path("05_analysis", "tab-primary.tex"))

# ── 8. tab-rate ─────────────────────────────────────────────────────────────
rs <- rate_slopes
ct_u <- coeftable(m_rate_u, vcov = NW(14))
ct_i <- coeftable(m_rate,   vcov = NW(14))
g_u  <- ct_u["log_crossing_attempts", 1]; g_se_u <- ct_u["log_crossing_attempts", 2]
g_i  <- ct_i["log_crossing_attempts", 1]; g_se_i <- ct_i["log_crossing_attempts", 2]
g_p_u <- 2 * pnorm(-abs(g_u  / g_se_u))
g_p_i <- 2 * pnorm(-abs(g_i  / g_se_i))

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[h!]")
add("\\centering")
add("\\small")
add("\\caption{Volume-controlled Poisson model with crossing-volume control.}")
add("\\label{tab:rate}")
add("\\begin{tabular}{lcc}")
add("\\hline")
add("                              & UNITED & IOM (comparison) \\\\")
add("\\hline")
add("$\\beta_1$: SWH$_{t-1:t-5}$            & ",
    fcoef(rs$b_pre[1],   rs$p_pre[1]),   " & ", fcoef(rs$b_pre[2],   rs$p_pre[2]),   " \\\\")
add("                                      & ",
    fse(rs$se_pre[1]),   " & ", fse(rs$se_pre[2]),   " \\\\")
add("$\\beta_3$: SWH$_{t-1:t-5}$ $\\times$ Post-MoU & ",
    fcoef(rs$b_shift[1], rs$p_shift[1]), " & ", fcoef(rs$b_shift[2], rs$p_shift[2]), " \\\\")
add("                                      & ",
    fse(rs$se_shift[1]), " & ", fse(rs$se_shift[2]), " \\\\")
add("$\\gamma$: $\\log C_t$                  & ",
    fcoef(g_u, g_p_u), " & ", fcoef(g_i, g_p_i), " \\\\")
add("                                      & ",
    fse(g_se_u), " & ", fse(g_se_i), " \\\\")
add("$\\beta_1+\\beta_3$: implied post-MoU slope    & ",
    fcoef(rs$b_post[1],  rs$p_post[1]),  " & ", fcoef(rs$b_post[2],  rs$p_post[2]),  " \\\\")
add("                                      & ",
    fse(rs$se_post[1]),  " & ", fse(rs$se_post[2]),  " \\\\")
add("\\hline")
add("Wald: $H_{0}: \\gamma = 1$, $z$        & ",
    sprintf("$%+.2f$", elast_tbl$z_vs_1[1]), " & ",
    sprintf("$%+.2f$", elast_tbl$z_vs_1[2]), " \\\\")
add("$p$-value                             & ", fp(elast_tbl$p_vs_1[1]),
    " & ", fp(elast_tbl$p_vs_1[2]), " \\\\")
add("Month--year fixed effects             & Yes            & Yes            \\\\")
add("Newey--West SEs (lag 14)              & Yes            & Yes            \\\\")
add("Observations (days, $C_t>0$)          & ",
    fint(rs$n_obs[1]), " & ", fint(rs$n_obs[2]), " \\\\")
add("\\hline")
add("\\multicolumn{3}{l}{\\footnotesize Stars denote two-sided $p$-values: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\multicolumn{3}{l}{\\footnotesize Newey--West standard errors in parentheses; delta-method SEs for $\\beta_1+\\beta_3$.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
writeLines(L, tbl_path("05_analysis", "tab-rate.tex"))

# ── 9. tab-appx-exposure ────────────────────────────────────────────────────
window_order <- c("lag 1d", "lag 1-3d", "lag 1-5d", "lag 1-7d",
                  "lead 1d", "lead 1-3d", "lead 1-5d", "lead 1-7d")
sens_wide <- sens_tbl |>
  dplyr::select(source, timing, window, family, b3, se3, p3) |>
  tidyr::pivot_wider(names_from = family,
                     values_from = c(b3, se3, p3)) |>
  dplyr::arrange(source, timing, match(window, window_order))

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[h!]")
add("\\centering")
add("\\small")
add("\\caption{SWH-window grid: SWH$\\times$Post-MoU shift coefficient ($\\beta_3$).}")
add("\\label{tab:appx-exposure}")
add("\\begin{tabular}{llcccc}")
add("\\hline")
add("Source & Window      & \\multicolumn{2}{c}{NegBin} & \\multicolumn{2}{c}{Poisson} \\\\")
add("       &             & $\\beta_3$ & SE    & $\\beta_3$ & SE \\\\")
add("\\hline")
write_block <- function(timing_label, header) {
  add(sprintf("\\multicolumn{6}{l}{\\textit{%s}} \\\\", header))
  rows <- sens_wide[sens_wide$timing == timing_label, ]
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    add(sprintf("%-7s & %-11s & %-15s & %.3f & %-15s & %.3f \\\\",
                as.character(r$source), r$window,
                fcoef(r$b3_NegBin,  r$p3_NegBin),  r$se3_NegBin,
                fcoef(r$b3_Poisson, r$p3_Poisson), r$se3_Poisson))
  }
}
write_block("past",           "Past (lagged) windows")
write_block("future_placebo", "Future (placebo) windows")
add("\\hline")
add("\\multicolumn{6}{l}{\\footnotesize SE columns report Newey--West standard errors. Stars: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add(sprintf("\\multicolumn{6}{l}{\\footnotesize Shared sample: %s days.}",
            fint(nrow(d_sens))))
add("\\end{tabular}")
add("\\end{table}")
writeLines(L, tbl_path("05_analysis", "tab-appx-exposure.tex"))
