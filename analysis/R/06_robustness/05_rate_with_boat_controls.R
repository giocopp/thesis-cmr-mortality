# Volume-controlled Poisson with boat-composition controls (additive + SWH×inflatable).
# Produces tab-appx-boat.tex.

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

# ── 1. Build boat-observable rate sample ────────────────────────────────────
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
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    unit           = 1L,
    month_year_fac = factor(month_year)
  )

d_rate_boat <- panel |>
  filter(!is.na(lc_lag14), !is.na(swh_prev5days),
         crossing_attempts > 0,
         frx_incidents > 0) |>
  mutate(log_crossing_attempts = log(crossing_attempts))

# ── 2. V1, V2, V3 for IOM and UNITED ────────────────────────────────────────
v1_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

v2_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               frx_inflatable_share + frx_wooden_share +
               log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

v3_iom <- fepois(
  n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou +
               frx_inflatable_share + frx_wooden_share +
               swh_prev5days:frx_inflatable_share +
               log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

v1_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

v2_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  frx_inflatable_share + frx_wooden_share +
                  log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

v3_united <- fepois(
  n_dead_united ~ swh_prev5days + swh_prev5days:post_mou +
                  frx_inflatable_share + frx_wooden_share +
                  swh_prev5days:frx_inflatable_share +
                  log_crossing_attempts | month_year_fac,
  data = d_rate_boat, vcov = NW(14), panel.id = ~unit + date)

# ── 3. LaTeX (tab-appx-boat) ────────────────────────────────────────────────
sig_stars <- function(p) {
  ifelse(p < 0.001, "^{***}",
  ifelse(p < 0.01,  "^{**}",
  ifelse(p < 0.05,  "^{*}", "")))
}
fcoef <- function(b, p) sprintf("$%+.3f%s$", b, sig_stars(p))
fse   <- function(se)   sprintf("(%.3f)", se)
fint  <- function(x)    formatC(round(x), format = "d", big.mark = ",")

get_coef <- function(m, name) {
  ct <- coeftable(m, vcov = NW(14))
  if (!(name %in% rownames(ct))) return(list(present = FALSE))
  r <- ct[name, ]
  list(present = TRUE, b = r[1], se = r[2],
       p = 2 * pnorm(-abs(r[1] / r[2])))
}

N_utd <- nobs(v1_united)
N_iom <- nobs(v1_iom)

row_triple <- function(label, models, name) {
  cs <- lapply(models, get_coef, name = name)
  bs <- sapply(cs, function(c) if (c$present) fcoef(c$b, c$p) else "")
  ss <- sapply(cs, function(c) if (c$present) fse(c$se)     else "")
  list(
    sprintf("%-31s & %-14s & %-14s & %-14s \\\\", label, bs[1], bs[2], bs[3]),
    sprintf("%-31s & %-14s & %-14s & %-14s \\\\", "",    ss[1], ss[2], ss[3])
  )
}

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[H]")
add("\\centering")
add("\\small")
add("\\caption{Volume-controlled Poisson with boat-composition controls. Boat-observable")
add(sprintf("sample: %s days (UNITED), %s days (IOM).}", fint(N_utd), fint(N_iom)))
add("\\label{tab:appx-boat}")
add("\\begin{tabular}{lccc}")
add("\\hline")
add("                                & V1: baseline & V2: + boat shares & V3: + boat $\\times$ SWH \\\\")
add("\\hline")
add("\\multicolumn{4}{l}{\\textit{UNITED}} \\\\")
models_u <- list(v1_united, v2_united, v3_united)
for (line in row_triple("SWH$_{t-1:t-5}$",                models_u, "swh_prev5days"))            add(line)
for (line in row_triple("SWH $\\times$ Post-MoU",         models_u, "swh_prev5days:post_mou"))   add(line)
for (line in row_triple("$\\log C_t$",                     models_u, "log_crossing_attempts"))    add(line)
for (line in row_triple("SWH $\\times$ inflatable share", models_u, "swh_prev5days:frx_inflatable_share")) add(line)
add("\\multicolumn{4}{l}{\\textit{IOM (comparison)}} \\\\")
models_i <- list(v1_iom, v2_iom, v3_iom)
for (line in row_triple("SWH$_{t-1:t-5}$",                models_i, "swh_prev5days"))            add(line)
for (line in row_triple("SWH $\\times$ Post-MoU",         models_i, "swh_prev5days:post_mou"))   add(line)
for (line in row_triple("$\\log C_t$",                     models_i, "log_crossing_attempts"))    add(line)
for (line in row_triple("SWH $\\times$ inflatable share", models_i, "swh_prev5days:frx_inflatable_share")) add(line)
add("\\hline")
add("Month-year FE                   & Yes            & Yes            & Yes            \\\\")
add("Newey-West SEs (lag 14)         & Yes            & Yes            & Yes            \\\\")
add("\\hline")
add("\\multicolumn{4}{l}{\\footnotesize Newey--West standard errors in parentheses. Stars: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
writeLines(L, tbl_path("06_robustness", "tab-appx-boat.tex"))
