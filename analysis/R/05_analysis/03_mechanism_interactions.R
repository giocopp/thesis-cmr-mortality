# Mechanism: SWH × weekly-lagged SAR moderator (share + log persons), NB + Poisson.
# Produces tab-mechanism.tex.

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

# ── 1. Build daily panel with SAR moderators ────────────────────────────────
panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

iom_daily    <- build_iom_daily()
united_daily <- build_united_daily()

sar_persons_cols <- c(
  "frx_persons_sar_ngo",   "frx_persons_sar_eu",    "frx_persons_sar_ita",
  "frx_persons_sar_comm",  "frx_persons_sar_cg",    "frx_persons_sar_land",
  "frx_persons_sar_noint", "frx_persons_sar_other", "frx_persons_sar_na"
)
panel$frx_persons_sar_total <- rowSums(panel[, sar_persons_cols], na.rm = TRUE)

panel <- panel |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  left_join(united_daily,                                       by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) |>
  arrange(date) |>
  mutate(
    sar_events_pw    = dplyr::lag(zoo::rollsumr(frx_n_sar,             k = 7, fill = NA), 1),
    sar_persons_pw   = dplyr::lag(zoo::rollsumr(frx_persons_sar_total, k = 7, fill = NA), 1),
    incidents_pw     = dplyr::lag(zoo::rollsumr(frx_incidents,         k = 7, fill = NA), 1),

    sar_share_pw     = ifelse(incidents_pw > 0, sar_events_pw / incidents_pw,
                                NA_real_),
    sar_share_pw_z   = (sar_share_pw - mean(sar_share_pw, na.rm = TRUE)) /
                         sd(sar_share_pw, na.rm = TRUE),

    log1p_sar_persons_pw   = log1p(sar_persons_pw),
    log1p_sar_persons_pw_z = (log1p_sar_persons_pw - mean(log1p_sar_persons_pw, na.rm = TRUE)) /
                                sd(log1p_sar_persons_pw, na.rm = TRUE),

    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag14         = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    unit           = 1L,
    month_year_fac = factor(month_year)
  )

d <- panel |> filter(!is.na(sar_share_pw), !is.na(swh_prev5days),
                       !is.na(lc_lag14))

# ── 2. Fit the four mechanism models per source ─────────────────────────────
fit_mech_set <- function(dep, d) {
  nb <- function(rhs) fenegbin(as.formula(sprintf("%s ~ %s", dep, rhs)),
                                data = d, vcov = NW(14), panel.id = ~unit + date)
  po <- function(rhs) fepois  (as.formula(sprintf("%s ~ %s", dep, rhs)),
                                data = d, vcov = NW(14), panel.id = ~unit + date)

  m1 <- "swh_prev5days + swh_prev5days:sar_share_pw_z + sar_share_pw_z | month_year_fac"
  p1 <- "swh_prev5days + swh_prev5days:log1p_sar_persons_pw_z + log1p_sar_persons_pw_z | month_year_fac"

  list(
    my_sar_nb = nb(m1),
    my_sar_po = po(m1),
    my_per_nb = nb(p1),
    my_per_po = po(p1)
  )
}

fits_iom    <- fit_mech_set("n_dead_iom",    d)
fits_united <- fit_mech_set("n_dead_united", d)

# ── 3. LaTeX (tab-mechanism) ────────────────────────────────────────────────
sig_stars <- function(p) {
  ifelse(p < 0.001, "^{***}",
  ifelse(p < 0.01,  "^{**}",
  ifelse(p < 0.05,  "^{*}", "")))
}
fcoef <- function(b, p) sprintf("$%+.3f%s$", b, sig_stars(p))
fse   <- function(se)   sprintf("(%.3f)", se)
fint  <- function(x)    formatC(round(x), format = "d", big.mark = ",")

get_int <- function(m, moderator) {
  ct  <- coeftable(m, vcov = NW(14))
  pat <- paste0("swh_prev5days:", moderator)
  r   <- grep(pat, rownames(ct), fixed = TRUE)
  list(b = ct[r, 1], se = ct[r, 2],
       p = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])),
       n = nobs(m))
}

share_u_nb <- get_int(fits_united$my_sar_nb, "sar_share_pw_z")
share_u_po <- get_int(fits_united$my_sar_po, "sar_share_pw_z")
share_i_nb <- get_int(fits_iom$my_sar_nb,    "sar_share_pw_z")
share_i_po <- get_int(fits_iom$my_sar_po,    "sar_share_pw_z")
pers_u_nb  <- get_int(fits_united$my_per_nb, "log1p_sar_persons_pw_z")
pers_u_po  <- get_int(fits_united$my_per_po, "log1p_sar_persons_pw_z")
pers_i_nb  <- get_int(fits_iom$my_per_nb,    "log1p_sar_persons_pw_z")
pers_i_po  <- get_int(fits_iom$my_per_po,    "log1p_sar_persons_pw_z")

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[H]")
add("\\centering")
add("\\small")
add("\\caption{SWH $\\times$ SAR-proxy moderator interaction. Coefficients are per one standard deviation of the standardized weekly SAR proxy.}")
add("\\label{tab:mechanism}")
add("\\begin{tabular}{lcccc}")
add("\\hline")
add("                                 & \\multicolumn{2}{c}{UNITED} & \\multicolumn{2}{c}{IOM (comparison)} \\\\")
add("                                 & NegBin & Poisson & NegBin & Poisson \\\\")
add("\\hline")
add("\\multicolumn{5}{l}{\\textit{(i) SAR-share moderator}} \\\\")
add("SWH $\\times$ SAR share           & ",
    fcoef(share_u_nb$b, share_u_nb$p), " & ", fcoef(share_u_po$b, share_u_po$p), " & ",
    fcoef(share_i_nb$b, share_i_nb$p), " & ", fcoef(share_i_po$b, share_i_po$p), " \\\\")
add("                                 & ",
    fse(share_u_nb$se), " & ", fse(share_u_po$se), " & ",
    fse(share_i_nb$se), " & ", fse(share_i_po$se), " \\\\")
add("\\multicolumn{5}{l}{\\textit{(ii) SAR-persons moderator}} \\\\")
add("SWH $\\times$ SAR persons         & ",
    fcoef(pers_u_nb$b,  pers_u_nb$p),  " & ", fcoef(pers_u_po$b,  pers_u_po$p),  " & ",
    fcoef(pers_i_nb$b,  pers_i_nb$p),  " & ", fcoef(pers_i_po$b,  pers_i_po$p),  " \\\\")
add("                                 & ",
    fse(pers_u_nb$se),  " & ", fse(pers_u_po$se),  " & ",
    fse(pers_i_nb$se),  " & ", fse(pers_i_po$se),  " \\\\")
add("\\hline")
add("Month--year fixed effects        & Yes            & Yes            & Yes           & Yes           \\\\")
add("Newey--West SEs (lag 14)         & Yes            & Yes            & Yes           & Yes           \\\\")
add("Observations (days)              & ",
    fint(share_u_nb$n), " & ", fint(share_u_po$n), " & ",
    fint(share_i_nb$n), " & ", fint(share_i_po$n), " \\\\")
add("\\hline")
add("\\multicolumn{5}{l}{\\footnotesize Stars denote two-sided $p$-values: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\multicolumn{5}{l}{\\footnotesize Newey--West standard errors in parentheses.} \\\\")
add("\\multicolumn{5}{l}{\\footnotesize SAR share is $\\mathrm{sar\\_share}_{t-7:t-1}$; SAR persons is $\\log(1+\\mathrm{sar\\_persons}_{t-7:t-1})$.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
writeLines(L, tbl_path("05_analysis", "tab-mechanism.tex"))
