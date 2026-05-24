# 3-period regime decomposition of SWH-mortality slope on the full ERA5 span
# (no crossing-volume control), clipped at the eve of the Piantedosi decree.
# Produces tab-period.tex.

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

D_12        <- PERIOD_END_1   # 2017-02-01
D_23        <- PERIOD_END_2   # 2020-10-20
SAMPLE_END  <- PERIOD_END_3   # 2023-01-01

# ── 1. Build wider source-specific samples on the full ERA5 span ────────────
iom_daily    <- build_iom_daily()
united_daily <- build_united_daily()

swh_daily_full <- readRDS(file.path(BASE_DIR, "data", "processed",
                                     "era5_swh_daily.RDS"))
dim(swh_daily_full$date) <- NULL

build_period_sample <- function(deaths_daily, count_col, date_floor) {
  swh_daily_full |>
    dplyr::select(date, swh_prev5days) |>
    dplyr::left_join(deaths_daily |> dplyr::select(date, dplyr::all_of(count_col)),
                     by = "date") |>
    tidyr::replace_na(setNames(list(0), count_col)) |>
    dplyr::filter(date >= date_floor, date <= SAMPLE_END,
                  !is.na(swh_prev5days)) |>
    dplyr::mutate(
      unit           = 1L,
      month_year_fac = factor(format(date, "%Y-%m")),
      period = factor(dplyr::case_when(
        date <= D_12 ~ "1. SAR + border control",
        date <= D_23 ~ "2. MoU + NGO containment",
        TRUE         ~ "3. Lamorgese rollback"
      ))
    )
}

d_utd_3p <- build_period_sample(united_daily, "n_dead_united",
                                 as.Date("2013-01-01"))
d_iom_3p <- build_period_sample(iom_daily |> rename(n_dead_iom = n_dead_missing),
                                 "n_dead_iom",
                                 as.Date("2014-01-01"))

# ── 2. NegBin: swh × period interaction ─────────────────────────────────────
f1_3p_iom_nb <- fenegbin(n_dead_iom ~ swh_prev5days:period | month_year_fac,
                         data = d_iom_3p, vcov = NW(14), panel.id = ~unit + date)
f1_3p_utd_nb <- fenegbin(n_dead_united ~ swh_prev5days:period | month_year_fac,
                         data = d_utd_3p, vcov = NW(14), panel.id = ~unit + date)

# ── 3. Wald test: all 3 period gradients equal ──────────────────────────────
wald_periods <- function(m) {
  co <- coef(m); V <- vcov(m, vcov = NW(14))
  idx <- grep("swh_prev5days:period", names(co))
  b   <- co[idx]; Vs <- V[idx, idx, drop = FALSE]; k <- length(b)
  R <- matrix(0, nrow = k - 1, ncol = k)
  for (j in seq_len(k - 1)) { R[j, 1] <- -1; R[j, j + 1] <- 1 }
  Rb  <- R %*% b
  RVR <- R %*% Vs %*% t(R)
  stat <- as.numeric(t(Rb) %*% solve(RVR) %*% Rb)
  p    <- pchisq(stat, df = k - 1, lower.tail = FALSE)
  list(chi2 = stat, df = k - 1, p = p)
}

w_u <- wald_periods(f1_3p_utd_nb)
w_i <- wald_periods(f1_3p_iom_nb)

# ── 4. LaTeX (tab-period) ───────────────────────────────────────────────────
sig_stars <- function(p) {
  ifelse(p < 0.001, "^{***}",
  ifelse(p < 0.01,  "^{**}",
  ifelse(p < 0.05,  "^{*}", "")))
}
fcoef <- function(b, p) sprintf("$%+.3f%s$", b, sig_stars(p))
fse   <- function(se)   sprintf("(%.3f)", se)
fint  <- function(x)    formatC(round(x), format = "d", big.mark = ",")
fp    <- function(p)    ifelse(p < 0.001, "$< 0.001$", sprintf("$%.3f$", p))

get_period <- function(m, p_label) {
  ct <- coeftable(m, vcov = NW(14))
  rn <- grep(paste0("swh_prev5days:period", p_label), rownames(ct), fixed = TRUE)
  list(b = ct[rn, 1], se = ct[rn, 2],
       p = 2 * pnorm(-abs(ct[rn, 1] / ct[rn, 2])))
}
periods_labels <- c("1. SAR + border control",
                    "2. MoU + NGO containment",
                    "3. Lamorgese rollback")
u <- lapply(periods_labels, function(p) get_period(f1_3p_utd_nb, p))
i <- lapply(periods_labels, function(p) get_period(f1_3p_iom_nb, p))

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[h!]")
add("\\centering")
add("\\small")
add("\\caption{Period-specific SWH--mortality gradient. NegBin with month--year FE and NW(14) SEs. UNITED is fit from 2013-01-01 to take advantage of its longer span; IOM is fit from 2014-01-01, when MMP coverage begins. Both samples end on the eve of the January 2023 Piantedosi decree.}")
add("\\label{tab:period}")
add("\\begin{tabular}{lcc}")
add("\\hline")
add(sprintf("SAR/policy regime                                & UNITED (%d--%d) & IOM (%d--%d) \\\\",
            year(min(d_utd_3p$date)), year(max(d_utd_3p$date)),
            year(min(d_iom_3p$date)), year(max(d_iom_3p$date))))
add("\\hline")
add("1. SAR + border control (through Feb 2017)       & ",
    fcoef(u[[1]]$b, u[[1]]$p), " & ", fcoef(i[[1]]$b, i[[1]]$p), " \\\\")
add("                                                 & ",
    fse(u[[1]]$se), " & ", fse(i[[1]]$se), " \\\\")
add("2. MoU + NGO containment (Feb 2017--Oct 2020)    & ",
    fcoef(u[[2]]$b, u[[2]]$p), " & ", fcoef(i[[2]]$b, i[[2]]$p), " \\\\")
add("                                                 & ",
    fse(u[[2]]$se), " & ", fse(i[[2]]$se), " \\\\")
add("3. Lamorgese rollback (Oct 2020--Jan 2023)       & ",
    fcoef(u[[3]]$b, u[[3]]$p), " & ", fcoef(i[[3]]$b, i[[3]]$p), " \\\\")
add("                                                 & ",
    fse(u[[3]]$se), " & ", fse(i[[3]]$se), " \\\\")
add("\\hline")
add("Wald $H_{0}$: all equal, $\\chi^{2}(2)$           & ",
    sprintf("$%.2f$", w_u$chi2), " & ", sprintf("$%.2f$", w_i$chi2), " \\\\")
add("$p$-value                                        & ",
    fp(w_u$p), " & ", fp(w_i$p), " \\\\")
add("Observations (days)                              & ",
    fint(nrow(d_utd_3p)), " & ", fint(nrow(d_iom_3p)), " \\\\")
add("\\hline")
add("\\multicolumn{3}{l}{\\footnotesize Stars denote two-sided $p$-values: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\multicolumn{3}{l}{\\footnotesize Newey--West standard errors in parentheses.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
writeLines(L, tbl_path("05_analysis", "tab-period.tex"))
