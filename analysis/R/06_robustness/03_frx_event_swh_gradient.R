# Frontex event-count SWH gradient (SAR vs Not-SAR placebo). Produces
# tab-appx-frx-events.tex.

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

# ── 1. Build daily SAR / non-SAR event counts ───────────────────────────────
panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

panel <- panel |>
  mutate(
    frx_n_notsar   = frx_incidents - frx_n_sar,
    unit           = 1L,
    month_year_fac = factor(month_year)
  )

d <- panel |> filter(!is.na(swh_prev5days))

# ── 2. NegBin + Poisson on SAR and Non-SAR counts ───────────────────────────
m_nb_sar   <- fenegbin(frx_n_sar    ~ swh_prev5days + swh_prev5days:post_mou |
                         month_year_fac, data = d, vcov = NW(14),
                       panel.id = ~unit + date)
m_pois_sar <- fepois  (frx_n_sar    ~ swh_prev5days + swh_prev5days:post_mou |
                         month_year_fac, data = d, vcov = NW(14),
                       panel.id = ~unit + date)
m_nb_not   <- fenegbin(frx_n_notsar ~ swh_prev5days + swh_prev5days:post_mou |
                         month_year_fac, data = d, vcov = NW(14),
                       panel.id = ~unit + date)
m_pois_not <- fepois  (frx_n_notsar ~ swh_prev5days + swh_prev5days:post_mou |
                         month_year_fac, data = d, vcov = NW(14),
                       panel.id = ~unit + date)

# ── 3. LaTeX (tab-appx-frx-events) ──────────────────────────────────────────
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
  r <- which(rownames(ct) == name)
  if (length(r) == 0) {
    r <- grep(name, rownames(ct), fixed = TRUE)
  }
  list(b = ct[r, 1], se = ct[r, 2],
       p = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

models <- list(m_nb_sar, m_pois_sar, m_nb_not, m_pois_not)
b1 <- lapply(models, get_coef, name = "swh_prev5days")
b3 <- lapply(models, get_coef, name = "swh_prev5days:post_mou")
Ns <- sapply(models, nobs)

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[H]")
add("\\centering")
add("\\small")
add("\\caption{SWH$\\times$Post-MoU shift on Frontex daily event counts.}")
add("\\label{tab:appx-frx-events}")
add("\\begin{tabular}{lcccc}")
add("\\hline")
add("                              & \\multicolumn{2}{c}{SAR events} & \\multicolumn{2}{c}{Non-SAR (placebo)} \\\\")
add("                              & NegBin & Poisson & NegBin & Poisson \\\\")
add("\\hline")
add("SWH$_{t-1:t-5}$               & ",
    fcoef(b1[[1]]$b, b1[[1]]$p), " & ", fcoef(b1[[2]]$b, b1[[2]]$p), " & ",
    fcoef(b1[[3]]$b, b1[[3]]$p), " & ", fcoef(b1[[4]]$b, b1[[4]]$p), " \\\\")
add("                              & ",
    fse(b1[[1]]$se), " & ", fse(b1[[2]]$se), " & ",
    fse(b1[[3]]$se), " & ", fse(b1[[4]]$se), " \\\\")
add("SWH $\\times$ Post-MoU         & ",
    fcoef(b3[[1]]$b, b3[[1]]$p), " & ", fcoef(b3[[2]]$b, b3[[2]]$p), " & ",
    fcoef(b3[[3]]$b, b3[[3]]$p), " & ", fcoef(b3[[4]]$b, b3[[4]]$p), " \\\\")
add("                              & ",
    fse(b3[[1]]$se), " & ", fse(b3[[2]]$se), " & ",
    fse(b3[[3]]$se), " & ", fse(b3[[4]]$se), " \\\\")
add("\\hline")
add("Month-year FE                 & Yes            & Yes             & Yes      & Yes      \\\\")
add("Newey-West SEs (lag 14)       & Yes            & Yes             & Yes      & Yes      \\\\")
add("Observations (days)           & ",
    fint(Ns[1]), " & ", fint(Ns[2]), " & ", fint(Ns[3]), " & ", fint(Ns[4]), " \\\\")
add("\\hline")
add("\\multicolumn{5}{l}{\\footnotesize Newey--West standard errors in parentheses. Stars: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
writeLines(L, tbl_path("06_robustness", "tab-appx-frx-events.tex"))
