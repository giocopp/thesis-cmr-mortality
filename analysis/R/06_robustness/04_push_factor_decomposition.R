# Push-factor robustness: triple interaction SWH × Post-MoU × Libya conflict.
# Produces tab-appx-pushfactor.tex.

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

# ── 1. Build sample with ACLED lag-1w moderator ─────────────────────────────
panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

iom_daily    <- build_iom_daily()
united_daily <- build_united_daily()

panel <- panel |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  left_join(united_daily,                                       by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0)) |>
  arrange(date) |>
  mutate(
    libya_conflict_lag1w = dplyr::lag(libya_conflict, 7),
    living_crossings     = frx_persons + lcg_tcg_pushbacks,
    lc_lag14 = dplyr::lag(
      zoo::rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8),
    unit           = 1L,
    month_year_fac = factor(month_year)
  )

d <- panel |>
  filter(!is.na(lc_lag14), !is.na(swh_prev5days),
         !is.na(libya_conflict_lag1w)) |>
  mutate(libya_conflict_z = (libya_conflict_lag1w - mean(libya_conflict_lag1w)) /
                              sd(libya_conflict_lag1w))

# ── 2. Triple-interaction count models ──────────────────────────────────────
mod_z       <- "libya_conflict_z"
triple_term <- sprintf("swh_prev5days:post_mou:%s", mod_z)

l1_terms <- c(
  "swh_prev5days",
  "swh_prev5days:post_mou",
  mod_z,
  sprintf("%s:post_mou",                mod_z),
  sprintf("swh_prev5days:%s",           mod_z),
  sprintf("swh_prev5days:post_mou:%s",  mod_z)
)
rhs_l1 <- paste(paste(l1_terms, collapse = " + "), "| month_year_fac")

fit_one <- function(dep, family) {
  fn  <- if (family == "NegBin") fenegbin else fepois
  fn(as.formula(sprintf("%s ~ %s", dep, rhs_l1)),
     data = d, vcov = NW(14), panel.id = ~unit + date)
}

l1 <- list(
  united_nb = fit_one("n_dead_united", "NegBin"),
  united_po = fit_one("n_dead_united", "Poisson"),
  iom_nb    = fit_one("n_dead_iom",    "NegBin"),
  iom_po    = fit_one("n_dead_iom",    "Poisson")
)

# ── 3. LaTeX (tab-appx-pushfactor) ──────────────────────────────────────────
extract_coef <- function(m, term) {
  ct <- coeftable(m, vcov = NW(14))
  rn <- rownames(ct)
  r  <- which(rn == term)
  if (length(r) == 0) {
    parts <- strsplit(term, ":")[[1]]
    if (length(parts) == 2) r <- which(rn == paste(parts[2], parts[1], sep = ":"))
  }
  if (length(r) == 0) return(list(b = NA_real_, se = NA_real_, p = NA_real_))
  list(b = ct[r, 1], se = ct[r, 2],
       p = 2 * pnorm(-abs(ct[r, 1] / ct[r, 2])))
}

sig_stars <- function(p) {
  ifelse(p < 0.001, "^{***}",
  ifelse(p < 0.01,  "^{**}",
  ifelse(p < 0.05,  "^{*}", "")))
}
fcoef <- function(b, p) sprintf("$%+.3f%s$", b, sig_stars(p))
fse   <- function(se)   sprintf("(%.3f)", se)
fint  <- function(x)    formatC(round(x), format = "d", big.mark = ",")

models <- list(l1$united_nb, l1$united_po, l1$iom_nb, l1$iom_po)
b1 <- lapply(models, extract_coef, term = "swh_prev5days")
b3 <- lapply(models, extract_coef, term = "swh_prev5days:post_mou")
bt <- lapply(models, extract_coef, term = triple_term)
Ns <- sapply(models, nobs)

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[H]")
add("\\centering")
add("\\small")
add("\\caption{Triple interaction SWH$\\times$Post-MoU$\\times$Libya conflict.}")
add("\\label{tab:appx-pushfactor}")
add("\\begin{tabular}{lcccc}")
add("\\hline")
add("                                          & \\multicolumn{2}{c}{UNITED}    & \\multicolumn{2}{c}{IOM} \\\\")
add("                                          & NegBin & Poisson & NegBin & Poisson \\\\")
add("\\hline")
add("SWH$_{t-1:t-5}$                           & ",
    fcoef(b1[[1]]$b, b1[[1]]$p), " & ", fcoef(b1[[2]]$b, b1[[2]]$p), " & ",
    fcoef(b1[[3]]$b, b1[[3]]$p), " & ", fcoef(b1[[4]]$b, b1[[4]]$p), " \\\\")
add("                                          & ",
    fse(b1[[1]]$se), " & ", fse(b1[[2]]$se), " & ",
    fse(b1[[3]]$se), " & ", fse(b1[[4]]$se), " \\\\")
add("SWH $\\times$ Post-MoU                     & ",
    fcoef(b3[[1]]$b, b3[[1]]$p), " & ", fcoef(b3[[2]]$b, b3[[2]]$p), " & ",
    fcoef(b3[[3]]$b, b3[[3]]$p), " & ", fcoef(b3[[4]]$b, b3[[4]]$p), " \\\\")
add("                                          & ",
    fse(b3[[1]]$se), " & ", fse(b3[[2]]$se), " & ",
    fse(b3[[3]]$se), " & ", fse(b3[[4]]$se), " \\\\")
add("SWH $\\times$ Post-MoU $\\times$ Libya conflict & ",
    fcoef(bt[[1]]$b, bt[[1]]$p), " & ", fcoef(bt[[2]]$b, bt[[2]]$p), " & ",
    fcoef(bt[[3]]$b, bt[[3]]$p), " & ", fcoef(bt[[4]]$b, bt[[4]]$p), " \\\\")
add("                                              & ",
    fse(bt[[1]]$se), " & ", fse(bt[[2]]$se), " & ",
    fse(bt[[3]]$se), " & ", fse(bt[[4]]$se), " \\\\")
add("\\hline")
add("Month-year FE                             & Yes            & Yes           & Yes            & Yes           \\\\")
add("Newey-West SEs (lag 14)                   & Yes            & Yes           & Yes            & Yes           \\\\")
add("Observations (days)                       & ",
    fint(Ns[1]), " & ", fint(Ns[2]), " & ", fint(Ns[3]), " & ", fint(Ns[4]), " \\\\")
add("\\hline")
add("\\multicolumn{5}{l}{\\footnotesize Libya conflict is the standardized")
add("one-week lag of ACLED battles + explosions + VAC.} \\\\")
add("\\multicolumn{5}{l}{\\footnotesize Newey--West standard errors in parentheses. Stars: $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.} \\\\")
add("\\end{tabular}")
add("\\end{table}")
writeLines(L, tbl_path("06_robustness", "tab-appx-pushfactor.tex"))
