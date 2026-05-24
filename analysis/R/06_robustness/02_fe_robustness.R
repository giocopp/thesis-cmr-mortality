# FE / variance / sample robustness of the primary reduced form (IOM NegBin).
# Produces tab-appx-fe.tex.

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

# ── 1. Build sample to match the primary model ──────────────────────────────
panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))
iom_daily <- build_iom_daily()

panel <- panel |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  replace_na(list(n_dead_iom = 0)) |>
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag14 = dplyr::lag(zoo::rollmeanr(living_crossings, k = 7,
                                           fill = NA, align = "right"), 8),
    unit           = 1L,
    year_fac       = factor(year),
    month_year_fac = factor(month_year),
    month_of_year  = factor(month(date)),
    quarter_year   = factor(paste0(year(date), "Q", quarter(date))),
    dow            = factor(wday(date, week_start = 1))
  )

d <- panel |> filter(!is.na(lc_lag14), !is.na(swh_prev5days))

# ── 2. SWH variance absorbed by each candidate FE ───────────────────────────
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

# ── 3. NegBin fits across FE specs ──────────────────────────────────────────
fit_fam <- function(fn, fe_rhs) {
  fml <- as.formula(paste(
    "n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou |", fe_rhs
  ))
  fn(fml, data = d, vcov = NW(14), panel.id = ~unit + date)
}

models_nb <- map(fe_specs, \(fe) fit_fam(fenegbin, fe))

extract_row <- function(m, nm) {
  ct <- coeftable(m, vcov = NW(14))
  rn <- rownames(ct)
  b1_row <- which(rn == "swh_prev5days")
  b3_row <- grep(":post_mou$", rn)
  get <- function(r) if (length(r) == 1) unname(ct[r, 1:2]) else c(NA_real_, NA_real_)
  v1 <- get(b1_row); v3 <- get(b3_row)
  tibble(
    fe_spec = nm, family = "NegBin",
    b1 = v1[1], b1_se = v1[2],
    b3 = v3[1], b3_se = v3[2]
  )
}

results <- imap_dfr(models_nb, \(m, nm) extract_row(m, nm)) |>
  left_join(absorbed |> select(fe_spec, abs_r2), by = "fe_spec") |>
  mutate(b3_p = 2 * pnorm(-abs(b3 / b3_se)))

# ── 4. SE variants on the primary FE ────────────────────────────────────────
m_nb_base <- fit_fam(fenegbin, "month_year_fac")

vcv_variants <- list(
  "NW(7)"                  = NW(7),
  "NW(14) [primary]"       = NW(14),
  "NW(21)"                 = NW(21),
  "cluster: month_year"    = ~ month_year_fac,
  "cluster: year"          = ~ year_fac
)

vcov_rows <- map_dfr(names(vcv_variants), function(nm) {
  ct <- coeftable(m_nb_base, vcov = vcv_variants[[nm]])
  rn <- rownames(ct)
  b3_r <- grep(":post_mou$", rn)
  tibble(family = "NegBin", variant = nm,
         b3 = ct[b3_r, 1], b3_se = ct[b3_r, 2],
         b3_p = 2 * pnorm(-abs(ct[b3_r, 1] / ct[b3_r, 2])))
})

# ── 5. Sample-restriction variants on the primary FE ────────────────────────
sample_variants <- list(
  "full"                 = d,
  "drop zero-death days" = d |> filter(n_dead_iom > 0),
  "cap deaths at 100"    = d |> filter(n_dead_iom <= 100)
)

samp_rows <- imap_dfr(sample_variants, function(dd, nm) {
  dd <- dd |> mutate(month_year_fac = droplevels(factor(month_year_fac)))
  m <- tryCatch(
    fenegbin(n_dead_iom ~ swh_prev5days + swh_prev5days:post_mou | month_year_fac,
             data = dd, vcov = NW(14), panel.id = ~unit + date),
    error = function(e) NULL)
  if (is.null(m)) return(tibble(variant = nm, N = nrow(dd),
                                b3 = NA, b3_se = NA, b3_p = NA))
  ct <- coeftable(m, vcov = NW(14))
  rn <- rownames(ct)
  b3_r <- grep(":post_mou$", rn)
  tibble(
    variant = nm, N = nobs(m),
    b3 = ct[b3_r, 1], b3_se = ct[b3_r, 2],
    b3_p = 2 * pnorm(-abs(ct[b3_r, 1] / ct[b3_r, 2]))
  )
})

# ── 6. LaTeX (tab-appx-fe) ──────────────────────────────────────────────────
fb  <- function(b) sprintf("$%+.3f$", b)
fse <- function(s) sprintf("%.3f", s)
fp  <- function(p) sprintf("%.3f", p)

fe_rows <- list(
  list(label = "year",                       script = "year"),
  list(label = "month-of-year",              script = "month-of-year"),
  list(label = "year + month",               script = "year + month"),
  list(label = "quarter-year",               script = "quarter-year"),
  list(label = "month-year (primary)",       script = "month-year (primary/20)"),
  list(label = "month-year + day-of-week",   script = "month-year + dow")
)
se_rows <- list(
  list(label = "Newey-West, lag 7",            script = "NW(7)"),
  list(label = "Newey-West, lag 14 (primary)", script = "NW(14) [primary]"),
  list(label = "Newey-West, lag 21",           script = "NW(21)"),
  list(label = "Cluster by month-year",        script = "cluster: month_year"),
  list(label = "Cluster by year",              script = "cluster: year")
)
samp_rows_use <- list(
  list(label_stub = "Full sample",          script = "full"),
  list(label_stub = "Drop zero-death days", script = "drop zero-death days"),
  list(label_stub = "Cap deaths at 100",    script = "cap deaths at 100")
)

L <- character(); add <- function(...) L <<- c(L, paste0(...))
add("\\begin{table}[h!]")
add("\\centering")
add("\\small")
add("\\caption{Fixed-effects and inference robustness, IOM NegBin")
add("$\\beta_3=\\mathrm{SWH}_{t-1:t-5}\\times\\mathrm{Post\\text{-}MoU}$.}")
add("\\label{tab:appx-fe}")
add("\\begin{tabular}{lccc}")
add("\\hline")
add("Variant & $\\beta_3$ & SE & $p$ \\\\")
add("\\hline")
add("\\multicolumn{4}{l}{\\textit{(i) FE specifications (NegBin, NW(14))}} \\\\")
for (r in fe_rows) {
  row <- results[results$fe_spec == r$script, ]
  add(sprintf("%-26s & %-9s & %s & %s \\\\",
              r$label, fb(row$b3), fse(row$b3_se), fp(row$b3_p)))
}
add("\\multicolumn{4}{l}{\\textit{(ii) SE variants on month-year FE}} \\\\")
for (r in se_rows) {
  row <- vcov_rows[vcov_rows$variant == r$script, ]
  add(sprintf("%-26s & %-9s & %s & %s \\\\",
              r$label, fb(row$b3), fse(row$b3_se), fp(row$b3_p)))
}
add("\\multicolumn{4}{l}{\\textit{(iii) Sample restrictions}} \\\\")
for (r in samp_rows_use) {
  row <- samp_rows[samp_rows$variant == r$script, ]
  add(sprintf("%s, $N=%d$  & %-9s & %s & %s \\\\",
              r$label_stub, row$N, fb(row$b3), fse(row$b3_se), fp(row$b3_p)))
}
add("\\hline")
add("\\multicolumn{4}{l}{\\footnotesize SE column reports the standard error from the variance estimator named in each row.}")
add("\\end{tabular}")
add("\\end{table}")
writeLines(L, tbl_path("06_robustness", "tab-appx-fe.tex"))
