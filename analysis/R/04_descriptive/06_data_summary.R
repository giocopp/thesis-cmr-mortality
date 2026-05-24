# Data and analytical-sample summary table for Methods & Data.
#
# Single descriptive table summarizing the 2014-2023 route series across
# three blocks: (a) sample windows, (b) mortality outcomes (UNITED + IOM,
# analytical filters), (c) constructed exposure components (full panel,
# broad IOM filter that feeds C_t). Numbers come from the same panel and
# shared builders the analytical models use.
#
# Output: a clean booktabs LaTeX table written directly (no gt). Italic
# group headers, indented rows, bolded total, no font-size acrobatics.
#
# Input:
#   analysis/data/daily_panel_complete.RDS
#   data/processed/{iom_mmp_incidents,united_incidents,core_corridor}.RDS
# Output:
#   output/tables/04_descriptive/06_data_summary.tex

library(tidyverse)
library(lubridate)
library(zoo)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("DATA SUMMARY TABLE (Methods & Data)\n")
cat("============================================================\n\n")

# ── 1. Load full panel and analytical death series ─────────
panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS"))

iom_daily    <- build_iom_daily()    # analytical: incident-only, drown+mixed, central
united_daily <- build_united_daily()

panel <- panel |>
  left_join(iom_daily |> rename(n_dead_iom = n_dead_missing), by = "date") |>
  left_join(united_daily, by = "date") |>
  replace_na(list(n_dead_iom = 0, n_dead_united = 0))

# Build the same lc_lag14 + swh_prev5days filter as the primary model so the
# regression-sample size matches the count estimator exactly.
panel <- panel |>
  arrange(date) |>
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag14 = dplyr::lag(
      rollmeanr(living_crossings, k = 7, fill = NA, align = "right"), 8)
  )

d_reg <- panel |> filter(!is.na(lc_lag14), !is.na(swh_prev5days))

cat(sprintf("  Full panel: %d days (%s to %s)\n",
            nrow(panel), min(panel$date), max(panel$date)))
cat(sprintf("  Regression sample: %d days\n", nrow(d_reg)))

# ── 2. Compute numbers ─────────────────────────────────────
window_str <- sprintf("%s to %s",
                       format(min(panel$date), "%Y-%m-%d"),
                       format(max(panel$date), "%Y-%m-%d"))
N_calendar <- nrow(panel)
N_regress  <- nrow(d_reg)
N_rate     <- sum(d_reg$crossing_attempts > 0)

# Mortality (analytical filters, full 2014-2023 panel)
N_united   <- sum(panel$n_dead_united)
N_iom      <- sum(panel$n_dead_iom)
dd_united  <- sum(panel$n_dead_united > 0)
dd_iom     <- sum(panel$n_dead_iom > 0)
max_united <- max(panel$n_dead_united)
max_iom    <- max(panel$n_dead_iom)

# Cross-source correlations
monthly <- panel |>
  mutate(ym = floor_date(date, "month")) |>
  group_by(ym) |>
  summarise(u = sum(n_dead_united), i = sum(n_dead_iom), .groups = "drop")
r_daily   <- cor(panel$n_dead_united, panel$n_dead_iom)
r_monthly <- cor(monthly$u, monthly$i)

# Constructed exposure components.
# The death term in C_t now uses UNITED (analytical filter, same series as
# the primary outcome), stored in the panel as n_dead_united_for_ct.
frx_total          <- sum(panel$frx_persons)
lcg_tcg_total      <- sum(panel$lcg_tcg_pushbacks)
united_in_ct_total <- sum(panel$n_dead_united_for_ct)
Ct_total           <- frx_total + lcg_tcg_total + united_in_ct_total

# Component shares within C_t
frx_share    <- 100 * frx_total          / Ct_total
lcg_share    <- 100 * lcg_tcg_total      / Ct_total
united_share <- 100 * united_in_ct_total / Ct_total

# Sanity check vs hand totals expected from earlier verification
stopifnot(
  N_calendar == 3438L,
  N_united   == 20368,
  N_iom      == 15285,
  dd_united  == 503L,
  dd_iom     == 616L
)

# ── 3. Build LaTeX directly (booktabs) ─────────────────────
fmt_int <- function(x) formatC(round(x), format = "d", big.mark = ",")
fmt_r   <- function(x) formatC(x, format = "f", digits = 2)
fmt_pct <- function(x) sprintf("%.1f\\%%", x)

# Two-column body for sample window and mortality, three-column body for
# exposure (component | count | share). We use a single tabular with three
# columns and let multicolumn handle the cells that only need two.
L <- character()
add <- function(...) L <<- c(L, paste0(...))

add("\\begin{table}[H]")
add("\\centering")
add("\\caption{Summary of data used in the analysis.}")
add("\\label{tab:data-summary}")
add("\\small")
add("\\begin{tabular}{@{}l r r@{}}")
add("\\toprule")

# Panel A: Sample window
add("\\multicolumn{3}{@{}l}{\\textit{Sample window}} \\\\")
add("\\addlinespace[2pt]")
add("\\quad Date range                                        & \\multicolumn{2}{r}{", window_str, "} \\\\")
add("\\quad Calendar days                                     & \\multicolumn{2}{r}{", fmt_int(N_calendar), "} \\\\")
add("\\quad Regression sample (lagged windows filled)         & \\multicolumn{2}{r}{", fmt_int(N_regress),  "} \\\\")
add("\\quad Volume-controlled sample (days with $C_t>0$)      & \\multicolumn{2}{r}{", fmt_int(N_rate),     "} \\\\")
add("\\midrule")

# Panel B: Mortality outcomes
add("\\multicolumn{3}{@{}l}{\\textit{Mortality outcomes (analysis sample)}} \\\\")
add("\\addlinespace[2pt]")
add("                                                         & \\emph{UNITED} & \\emph{IOM} \\\\")
add("\\cmidrule(lr){2-3}")
add("\\quad Deaths (total)                                    & ", fmt_int(N_united),  " & ", fmt_int(N_iom),  " \\\\")
add("\\quad Death-days                                        & ", fmt_int(dd_united), " & ", fmt_int(dd_iom), " \\\\")
add("\\quad Maximum daily count                               & ", fmt_int(max_united)," & ", fmt_int(max_iom)," \\\\")
add("\\addlinespace[2pt]")
add("\\quad Cross-source correlation, daily                   & \\multicolumn{2}{r}{", fmt_r(r_daily),   "} \\\\")
add("\\quad Cross-source correlation, monthly                 & \\multicolumn{2}{r}{", fmt_r(r_monthly), "} \\\\")
add("\\midrule")

# Panel C: Observed crossing attempts
add("\\multicolumn{3}{@{}l}{\\textit{Lower-bound crossing-attempt proxy $C_t$}} \\\\")
add("\\addlinespace[2pt]")
add("                                                         & \\emph{Count} & \\emph{Share} \\\\")
add("\\cmidrule(lr){2-3}")
add("\\quad Persons in Frontex-recorded crossing events       & ", fmt_int(frx_total),          " & ", fmt_pct(frx_share),    " \\\\")
add("\\quad Persons pulled back by Libyan/Tunisian coast guards & ", fmt_int(lcg_tcg_total),      " & ", fmt_pct(lcg_share),    " \\\\")
add("\\quad Deaths recorded by UNITED                         & ", fmt_int(united_in_ct_total), " & ", fmt_pct(united_share), " \\\\")
add("\\addlinespace[2pt]")
add("\\quad \\textbf{Total lower-bound attempts $C_t$}         & \\textbf{", fmt_int(Ct_total), "} & \\textbf{100.0\\%} \\\\")
add("\\bottomrule")
add("\\end{tabular}")
add("")
add("\\vspace{0.4em}")
add("\\begin{minipage}{0.95\\linewidth}")
add("\\footnotesize\\itshape\\raggedright")
add("Notes: Mortality data are filtered by geographic polygonal area and by cause of death (drowning, mixed, or unknown); IOM figures exclude split incidents.")
add("\\end{minipage}")
add("\\end{table}")

writeLines(L, tbl_path("04_descriptive", "06_data_summary.tex"))

cat(sprintf("\nSaved: %s\n", tbl_path("04_descriptive", "06_data_summary.tex")))

# Also print the raw values for quick visual check
cat("\n=== Numbers ===\n")
cat(sprintf("  Window: %s\n", window_str))
cat(sprintf("  Calendar days: %s  Regression: %s  Rate: %s\n",
            fmt_int(N_calendar), fmt_int(N_regress), fmt_int(N_rate)))
cat(sprintf("  UNITED: %s deaths over %s days (max %s)\n",
            fmt_int(N_united), fmt_int(dd_united), fmt_int(max_united)))
cat(sprintf("  IOM:    %s dead/missing over %s days (max %s)\n",
            fmt_int(N_iom), fmt_int(dd_iom), fmt_int(max_iom)))
cat(sprintf("  Corr daily=%s monthly=%s\n", fmt_r(r_daily), fmt_r(r_monthly)))
cat(sprintf("  C_t = %s (Frontex %s | LCG/TCG %s | UNITED %s)\n",
            fmt_int(Ct_total), fmt_int(frx_total),
            fmt_int(lcg_tcg_total), fmt_int(united_in_ct_total)))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
