# Deiana, Maheshri, Mastrobuoni (2024) — Table 3
# "Elasticities of Crossing Attempts to Crossing Conditions"
# (AEJ: Economic Policy 16(2), 335-365)
#
# Poisson GLM with week-year fixed effects; Newey-West HAC SEs (lag 28).
# The paper interprets coefficients on continuous regressors as
# semielasticities (100 * [exp(beta) - 1] %).
#
# R replication of the Stata spec in do/dofile.do (lines 274-366):
#   glm totacross onda onda_frac onda_post onda_post_frac i.weekanno,
#       family(poisson) vce(hac nwest 28) t(data)

suppressPackageStartupMessages({
  library(haven)
  library(sandwich)
  library(lmtest)
  library(dplyr)
})

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  if (!is.null(sys.frames()[[1]]$ofile))
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
})
# R-rep/ lives alongside Replication/ (Deiana's original files)
dei_data <- normalizePath(file.path(script_dir, "..", "Replication", "data"))
out_dir  <- file.path(script_dir, "out")

d <- read_dta(file.path(dei_data, "data_tab_main.dta")) |>
  arrange(data)

fr_defs <- list(
  `(1) Inflatable`                      = "fr_across3",
  `(2) Inflatable + Unknown`            = "fr_across1_3",
  `(3) Inflatable + Unknown + Other`    = "fr_across_123"
)

fit_col <- function(fr) {
  d$onda           <- d$swh_lyb_3
  d$onda_frac      <- d$swh_lyb_3 * d[[fr]]
  d$onda_post      <- d$swh_lyb_3 * d$postM
  d$onda_post_frac <- d$swh_lyb_3 * d$postM * d[[fr]]

  m <- glm(
    totacross ~ onda + onda_frac + onda_post + onda_post_frac + factor(weekanno),
    family = poisson(link = "log"),
    data   = d
  )

  # Stata vce(hac nwest 28) on glm: Bartlett kernel, max lag 28,
  # no prewhitening, no n/(n-k) small-sample correction.
  v  <- NeweyWest(m, lag = 28, prewhite = FALSE, adjust = FALSE)
  ct <- coeftest(m, vcov. = v)

  keep <- c("onda_post_frac", "onda", "onda_frac", "onda_post")
  list(model = m, vcov = v, coef = ct[keep, , drop = FALSE], n = nobs(m))
}

results <- lapply(fr_defs, fit_col)

pre <- d[d$postM == 0, ]
pre_stats <- c(
  mean_totacross    = round(mean(pre$totacross), 0),
  mean_wave_height  = round(mean(pre$swh_lyb_3), 2),
  mean_fr_across3   = round(mean(pre$fr_across3), 2),
  mean_fr_across1_3 = round(mean(pre$fr_across1_3), 2),
  mean_fr_across_123= round(mean(pre$fr_across_123), 2)
)

fmt <- function(x) formatC(x, digits = 3, format = "f")
row_line <- function(label, idx) {
  vals <- sapply(results, function(r) sprintf("%s", fmt(r$coef[idx, "Estimate"])))
  ses  <- sapply(results, function(r) sprintf("(%s)", fmt(r$coef[idx, "Std. Error"])))
  paste0(
    sprintf("%-38s %12s %12s %12s\n", label, vals[1], vals[2], vals[3]),
    sprintf("%-38s %12s %12s %12s\n", "",    ses[1],  ses[2],  ses[3])
  )
}

out <- c(
  "Table 3 — Elasticities of Crossing Attempts to Crossing Conditions\n",
  "Outcome: totacross  |  Poisson GLM, HAC Newey-West (lag 28)\n",
  strrep("-", 78), "\n",
  sprintf("%-38s %12s %12s %12s\n", "",
          "(1) Infl.", "(2) +Unkn.", "(3) +Other"),
  strrep("-", 78), "\n",
  row_line("Wave Height * Post SAR * Fr. Boat", "onda_post_frac"),
  row_line("Wave Height",                       "onda"),
  row_line("Wave Height * Fr. Boat",            "onda_frac"),
  row_line("Wave Height * Post SAR",            "onda_post"),
  strrep("-", 78), "\n",
  "Week-Year FE                           ", "yes          yes          yes\n",
  sprintf("%-38s %12d %12d %12d\n", "Observations",
          results[[1]]$n, results[[2]]$n, results[[3]]$n),
  strrep("-", 78), "\n",
  "Pre-SAR means:\n",
  sprintf("  Mean total attempts : %d\n",   pre_stats["mean_totacross"]),
  sprintf("  Mean wave height    : %.2f\n", pre_stats["mean_wave_height"]),
  sprintf("  Mean frac. unsafe   : %.2f / %.2f / %.2f\n",
          pre_stats["mean_fr_across3"],
          pre_stats["mean_fr_across1_3"],
          pre_stats["mean_fr_across_123"])
)

cat(out, sep = "")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(paste0(out, collapse = ""), file.path(out_dir, "tab3_elasticities.txt"))
saveRDS(results, file.path(out_dir, "tab3_elasticities.rds"))
