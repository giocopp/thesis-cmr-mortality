# 03b3_placebo_lag1.R
# ===================
# Placebo treatment dates for the SWH lag-1 specification.
# Re-estimates beta_3 at quarterly placebo break dates with both
# week-year and month-year FE. If the MoU date is special, its
# beta_3 should stand out relative to other dates.
#
# Input:  data/processed/cmr_daily_weather_panel.RDS
# Output: output/figures/placebo_lag1.pdf / .png
#         output/tables/placebo_lag1_results.csv

library(fixest)
library(data.table)
library(ggplot2)

BASE_DIR <- here::here()
d <- as.data.table(readRDS(file.path(BASE_DIR, "data", "processed",
                                      "cmr_daily_weather_panel.RDS")))

MOU_DATE <- as.Date("2017-07-01")

cat("============================================================\n")
cat("PLACEBO TEST: SWH lag-1 x Post at quarterly break dates\n")
cat("============================================================\n\n")

# Candidate dates: quarterly from 2015-Q1 to 2021-Q1
placebo_dates <- as.Date(c(
  "2015-01-01", "2015-04-01", "2015-07-01", "2015-10-01",
  "2016-01-01", "2016-04-01", "2016-07-01", "2016-10-01",
  "2017-04-01",
  "2017-07-01",  # actual MoU
  "2017-10-01",
  "2018-01-01", "2018-04-01", "2018-07-01",
  "2019-01-01", "2019-07-01",
  "2020-01-01", "2020-07-01",
  "2021-01-01"
))

estimate_placebo <- function(pd, fe_type = "weekly") {
  d[, post_placebo := as.integer(date >= pd)]

  n_pre  <- sum(d$post_placebo == 0)
  n_post <- sum(d$post_placebo == 1)
  if (n_pre < 60 || n_post < 60) return(NULL)

  if (fe_type == "weekly") {
    fe_var <- "week_year_fac"
  } else {
    fe_var <- "month_year_fac"
  }

  warn_msg <- NULL
  fml <- as.formula(paste0(
    "n_dead_missing ~ swh_core_lag1 + swh_core_lag1:post_placebo | ", fe_var))

  m <- tryCatch(
    withCallingHandlers(
      fenegbin(fml, data = d[!is.na(swh_core_lag1)], vcov = "hetero"),
      warning = function(w) { warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
    ),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)

  converged <- is.null(warn_msg) || !grepl("did not converge|singular", warn_msg)

  ct <- summary(m, vcov = "hetero")$coeftable
  int_row <- grep("swh_core_lag1:post_placebo|post_placebo:swh_core_lag1",
                  rownames(ct), fixed = FALSE)
  if (length(int_row) == 0) return(NULL)

  data.table(
    date      = pd,
    is_mou    = pd == MOU_DATE,
    fe        = ifelse(fe_type == "weekly", "Week-year FE", "Month-year FE"),
    beta      = ct[int_row, 1],
    se        = if (converged) ct[int_row, 2] else NA_real_,
    p         = if (converged) ct[int_row, 4] else NA_real_,
    irr       = exp(ct[int_row, 1]),
    n_obs     = nobs(m),
    converged = converged
  )
}

# --- Run for both FE structures ---
cat("--- Week-year FE ---\n")
res_weekly <- rbindlist(lapply(placebo_dates, function(pd) {
  estimate_placebo(pd, "weekly")
}))

cat("--- Month-year FE ---\n")
res_monthly <- rbindlist(lapply(placebo_dates, function(pd) {
  estimate_placebo(pd, "monthly")
}))

results <- rbind(res_weekly, res_monthly)
results[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]

# --- Print table ---
for (fe in c("Week-year FE", "Month-year FE")) {
  cat(sprintf("\n%s:\n", fe))
  cat(sprintf("  %-12s %+8s %7s %7s %8s\n", "Date", "Beta", "SE", "p", "IRR"))
  cat(paste0("  ", paste(rep("-", 55), collapse = "")), "\n")
  sub <- results[fe == get("fe")]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i]
    marker <- ifelse(r$is_mou, " <-- MoU", "")
    if (!r$converged) {
      cat(sprintf("  %-12s %+8.4f %7s %7s %8.4f [NC]%s\n",
          as.character(r$date), r$beta, "NA", "NA", r$irr, marker))
    } else {
      stars <- ifelse(r$p < 0.01, "***",
               ifelse(r$p < 0.05, "**",
               ifelse(r$p < 0.1, "*", "")))
      cat(sprintf("  %-12s %+8.4f %7.4f %7.4f %8.4f %3s%s\n",
          as.character(r$date), r$beta, r$se, r$p, r$irr, stars, marker))
    }
  }
}

# --- Save table ---
fwrite(results, file.path(BASE_DIR, "output", "tables", "placebo_lag1_results.csv"))
cat("\nSaved: output/tables/placebo_lag1_results.csv\n")

# --- Plot ---
cat("\n--- Generating plot ---\n")

plot_dt <- results[converged == TRUE]

p <- ggplot(plot_dt, aes(x = date, y = beta, colour = fe, fill = fe)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dashed", colour = "red",
             linewidth = 0.6) +
  annotate("text", x = MOU_DATE, y = max(plot_dt$ci_hi, na.rm = TRUE) * 0.95,
           label = "MoU (Jul 2017)", colour = "red", hjust = -0.1, size = 3) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.5) +
  geom_point(aes(shape = is_mou, size = is_mou)) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4), guide = "none") +
  scale_colour_manual(values = c("Week-year FE" = "steelblue",
                                  "Month-year FE" = "darkorange")) +
  scale_fill_manual(values = c("Week-year FE" = "steelblue",
                                "Month-year FE" = "darkorange")) +
  labs(
    title = expression("Placebo treatment dates: " * hat(beta)[3] *
                        " (SWH lag-1 × Post)"),
    subtitle = "NegBin, drowning/suspected drowning, core corridor.\nDiamond = actual MoU. Red line = Jul 2017.",
    x = "Placebo treatment date",
    y = expression(hat(beta)[3]),
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(BASE_DIR, "output", "figures", "placebo_lag1.pdf"),
       p, width = 10, height = 6)
ggsave(file.path(BASE_DIR, "output", "figures", "placebo_lag1.png"),
       p, width = 10, height = 6, dpi = 200)
cat("Saved: output/figures/placebo_lag1.pdf + .png\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
