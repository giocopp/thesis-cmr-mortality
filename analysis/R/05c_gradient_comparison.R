# 05c_gradient_comparison.R
# =========================
# Year-by-year SWH-mortality gradient, comparing IOM vs UNITED death sources
# across multiple time windows.
#
# Model: n_dead_missing ~ swh_prevweek:year_fac | month_year
# NegBin, Newey-West(14) SEs.
#
# Windows:
#   IOM:    2014-2023, 2014-2025
#   UNITED: 2014-2023, 2014-2025, 2009-2025
#
# Builds an extended panel from ERA5 SWH + death series directly (no Frontex
# needed), so windows can extend beyond the Frontex-bounded daily panel.
#
# In:  data/processed/era5_swh_daily.RDS
#      data/processed/iom_mmp_incidents.RDS
#      data/processed/united_incidents.RDS
#      data/processed/core_corridor.RDS
# Out: output/figures/05c_gradient_comparison.png
#      output/tables/05c_gradient_comparison.txt

library(tidyverse)
library(lubridate)
library(fixest)
library(patchwork)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

CMR_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia", "Mediterranean")

cat("============================================================\n")
cat("05c  YEARLY GRADIENT COMPARISON: IOM vs UNITED\n")
cat("============================================================\n\n")

# ── 1. Build extended SWH panel ────────────────────────────
# Use the daily_panel_complete.RDS as the base (its SWH works with fixest),
# then extend with raw ERA5 for dates outside the panel window.
cat("--- 1. Loading SWH ---\n")

panel_base <- readRDS(file.path(BASE_DIR, "analysis", "data",
                                 "daily_panel_complete.RDS")) %>%
  select(date, swh, swh_prevweek)

era5_extra <- readRDS(file.path(BASE_DIR, "data", "processed", "era5_swh_daily.RDS")) %>%
  select(date, swh, swh_prevweek) %>%
  filter(!date %in% panel_base$date)

swh <- bind_rows(panel_base, era5_extra) %>%
  filter(!is.na(swh_prevweek)) %>%
  arrange(date)

cat(sprintf("  SWH panel: %s to %s (%d days, %d from panel + %d extended)\n",
            min(swh$date), max(swh$date), nrow(swh),
            nrow(panel_base), nrow(era5_extra)))

# ── 2. Build daily death series ────────────────────────────
cat("\n--- 2. Building death series ---\n")

# IOM: central corridor, sea causes (default helper)
iom_daily <- build_iom_daily()
cat(sprintf("  IOM: %d days, %.0f deaths\n",
            nrow(iom_daily), sum(iom_daily$n_dead_missing)))

# UNITED: same 5 CMR countries + Mediterranean, sea deaths
united_raw <- readRDS(file.path(BASE_DIR, "data", "processed", "united_incidents.RDS"))
united_daily <- united_raw %>%
  filter(country_of_death %in% CMR_COUNTRIES,
         (manner_of_death == "drowned" & !is.na(manner_of_death)) |
         (transport_means == "boat_ship_ferry" & !is.na(transport_means))) %>%
  group_by(date = incident_date_clean) %>%
  summarise(n_dead_missing = sum(n_deaths, na.rm = TRUE), .groups = "drop")
cat(sprintf("  UNITED: %d days, %.0f deaths\n",
            nrow(united_daily), sum(united_daily$n_dead_missing)))

# ── 3. Assemble source-specific panels ─────────────────────
# NOTE: fixest 0.14.0 has an internal scoping issue with tibbles built via
# dplyr pipes on ERA5 data. Workaround: construct with base R data.frame().
cat("\n--- 3. Assembling panels ---\n")

base_panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                                 "daily_panel_complete.RDS"))
era5_raw   <- readRDS(file.path(BASE_DIR, "data", "processed",
                                 "era5_swh_daily.RDS"))

# Combine dates and SWH from panel + ERA5 extension
all_dates <- sort(unique(c(base_panel$date,
                            era5_raw$date[!era5_raw$date %in% base_panel$date])))
all_swh   <- rep(NA_real_, length(all_dates))
all_swh[match(base_panel$date, all_dates)] <- base_panel$swh_prevweek
extra <- era5_raw[!era5_raw$date %in% base_panel$date, ]
all_swh[match(extra$date, all_dates)] <- extra$swh_prevweek

make_panel <- function(death_daily) {
  deaths <- rep(0, length(all_dates))
  idx <- match(death_daily$date, all_dates)
  idx <- idx[!is.na(idx)]
  deaths[idx] <- death_daily$n_dead_missing[match(all_dates[idx], death_daily$date)]

  d <- data.frame(date = all_dates, swh_prevweek = all_swh,
                  n_dead_missing = deaths,
                  year = as.integer(format(all_dates, "%Y")),
                  unit = 1L, stringsAsFactors = FALSE)
  d[!is.na(d$swh_prevweek), ]
}

panel_iom    <- make_panel(iom_daily)
panel_united <- make_panel(united_daily)

cat(sprintf("  IOM panel: %d days (%s to %s)\n",
            nrow(panel_iom), min(panel_iom$date), max(panel_iom$date)))
cat(sprintf("  UNITED panel: %d days (%s to %s)\n",
            nrow(panel_united), min(panel_united$date), max(panel_united$date)))

# ── 4. Yearly gradient estimation + plot ───────────────────
cat("\n--- 4. Estimating yearly gradients ---\n")

make_gradient_plot <- function(yr, label) {
  ggplot(yr, aes(year, beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 2017.5, linetype = "dotted",
               colour = "#D32F2F", linewidth = 0.5) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "grey40") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2) +
    annotate("text", x = 2017.7, y = max(yr$ci_hi) * 0.95,
             label = "MoU", colour = "#D32F2F", size = 3, hjust = 0) +
    scale_x_continuous(breaks = seq(min(yr$year), max(yr$year), by = 2)) +
    labs(title = label,
         x = NULL, y = expression(beta[SWH])) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank())
}

# fixest 0.14.0 has an env-scoping bug when fenegbin is called inside a
# function. Workaround: prepare each dataset at the top level, estimate
# inline, and collect results.

run_gradient <- function(panel_data, start_date, end_date, label) {
  sub <- panel_data[panel_data$date >= start_date & panel_data$date <= end_date &
                     !is.na(panel_data$swh_prevweek), ]
  sub$year_fac   <- factor(sub$year)
  sub$month_year <- factor(format(sub$date, "%Y-%m"))
  .d <<- sub

  cat(sprintf("  %s: N=%d, deaths=%.0f, years=%d-%d\n",
              label, nrow(.d), sum(.d$n_dead_missing),
              min(.d$year), max(.d$year)))
  label
}

extract_results <- function(m, label) {
  co <- coef(m)
  V  <- vcov(m, vcov = NW(14))[seq_along(co), seq_along(co)]
  tibble(year  = parse_number(names(co)),
         beta  = co,
         se    = sqrt(diag(V)),
         ci_lo = beta - 1.96 * se,
         ci_hi = beta + 1.96 * se,
         label = label)
}

# Run all 5 combinations — fenegbin at top level to avoid scoping bug
specs <- list(
  list(panel_iom,    as.Date("2014-01-01"), as.Date("2023-05-31"), "IOM (2014-2023)"),
  list(panel_iom,    as.Date("2014-01-01"), as.Date("2025-12-31"), "IOM (2014-2025)"),
  list(panel_united, as.Date("2014-01-01"), as.Date("2023-05-31"), "UNITED (2014-2023)"),
  list(panel_united, as.Date("2014-01-01"), as.Date("2025-12-31"), "UNITED (2014-2025)"),
  list(panel_united, as.Date("2009-01-01"), as.Date("2025-12-31"), "UNITED (2009-2025)")
)

results_list <- list()
plots_list   <- list()

for (s in specs) {
  lbl <- run_gradient(s[[1]], s[[2]], s[[3]], s[[4]])
  m <- fenegbin(n_dead_missing ~ swh_prevweek:year_fac | month_year,
                data = .d, vcov = NW(14), panel.id = ~unit + date)
  yr <- extract_results(m, lbl)
  results_list[[lbl]] <- yr
  plots_list[[lbl]]   <- make_gradient_plot(yr, lbl)
}

r1 <- results_list[[1]]; r2 <- results_list[[2]]; r3 <- results_list[[3]]
r4 <- results_list[[4]]; r5 <- results_list[[5]]

# ── 5. Combine plots ──────────────────────────────────────
cat("\n--- 5. Saving figures ---\n")

panel_fig <- (plots_list[[1]] | plots_list[[3]]) /
             (plots_list[[2]] | plots_list[[4]]) /
             (plot_spacer()   | plots_list[[5]])
ggsave(file.path(BASE_DIR, "output", "figures", "05c_gradient_comparison.png"),
       panel_fig, width = 12, height = 14, dpi = 200)
cat("Saved: output/figures/05c_gradient_comparison.png\n")

# ── 6. Text results ───────────────────────────────────────
cat("\n--- 6. Saving results table ---\n")

all_results <- bind_rows(results_list)

sink(file.path(BASE_DIR, "output", "tables", "05c_gradient_comparison.txt"))
cat("YEARLY SWH-MORTALITY GRADIENT: IOM vs UNITED\n")
cat("NegBin | swh_prevweek:year_fac | month-year FE | NW(14) SEs\n\n")

all_results %>%
  group_split(label) %>%
  walk(\(x) {
    cat(sprintf("=== %s ===\n", unique(x$label)))
    walk(seq_len(nrow(x)), \(i) {
      r <- x[i, ]
      cat(sprintf("  %d: %+.3f (SE=%.3f)%s\n",
                  r$year, r$beta, r$se,
                  if (abs(r$beta / r$se) > 1.96) " *" else ""))
    })
    cat("\n")
  })
sink()
cat("Saved: output/tables/05c_gradient_comparison.txt\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
