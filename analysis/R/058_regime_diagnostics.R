# 058_regime_diagnostics.R
# ========================
# Diagnostic companion to 05_reduced_form_primary.R — does NOT replace the
# binary post_mou model. Motivated by the rolling-beta pattern in
# 053_rolling_beta.R, which shows the SWH-deaths gradient drifting gradually
# across 2014-2022 rather than jumping at MoU. The question here is simple:
# if we let the data pick the break date instead of pre-committing to
# 2017-07-01, where does it sit?
#
# Method: Andrews sup-Wald grid sweep. For each candidate break date tau in
# a weekly grid, fit the primary model with tau replacing MoU_date:
#
#     n_dead_missing ~ swh_prevweek + swh_prevweek:1[date >= tau]
#                      | month_year
#
# with fenegbin + NW(14) — exactly the specification of 05_reduced_form_primary.R,
# only the break date moves. Record the z-stat on the interaction and find
# tau* = argmax |z|. Compare tau* to the institutional dates (Mare Nostrum
# end, MoU+Minniti, Salvini crackdown).
#
# The test statistic is non-standard under the null because tau is
# estimated (Andrews 1993 / Andrews-Ploberger 1994), so we do NOT read it
# as a p-value — it is a *visual* diagnostic. What we want to know is:
# where does the peak of the z-profile sit relative to the institutional
# dates, and how sharp/broad is it.
#
# In:  analysis/data/daily_panel_complete.RDS
# Out: output/tables/058_regime_diagnostics.txt
#      output/figures/058_break_sweep.png

library(tidyverse)
library(lubridate)
library(fixest)

BASE_DIR    <- here::here()
MOU_DATE    <- as.Date("2017-07-01")
START_DATE  <- as.Date("2014-01-01")

# Institutional reference dates used for visual comparison. Not inputs to
# the statistical tests — only overlaid on plots.
INSTITUTIONAL <- tibble(
  date  = as.Date(c("2014-11-01", "2017-07-01", "2018-06-01")),
  label = c("Mare Nostrum ends", "MoU + Minniti code", "Salvini crackdown"),
  col   = c("#1B7837",            "#D6604D",             "#762A83")
)

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("058  REGIME DIAGNOSTICS: BREAK-POINT SWEEP + SAR COVERAGE\n")
cat("============================================================\n\n")

# ── 1. Load panel and rebuild primary death series ─────────
cat("--- 1. Loading panel ---\n")

panel <- readRDS(file.path(BASE_DIR, "analysis", "data",
                            "daily_panel_complete.RDS")) %>%
  select(-n_dead_missing) %>%
  arrange(date) %>%
  left_join(build_iom_daily(), by = "date") %>%
  replace_na(list(n_dead_missing = 0)) %>%
  mutate(
    unit            = 1L,
    year            = year(date),
    month_year      = factor(format(date, "%Y-%m")),
    year_fac        = factor(year),
    month_of_year   = factor(month(date))
  )

PANEL_END <- max(panel$date)
cat(sprintf("  N = %d days, %s to %s, %.0f deaths\n",
            nrow(panel), min(panel$date), PANEL_END,
            sum(panel$n_dead_missing)))

# ────────────────────────────────────────────────────────────
# OPTION B — Andrews sup-Wald break-point sweep
# ────────────────────────────────────────────────────────────
cat("\n============================================================\n")
cat("OPTION B — SUP-WALD BREAK-POINT SWEEP\n")
cat("============================================================\n\n")

# Candidate break dates on a weekly grid. Leave a 1-year buffer on each
# side so every fit has enough pre and post observations for the NW(14)
# bandwidth and the month_year FE to stay identified.
buffer_days <- 365
tau_grid <- seq(min(panel$date) + buffer_days,
                PANEL_END      - buffer_days,
                by = "7 days")
cat(sprintf("  Candidate grid: %d weekly dates from %s to %s\n",
            length(tau_grid), min(tau_grid), max(tau_grid)))

# One fit per candidate tau, fenegbin + NW(14) (same family + SE as the
# primary model), record z on the interaction.
fit_tau <- function(tau) {
  d_tau <- panel %>% mutate(post_tau = as.integer(date >= tau))
  m <- tryCatch(
    fenegbin(n_dead_missing ~ swh_prevweek + swh_prevweek:post_tau |
               month_year,
             data = d_tau, vcov = NW(14), panel.id = ~ unit + date),
    error = function(e) NULL
  )
  if (is.null(m)) return(tibble(tau = tau, b_int = NA_real_, se_int = NA_real_,
                                 z_int = NA_real_, b_main = NA_real_))
  ct <- coeftable(m)
  rn <- rownames(ct)
  r_int  <- grep("swh_prevweek:post_tau|post_tau:swh_prevweek", rn)
  r_main <- which(rn == "swh_prevweek")
  if (length(r_int) != 1) return(tibble(tau = tau, b_int = NA_real_,
                                         se_int = NA_real_, z_int = NA_real_,
                                         b_main = NA_real_))
  tibble(
    tau    = tau,
    b_main = ct[r_main, 1],
    b_int  = ct[r_int, 1],
    se_int = ct[r_int, 2],
    z_int  = ct[r_int, 1] / ct[r_int, 2]
  )
}

t0 <- Sys.time()
sup_wald <- map_dfr(tau_grid, fit_tau)
cat(sprintf("  Swept %d candidate dates in %.1f s\n",
            nrow(sup_wald), as.numeric(Sys.time() - t0, units = "secs")))

# Peak: tau* = argmax |z| on the interaction. This is the single-break
# Andrews estimate. It is not unique if several nearby dates give similar
# |z|, so we also report the top-10 dates as a region.
sup_wald <- sup_wald %>%
  mutate(abs_z = abs(z_int)) %>%
  arrange(desc(abs_z))

tau_star <- sup_wald$tau[1]
cat(sprintf("\n  tau* (argmax |z|) = %s   (b_int = %+.3f, z = %+.2f)\n",
            tau_star, sup_wald$b_int[1], sup_wald$z_int[1]))

cat("\n  Top 10 candidate break dates by |z|:\n")
print(sup_wald %>% head(10) %>%
  transmute(tau, b_main = round(b_main, 3),
            b_int = round(b_int, 3), z = round(z_int, 2)))

# Distance from tau* to each institutional date
inst_cmp <- INSTITUTIONAL %>%
  mutate(days_from_tau_star = as.numeric(date - tau_star))
cat("\n  Distance (days) from tau* to institutional dates:\n")
print(inst_cmp)

# Plot: z-stat profile with the institutional dates as vertical markers.
sup_wald_plot <- sup_wald %>% arrange(tau)

p_B <- ggplot(sup_wald_plot, aes(tau, z_int)) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed",
             colour = "grey60") +
  geom_hline(yintercept = 0, colour = "grey30") +
  geom_line(linewidth = 0.6, colour = "#2166AC") +
  geom_point(data = filter(sup_wald_plot, tau == tau_star),
             aes(tau, z_int), size = 3.5, colour = "#D6604D") +
  geom_vline(data = INSTITUTIONAL,
             aes(xintercept = date, colour = label),
             linetype = "dotted", linewidth = 0.6,
             show.legend = TRUE) +
  scale_colour_manual(values = setNames(INSTITUTIONAL$col, INSTITUTIONAL$label),
                      name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title    = "(B) Andrews sup-Wald: z-stat on SWH x 1[t>=tau] across candidate break dates",
    subtitle = sprintf(
      paste("fenegbin, NW(14), month-year FE. Grid = weekly, %d-day buffer.",
            "Red dot = argmax |z| at %s (z=%+.2f). Dotted lines = policy dates."),
      buffer_days, format(tau_star, "%Y-%m-%d"), sup_wald$z_int[1]),
    x = NULL,
    y = expression(z-statistic ~ on ~ beta[SWH %*% "post"])
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

ggsave(file.path(BASE_DIR, "output", "figures", "058_break_sweep.png"),
       p_B, width = 10, height = 5.5, dpi = 200)
cat("\n  Saved: output/figures/058_break_sweep.png\n")

# ── Text output ─────────────────────────────────────────────
cat("\n--- Saving text output ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables", "058_regime_diagnostics.txt")
sink(sink_file)
old_opts <- options(tibble.width = Inf, tibble.print_max = Inf)
on.exit(options(old_opts), add = TRUE)

cat("058  REGIME DIAGNOSTICS — structural break sweep\n")
cat("================================================\n")
cat(sprintf("Sample: %s to %s (N = %d days, %.0f deaths)\n",
            min(panel$date), PANEL_END, nrow(panel),
            sum(panel$n_dead_missing)))
cat("\n")

cat("Andrews sup-Wald break-point sweep\n")
cat("----------------------------------\n")
cat("Model: fenegbin(n_dead_missing ~ swh_prevweek + swh_prevweek:1[date>=tau] | month_year)\n")
cat("SE: NW(14). Candidate tau grid: weekly, 1-year buffer on each side.\n")
cat("(Same family + SE as 05_reduced_form_primary.R; only the break date moves.)\n\n")
cat(sprintf("tau* (argmax |z|): %s   b_int = %+.3f   z = %+.2f\n",
            tau_star, sup_wald$b_int[1], sup_wald$z_int[1]))
cat("\nTop 10 candidate break dates by |z|:\n")
print(sup_wald %>% head(10) %>%
  transmute(tau,
            b_main = round(b_main, 3),
            b_int  = round(b_int,  3),
            z_int  = round(z_int,  2)))
cat("\nDistance (days) from tau* to institutional dates:\n")
print(inst_cmp)

sink()
cat(sprintf("Saved: %s\n", sink_file))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
