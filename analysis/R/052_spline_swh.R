# 052_spline_swh.R
# ================
# Enhancement #3: NON-PARAMETRIC functional form via natural cubic spline.
#
# Replaces linear swh_prevweek (metres) with a natural cubic spline (df=3) and
# interacts the spline basis with post_mou. Answers: is the linear
# assumption hiding a threshold effect? Does the curve shape change post-MoU?
#
# Spec (schematic):
#   deaths ~ ns(swh_prevweek, df=3) + ns(swh_prevweek, df=3):post_mou | FE
#
# Three flavors:
#   (A) daily-agg:  FE=month_year,           SE=NW(14)
#   (B) 2-bloc:     FE=month_year + sar_bloc, SE=NW(14)
#   (C) 4-country:  FE=month_year + country,  SE=NW(14)
#
# The individual spline basis coefficients are not directly interpretable;
# we report the JOINT Wald test that all 3 spline-by-post interactions = 0
# (i.e., "does the curve shape change post-MoU?") and save the predicted
# curves figure as the substantive output.
#
# Output: output/tables/052_spline_swh.txt
#         output/figures/052_spline_shape.png

library(tidyverse)
library(fixest)
library(lubridate)
library(splines)

BASE_DIR <- here::here()
MOU_DATE <- as.Date("2017-07-01")
YEAR_START <- 2014
PERIODS <- c(2020, 2023)
SPLINE_DF <- 3

source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

cat("============================================================\n")
cat("052  SPLINE FUNCTIONAL FORM (ns df=3)\n")
cat("============================================================\n\n")

# ── 1. Load data ────────────────────────────────────────────
cat("--- 1. Loading panels ---\n")

# Drop the panel's broad n_dead_missing and replace with the analytical
# series via the shared helper. Default = incident-only, core corridor,
# all causes. Change the call to test sensitivity variants.
da <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_complete.RDS")) %>%
  select(-n_dead_missing) %>%
  left_join(build_iom_daily(), by = "date") %>%
  replace_na(list(n_dead_missing = 0))

zp <- readRDS(file.path(BASE_DIR, "analysis", "data", "daily_panel_zone.RDS"))

# Collapse zone -> 2 blocs: sum deaths, mean SWH (SWH is identical across
# zones on any given day under A2).
bloc <- zp %>%
  group_by(date, sar_bloc) %>%
  summarise(
    n_dead_missing = sum(n_dead_missing),
    swh            = mean(swh, na.rm = TRUE),
    swh_prevweek   = mean(swh_prevweek, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    post_mou   = as.integer(date >= MOU_DATE),
    year       = year(date),
    month_year = factor(format(date, "%Y-%m")),
    iso_week   = paste0(isoyear(date), "_w", sprintf("%02d", isoweek(date)))
  )
dim(bloc$date) <- NULL

cat(sprintf("  daily-agg: %d rows\n", nrow(da)))
cat(sprintf("  zone 4-country panel: %d rows\n", nrow(zp)))
cat(sprintf("  zone 2-bloc panel: %d rows\n", nrow(bloc)))

# ── 2. Helper: build spline basis columns on a data frame ───
# Captures the internal knot positions chosen by ns() so that predict_curve
# below can reuse them. Without this, ns(grid, df = 3) places knots at
# quantiles of the uniform grid — different from the quantiles of the
# fitted data — and the coefficient-times-basis product gives wrong
# predictions (coefficient j then multiplies a different basis function
# than it was fit on).
add_spline_cols <- function(d, var = "swh_prevweek", df = SPLINE_DF) {
  bnd <- range(d[[var]], na.rm = TRUE)
  B <- ns(d[[var]], df = df, Boundary.knots = bnd)
  knots <- attr(B, "knots")
  colnames(B) <- paste0("spl", 1:df)
  d <- cbind(d, B)
  for (j in 1:df) {
    d[[paste0("spl", j, "_post")]] <- d[[paste0("spl", j)]] * d$post_mou
  }
  attr(d, "spline_bounds") <- bnd
  attr(d, "spline_knots")  <- knots
  d
}

spl_cols  <- paste0("spl", 1:SPLINE_DF)
intr_cols <- paste0("spl", 1:SPLINE_DF, "_post")
spl_rhs <- paste(c(spl_cols, intr_cols), collapse = " + ")

# ── 3. Estimation ───────────────────────────────────────────
cat("\n--- 2. Estimation ---\n")

sink_file <- file.path(BASE_DIR, "output", "tables", "052_spline_swh.txt")
sink(sink_file)

cat("052  SPLINE SWH (natural cubic spline, df=3)\n")
cat("============================================\n")
cat("Functional form: ns(swh_prevweek, df=3), interacted with post_mou.\n")
cat("Joint Wald test: H0 = all 3 spline-by-post interactions are 0\n")
cat("                 (i.e., the SWH->deaths curve shape did NOT change post-MoU).\n\n")
cat("(A) daily-agg:  FE=month_year,           NW(14)\n")
cat("(B) 2-bloc:     FE=month_year+sar_bloc,  NW(14), panel=sar_bloc+date\n")
cat("(C) 4-country:  FE=month_year+country,   NW(14), panel=country+date\n\n")

all_models <- list()

for (ye in PERIODS) {
  label <- sprintf("%d-%d", YEAR_START, ye)
  cat(sprintf("\n=== %s ===\n", label))

  # --- (A) daily-agg ---
  d_da <- da %>%
    filter(year(date) >= YEAR_START, year(date) <= ye,
           !is.na(swh_prevweek)) %>%
    mutate(unit = 1L) %>%
    add_spline_cols()

  f_da <- as.formula(paste("n_dead_missing ~", spl_rhs, "| month_year"))
  m_da <- fenegbin(f_da, data = d_da,
                    vcov = NW(14), panel.id = ~unit + date)

  # --- (B) 2-bloc ---
  d_bl <- bloc %>%
    filter(year >= YEAR_START, year <= ye, !is.na(swh_prevweek)) %>%
    add_spline_cols()

  f_bl <- as.formula(paste("n_dead_missing ~", spl_rhs,
                            "| month_year + sar_bloc"))
  m_bl <- fenegbin(f_bl, data = d_bl,
                    vcov = NW(14), panel.id = ~sar_bloc + date)

  # --- (C) 4-country ---
  d_zp <- zp %>%
    filter(year >= YEAR_START, year <= ye, !is.na(swh_prevweek)) %>%
    add_spline_cols()

  f_zp <- as.formula(paste("n_dead_missing ~", spl_rhs,
                            "| month_year + country"))
  m_zp <- fenegbin(f_zp, data = d_zp,
                    vcov = NW(14), panel.id = ~country + date)

  cat(sprintf("  [A] daily-agg  N = %d\n", nrow(d_da)))
  cat(sprintf("  [B] 2-bloc     N = %d\n", nrow(d_bl)))
  cat(sprintf("  [C] 4-country  N = %d\n", nrow(d_zp)))

  # Joint Wald test (null: all 3 spline-by-post interactions = 0)
  cat("\n  Joint Wald test H0: all spline-by-post interactions = 0\n")
  joint_test <- function(m, vcov_type = NW(14)) {
    wt <- tryCatch(
      wald(m, keep = "_post$", vcov = vcov_type, print = FALSE),
      error = function(e) NULL
    )
    if (!is.null(wt)) {
      cat(sprintf("    stat = %.3f, p = %.4f\n", wt$stat, wt$p))
      invisible(wt)
    } else {
      cat("    (wald test failed)\n")
      invisible(NULL)
    }
  }

  cat("    [A] daily-agg:\n")
  joint_test(m_da)
  cat("    [B] 2-bloc:\n")
  joint_test(m_bl)
  cat("    [C] 4-country:\n")
  joint_test(m_zp)

  all_models[[label]] <- list(
    da = list(model = m_da, data = d_da,
              bnd   = attr(d_da, "spline_bounds"),
              knots = attr(d_da, "spline_knots")),
    bl = list(model = m_bl, data = d_bl,
              bnd   = attr(d_bl, "spline_bounds"),
              knots = attr(d_bl, "spline_knots")),
    zp = list(model = m_zp, data = d_zp,
              bnd   = attr(d_zp, "spline_bounds"),
              knots = attr(d_zp, "spline_knots"))
  )
}

sink()
cat(sprintf("\nSaved: %s\n", sink_file))

# ── 4. Predicted curves ─────────────────────────────────────
cat("\n--- 3. Predicted-curve plot ---\n")

# For each flavor and period, predict the relative IRR across a SWH grid
# (relative to the median SWH), comparing pre-MoU vs post-MoU. FE contribution
# cancels when we take ratios against the median.

predict_curve <- function(m, grid, bnd, knots) {
  # Reuse the internal knots from the fitted basis so B_grid matches the
  # basis functions the coefficients were estimated on.
  B_grid <- ns(grid, knots = knots, Boundary.knots = bnd)
  colnames(B_grid) <- paste0("spl", 1:SPLINE_DF)
  co <- coef(m)

  pre_main <- rep(0, length(grid))
  for (j in 1:SPLINE_DF) {
    nm <- paste0("spl", j)
    if (nm %in% names(co)) pre_main <- pre_main + co[[nm]] * B_grid[, j]
  }

  post_main <- pre_main
  for (j in 1:SPLINE_DF) {
    nm <- paste0("spl", j, "_post")
    if (nm %in% names(co)) post_main <- post_main + co[[nm]] * B_grid[, j]
  }

  med_idx <- which.min(abs(grid - median(grid)))
  pre_rel <- exp(pre_main - pre_main[med_idx])
  post_rel <- exp(post_main - post_main[med_idx])
  tibble(swh = grid, pre = pre_rel, post = post_rel)
}

plot_list <- list()
for (label in names(all_models)) {
  am <- all_models[[label]]
  grid_da <- seq(am$da$bnd[1], am$da$bnd[2], length.out = 100)
  grid_bl <- seq(am$bl$bnd[1], am$bl$bnd[2], length.out = 100)
  grid_zp <- seq(am$zp$bnd[1], am$zp$bnd[2], length.out = 100)

  df_da <- predict_curve(am$da$model, grid_da, am$da$bnd, am$da$knots) %>%
    pivot_longer(c(pre, post), names_to = "era", values_to = "IRR") %>%
    mutate(flavor = "(A) daily-agg", period = label)
  df_bl <- predict_curve(am$bl$model, grid_bl, am$bl$bnd, am$bl$knots) %>%
    pivot_longer(c(pre, post), names_to = "era", values_to = "IRR") %>%
    mutate(flavor = "(B) 2-bloc", period = label)
  df_zp <- predict_curve(am$zp$model, grid_zp, am$zp$bnd, am$zp$knots) %>%
    pivot_longer(c(pre, post), names_to = "era", values_to = "IRR") %>%
    mutate(flavor = "(C) 4-country", period = label)

  plot_list[[label]] <- bind_rows(df_da, df_bl, df_zp)
}

plot_df <- bind_rows(plot_list) %>%
  mutate(era = factor(era, levels = c("pre", "post"),
                       labels = c("Pre-MoU", "Post-MoU")),
         flavor = factor(flavor,
                          levels = c("(A) daily-agg", "(B) 2-bloc",
                                     "(C) 4-country")))

p_spline <- ggplot(plot_df, aes(x = swh, y = IRR, colour = era, linetype = era)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Pre-MoU" = "#2166AC", "Post-MoU" = "#B2182B")) +
  scale_y_log10() +
  facet_grid(flavor ~ period, scales = "free_y") +
  labs(
    title = "Spline SWH: predicted IRR relative to median SWH, pre vs post MoU",
    subtitle = "Natural cubic spline df=3. y-axis is log IRR (death count rate).",
    x = "Previous-week mean SWH (m)",
    y = "IRR relative to median SWH",
    colour = "Era", linetype = "Era"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

fig_path <- file.path(BASE_DIR, "output", "figures", "052_spline_shape.png")
ggsave(fig_path, p_spline, width = 10, height = 9, dpi = 200)
cat(sprintf("Saved: %s\n", fig_path))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
