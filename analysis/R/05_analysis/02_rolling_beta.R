# Rolling-window β(SWH): 2-year Poisson fits stepped weekly, on UNITED deaths,
# corridor-wide + AFR/EU SAR-bloc panels. Produces fig-rolling-beta.png.

library(tidyverse)
library(fixest)
library(lubridate)
library(patchwork)

BASE_DIR <- here::here()
source(file.path(BASE_DIR, "analysis", "R", "_helpers.R"))

STEP_DAYS      <- 7
MIN_DEATH_DAYS <- 30
MIN_ROWS       <- 200
WINDOW_PRIMARY <- 730

# ── 1. Load panels and replicate primary-sample filter ──────────────────────
panel_daily <- readRDS(file.path(BASE_DIR, "analysis", "data",
                                   "daily_panel_complete.RDS"))
united_daily <- build_united_daily()

panel_daily <- panel_daily |>
  left_join(united_daily, by = "date") |>
  replace_na(list(n_dead_united = 0)) |>
  mutate(
    living_crossings = frx_persons + lcg_tcg_pushbacks,
    lc_lag14 = dplyr::lag(zoo::rollmeanr(living_crossings, k = 7,
                                           fill = NA, align = "right"), 8),
    unit           = 1L,
    month_year_fac = factor(month_year)
  )

sample_dates <- panel_daily |>
  filter(!is.na(lc_lag14), !is.na(swh_prev5days)) |>
  pull(date)

da <- panel_daily |>
  filter(date %in% sample_dates) |>
  arrange(date)

zp <- readRDS(file.path(BASE_DIR, "analysis", "data",
                          "daily_panel_zone.RDS")) |>
  filter(date %in% sample_dates)

bloc <- zp |>
  group_by(date, sar_bloc) |>
  summarise(n_dead_united = sum(n_dead_united),
            swh_prev5days = first(swh_prev5days),
            .groups = "drop") |>
  mutate(year = year(date))
dim(bloc$date) <- NULL

# ── 2. Rolling estimator ────────────────────────────────────────────────────
run_rolling <- function(df, wind = WINDOW_PRIMARY, step = STEP_DAYS) {
  fit_one <- function(sub, mid_date) {
    if (nrow(sub) < MIN_ROWS) return(NULL)
    if (sum(sub$n_dead_united > 0) < MIN_DEATH_DAYS) return(NULL)
    sub$unit <- 1L
    if (length(unique(sub$month_year_fac)) < 2) return(NULL)

    m <- tryCatch(
      fepois(n_dead_united ~ swh_prev5days | month_year_fac,
             data = sub, vcov = NW(14), panel.id = ~unit + date),
      error = function(e) NULL
    )
    if (is.null(m) || !"swh_prev5days" %in% names(coef(m))) return(NULL)

    ct <- tryCatch(coeftable(m, vcov = NW(14)), error = function(e) NULL)
    if (is.null(ct)) return(NULL)
    row <- which(rownames(ct) == "swh_prev5days")
    tibble(date_mid     = mid_date,
           beta         = ct[row, 1],
           se           = ct[row, 2],
           n_obs        = nrow(sub),
           n_death_days = sum(sub$n_dead_united > 0))
  }

  dates_all  <- sort(unique(df$date))
  start_pool <- seq(min(dates_all) + wind / 2,
                     max(dates_all) - wind / 2,
                     by = step)

  map_dfr(start_pool, function(mid) {
    lo <- mid - wind / 2
    hi <- mid + wind / 2
    sub <- df |> filter(date >= lo, date <= hi) |>
      mutate(month_year_fac = factor(format(date, "%Y-%m")))
    fit_one(sub, mid)
  })
}

# ── 3. Run rolling per panel ────────────────────────────────────────────────
res_da_2y <- run_rolling(da, wind = WINDOW_PRIMARY) |>
  mutate(flavor = "daily-agg", label = "daily-agg", window = "2-year")

res_bloc_list <- list()
for (bl in c("AFR", "EU")) {
  sub_bl <- bloc |> filter(sar_bloc == bl)
  res_bloc_list[[bl]] <- run_rolling(sub_bl, wind = WINDOW_PRIMARY) |>
    mutate(flavor = "2-bloc", label = bl, window = "2-year")
}
res_bloc <- bind_rows(res_bloc_list)

all_res <- bind_rows(res_da_2y, res_bloc)

# ── 4. Pre / post-MoU means (descriptive segments on panel a) ───────────────
HALF_WIND <- WINDOW_PRIMARY / 2

per_series <- all_res |>
  group_by(flavor, label, window) |>
  group_modify(~ {
    fully_pre  <- .x |> filter(date_mid <= MOU_DATE - HALF_WIND)
    fully_post <- .x |> filter(date_mid >= MOU_DATE + HALF_WIND)
    tibble(
      beta_pre_mean  = if (nrow(fully_pre)  > 0) mean(fully_pre$beta)  else NA_real_,
      beta_post_mean = if (nrow(fully_post) > 0) mean(fully_post$beta) else NA_real_
    )
  }) |>
  ungroup()

# ── 5. Plot ─────────────────────────────────────────────────────────────────
add_segments <- function(df, gap_threshold = 14) {
  df |>
    group_by(flavor, label, window) |>
    arrange(date_mid) |>
    mutate(gap = c(0, as.numeric(diff(date_mid))),
           segment = cumsum(gap > gap_threshold)) |>
    ungroup() |>
    mutate(series_id = paste(flavor, label, window, segment, sep = "_"))
}

plot_df <- all_res |>
  mutate(ci_lo = beta - 1.96 * se,
         ci_hi = beta + 1.96 * se) |>
  add_segments()

x_min  <- min(plot_df$date_mid, na.rm = TRUE)
x_max  <- max(plot_df$date_mid, na.rm = TRUE)
x_lims <- as.Date(c(x_min, x_max))

x_min_da <- min(plot_df$date_mid[plot_df$flavor == "daily-agg"], na.rm = TRUE)
x_max_da <- max(plot_df$date_mid[plot_df$flavor == "daily-agg"], na.rm = TRUE)

da_summary <- per_series |> filter(flavor == "daily-agg")

mean_segs <- bind_rows(
  da_summary |>
    filter(!is.na(beta_pre_mean)) |>
    transmute(x    = x_min_da,
              xend = MOU_DATE - HALF_WIND,
              y    = beta_pre_mean,
              yend = beta_pre_mean,
              side = "pre"),
  da_summary |>
    filter(!is.na(beta_post_mean)) |>
    transmute(x    = MOU_DATE + HALF_WIND,
              xend = x_max_da,
              y    = beta_post_mean,
              yend = beta_post_mean,
              side = "post")
)

boxed_legend_theme <- theme(
  legend.position       = "right",
  legend.background     = element_blank(),
  legend.box.background = element_rect(fill = "grey97", colour = "grey80",
                                       linewidth = 0.5),
  legend.box.margin     = margin(3, 5, 3, 5),
  legend.box            = "vertical",
  legend.key            = element_blank(),
  legend.margin         = margin(0, 0, 0, 0),
  legend.text           = element_text(size = 9.5, lineheight = 0.9),
  legend.key.size       = unit(0.35, "cm"),
  legend.key.width      = grid::unit(1, "lines"),
  legend.title          = element_blank(),
  legend.spacing        = unit(0, "cm"),
  legend.spacing.y      = unit(0.15, "cm")
)

panel1_col   <- "#6A3D9A"
panel1_label <- "Mean of all incidents\n(full area)"

p1 <- plot_df |>
  filter(flavor == "daily-agg") |>
  ggplot(aes(x = date_mid, y = beta, group = series_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dotted",
             colour = "#D32F2F", linewidth = 0.6) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.15, fill = panel1_col, colour = NA) +
  geom_line(aes(colour = panel1_label), linewidth = 0.7) +
  geom_segment(data = mean_segs,
               aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 1.1, linetype = "longdash", colour = panel1_col,
               inherit.aes = FALSE) +
  scale_colour_manual(values = setNames(panel1_col, panel1_label)) +
  scale_x_date(limits = x_lims, date_breaks = "2 years",
               date_labels = "%Y") +
  labs(
    title  = "(a) Daily-aggregate rolling β",
    x      = NULL, y = expression(beta(SWH[prevweek])),
    colour = NULL
  ) +
  theme_minimal(base_size = 10) +
  boxed_legend_theme +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 11),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.title.y     = element_text(margin = margin(r = 5)),
    plot.margin      = margin(3, 4, 1, 4)
  )

bloc_levels <- c("AFR", "EU")
bloc_labels <- c("AFR" = "African SAR zone\n(Libya + Tunisia)",
                 "EU"  = "European SAR zone\n(Italy + Malta)")
bloc_colors <- c("AFR" = "#D6604D", "EU" = "#2166AC")

p2 <- plot_df |>
  filter(flavor == "2-bloc") |>
  mutate(label = factor(label, levels = bloc_levels)) |>
  ggplot(aes(x = date_mid, y = beta, colour = label, fill = label,
             group = series_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = MOU_DATE, linetype = "dotted",
             colour = "#D32F2F", linewidth = 0.6) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = bloc_colors, labels = bloc_labels) +
  scale_fill_manual  (values = bloc_colors, labels = bloc_labels) +
  scale_x_date(limits = x_lims, date_breaks = "2 years",
               date_labels = "%Y") +
  labs(
    title  = "(b) By SAR responsibility zone",
    x      = NULL, y = expression(beta(SWH[prevweek])),
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  boxed_legend_theme +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 11),
    axis.title.y     = element_text(margin = margin(r = 5)),
    plot.margin      = margin(1, 4, 3, 4)
  )

combined_plot <- p1 / p2

combined_plot_framed <- cowplot::ggdraw(combined_plot) +
  cowplot::draw_grob(grid::rectGrob(
    gp = grid::gpar(col = "black", fill = NA, lwd = 2)
  ))

ggsave(fig_path("05_analysis", "fig-rolling-beta.png"),
       combined_plot_framed, width = 10, height = 6.0, dpi = 200)
