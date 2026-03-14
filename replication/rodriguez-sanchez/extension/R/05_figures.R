# 05_figures.R
# ============
# All plotting functions for the mortality analysis.
# Extracted from run_new_model.R lines 254-920 and 1268-1826.


# --- Constants ---
EPSILON <- 0.01
Y_AXIS_CLIP_QUANTILES <- c(0.02, 0.98)
BACKTRANSFORM_LOG_CAP <- 8


# --- Helper functions ---

moving_average <- function(x, k = 6) {
  n <- length(x)
  half <- floor(k / 2)
  vapply(seq_len(n), function(i) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    mean(x[lo:hi], na.rm = TRUE)
  }, numeric(1))
}

panel_ylim <- function(df, cols,
                       probs = Y_AXIS_CLIP_QUANTILES, pad = 0.08) {
  vals <- unlist(df[, cols], use.names = FALSE)
  vals <- vals[is.finite(vals)]
  if (length(vals) < 5) return(NULL)
  qv <- as.numeric(quantile(vals, probs = probs, na.rm = TRUE))
  if (qv[1] == qv[2]) {
    span <- max(1, abs(qv[1])) * 0.1
    return(c(qv[1] - span, qv[2] + span))
  }
  span <- qv[2] - qv[1]
  c(qv[1] - pad * span, qv[2] + pad * span)
}

add_natural_scale <- function(df, outcome_type = c("rate", "deaths")) {
  outcome_type <- match.arg(outcome_type)
  cap_log <- function(x) pmin(x, BACKTRANSFORM_LOG_CAP)
  inv_fun <- if (outcome_type == "rate") {
    function(x) pmax(exp(cap_log(x)) - EPSILON, 0)
  } else {
    function(x) pmax(exp(cap_log(x)) - 1, 0)
  }

  df %>%
    dplyr::mutate(
      original_nat               = inv_fun(original),
      original_smooth_nat        = moving_average(original_nat, k = 6),
      prediction_nat             = inv_fun(prediction),
      prediction_lower_nat       = inv_fun(prediction_lower),
      prediction_upper_nat       = inv_fun(prediction_upper),
      pointwise_effect_nat       = original_nat - prediction_nat,
      pointwise_effect_lower_nat = original_nat - prediction_upper_nat,
      pointwise_effect_upper_nat = original_nat - prediction_lower_nat
    )
}


# --- Theme functions ---

theme_paper <- function() {
  ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 16, hjust = 0),
      axis.title  = ggplot2::element_text(size = 14, color = "black"),
      axis.text   = ggplot2::element_text(size = 12, color = "black"),
      axis.line   = ggplot2::element_line(color = "black", linewidth = 0.8),
      axis.ticks  = ggplot2::element_line(color = "black", linewidth = 0.8),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
    )
}

theme_panel <- function() {
  ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
      axis.title  = ggplot2::element_text(size = 14, color = "black"),
      axis.text   = ggplot2::element_text(size = 12, color = "black"),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.line   = ggplot2::element_line(color = "black", linewidth = 0.8),
      axis.ticks  = ggplot2::element_line(color = "black", linewidth = 0.8),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
    )
}

theme_diag <- function() {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      axis.title = ggplot2::element_text(size = 11),
      axis.text  = ggplot2::element_text(size = 10)
    )
}


# --- Color palette ---
fill_pink  <- "#FDE6E6"
fill_grey  <- "#D9D9D9"
fill_ci    <- "#B7C9E2"
col_cf     <- "#6A6FB0"
col_smooth <- "grey40"
col_obs    <- "black"


#' Save a figure as PNG + PDF, return file paths
#'
#' @param plot ggplot or grob object
#' @param name File name without extension
#' @param output_dir Base output directory (will create figures/ sub-dirs)
#' @param width Width in inches
#' @param height Height in inches
#' @return character vector of saved file paths
save_figure <- function(plot, name, output_dir = "output/figures",
                        width = 8, height = 10) {
  png_dir <- file.path(output_dir, "png")
  pdf_dir <- file.path(output_dir, "pdf")
  dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

  png_path <- file.path(png_dir, paste0(name, ".png"))
  pdf_path <- file.path(pdf_dir, paste0(name, ".pdf"))

  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 300)
  ggplot2::ggsave(pdf_path, plot, width = width, height = height)

  c(png_path, pdf_path)
}


#' Plot Figure 1: Descriptive time series
#'
#' @param df_model Output of prepare_model_data() (list)
#' @param output_dir Output directory
#' @return character vector of saved file paths
plot_descriptive <- function(df_model, output_dir = "output/figures") {
  df <- df_model$df_full

  intervention_dates <- list(
    mare_nostrum_start = lubridate::ymd("2013-10-01"),
    mare_nostrum_end   = lubridate::ymd("2014-10-01"),
    ngo_sar_start      = lubridate::ymd("2014-11-01"),
    eu_libya_start     = lubridate::ymd("2017-02-01"),
    lcg_sar_zone       = lubridate::ymd("2017-08-01")
  )

  x_breaks <- lubridate::ymd(c("2010-01-01", "2015-01-01", "2020-01-01"))

  fig1a <- ggplot2::ggplot(df, ggplot2::aes(x = date,
                                             y = dead_and_missing_Central_Mediterranean)) +
    ggplot2::annotate("rect",
      xmin = intervention_dates$mare_nostrum_start,
      xmax = intervention_dates$mare_nostrum_end,
      ymin = -Inf, ymax = Inf, fill = "#B8D4E8", alpha = 0.35) +
    ggplot2::annotate("rect",
      xmin = intervention_dates$ngo_sar_start,
      xmax = intervention_dates$eu_libya_start,
      ymin = -Inf, ymax = Inf, fill = "#FFCCCC", alpha = 0.30) +
    ggplot2::annotate("rect",
      xmin = intervention_dates$eu_libya_start,
      xmax = intervention_dates$lcg_sar_zone,
      ymin = -Inf, ymax = Inf, fill = "#FFE4B5", alpha = 0.30) +
    ggplot2::annotate("rect",
      xmin = intervention_dates$lcg_sar_zone,
      xmax = max(df$date, na.rm = TRUE),
      ymin = -Inf, ymax = Inf, fill = "#F5DEB3", alpha = 0.30) +
    ggplot2::annotate("text", x = lubridate::ymd("2014-04-01"), y = 1350,
                      label = "EU\nMare Nostrum", size = 3.2, color = "blue") +
    ggplot2::annotate("text", x = lubridate::ymd("2016-01-01"), y = 1350,
                      label = "NGOs\nSAR", size = 3.2, color = "red") +
    ggplot2::annotate("text", x = lubridate::ymd("2017-05-01"), y = 1350,
                      label = "EU-LCG\nDeal", size = 3.2, color = "darkorange") +
    ggplot2::annotate("text", x = lubridate::ymd("2019-06-01"), y = 1350,
                      label = "Expansion\nLCG SAR-Zone",
                      size = 3.2, color = "darkgreen") +
    ggplot2::geom_line(linewidth = 0.6, color = "black") +
    ggplot2::scale_x_date(breaks = x_breaks, date_labels = "%Y",
                          expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(title = "A. Dead and missing in the Central Mediterranean",
                  x = "Date", y = "Dead and missing (monthly count)") +
    theme_paper()

  fig1b <- ggplot2::ggplot(df, ggplot2::aes(x = date,
                                             y = mortality_rate_100)) +
    ggplot2::annotate("rect",
      xmin = intervention_dates$mare_nostrum_start,
      xmax = intervention_dates$eu_libya_start,
      ymin = -Inf, ymax = Inf, fill = "#E6E6FA", alpha = 0.25) +
    ggplot2::annotate("text", x = lubridate::ymd("2015-06-01"), y = 30,
                      label = "Search-and-rescue", size = 4.5, color = "purple") +
    ggplot2::geom_line(linewidth = 0.6, color = "black") +
    ggplot2::scale_x_date(breaks = x_breaks, date_labels = "%Y",
                          expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(limits = c(0, 50),
                                breaks = seq(0, 50, by = 10)) +
    ggplot2::labs(title = "B. Mortality rate per 100 attempted crossings",
                  x = "Date", y = "Deaths per 100 crossings") +
    theme_paper()

  figure1 <- gridExtra::grid.arrange(
    fig1a, fig1b, ncol = 1, heights = c(1, 1),
    bottom = grid::textGrob("Note: own calculations.",
                            x = 0.5, hjust = 0.5,
                            gp = grid::gpar(fontsize = 9))
  )

  save_figure(figure1, "Figure1_mortality_descriptive",
              output_dir, width = 8.5, height = 8.5)
}


#' Plot counterfactual figures (Figures 2, 6, 8, 9)
#'
#' @param model_a Model A fit (from fit_causal_impact)
#' @param model_b Model B fit
#' @param model_c_rate Model C mortality rate fit
#' @param model_c_deaths Model C death count fit
#' @param output_dir Output directory
#' @return character vector of all saved file paths
plot_counterfactual <- function(model_a, model_b, model_c_rate,
                                model_a_deaths, model_b_deaths,
                                model_c_deaths, output_dir = "output/figures") {
  all_paths <- character(0)

  prep <- function(m) {
    m$results %>%
      dplyr::mutate(date = as.Date(date),
                    original_smooth = moving_average(original, k = 6))
  }

  # Prepare results with smoothing
  results_A <- prep(model_a)
  results_B <- prep(model_b)
  results_C <- prep(model_c_rate)
  results_A_deaths <- prep(model_a_deaths)
  results_B_deaths <- prep(model_b_deaths)
  results_C_deaths <- prep(model_c_deaths)

  date_min <- min(results_A$date, na.rm = TRUE)
  date_max <- max(results_A$date, na.rm = TRUE)
  x_breaks_2y <- seq(as.Date("2012-01-01"), date_max, by = "2 years")

  mn_start  <- as.Date("2013-10-01")
  mn_end    <- as.Date("2014-10-01")
  ngo_start <- as.Date("2014-11-01")
  ngo_end   <- as.Date("2017-02-01")
  eu_start  <- as.Date("2017-02-01")

  # --- Figure 2: Log mortality rate counterfactual ---
  make_cf_panel <- function(results, post_start, post_end, title, ylab, has_grey = FALSE) {
    ylims <- panel_ylim(results, c("original", "original_smooth", "prediction",
                                   "prediction_lower", "prediction_upper"))
    p <- ggplot2::ggplot(results, ggplot2::aes(x = date))

    if (has_grey) {
      p <- p +
        ggplot2::annotate("rect", xmin = post_start, xmax = post_end,
                          ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60) +
        ggplot2::annotate("rect", xmin = post_end, xmax = date_max,
                          ymin = -Inf, ymax = Inf, fill = fill_grey, alpha = 0.60)
    } else {
      p <- p +
        ggplot2::annotate("rect", xmin = post_start, xmax = date_max,
                          ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60)
    }

    p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = prediction_lower,
                                         ymax = prediction_upper),
                           fill = fill_ci, alpha = 0.35) +
      ggplot2::geom_line(ggplot2::aes(y = original_smooth),
                         color = col_smooth, linewidth = 0.8) +
      ggplot2::geom_line(ggplot2::aes(y = original),
                         color = col_obs, linewidth = 0.8) +
      ggplot2::geom_line(ggplot2::aes(y = prediction),
                         linetype = 2, color = col_cf, linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = post_start, color = "darkred",
                          linetype = 2, linewidth = 0.8) +
      ggplot2::coord_cartesian(ylim = ylims) +
      ggplot2::scale_x_date(breaks = x_breaks_2y, date_labels = "%y",
                            limits = c(date_min, date_max),
                            expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
      ggplot2::labs(title = title, x = "", y = ylab) +
      theme_panel()
  }

  fig2a <- make_cf_panel(results_A, mn_start, mn_end,
    "A. Mare Nostrum (state-led SAR and anti-smuggler operations)",
    "Log mortality rate (deaths per 100 crossings)", has_grey = TRUE)
  fig2b <- make_cf_panel(results_B, ngo_start, ngo_end,
    "B. NGOs (private-led SAR by various actors)",
    "Log mortality rate (deaths per 100 crossings)", has_grey = TRUE)
  fig2c <- make_cf_panel(results_C, eu_start, date_max,
    "C. EU-Libya cooperation / MoU (pushbacks and LCG SAR-zone extension)",
    "Log mortality rate (deaths per 100 crossings)") +
    ggplot2::labs(x = "Date")

  note_fig2 <- paste(
    "Note: own calculations.",
    "Corrected specification: ERA5 sea conditions (lag 0) + weather + oil price + month dummies.",
    "Y-axis clipped to central quantiles to keep trajectories readable."
  )

  figure2 <- gridExtra::grid.arrange(
    fig2a, fig2b, fig2c, ncol = 1,
    bottom = grid::textGrob(note_fig2, x = 0.01, hjust = 0,
                            gp = grid::gpar(fontsize = 9))
  )

  all_paths <- c(all_paths,
    save_figure(figure2, "Figure2_mortality_counterfactual",
                output_dir, width = 8, height = 10))

  # --- Figure 6: Death count counterfactual (3 panels) ---
  fig6a <- make_cf_panel(results_A_deaths, mn_start, mn_end,
    "A. Mare Nostrum (death count)",
    "Log deaths (dead and missing + 1)", has_grey = TRUE)
  fig6b <- make_cf_panel(results_B_deaths, ngo_start, ngo_end,
    "B. NGO search-and-rescue (death count)",
    "Log deaths (dead and missing + 1)", has_grey = TRUE)
  fig6c <- make_cf_panel(results_C_deaths, eu_start, date_max,
    "C. EU-Libya cooperation / MoU (death count)",
    "Log deaths (dead and missing + 1)") +
    ggplot2::labs(x = "Date")

  figure6 <- gridExtra::grid.arrange(
    fig6a, fig6b, fig6c, ncol = 1,
    bottom = grid::textGrob(
      "Note: own calculations. Death-count outcome: log(dead and missing + 1).",
      x = 0.01, hjust = 0, gp = grid::gpar(fontsize = 9))
  )

  all_paths <- c(all_paths,
    save_figure(figure6, "Figure6_deaths_counterfactual",
                output_dir, width = 8, height = 10))

  # --- Natural-scale figures ---
  results_A_nat <- add_natural_scale(results_A, "rate")
  results_B_nat <- add_natural_scale(results_B, "rate")
  results_C_nat <- add_natural_scale(results_C, "rate")
  results_A_deaths_nat <- add_natural_scale(results_A_deaths, "deaths")
  results_B_deaths_nat <- add_natural_scale(results_B_deaths, "deaths")
  results_C_deaths_nat <- add_natural_scale(results_C_deaths, "deaths")

  make_nat_panel <- function(results, post_start, post_end, title, ylab,
                              has_grey = FALSE) {
    ylims <- panel_ylim(results, c("original_nat", "original_smooth_nat",
                                   "prediction_nat", "prediction_lower_nat",
                                   "prediction_upper_nat"))
    p <- ggplot2::ggplot(results, ggplot2::aes(x = date))
    if (has_grey) {
      p <- p +
        ggplot2::annotate("rect", xmin = post_start, xmax = post_end,
                          ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60) +
        ggplot2::annotate("rect", xmin = post_end, xmax = date_max,
                          ymin = -Inf, ymax = Inf, fill = fill_grey, alpha = 0.60)
    } else {
      p <- p +
        ggplot2::annotate("rect", xmin = post_start, xmax = date_max,
                          ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60)
    }
    p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = prediction_lower_nat,
                                         ymax = prediction_upper_nat),
                           fill = fill_ci, alpha = 0.35) +
      ggplot2::geom_line(ggplot2::aes(y = original_smooth_nat),
                         color = col_smooth, linewidth = 0.8) +
      ggplot2::geom_line(ggplot2::aes(y = original_nat),
                         color = col_obs, linewidth = 0.8) +
      ggplot2::geom_line(ggplot2::aes(y = prediction_nat),
                         linetype = 2, color = col_cf, linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = post_start, color = "darkred",
                          linetype = 2, linewidth = 0.8) +
      ggplot2::coord_cartesian(ylim = ylims) +
      ggplot2::scale_x_date(breaks = x_breaks_2y, date_labels = "%y",
                            limits = c(date_min, date_max),
                            expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
      ggplot2::labs(title = title, x = "", y = ylab) +
      theme_panel()
  }

  fig8a <- make_nat_panel(results_A_nat, mn_start, mn_end,
    "A. Mare Nostrum (natural units)",
    "Mortality rate (deaths per 100 crossings)", has_grey = TRUE)
  fig8b <- make_nat_panel(results_B_nat, ngo_start, ngo_end,
    "B. NGO search-and-rescue (natural units)",
    "Mortality rate (deaths per 100 crossings)", has_grey = TRUE)
  fig8c <- make_nat_panel(results_C_nat, eu_start, date_max,
    "C. EU-Libya cooperation / MoU (natural units)",
    "Mortality rate (deaths per 100 crossings)") +
    ggplot2::labs(x = "Date")

  figure8 <- gridExtra::grid.arrange(
    fig8a, fig8b, fig8c, ncol = 1,
    bottom = grid::textGrob(
      paste0("Note: own calculations. Back-transformed from log scale. ",
             "Back-transform capped at log=", BACKTRANSFORM_LOG_CAP,
             " for visual stability."),
      x = 0.01, hjust = 0, gp = grid::gpar(fontsize = 9))
  )

  all_paths <- c(all_paths,
    save_figure(figure8, "Figure8_mortality_counterfactual_natural",
                output_dir, width = 8, height = 10))

  fig9a <- make_nat_panel(results_A_deaths_nat, mn_start, mn_end,
    "A. Mare Nostrum (death count, natural units)",
    "Deaths (monthly count)", has_grey = TRUE)
  fig9b <- make_nat_panel(results_B_deaths_nat, ngo_start, ngo_end,
    "B. NGO search-and-rescue (death count, natural units)",
    "Deaths (monthly count)", has_grey = TRUE)
  fig9c <- make_nat_panel(results_C_deaths_nat, eu_start, date_max,
    "C. EU-Libya cooperation / MoU (death count, natural units)",
    "Deaths (monthly count)") +
    ggplot2::labs(x = "Date")

  figure9 <- gridExtra::grid.arrange(
    fig9a, fig9b, fig9c, ncol = 1,
    bottom = grid::textGrob(
      paste0("Note: own calculations. Back-transformed from log(deaths + 1). ",
             "Back-transform capped at log=", BACKTRANSFORM_LOG_CAP,
             " for visual stability."),
      x = 0.01, hjust = 0, gp = grid::gpar(fontsize = 9))
  )

  all_paths <- c(all_paths,
    save_figure(figure9, "Figure9_deaths_counterfactual_natural",
                output_dir, width = 8, height = 10))

  all_paths
}


#' Plot pointwise effects (Figures 3, 7)
#'
#' @param model_a Model A fit
#' @param model_b Model B fit
#' @param model_c_rate Model C mortality rate fit
#' @param model_c_deaths Model C death count fit
#' @param output_dir Output directory
#' @return character vector of saved file paths
plot_pointwise <- function(model_a, model_b, model_c_rate,
                           model_a_deaths, model_b_deaths,
                           model_c_deaths, output_dir = "output/figures") {
  all_paths <- character(0)

  prep <- function(m) m$results %>% dplyr::mutate(date = as.Date(date))

  results_A <- prep(model_a)
  results_B <- prep(model_b)
  results_C <- prep(model_c_rate)
  results_A_deaths <- prep(model_a_deaths)
  results_B_deaths <- prep(model_b_deaths)
  results_C_deaths <- prep(model_c_deaths)

  date_min <- min(results_A$date, na.rm = TRUE)
  date_max <- max(results_A$date, na.rm = TRUE)
  x_breaks_2y <- seq(as.Date("2012-01-01"), date_max, by = "2 years")

  mn_start  <- as.Date("2013-10-01")
  mn_end    <- as.Date("2014-10-01")
  ngo_start <- as.Date("2014-11-01")
  ngo_end   <- as.Date("2017-02-01")
  eu_start  <- as.Date("2017-02-01")

  make_pw_panel <- function(results, post_start, post_end, title, ylab,
                            has_grey = FALSE) {
    ylims <- panel_ylim(results, c("pointwise_effect",
                                   "pointwise_effect_lower",
                                   "pointwise_effect_upper"))
    p <- ggplot2::ggplot(results, ggplot2::aes(x = date))
    if (has_grey) {
      p <- p +
        ggplot2::annotate("rect", xmin = post_start, xmax = post_end,
                          ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60) +
        ggplot2::annotate("rect", xmin = post_end, xmax = date_max,
                          ymin = -Inf, ymax = Inf, fill = fill_grey, alpha = 0.60)
    } else {
      p <- p +
        ggplot2::annotate("rect", xmin = post_start, xmax = date_max,
                          ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60)
    }
    p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = pointwise_effect_lower,
                                         ymax = pointwise_effect_upper),
                           fill = fill_ci, alpha = 0.35) +
      ggplot2::geom_line(ggplot2::aes(y = pointwise_effect),
                         linetype = 2, color = col_cf, linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = 0, color = "black",
                          linetype = 2, linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = post_start, color = "darkred",
                          linetype = 2, linewidth = 0.8) +
      ggplot2::coord_cartesian(ylim = ylims) +
      ggplot2::scale_x_date(breaks = x_breaks_2y, date_labels = "%y",
                            limits = c(date_min, date_max),
                            expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
      ggplot2::labs(title = title, x = "", y = ylab) +
      theme_panel()
  }

  fig3a <- make_pw_panel(results_A, mn_start, mn_end,
    "A. Mare Nostrum", "Pointwise effect (log mortality-rate units)",
    has_grey = TRUE)
  fig3b <- make_pw_panel(results_B, ngo_start, ngo_end,
    "B. NGO search-and-rescue", "Pointwise effect (log mortality-rate units)",
    has_grey = TRUE)
  fig3c <- make_pw_panel(results_C, eu_start, date_max,
    "C. EU-Libya cooperation / MoU",
    "Pointwise effect (log mortality-rate units)") +
    ggplot2::labs(x = "Date")

  note_fig3 <- paste(
    "Note: own calculations.",
    "Positive values = observed mortality HIGHER than predicted (route more dangerous).",
    "Y-axis clipped to central quantiles to keep trajectories readable."
  )

  figure3 <- gridExtra::grid.arrange(
    fig3a, fig3b, fig3c, ncol = 1,
    bottom = grid::textGrob(note_fig3, x = 0.01, hjust = 0,
                            gp = grid::gpar(fontsize = 9))
  )

  all_paths <- c(all_paths,
    save_figure(figure3, "Figure3_mortality_pointwise",
                output_dir, width = 8, height = 10))

  # --- Figure 7: Death count pointwise (3 panels) ---
  fig7a <- make_pw_panel(results_A_deaths, mn_start, mn_end,
    "A. Mare Nostrum (death count)",
    "Pointwise effect (log death-count units)", has_grey = TRUE)
  fig7b <- make_pw_panel(results_B_deaths, ngo_start, ngo_end,
    "B. NGO search-and-rescue (death count)",
    "Pointwise effect (log death-count units)", has_grey = TRUE)
  fig7c <- make_pw_panel(results_C_deaths, eu_start, date_max,
    "C. EU-Libya cooperation / MoU (death count)",
    "Pointwise effect (log death-count units)") +
    ggplot2::labs(x = "Date")

  figure7 <- gridExtra::grid.arrange(
    fig7a, fig7b, fig7c, ncol = 1,
    bottom = grid::textGrob(
      "Note: own calculations. Positive values = observed deaths higher than predicted.",
      x = 0.01, hjust = 0, gp = grid::gpar(fontsize = 9))
  )

  all_paths <- c(all_paths,
    save_figure(figure7, "Figure7_deaths_pointwise",
                output_dir, width = 8, height = 10))

  all_paths
}


#' Plot placebo diagnostic distributions (Figure 4)
#'
#' @param placebos_a Placebo results for Model A
#' @param placebos_b Placebo results for Model B
#' @param placebos_c_rate Placebo results for Model C (mortality rate)
#' @param placebos_c_deaths Placebo results for Model C (death count)
#' @param model_a Model A fit
#' @param model_b Model B fit
#' @param model_c_rate Model C rate fit
#' @param model_c_deaths Model C deaths fit
#' @param output_dir Output directory
#' @return character vector of saved file paths
plot_placebo_distributions <- function(placebos_a, placebos_b,
                                       placebos_c_rate,
                                       placebos_a_deaths, placebos_b_deaths,
                                       placebos_c_deaths,
                                       model_a, model_b, model_c_rate,
                                       model_a_deaths, model_b_deaths,
                                       model_c_deaths,
                                       output_dir = "output/figures") {
  all_placebos <- dplyr::bind_rows(placebos_a, placebos_b, placebos_c_rate,
                                   placebos_a_deaths, placebos_b_deaths,
                                   placebos_c_deaths)

  all_models <- list(model_a, model_b, model_c_rate,
                     model_a_deaths, model_b_deaths, model_c_deaths)

  actual_results <- tibble::tibble(
    model    = sapply(all_models, function(m) m$label),
    actual_p = sapply(all_models, function(m) {
      as.numeric(m$impact$summary["Cumulative", "p"])
    })
  )

  model_labels <- c(
    "A_mortality" = "A. Mare Nostrum\n(mortality rate)",
    "B_mortality" = "B. NGO SAR\n(mortality rate)",
    "C_mortality" = "C. MoU\n(mortality rate)",
    "A_deaths"    = "A. Mare Nostrum\n(death count)",
    "B_deaths"    = "B. NGO SAR\n(death count)",
    "C_deaths"    = "C. MoU\n(death count)"
  )

  plots <- list()
  for (mod in names(model_labels)) {
    plac   <- all_placebos %>% dplyr::filter(model == mod)
    actual <- actual_results %>% dplyr::filter(model == mod)
    if (nrow(plac) == 0) next

    p <- ggplot2::ggplot(plac, ggplot2::aes(x = p_value)) +
      ggplot2::geom_histogram(binwidth = 0.05, fill = "grey70",
                              color = "grey40", boundary = 0) +
      ggplot2::geom_vline(xintercept = actual$actual_p, color = "red",
                          linetype = "dashed", linewidth = 1) +
      ggplot2::geom_vline(xintercept = 0.05, color = "blue",
                          linetype = "dotted", linewidth = 0.8) +
      ggplot2::annotate("text", x = actual$actual_p + 0.03, y = Inf,
                        label = paste0("Actual\np=", round(actual$actual_p, 3)),
                        color = "red", vjust = 1.5, size = 3.5) +
      ggplot2::scale_x_continuous(limits = c(0, 1),
                                  breaks = seq(0, 1, 0.2)) +
      ggplot2::labs(title = model_labels[mod],
                    x = "Placebo p-value", y = "Count") +
      theme_diag()

    plots[[mod]] <- p
  }

  if (length(plots) == 0) return(character(0))

  # Order: rate models first column, death models second column
  plot_order <- c("A_mortality", "A_deaths",
                  "B_mortality", "B_deaths",
                  "C_mortality", "C_deaths")
  ordered_plots <- plots[intersect(plot_order, names(plots))]

  fig_placebo <- do.call(gridExtra::grid.arrange, c(ordered_plots, list(
    ncol = 2,
    top = grid::textGrob(
      "Placebo Test Results: Distribution of p-values under no intervention",
      gp = grid::gpar(fontsize = 14, fontface = "bold")),
    bottom = grid::textGrob(
      "Grey: placebo p-values | Red dashed: actual intervention p-value | Blue dotted: 0.05 threshold",
      gp = grid::gpar(fontsize = 9))
  )))

  save_figure(fig_placebo, "Figure4_placebo_diagnostics",
              output_dir, width = 10, height = 12)
}


#' Plot truncation test results (Figure 5)
#'
#' @param truncation_results Output of run_truncation_test()
#' @param output_dir Output directory
#' @return character vector of saved file paths
plot_truncation <- function(truncation_results,
                            output_dir = "output/figures") {
  if (nrow(truncation_results) == 0) return(character(0))

  fig_trunc <- ggplot2::ggplot(truncation_results,
                               ggplot2::aes(x = post_months, y = p_value)) +
    ggplot2::geom_point(size = 4, color = "darkred") +
    ggplot2::geom_line(color = "darkred", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0.05, color = "blue",
                        linetype = "dotted", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -1.2, size = 3.5) +
    ggplot2::scale_y_continuous(
      limits = c(0, max(truncation_results$p_value) * 1.3),
      breaks = seq(0, 1, 0.1)) +
    ggplot2::labs(
      title = "Model A: Effect of Post-Period Length on Significance",
      subtitle = "Does the Mare Nostrum result hold when we exclude later interventions?",
      x = "Post-period length (months)",
      y = "p-value"
    ) +
    theme_diag()

  save_figure(fig_trunc, "Figure5_model_A_truncation",
              output_dir, width = 8, height = 5)
}


#' Plot forecast diagnostic figures (Approach A and B)
#'
#' @param forecast_diag Output of run_forecast_diagnostic()
#' @param cov_cols_diag Character vector of covariate column names
#' @param output_dir Output directory
#' @return character vector of saved file paths
plot_forecast_diagnostic <- function(forecast_diag,
                                     output_dir = "output/figures") {
  all_paths <- character(0)
  diag_results_list <- forecast_diag$fold_details
  n_cov <- forecast_diag$settings$n_cov

  # Approach A
  plot_data_a <- dplyr::bind_rows(lapply(diag_results_list, function(r) {
    series <- c("Actual", "M0a (month mean)", "M0c (state-only)",
                "M1-A1 (LL+dummies)", "M1-A2 (LLT+dummies)")
    vals <- c(r$actuals, r$preds[["M0a (month mean)"]],
              r$preds[["M0c (state-only)"]],
              r$preds[["M1-A1 (LL+dummies)"]],
              r$preds[["M1-A2 (LLT+dummies)"]])
    data.frame(
      date   = rep(r$dates, length(series)),
      value  = vals,
      series = rep(series, each = length(r$dates)),
      fold   = r$fold,
      stringsAsFactors = FALSE
    )
  }))

  fig_diag_a <- ggplot2::ggplot(plot_data_a,
    ggplot2::aes(x = date, y = value, color = series, linetype = series)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~fold, scales = "free_x", ncol = 1) +
    ggplot2::scale_color_manual(values = c(
      "Actual"               = "black",
      "M0a (month mean)"     = "steelblue",
      "M0c (state-only)"     = "darkorange",
      "M1-A1 (LL+dummies)"   = "red3",
      "M1-A2 (LLT+dummies)"  = "purple"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "Actual"               = "solid",
      "M0a (month mean)"     = "dashed",
      "M0c (state-only)"     = "dotdash",
      "M1-A1 (LL+dummies)"   = "dashed",
      "M1-A2 (LLT+dummies)"  = "dotted"
    )) +
    ggplot2::labs(
      title    = "Forecasting Diagnostic v2 -- Approach A (month factor dummies)",
      subtitle = paste0("log(mortality_rate + 0.01) | ", n_cov,
                        " exogenous + 11 month dummies"),
      x = "Date", y = "log(mortality rate + 0.01)",
      color = "Series", linetype = "Series"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(face = "bold"),
                   strip.text = ggplot2::element_text(face = "bold", size = 11))

  all_paths <- c(all_paths,
    save_figure(fig_diag_a, "diagnostic_v2_approach_A",
                output_dir, width = 10, height = 9))

  # Approach B
  plot_data_b <- dplyr::bind_rows(lapply(diag_results_list, function(r) {
    series <- c("Actual", "M0a (month mean)", "M0c (state-only)",
                "M1-B (LL+Seasonal)")
    vals <- c(r$actuals, r$preds[["M0a (month mean)"]],
              r$preds[["M0c (state-only)"]],
              r$preds[["M1-B (LL+Seasonal)"]])
    data.frame(
      date   = rep(r$dates, length(series)),
      value  = vals,
      series = rep(series, each = length(r$dates)),
      fold   = r$fold,
      stringsAsFactors = FALSE
    )
  }))

  fig_diag_b <- ggplot2::ggplot(plot_data_b,
    ggplot2::aes(x = date, y = value, color = series, linetype = series)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~fold, scales = "free_x", ncol = 1) +
    ggplot2::scale_color_manual(values = c(
      "Actual"              = "black",
      "M0a (month mean)"    = "steelblue",
      "M0c (state-only)"    = "darkorange",
      "M1-B (LL+Seasonal)"  = "forestgreen"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "Actual"              = "solid",
      "M0a (month mean)"    = "dashed",
      "M0c (state-only)"    = "dotdash",
      "M1-B (LL+Seasonal)"  = "solid"
    )) +
    ggplot2::labs(
      title    = "Forecasting Diagnostic v2 -- Approach B (seasonal state component)",
      subtitle = paste0("log(mortality_rate + 0.01) | AddLocalLevel + AddSeasonal(12) + ",
                        n_cov, " exogenous covariates"),
      x = "Date", y = "log(mortality rate + 0.01)",
      color = "Series", linetype = "Series"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(face = "bold"),
                   strip.text = ggplot2::element_text(face = "bold", size = 11))

  all_paths <- c(all_paths,
    save_figure(fig_diag_b, "diagnostic_v2_approach_B",
                output_dir, width = 10, height = 9))

  all_paths
}
