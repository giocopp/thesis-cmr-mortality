library(tidyverse)
library(lubridate)
library(gridExtra)
library(scales)
library(grid)

# Output directories
dir.create("replicated_results/figures/pdf", recursive = TRUE, showWarnings = FALSE)
dir.create("replicated_results/figures/png", recursive = TRUE, showWarnings = FALSE)

# Load raw data
df <- readRDS(file = "../Original code and data/df.RDS")
if (file.exists("../Original code and data/df_y_rec.RDS")) {
  df_y_rec <- readRDS(file = "../Original code and data/df_y_rec.RDS")
  df <- left_join(df, df_y_rec, by = "date")
}
if (!"y_rec" %in% names(df)) df$y_rec <- 0

df <- df %>%
  mutate(
    LCG_pushbacks_count = as.numeric(ifelse(is.na(LCG_pushbacks_count), 0, LCG_pushbacks_count)),
    TCG_pushbacks_count = as.numeric(ifelse(is.na(TCG_pushbacks_count), 0, TCG_pushbacks_count)),
    dead_and_missing_Central_Mediterranean = as.numeric(ifelse(is.na(dead_and_missing_Central_Mediterranean), 0, dead_and_missing_Central_Mediterranean)),
    crossings_CMR = arrivals_CMR + LCG_pushbacks_count + TCG_pushbacks_count + dead_and_missing_Central_Mediterranean,
    mortality_rate_100 = (dead_and_missing_Central_Mediterranean / crossings_CMR) * 100
  )

# Load saved CausalImpact model results
model_df_results_A <- readRDS("replicated_results/data/model_df_results_A_all.RDS")
model_df_results_B <- readRDS("replicated_results/data/model_df_results_B_all.RDS")
model_df_results_C <- readRDS("replicated_results/data/model_df_results_C_all.RDS")

###############################################################################
# FIGURE 1: Time-series of crossing attempts, mortality rate, and
#            main intervention periods, 2009-2021
###############################################################################

# Intervention period dates
mn_start   <- as.Date("2013-10-01")
mn_end     <- as.Date("2014-10-31")
ngo_start  <- as.Date("2014-11-01")
ngo_end    <- as.Date("2017-01-31")
eulcg_date <- as.Date("2017-02-01")
exp_start  <- as.Date("2017-08-01")
exp_end    <- max(df$date, na.rm = TRUE)

# Panel A: Attempted Crossings with intervention period shading
fig1A <- ggplot(df, aes(x = date, y = crossings_CMR)) +
  annotate("rect", xmin = mn_start, xmax = mn_end,
           ymin = -Inf, ymax = Inf, fill = "#C5CAE9", alpha = 0.5) +
  annotate("rect", xmin = ngo_start, xmax = ngo_end,
           ymin = -Inf, ymax = Inf, fill = "#FFCCBC", alpha = 0.45) +
  annotate("rect", xmin = exp_start, xmax = exp_end,
           ymin = -Inf, ymax = Inf, fill = "#D7CCC8", alpha = 0.5) +
  # Labels
  annotate("label", x = mn_start + days(180), y = 29000,
           label = "EU\nMare Nostrum", color = "darkgreen", size = 2.8,
           fontface = "bold", fill = NA, label.size = 0.3) +
  annotate("label", x = ngo_start + days(400), y = 29000,
           label = "NGOs\nsearch-and-rescue", color = "darkorange3", size = 2.8,
           fontface = "bold", fill = NA, label.size = 0.3) +
  annotate("label", x = eulcg_date + days(40), y = 29000,
           label = "EU-LCG\nDeal", color = "red3", size = 2.8,
           fontface = "bold", fill = NA, label.size = 0.3) +
  annotate("label", x = exp_start + days(550), y = 29000,
           label = "Expansion\nLCG SAR-Zone", color = "brown", size = 2.8,
           fontface = "bold", fill = NA, label.size = 0.3) +
  geom_line(linewidth = 0.4) +
  scale_x_date(breaks = seq(as.Date("2010-01-01"), as.Date("2020-01-01"), by = "5 years"),
               date_labels = "%Y") +
  scale_y_continuous(labels = scales::comma, breaks = seq(0, 30000, 10000)) +
  coord_cartesian(ylim = c(0, 31000)) +
  labs(
    title = "A. Attempted Crossings (arrivals, deaths, and pushbacks),\n    intervention periods, and structural breakpoints",
    x = "", y = "Attempted crossings"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    plot.margin = margin(5, 10, 5, 5)
  )

# Panel B: Mortality rate per 100 attempted crossings
fig1B <- ggplot(df, aes(x = date, y = mortality_rate_100)) +
  annotate("rect", xmin = mn_start, xmax = ngo_end,
           ymin = -Inf, ymax = Inf, fill = "#E1BEE7", alpha = 0.35) +
  annotate("text", x = as.Date("2015-06-01"), y = 30,
           label = "Search-and-rescue", color = "purple3", size = 5,
           family = "mono") +
  geom_line(linewidth = 0.4) +
  scale_x_date(breaks = seq(as.Date("2010-01-01"), as.Date("2020-01-01"), by = "5 years"),
               date_labels = "%Y") +
  scale_y_continuous(breaks = seq(0, 50, 10)) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(
    title = "B. Mortality rate per 100 attempted crossings",
    x = "Date",
    y = expression(paste("Rate per 100 crossings", D[t]/Y[t], "*100"))
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    plot.margin = margin(5, 10, 5, 5)
  )

figure1 <- arrangeGrob(
  fig1A, fig1B, ncol = 1, heights = c(1, 1),
  bottom = textGrob("Note: own calculations.", x = 0.02, hjust = 0,
                    gp = gpar(fontsize = 9, fontface = "italic"))
)

ggsave("replicated_results/figures/pdf/Figure1.pdf", figure1, width = 10, height = 10)
ggsave("replicated_results/figures/png/Figure1.png", figure1, width = 10, height = 10, dpi = 300)
cat("Figure 1 saved.\n")

###############################################################################
# FIGURE 2: Intervention periods and the predicted counterfactual and
#            observed time-series (log scale), 2011-2020
###############################################################################

date_min <- min(model_df_results_A$date, na.rm = TRUE)
date_max <- max(model_df_results_A$date, na.rm = TRUE)
x_breaks <- seq(as.Date("2012-01-01"), date_max, by = "2 years")

fill_pink <- "#FDE6E6"
fill_grey <- "#D9D9D9"
fill_ci   <- "#B7C9E2"
col_cf    <- "#4A4FB0"

theme_fig <- function() {
  theme_classic(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", size = 12, hjust = 0),
      axis.title    = element_text(size = 12, color = "black"),
      axis.text     = element_text(size = 10, color = "black"),
      axis.text.x   = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.line     = element_line(color = "black"),
      plot.margin   = margin(6, 10, 6, 8)
    )
}

make_cf_panel <- function(results, post_start, post_end, title,
                          ylab = "", has_grey = FALSE, show_xlabel = FALSE) {
  p <- ggplot(results, aes(x = date))

  if (has_grey) {
    p <- p +
      annotate("rect", xmin = post_start, xmax = post_end,
               ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.55) +
      annotate("rect", xmin = post_end, xmax = date_max,
               ymin = -Inf, ymax = Inf, fill = fill_grey, alpha = 0.50)
  } else {
    p <- p +
      annotate("rect", xmin = post_start, xmax = date_max,
               ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.55)
  }

  p +
    geom_ribbon(aes(ymin = prediction_lower, ymax = prediction_upper),
                fill = fill_ci, alpha = 0.4) +
    geom_line(aes(y = original), color = "black", linewidth = 0.6) +
    geom_line(aes(y = prediction), linetype = "dashed", color = col_cf, linewidth = 0.6) +
    geom_vline(xintercept = post_start, color = "darkred", linetype = "dashed", linewidth = 0.8) +
    coord_cartesian(ylim = c(-25, 45)) +
    scale_x_date(
      breaks = x_breaks, date_labels = "%Y",
      limits = c(date_min, date_max),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = title,
      x = if (show_xlabel) "Date" else "",
      y = ylab
    ) +
    theme_fig()
}

plotA <- make_cf_panel(
  model_df_results_A,
  as.Date("2013-10-01"), as.Date("2014-10-01"),
  "A. Mare Nostrum (state-led search-and-rescue and\n    anti-smuggler operations)",
  has_grey = TRUE
)

plotB <- make_cf_panel(
  model_df_results_B,
  as.Date("2014-11-01"), as.Date("2017-02-01"),
  "B. NGOs (private-led search-and-rescue\n    by various actors)",
  ylab = "Log of attempted crossings",
  has_grey = TRUE
)

plotC <- make_cf_panel(
  model_df_results_C,
  as.Date("2017-02-01"), date_max,
  "C. EU and Libya cooperation (pushbacks and Libyan\n    search-and-rescue zone extension)",
  show_xlabel = TRUE
)

note_fig2 <- "Note: own calculations. A. Mare nostrum = not significant; B. NGO-led search-and-rescue = not significant; and C. Pushbacks = significant."

figure2 <- arrangeGrob(
  plotA, plotB, plotC, ncol = 1,
  bottom = textGrob(note_fig2, x = 0.01, hjust = 0,
                    gp = gpar(fontsize = 9, fontface = "italic"))
)

ggsave("replicated_results/figures/pdf/Figure2.pdf", figure2, width = 10, height = 13, bg = "white")
ggsave("replicated_results/figures/png/Figure2.png", figure2, width = 10, height = 13, dpi = 300, bg = "white")
cat("Figure 2 saved.\n")

###############################################################################
# FIGURE 3: Pointwise effects (difference between observed and predicted)
#            for the three intervention periods, 2011-2020
###############################################################################

make_pw_panel <- function(results, post_start, post_end, title,
                          ylab = "", has_grey = FALSE, show_xlabel = FALSE) {
  p <- ggplot(results, aes(x = date))

  if (has_grey) {
    p <- p +
      annotate("rect", xmin = post_start, xmax = post_end,
               ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.55) +
      annotate("rect", xmin = post_end, xmax = date_max,
               ymin = -Inf, ymax = Inf, fill = fill_grey, alpha = 0.50)
  } else {
    p <- p +
      annotate("rect", xmin = post_start, xmax = date_max,
               ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.55)
  }

  p +
    geom_ribbon(aes(ymin = pointwise_effect_lower, ymax = pointwise_effect_upper),
                fill = fill_ci, alpha = 0.4) +
    geom_line(aes(y = pointwise_effect), linetype = "dashed", color = col_cf, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_vline(xintercept = post_start, color = "darkred", linetype = "dashed", linewidth = 0.8) +
    coord_cartesian(ylim = c(-35, 25)) +
    scale_x_date(
      breaks = x_breaks, date_labels = "%Y",
      limits = c(date_min, date_max),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = title,
      x = if (show_xlabel) "Date" else "",
      y = ylab
    ) +
    theme_fig()
}

pwA <- make_pw_panel(
  model_df_results_A,
  as.Date("2013-10-01"), as.Date("2014-10-01"),
  "A. Mare Nostrum (state-led search-and-rescue and\n    anti-smuggler operations)",
  has_grey = TRUE
)

pwB <- make_pw_panel(
  model_df_results_B,
  as.Date("2014-11-01"), as.Date("2017-02-01"),
  "B. NGOs (private-led search-and-rescue\n    by various actors)",
  ylab = "Diff. between observed and predicted",
  has_grey = TRUE
)

pwC <- make_pw_panel(
  model_df_results_C,
  as.Date("2017-02-01"), date_max,
  "C. EU and Libya cooperation (pushbacks and Libyan\n    search-and-rescue zone extension)",
  show_xlabel = TRUE
)

note_fig3 <- "Note: own calculations. A. Mare nostrum = not significant; B. NGO-led search-and-rescue = not significant; and C. Pushbacks = significant."

figure3 <- arrangeGrob(
  pwA, pwB, pwC, ncol = 1,
  bottom = textGrob(note_fig3, x = 0.01, hjust = 0,
                    gp = gpar(fontsize = 9, fontface = "italic"))
)

ggsave("replicated_results/figures/pdf/Figure3.pdf", figure3, width = 10, height = 13, bg = "white")
ggsave("replicated_results/figures/png/Figure3.png", figure3, width = 10, height = 13, dpi = 300, bg = "white")
cat("Figure 3 saved.\n")

cat("\nAll 3 figures saved to replicated_results/figures/\n")
