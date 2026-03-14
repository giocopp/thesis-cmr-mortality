# plot_correlation_stability_panels.R
# ==================================
# Creates pre/post MoU correlation dot plots for:
#   (a) mortality rate  — matches existing 17_correlation_stability.png
#   (b) death count
#   (c) combined 2-panel version for the presentation slide

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --- Load saved data ---
cor_stability <- read.csv("output/tables/eda/correlations_pre_post_mou.csv",
                          stringsAsFactors = FALSE)

# --- Label ordering (match the original) ---
label_order <- c(
  "Wave height", "Days >2m waves", "Wave max", "Wave period",
  "Wave SD", "Precip. (coast)", "Cloud cover", "Cloud (coast)",
  "Low cloud", "Wind speed", "Current speed", "Wave direction",
  "Wind (coast)", "Current (opposing)", "Dewpoint dep.",
  "SST anomaly", "SST", "Air temp.", "Temp. (coast)"
)

theme_eda <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40"))

# --- Helper function ---
make_stability_plot <- function(data, r_col, x_label, title, subtitle) {
  df_wide <- data %>%
    select(label, period, r = all_of(r_col)) %>%
    pivot_wider(names_from = period, values_from = r) %>%
    mutate(label = factor(label, levels = label_order))

  ggplot(df_wide, aes(y = label)) +
    geom_segment(aes(x = `Pre-MoU`, xend = `Post-MoU`, yend = label),
                 color = "grey70", linewidth = 0.4) +
    geom_point(aes(x = `Pre-MoU`), color = "steelblue", size = 2.5) +
    geom_point(aes(x = `Post-MoU`), color = "firebrick", size = 2.5) +
    geom_vline(xintercept = 0, color = "grey50") +
    labs(x = x_label, y = NULL, title = title, subtitle = subtitle) +
    theme_eda +
    theme(axis.text.y = element_text(size = 8))
}

# --- (a) Mortality rate (standalone) ---
p_rate <- make_stability_plot(
  cor_stability, "r_rate",
  "Pearson r with log(mortality rate)",
  "Mortality Rate",
  "Blue = pre-MoU, red = post-MoU"
)

# --- (b) Death count (standalone) ---
p_deaths <- make_stability_plot(
  cor_stability, "r_deaths",
  "Pearson r with log(deaths + 1)",
  "Death Count",
  "Blue = pre-MoU, red = post-MoU"
)

# Save standalone death count figure
OUT_DIR <- "output/figures/eda"
ggsave(file.path(OUT_DIR, "17b_correlation_stability_deaths.png"), p_deaths,
       width = 9, height = 7, dpi = 200, bg = "white")
cat("Saved 17b_correlation_stability_deaths.png\n")

# --- (c) Combined panel ---
p_combined <- (p_rate + p_deaths) +
  plot_annotation(
    title = "Correlation Stability: Pre-MoU vs. Post-MoU",
    subtitle = "Blue = pre-MoU, red = post-MoU. 18/19 sign flips for rate; deaths retain pre-MoU signs.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave(file.path(OUT_DIR, "17c_correlation_stability_combined.png"), p_combined,
       width = 16, height = 7, dpi = 200, bg = "white")
cat("Saved 17c_correlation_stability_combined.png\n")
