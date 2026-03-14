# plot_conditioning_results.R
# Generate figure for conditioning on crossings sensitivity analysis

library(dplyr)
library(ggplot2)
library(targets)
library(CausalImpact)
library(gridExtra)
library(grid)

PROJECT_DIR <- "/Users/giocopp/Desktop/Uni/Hertie School/6th Semester/Thesis-MDS/Rodriguez-Sanchez-paper-replication/Extension-2-new-data"

# Load conditioned results
load(file.path(PROJECT_DIR, "output", "conditioning_results.RData"))

# Load original results for comparison
withr::with_dir(PROJECT_DIR, {
  orig_a <- tar_read(model_a_deaths)
  orig_b <- tar_read(model_b_deaths)
  orig_c <- tar_read(model_c_deaths)
})

# --- Helper: extract plot data from CausalImpact object ---
extract_plot_data <- function(impact_obj, dates, label) {
  ser <- impact_obj$series
  # Find pre/post boundary
  pre_idx <- which(ser$cum.effect == 0)
  pre_n <- max(pre_idx)

  data.frame(
    date = dates,
    observed = as.numeric(ser$response),
    predicted = as.numeric(ser$point.pred),
    pred_lower = as.numeric(ser$point.pred.lower),
    pred_upper = as.numeric(ser$point.pred.upper),
    is_post = seq_along(dates) > pre_n,
    label = label
  )
}

# --- Build comparison figure: C deaths original vs conditioned ---

# Get dates from the conditioned model
withr::with_dir(PROJECT_DIR, {
  df_model <- tar_read(df_model)
})
dates_deaths <- df_model$df_deaths$date

# Extract C deaths original
d_orig <- extract_plot_data(
  orig_c$impact, dates_deaths, "Original (no crossings covariate)"
)

# Extract C deaths conditioned
d_cond <- extract_plot_data(
  results_list[["C_deaths"]]$impact, dates_deaths,
  "Conditioned on log(crossings)"
)

# MoU date
mou_date <- as.Date("2017-02-01")

# Theme
theme_clean <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 10, 5, 10)
  )

# --- Panel 1: C deaths ORIGINAL ---
p1 <- ggplot(d_orig, aes(x = date)) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
              fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = predicted), color = "steelblue",
            linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = observed), color = "black", linewidth = 0.5) +
  geom_vline(xintercept = mou_date, linetype = "dashed",
             color = "red", linewidth = 0.5) +
  labs(
    title = "C. Deaths — Original (no crossings covariate)",
    subtitle = "p = 0.197 | Effect: +8.3% | Not significant",
    y = "log(deaths + 1)", x = NULL
  ) +
  theme_clean

# --- Panel 2: C deaths CONDITIONED ---
p2 <- ggplot(d_cond, aes(x = date)) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
              fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = predicted), color = "steelblue",
            linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = observed), color = "black", linewidth = 0.5) +
  geom_vline(xintercept = mou_date, linetype = "dashed",
             color = "red", linewidth = 0.5) +
  labs(
    title = "C. Deaths — Conditioned on log(crossings)",
    subtitle = "p = 0.021 | Effect: +16.4% | Significant",
    y = "log(deaths + 1)", x = NULL
  ) +
  theme_clean

# --- Combine ---
fig <- arrangeGrob(
  p1, p2,
  ncol = 2,
  top = textGrob(
    "Sensitivity: Conditioning on Crossing Attempts (Model C, Death Counts)",
    gp = gpar(fontsize = 13, fontface = "bold")
  ),
  bottom = textGrob(
    "Red dashed line = MoU (Feb 2017). Black = observed, blue dashed = counterfactual, shaded = 95% CI.",
    gp = gpar(fontsize = 8, col = "grey40")
  )
)

# Save
out_dir <- file.path(PROJECT_DIR, "output", "figures", "png")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "conditioning_on_crossings.png")

ggsave(out_path, fig, width = 12, height = 5, dpi = 200, bg = "white")
message("Saved: ", out_path)

# Also save PDF
out_pdf <- file.path(PROJECT_DIR, "output", "figures", "pdf",
                     "conditioning_on_crossings.pdf")
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(out_pdf, fig, width = 12, height = 5, bg = "white")
message("Saved: ", out_pdf)
