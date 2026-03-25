## 03a — Weather common-support diagnostic: pre vs post MoU
##
## Purpose: Compare the distributions of weather conditions at incident
## locations before and after the MoU (2017-02-02).  This is a prerequisite
## for Model A: if post-MoU incidents only happen in calm weather while
## pre-MoU incidents span the full range, the interaction β₃ would be
## identified off extrapolation rather than common support.
##
## Outputs:
##   output/figures/weather_overlap_densities.pdf
##   output/tables/weather_overlap_summary.csv   (descriptive stats + KS test)

library(data.table)
library(ggplot2)
library(patchwork)

# ---- 1. Load data -----------------------------------------------------------
d <- as.data.table(readRDS("data/processed/cmr_events_with_weather.RDS"))

cat("Total incidents:", nrow(d), "\n")
cat("Pre-MoU:", sum(d$post_mou == 0), " Post-MoU:", sum(d$post_mou == 1), "\n")

d[, period := fifelse(post_mou == 1, "Post-MoU (2017–2025)", "Pre-MoU (2014–2017)")]

# ---- 2. Define variables to compare -----------------------------------------
var_info <- data.table(
  var   = c("swh_day0", "wind_day0", "i10fg_day0", "mwp_lag0", "sst_day0"),
  label = c("Significant wave height (m)",
            "Wind speed (m/s)",
            "10-m wind gust (m/s)",
            "Mean wave period (s)",
            "Sea surface temp (°C)")
)

# ---- 3. Summary statistics table + KS tests ---------------------------------
stats_list <- lapply(seq_len(nrow(var_info)), function(i) {
  v <- var_info$var[i]
  lbl <- var_info$label[i]

  pre  <- d[post_mou == 0 & !is.na(get(v)), get(v)]
  post <- d[post_mou == 1 & !is.na(get(v)), get(v)]

  ks <- ks.test(pre, post)

  data.table(
    variable    = lbl,
    n_pre       = length(pre),
    mean_pre    = round(mean(pre), 3),
    sd_pre      = round(sd(pre), 3),
    median_pre  = round(median(pre), 3),
    q10_pre     = round(quantile(pre, 0.10), 3),
    q90_pre     = round(quantile(pre, 0.90), 3),
    n_post      = length(post),
    mean_post   = round(mean(post), 3),
    sd_post     = round(sd(post), 3),
    median_post = round(median(post), 3),
    q10_post    = round(quantile(post, 0.10), 3),
    q90_post    = round(quantile(post, 0.90), 3),
    ks_stat     = round(ks$statistic, 3),
    ks_pvalue   = signif(ks$p.value, 3)
  )
})

stats_dt <- rbindlist(stats_list)
fwrite(stats_dt, "output/tables/weather_overlap_summary.csv")
cat("\nSummary statistics written to output/tables/weather_overlap_summary.csv\n\n")
print(stats_dt[, .(variable, mean_pre, mean_post, median_pre, median_post, ks_stat, ks_pvalue)])

# ---- 4. Density plots -------------------------------------------------------
plot_list <- lapply(seq_len(nrow(var_info)), function(i) {
  v   <- var_info$var[i]
  lbl <- var_info$label[i]

  dd <- d[!is.na(get(v)), .(value = get(v), period)]

  ggplot(dd, aes(x = value, fill = period, colour = period)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    labs(x = lbl, y = "Density") +
    scale_fill_manual(values = c("Pre-MoU (2014–2017)" = "#2166AC",
                                 "Post-MoU (2017–2025)" = "#B2182B")) +
    scale_colour_manual(values = c("Pre-MoU (2014–2017)" = "#2166AC",
                                   "Post-MoU (2017–2025)" = "#B2182B")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank())
})

# Combine with patchwork
combined <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
            (plot_list[[4]] | plot_list[[5]] | plot_spacer()) +
  plot_annotation(
    title    = "Weather conditions at incident locations: Pre vs Post MoU",
    subtitle = "Kernel density estimates, day of incident",
    caption  = "Blue = Pre-MoU (2014–2017), Red = Post-MoU (2017–2025).\nSource: IOM MMP incidents matched to ERA5."
  ) &
  theme(plot.title = element_text(size = 13, face = "bold"))

ggsave("output/figures/weather_overlap_densities.pdf",
       combined, width = 12, height = 7.5)
cat("Density plots saved to output/figures/weather_overlap_densities.pdf\n")

# ---- 5. Print key takeaway --------------------------------------------------
cat("\n--- Overlap assessment ---\n")
for (i in seq_len(nrow(var_info))) {
  v <- var_info$var[i]
  pre_range  <- range(d[post_mou == 0 & !is.na(get(v)), get(v)])
  post_range <- range(d[post_mou == 1 & !is.na(get(v)), get(v)])
  overlap_lo <- max(pre_range[1], post_range[1])
  overlap_hi <- min(pre_range[2], post_range[2])
  cat(sprintf("%s: pre [%.2f, %.2f], post [%.2f, %.2f], overlap [%.2f, %.2f]\n",
              var_info$label[i],
              pre_range[1], pre_range[2],
              post_range[1], post_range[2],
              overlap_lo, overlap_hi))
}
