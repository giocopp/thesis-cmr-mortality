# Poster figure: a cleaner 2-panel version of the crossings figure, using only
# panel (a) "persons attempting crossing, by outcome" and panel (d) "share of
# persons intercepted, by operation type". Reuses the exact panel objects from
# the thesis script analysis/R/04_descriptive/01_panel.R (no cropping).

suppressPackageStartupMessages({
  library(ggplot2); library(cowplot); library(grid)
})

BASE_DIR <- here::here()
# Sourcing builds p_cross_united_display (panel a) and p_share_lines_display (panel d)
source(file.path(BASE_DIR, "analysis", "R", "04_descriptive", "01_panel.R"))

poster_panel_theme <- theme(
  plot.title = element_text(face = "bold", size = 21, lineheight = 0.98,
                            hjust = 0, margin = margin(b = 1)),
  plot.subtitle = element_text(size = 12.5, colour = "red3",
                               hjust = 0, margin = margin(b = 3)),
  axis.title = element_text(size = 15),
  axis.text = element_text(size = 13),
  legend.text = element_text(size = 14),
  legend.key.size = unit(0.54, "cm"),
  legend.box.margin = margin(3, 5, 3, 5),
  legend.spacing.x = unit(0.10, "cm"),
  legend.spacing.y = unit(0.03, "cm"),
  plot.margin = margin(t = 9, r = 8, b = 3, l = 8)
)

poster_cross_colours <- c(
  "Intercepted persons" = "#1F78B4",
  "Recorded Dead and Missing Migrants" = "#D32F2F",
  "Libyan CG operations" = "#333333",
  "Tunisian CG operations" = "#AAAAAA"
)
poster_cross_labels <- c(
  "Intercepted persons" = "Intercepted",
  "Recorded Dead and Missing Migrants" = "Dead/missing",
  "Libyan CG operations" = "Libyan CG",
  "Tunisian CG operations" = "Tunisian CG"
)
poster_share_colours <- c(
  "SAR operations recorded by Frontex" = "#1F78B4",
  "Non SAR operations recorded by Frontex" = "#F16913",
  "Libyan CG operations" = "#252525",
  "Tunisian CG operations" = "#969696"
)
poster_share_labels <- c(
  "SAR operations recorded by Frontex" = "SAR operations recorded by Frontex",
  "Non SAR operations recorded by Frontex" = "Non SAR operations recorded by Frontex",
  "Libyan CG operations" = "Libyan CG",
  "Tunisian CG operations" = "Tunisian CG"
)

pa <- p_cross_united_display +
  labs(
    title = "(a) Crossing attempts, by outcome",
    subtitle = "Red dashed line: 2 February 2017 Italy\u2013Libya MoU"
  ) +
  scale_fill_manual(
    values = poster_cross_colours,
    breaks = names(poster_cross_labels),
    labels = poster_cross_labels,
    name = NULL
  ) +
  poster_panel_theme +
  guides(fill = guide_legend(nrow = 2, byrow = FALSE, title = NULL))
pb <- p_share_lines_display +
  labs(
    title = "(b) Interception share",
    subtitle = "Red dashed line: 2 February 2017 Italy\u2013Libya MoU"
  ) +
  scale_colour_manual(
    values = poster_share_colours,
    breaks = names(poster_share_labels),
    labels = poster_share_labels,
    name = NULL
  ) +
  poster_panel_theme +
  guides(colour = guide_legend(nrow = 2, byrow = FALSE, title = NULL))

ad <- cowplot::plot_grid(pa, pb, ncol = 2, align = "hv", axis = "tblr")
ad_framed <- cowplot::ggdraw(ad) +
  cowplot::draw_grob(grid::rectGrob(gp = grid::gpar(col = "black", fill = NA, lwd = 2)))

out <- file.path(BASE_DIR, "paper", "poster", "fig-crossings-ad.png")
ggsave(out, ad_framed, width = 15.8, height = 7.55, dpi = 240)
cat("saved:", out, "\n")
