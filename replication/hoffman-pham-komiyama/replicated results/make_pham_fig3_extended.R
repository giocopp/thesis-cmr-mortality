#!/usr/bin/env Rscript

# Reproduce Hoffmann Pham & Komiyama (2024) Fig. 3 style and extend to latest IOM data.

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath("."))
  }
  raw_path <- sub("^--file=", "", file_arg[1])
  if (!file.exists(raw_path)) {
    raw_path <- gsub("~\\+~", " ", raw_path, fixed = FALSE)
  }
  dirname(normalizePath(raw_path))
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

centered_ma6 <- function(x) {
  as.numeric(stats::filter(x, rep(1 / 6, 6), sides = 2))
}

read_iom_crossings_from_excel <- function(excel_path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to read the IOM Excel file.")
  }

  raw <- readxl::read_excel(
    excel_path,
    sheet = "Crossings 2014-24 ALL ROUTE",
    col_names = FALSE,
    .name_repair = "minimal"
  )

  # Same row-boundary fix as in the cleaning pipeline (2024/2025 have no year marker).
  bounds <- list(
    `2014` = c(4, 15),
    `2015` = c(17, 28),
    `2016` = c(30, 41),
    `2017` = c(43, 54),
    `2018` = c(56, 67),
    `2019` = c(69, 80),
    `2020` = c(82, 93),
    `2021` = c(95, 106),
    `2022` = c(108, 119),
    `2023` = c(121, 132),
    `2024` = c(134, 145),
    `2025` = c(147, 158)
  )

  month_map <- c(
    january = 1, february = 2, march = 3, april = 4,
    may = 5, june = 6, july = 7, august = 8,
    september = 9, october = 10, november = 11, december = 12
  )

  rows <- list()
  j <- 1L

  for (yr in names(bounds)) {
    start_row <- bounds[[yr]][1]
    end_row <- bounds[[yr]][2]

    for (r in seq(start_row, end_row)) {
      month_raw <- raw[[2]][r]
      month_chr <- tolower(trimws(as.character(month_raw)))

      if (is.na(month_raw) || month_chr == "" ||
          grepl("total|as of|source", month_chr, ignore.case = TRUE)) {
        next
      }

      month_num <- month_map[[month_chr]]
      if (is.null(month_num) || is.na(month_num)) {
        next
      }

      rows[[j]] <- data.frame(
        date = as.Date(sprintf("%s-%02d-01", yr, month_num)),
        arrivals_italy = safe_num(raw[[7]][r]),
        arrivals_malta = safe_num(raw[[8]][r]),
        inter_libya = safe_num(raw[[16]][r]),
        inter_tunisia = safe_num(raw[[17]][r]),
        deaths_cmr = safe_num(raw[[21]][r]),
        stringsAsFactors = FALSE
      )
      j <- j + 1L
    }
  }

  out <- do.call(rbind, rows)
  out[order(out$date), ]
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))

input_path <- file.path(project_root, "Data", "IOM Data", "Raw", "ALL MED DATA 2010-2025_12.08.2025.xlsx")
output_dir <- file.path(script_dir, "figures")
output_path <- file.path(output_dir, "pham_komiyama_fig3_extended.png")

if (!file.exists(input_path)) {
  stop("Input Excel file not found: ", input_path)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_iom_crossings_from_excel(input_path)

# Match the paper-style crossing accounting:
# rescue (arrivals Italy+Malta) + interceptions (Libya+Tunisia) + deaths (CMR).
df$arrivals <- rowSums(cbind(df$arrivals_italy, df$arrivals_malta), na.rm = TRUE)
df$arrivals[is.na(df$arrivals_italy) & is.na(df$arrivals_malta)] <- NA_real_

df$interceptions <- ifelse(is.na(df$inter_libya), 0, df$inter_libya) +
  ifelse(is.na(df$inter_tunisia), 0, df$inter_tunisia)
df$interceptions[is.na(df$inter_libya) & is.na(df$inter_tunisia)] <- NA_real_

df$crossings <- df$arrivals + df$interceptions + df$deaths_cmr

plot_df <- subset(
  df,
  date >= as.Date("2016-01-01") &
    !is.na(crossings) & crossings > 0 &
    !is.na(arrivals) & !is.na(interceptions) & !is.na(deaths_cmr)
)

if (nrow(plot_df) == 0) {
  stop("No valid observations after filtering.")
}

plot_df$p_rescue <- plot_df$arrivals / plot_df$crossings
plot_df$p_interception <- plot_df$interceptions / plot_df$crossings
plot_df$p_death <- plot_df$deaths_cmr / plot_df$crossings

plot_df$ma_crossings <- centered_ma6(plot_df$crossings)
plot_df$ma_rescue <- centered_ma6(plot_df$p_rescue)
plot_df$ma_interception <- centered_ma6(plot_df$p_interception)
plot_df$ma_death <- centered_ma6(plot_df$p_death)

x_min <- min(plot_df$date, na.rm = TRUE)
x_max <- max(plot_df$date, na.rm = TRUE)

monthly_ticks <- seq(from = as.Date(format(x_min, "%Y-%m-01")), to = x_max, by = "month")
show_half_year_labels <- length(monthly_ticks) <= 72
if (show_half_year_labels) {
  label_ticks <- monthly_ticks[format(monthly_ticks, "%m") %in% c("01", "07")]
  label_text <- ifelse(format(label_ticks, "%m") == "01",
                       paste0("Jan\n", format(label_ticks, "%Y")),
                       "Jul")
} else {
  label_ticks <- monthly_ticks[format(monthly_ticks, "%m") == "01"]
  label_text <- paste0("Jan\n", format(label_ticks, "%Y"))
}

col_bg <- "#D9D9D9"
col_purple <- "#7F1484"
col_rescue <- "#0B7D1A"
col_interception <- "#1A27D8"
col_death <- "#E10C0C"

png(output_path, width = 1800, height = 2200, res = 220)
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = TRUE)

par(mfrow = c(2, 1), mar = c(5.6, 5.5, 3.4, 1.2), mgp = c(2.8, 0.9, 0), cex = 1.0)

# Panel 1: Number of crossings
ylim1 <- c(0, max(plot_df$crossings, na.rm = TRUE) * 1.05)
plot(plot_df$date, plot_df$crossings,
     type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = ylim1)
usr <- par("usr")
rect(usr[1], usr[3], usr[2], usr[4], col = col_bg, border = NA)
axis.Date(1, at = label_ticks, labels = label_text, tck = -0.03, cex.axis = 1.0)
axis(1, at = monthly_ticks, labels = FALSE, tck = -0.015)
axis(2, las = 1, tck = -0.02, cex.axis = 1.0)
box(lwd = 1.6)
points(plot_df$date, plot_df$crossings, pch = 16, cex = 0.72, col = col_purple)
lines(plot_df$date, plot_df$ma_crossings, lwd = 2.2, col = col_purple)
title(main = "Number of crossings", cex.main = 1.35, line = 0.75)
mtext("N crossing", side = 2, line = 1.95, cex = 1.25)

# Panel 2: Probabilities
plot(plot_df$date, plot_df$p_rescue,
     type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(-0.05, 1.03))
usr <- par("usr")
rect(usr[1], usr[3], usr[2], usr[4], col = col_bg, border = NA)
axis.Date(1, at = label_ticks, labels = label_text, tck = -0.03, cex.axis = 1.0)
axis(1, at = monthly_ticks, labels = FALSE, tck = -0.015)
axis(2, at = seq(0, 1, by = 0.2), las = 1, tck = -0.02, cex.axis = 1.0)
box(lwd = 1.6)

points(plot_df$date, plot_df$p_rescue, pch = 16, cex = 0.72,
       col = grDevices::adjustcolor(col_rescue, alpha.f = 0.22))
points(plot_df$date, plot_df$p_interception, pch = 16, cex = 0.72,
       col = grDevices::adjustcolor(col_interception, alpha.f = 0.22))
points(plot_df$date, plot_df$p_death, pch = 16, cex = 0.72,
       col = grDevices::adjustcolor(col_death, alpha.f = 0.22))

lines(plot_df$date, plot_df$ma_rescue, lwd = 2.2, col = col_rescue)
lines(plot_df$date, plot_df$ma_interception, lwd = 2.2, col = col_interception)
lines(plot_df$date, plot_df$ma_death, lwd = 2.2, col = col_death)

legend("topleft",
       legend = c("Rescue", "Interception", "Death"),
       col = c(col_rescue, col_interception, col_death),
       lty = 1, lwd = 2.2, bty = "o", cex = 1.05,
       pt.bg = "white", inset = c(0.02, 0.02), bg = "#DDDDDD")

title(main = "Probability of rescue, interception, and death", cex.main = 1.35, line = 0.75)
mtext("Probability", side = 2, line = 1.95, cex = 1.25)

message("Saved: ", output_path)
message("Coverage: ", format(min(plot_df$date), "%Y-%m"), " to ", format(max(plot_df$date), "%Y-%m"))
