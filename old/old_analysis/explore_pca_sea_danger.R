# explore_pca_sea_danger.R
# =======================
# Exploratory PCA on weather variables to construct a sea danger index.
# Uses the existing df_extended.RDS (base build) + daily wave extremes.
# Currents are deferred (large NetCDF, slow extraction).

library(ncdf4)
library(dplyr)
library(lubridate)

BASE_DIR <- file.path(
  "replication", "rodriguez-sanchez", "extension"
)

# ============================================================
# 1. Load base dataset
# ============================================================
df <- readRDS(file.path(BASE_DIR, "data", "df_extended.RDS"))
cat("Base dataset:", nrow(df), "rows x", ncol(df), "cols\n")

# ============================================================
# 2. Add wave extreme statistics from daily SWH
# ============================================================
cat("\n--- Computing wave extreme statistics from daily NetCDF ---\n")
nc <- nc_open(file.path(BASE_DIR, "data",
                        "era5_daily_waves_central_med.nc"))
swh_raw <- ncvar_get(nc, "swh")

time_name <- intersect(c("valid_time", "time"), names(nc$dim))[1]
time_vals <- ncvar_get(nc, time_name)
time_units <- ncatt_get(nc, time_name, "units")$value

if (grepl("seconds since 1970", time_units)) {
  dates_w <- as.Date(as.POSIXct(time_vals, origin = "1970-01-01", tz = "UTC"))
} else if (grepl("hours since 1900", time_units)) {
  dates_w <- as.Date(as.POSIXct("1900-01-01", tz = "UTC") + time_vals * 3600)
} else {
  ref <- sub("seconds since ", "", time_units)
  dates_w <- as.Date(as.POSIXct(ref, tz = "UTC") + time_vals)
}
nc_close(nc)

# Vectorized spatial mean (much faster than loop)
ndim <- length(dim(swh_raw))
n_time <- length(dates_w)
cat("  Daily SWH time steps:", n_time, ", dims:", paste(dim(swh_raw), collapse="x"), "\n")

if (ndim == 3) {
  # Collapse spatial dims: reshape to (space, time) then colMeans
  d <- dim(swh_raw)
  mat <- matrix(swh_raw, nrow = d[1] * d[2], ncol = d[3])
  swh_spatial <- colMeans(mat, na.rm = TRUE)
} else if (ndim == 2) {
  swh_spatial <- colMeans(swh_raw, na.rm = TRUE)
} else {
  stop("Unexpected SWH dimensions: ", ndim)
}

df_daily <- data.frame(
  month_date = floor_date(dates_w, "month"),
  swh_sm = swh_spatial
)

df_wave_ext <- df_daily %>%
  group_by(month_date) %>%
  summarise(
    wave_max_central_med = max(swh_sm, na.rm = TRUE),
    wave_sd_central_med  = sd(swh_sm, na.rm = TRUE),
    wave_days_above_2m   = sum(swh_sm > 2.0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(date = month_date)

df <- left_join(df, df_wave_ext, by = "date")
cat("  Wave extremes added. NAs:", sum(is.na(df$wave_max_central_med)), "\n")

# ============================================================
# 3. Prepare analysis sample
# ============================================================
PRE_START <- as.Date("2011-02-01")
END_DATE  <- as.Date("2021-09-01")
MOU_DATE  <- as.Date("2017-07-01")

df_anal <- df %>%
  filter(date >= PRE_START & date <= END_DATE) %>%
  mutate(
    neg_wave_period = -wave_period_central_med,
    post_mou = as.integer(date >= MOU_DATE)
  )

cat("\nAnalysis sample:", nrow(df_anal), "months\n")

# ============================================================
# 4. Summary statistics
# ============================================================
cat("\n============================================================\n")
cat("SUMMARY STATISTICS FOR PCA CANDIDATES\n")
cat("============================================================\n\n")

pca_candidates <- c(
  "wave_height_central_med",   # SWH
  "wind_speed_central_med",    # Wind
  "neg_wave_period",           # -MWP (higher = shorter period = more dangerous)
  "wave_max_central_med",      # Monthly max SWH
  "wave_sd_central_med",       # Monthly SD of SWH
  "wave_days_above_2m"         # Count of days with SWH > 2m
)

short_labels <- c("SWH", "Wind", "-MWP", "WMax", "WSD", "D>2m")

cat(sprintf("%-32s %8s %8s %8s %8s %4s\n",
    "Variable", "Mean", "SD", "Min", "Max", "NAs"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (i in seq_along(pca_candidates)) {
  v <- pca_candidates[i]
  x <- df_anal[[v]]
  cat(sprintf("%-32s %8.4f %8.4f %8.4f %8.4f %4d\n",
      v, mean(x, na.rm = TRUE), sd(x, na.rm = TRUE),
      min(x, na.rm = TRUE), max(x, na.rm = TRUE), sum(is.na(x))))
}

# ============================================================
# 5. Correlation matrix
# ============================================================
cat("\n============================================================\n")
cat("CORRELATION MATRIX\n")
cat("============================================================\n\n")

cor_mat <- cor(df_anal[, pca_candidates], use = "complete.obs")
colnames(cor_mat) <- rownames(cor_mat) <- short_labels
print(round(cor_mat, 3))

# ============================================================
# 6. PCA: 3 core variables (SWH, Wind, -MWP)
# ============================================================
cat("\n============================================================\n")
cat("PCA A: 3 CORE VARIABLES (SWH, Wind, -MWP)\n")
cat("============================================================\n")

core3 <- c("wave_height_central_med", "wind_speed_central_med",
           "neg_wave_period")

df_pca3 <- df_anal[complete.cases(df_anal[, core3]), ]
pca3 <- prcomp(df_pca3[, core3], center = TRUE, scale. = TRUE)

cat("\nLoadings:\n")
ld3 <- pca3$rotation
rownames(ld3) <- c("SWH", "Wind", "-MWP")
print(round(ld3, 3))
cat("\nVariance explained:\n")
print(round(summary(pca3)$importance, 3))

# ============================================================
# 7. PCA: 6 variables (add wave extremes)
# ============================================================
cat("\n============================================================\n")
cat("PCA B: 6 VARIABLES (SWH, Wind, -MWP, WMax, WSD, D>2m)\n")
cat("============================================================\n")

df_pca6 <- df_anal[complete.cases(df_anal[, pca_candidates]), ]
pca6 <- prcomp(df_pca6[, pca_candidates], center = TRUE, scale. = TRUE)

cat("\nLoadings:\n")
ld6 <- pca6$rotation
rownames(ld6) <- short_labels
print(round(ld6, 3))
cat("\nVariance explained:\n")
print(round(summary(pca6)$importance, 3))

# ============================================================
# 8. PCA: reduced set excluding redundant wave derivatives
# ============================================================
# wave_max and wave_days_above_2m are derived from the same daily SWH
# as wave_height_central_med. They may be mechanically redundant.
# Try: SWH, Wind, -MWP, WSD (variability is less redundant than max/count)

cat("\n============================================================\n")
cat("PCA C: 4 VARIABLES (SWH, Wind, -MWP, WSD)\n")
cat("============================================================\n")

core4b <- c("wave_height_central_med", "wind_speed_central_med",
            "neg_wave_period", "wave_sd_central_med")

df_pca4b <- df_anal[complete.cases(df_anal[, core4b]), ]
pca4b <- prcomp(df_pca4b[, core4b], center = TRUE, scale. = TRUE)

cat("\nLoadings:\n")
ld4b <- pca4b$rotation
rownames(ld4b) <- c("SWH", "Wind", "-MWP", "WSD")
print(round(ld4b, 3))
cat("\nVariance explained:\n")
print(round(summary(pca4b)$importance, 3))

# ============================================================
# 9. Validation: Beaufort + WMO + volume prediction
# ============================================================
cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n")

# Use PCA A (3-var) for validation
pc1_3 <- predict(pca3, newdata = df_pca3[, core3])[, 1]
if (pca3$rotation["wave_height_central_med", "PC1"] < 0) pc1_3 <- -pc1_3
df_pca3$pc1_3var <- pc1_3

# Use PCA B (6-var) for comparison
pc1_6 <- predict(pca6, newdata = df_pca6[, pca_candidates])[, 1]
if (pca6$rotation["wave_height_central_med", "PC1"] < 0) pc1_6 <- -pc1_6
df_pca6$pc1_6var <- pc1_6

# Beaufort
bf_breaks <- c(0, 0.3, 1.6, 3.4, 5.5, 8.0, 10.8, 13.9, 17.2, Inf)
bf_labels <- 0:8  # monthly means won't exceed ~8

df_pca3$beaufort <- as.integer(as.character(
  cut(df_pca3$wind_speed_central_med, breaks = bf_breaks,
      labels = bf_labels, right = FALSE, include.lowest = TRUE)))

# WMO
wmo_breaks <- c(0, 0.001, 0.1, 0.5, 1.25, 2.5, 4.0, 6.0, Inf)
wmo_labels <- 0:7

df_pca3$wmo <- as.integer(as.character(
  cut(df_pca3$wave_height_central_med, breaks = wmo_breaks,
      labels = wmo_labels, right = FALSE, include.lowest = TRUE)))

cat("\nBeaufort distribution (monthly means):\n")
print(table(df_pca3$beaufort))

cat("\nWMO sea state distribution (monthly means):\n")
print(table(df_pca3$wmo))

# Danger flag
df_pca3$danger <- as.integer(df_pca3$beaufort >= 6 | df_pca3$wmo >= 5)
cat("\nDanger months (Beaufort>=6 OR WMO>=5):",
    sum(df_pca3$danger, na.rm = TRUE), "out of", nrow(df_pca3), "\n")

# Correlations
cat("\n--- Rank correlations with PC1 (3-var) ---\n")
cat("  Spearman PC1 vs Beaufort:",
    round(cor(df_pca3$pc1_3var, df_pca3$beaufort,
              method = "spearman", use = "complete"), 3), "\n")
cat("  Spearman PC1 vs WMO:",
    round(cor(df_pca3$pc1_3var, df_pca3$wmo,
              method = "spearman", use = "complete"), 3), "\n")

# Cross-classification
top_q <- quantile(df_pca3$pc1_3var, 0.75, na.rm = TRUE)
if (sum(df_pca3$danger, na.rm = TRUE) > 0) {
  frac <- mean(df_pca3$pc1_3var[df_pca3$danger == 1] >= top_q, na.rm = TRUE)
  cat("  Danger months in top PC1 quartile:", round(frac, 3), "\n")
}

# Volume prediction
cat("\n--- Volume prediction ---\n")
df_pca3$log_cross <- log(df_pca3$crossings_CMR + 1)
df_pca3$month_fac <- factor(month(df_pca3$date))

m_vol_pc1 <- lm(log_cross ~ pc1_3var + month_fac, data = df_pca3)
cat("  PC1 coef:", round(coef(m_vol_pc1)["pc1_3var"], 4),
    ", p =", round(summary(m_vol_pc1)$coefficients["pc1_3var", 4], 4), "\n")

m_vol_swh <- lm(log_cross ~ wave_height_central_med + month_fac, data = df_pca3)
cat("  SWH coef:", round(coef(m_vol_swh)["wave_height_central_med"], 4),
    ", p =", round(summary(m_vol_swh)$coefficients["wave_height_central_med", 4], 4), "\n")

# ============================================================
# 10. Compare PC1 vs SWH alone in interaction regression
# ============================================================
cat("\n============================================================\n")
cat("INTERACTION REGRESSION: PC1 vs SWH alone\n")
cat("============================================================\n")

df_pca3$log_rate <- log(df_pca3$mortality_rate_100 + 0.01)

# PC1 interaction
m_pc1 <- lm(log_rate ~ pc1_3var * post_mou + month_fac, data = df_pca3)
s_pc1 <- summary(m_pc1)$coefficients

cat("\n--- PC1 x PostMoU ---\n")
cat("  PC1 main:     beta =", round(s_pc1["pc1_3var", 1], 4),
    ", p =", round(s_pc1["pc1_3var", 4], 4), "\n")
cat("  PC1 interact: beta =", round(s_pc1["pc1_3var:post_mou", 1], 4),
    ", p =", round(s_pc1["pc1_3var:post_mou", 4], 4), "\n")
cat("  R-squared:", round(summary(m_pc1)$r.squared, 4), "\n")

# SWH interaction
m_swh <- lm(log_rate ~ wave_height_central_med * post_mou + month_fac,
            data = df_pca3)
s_swh <- summary(m_swh)$coefficients

cat("\n--- SWH x PostMoU ---\n")
cat("  SWH main:     beta =", round(s_swh["wave_height_central_med", 1], 4),
    ", p =", round(s_swh["wave_height_central_med", 4], 4), "\n")
cat("  SWH interact: beta =", round(s_swh["wave_height_central_med:post_mou", 1], 4),
    ", p =", round(s_swh["wave_height_central_med:post_mou", 4], 4), "\n")
cat("  R-squared:", round(summary(m_swh)$r.squared, 4), "\n")

# Wind interaction
m_wind <- lm(log_rate ~ wind_speed_central_med * post_mou + month_fac,
             data = df_pca3)
s_wind <- summary(m_wind)$coefficients

cat("\n--- Wind x PostMoU ---\n")
cat("  Wind main:     beta =", round(s_wind["wind_speed_central_med", 1], 4),
    ", p =", round(s_wind["wind_speed_central_med", 4], 4), "\n")
cat("  Wind interact: beta =", round(s_wind["wind_speed_central_med:post_mou", 1], 4),
    ", p =", round(s_wind["wind_speed_central_med:post_mou", 4], 4), "\n")
cat("  R-squared:", round(summary(m_wind)$r.squared, 4), "\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
