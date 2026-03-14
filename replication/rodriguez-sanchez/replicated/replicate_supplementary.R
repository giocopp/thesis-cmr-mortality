################################################################################
# SUPPLEMENTARY MATERIALS FIGURE REPRODUCTION SCRIPT
# "Search-and-rescue in the Central Mediterranean Route does not induce migration"
# Rodriguez Sanchez et al. (2023) - Nature Scientific Reports
#
# This script reproduces ALL supplementary figures (S1-S10):
#   - S1: NGO operations timeline (daily resolution)
#   - S2: Time series decomposition
#   - S3: Death ratio over time
#   - S4: Structural breakpoints (deaths)
#   - S5: Structural breakpoints (crossings)
#   - S6: Cumulative prediction errors (8 model specifications)
#   - S7: Effect size comparison (8 model specifications)
#   - S8: ACF of residuals (using bsts::AcfDist)
#   - S9: Q-Q plots (using bsts::qqdist)
#   - S10: Coefficient inclusion probabilities
#
# NOTE: S11 (sensitivity analysis) requires re-running models and is not included
################################################################################

cat("=======================================================================\n")
cat("SUPPLEMENTARY FIGURES - Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=======================================================================\n\n")

# ==============================================================================
# 0. SETUP AND PACKAGE LOADING
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Loading packages...\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(gridExtra)
  library(CausalImpact)
  library(bsts)
  library(strucchange)
  library(forecast)
})

set.seed(270488)

# Create output directories
dir.create("replicated_results/figures/supplementary/png", recursive = TRUE, showWarnings = FALSE)
dir.create("replicated_results/figures/supplementary/pdf", recursive = TRUE, showWarnings = FALSE)
dir.create("replicated_results/validation", recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time(), "%H:%M:%S"), "] Packages loaded.\n\n")

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Loading and preparing data...\n")

df <- readRDS("../Original code and data/df.RDS")
df_y_rec_path <- "../Original code and data/df_y_rec.RDS"
if (file.exists(df_y_rec_path)) {
  df_y_rec <- readRDS(df_y_rec_path)
  df <- left_join(df, df_y_rec, by = "date")
}

df <- df %>%
  mutate(
    LCG_pushbacks_count = as.numeric(ifelse(is.na(LCG_pushbacks_count), 0, LCG_pushbacks_count)),
    TCG_pushbacks_count = as.numeric(ifelse(is.na(TCG_pushbacks_count), 0, TCG_pushbacks_count)),
    dead_and_missing_Central_Mediterranean = as.numeric(
      ifelse(is.na(dead_and_missing_Central_Mediterranean), 0, dead_and_missing_Central_Mediterranean)
    ),
    crossings_CMR = arrivals_CMR + LCG_pushbacks_count + TCG_pushbacks_count +
      dead_and_missing_Central_Mediterranean,
    mortality_rate = (dead_and_missing_Central_Mediterranean / crossings_CMR) * 1000
  )

# Prepare reduced dataset for model validation (S6-S9)
df_reduced <- df %>%
  filter(date >= "2011-02-01" & date < "2021-10-01") %>%
  dplyr::select(-c(
    contains("lag_24", ignore.case = TRUE),
    contains("lag_23", ignore.case = TRUE),
    contains("lag_22", ignore.case = TRUE),
    contains("lag_21", ignore.case = TRUE),
    contains("lag_20", ignore.case = TRUE),
    contains("lag_19", ignore.case = TRUE),
    contains("lag_18", ignore.case = TRUE),
    contains("lag_17", ignore.case = TRUE),
    contains("lag_16", ignore.case = TRUE),
    contains("lag_15", ignore.case = TRUE),
    contains("lag_14", ignore.case = TRUE),
    contains("lag_13", ignore.case = TRUE),
    contains("lag_12", ignore.case = TRUE),
    contains("lag_11", ignore.case = TRUE),
    contains("lag_10", ignore.case = TRUE),
    contains("lag_09", ignore.case = TRUE),
    contains("lag_08", ignore.case = TRUE),
    contains("lag_07", ignore.case = TRUE),
    starts_with("airflow_Palestinian.Territories"),
    starts_with("asylum"),
    "arrivals_BSR", "arrivals_CMR", "arrivals_CRAG", "arrivals_EBR",
    "arrivals_EMR", "arrivals_OR", "arrivals_WAR", "arrivals_WBR", "arrivals_WMR",
    "dead_and_missing_Eastern_Mediterranean", "dead_and_missing_Central_Mediterranean",
    "dead_and_missing_Western_Mediterranean",
    "sd_lat__Eastern_Mediterranean", "sd_lat__Central_Mediterranean",
    "sd_lat__Western_Mediterranean", "sd_lon__Eastern_Mediterranean",
    "sd_lon__Central_Mediterranean", "sd_lon__Western_Mediterranean",
    "frac_index_2_to_10_deads_Eastern_Mediterranean",
    "frac_index_2_to_10_deads_Central_Mediterranean",
    "frac_index_2_to_10_deads_Western_Mediterranean",
    "frac_index_less_than_1_dead_Eastern_Mediterranean",
    "frac_index_less_than_1_dead_Central_Mediterranean",
    "frac_index_less_than_1_dead_Western_Mediterranean",
    "frac_index_more_than_10_deads_Eastern_Mediterranean",
    "frac_index_more_than_10_deads_Central_Mediterranean",
    "frac_index_more_than_10_deads_Western_Mediterranean",
    "LCG_pushbacks_count", "TCG_pushbacks_count",
    "mortality_rate",
    any_of("y_rec")
  ))

df_reduced_A <- data.frame(
  date = df_reduced$date,
  crossings_CMR = log(df_reduced$crossings_CMR),
  dplyr::select(df_reduced, -c(date, crossings_CMR))
)

df_reduced_A <- df_reduced_A %>%
  mutate(
    month = month(date),
    semester = semester(date),
    quarter = quarter(date)
  )

df_min_A <- df_reduced_A %>% na.omit() %>% arrange(date)

cat("[", format(Sys.time(), "%H:%M:%S"), "] Data prepared. Rows:", nrow(df_min_A), "\n\n")

# ==============================================================================
# FIGURE S1: Number of NGO-led SAR operations (DAILY resolution)
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure S1...\n")

# Create daily date sequence
study_period <- seq(mdy("JAN 01, 2009"), mdy("Oct 31, 2021"), by = "days")
df_dates <- data.frame(date = study_period)

# Helper function for interval checking
in_interval <- function(dates, start, end) {
  dates >= mdy(start) & dates <= mdy(end)
}

# Initialize all SAR columns to 0
ngo_names <- c("SAR_SeaWatch123", "SAR_SeaWatch4", "SAR_Lifeline",
               "SAR_SeaEye_TheSeaEye", "SAR_SeaEye_TheSeefuchs", "SAR_SeaEye_AlanKurdi",
               "SAR_OpenArms_Astral", "SAR_OpenArms_GolfoAzzuro", "SAR_OpenArms_OpenArms",
               "SAR_MareLiberum", "SAR_Mediterranea", "SAR_SMH", "SAR_LouiseMichel",
               "SAR_RefugeeRescue", "SAR_MOAS", "SAR_JugendRettet", "SAR_MSFandSOS",
               "SAR_SavetheChildren", "SAR_MSF_BourbonArgos", "SAR_MSF_Dignity1",
               "SAR_MSF_VosPrudence", "SAR_Resqship", "SAR_Lifeboat")

for (ngo in ngo_names) {
  df_dates[[ngo]] <- 0
}

d <- df_dates$date

# SeaWatch 1/2/3 - multiple active periods
df_dates$SAR_SeaWatch123[in_interval(d, "Jun 20, 2015", "Jul 01, 2018")] <- 1
df_dates$SAR_SeaWatch123[in_interval(d, "Oct 21, 2018", "Jan 30, 2019")] <- 1
df_dates$SAR_SeaWatch123[in_interval(d, "Feb 23, 2019", "May 17, 2019")] <- 1
df_dates$SAR_SeaWatch123[in_interval(d, "Jun 02, 2019", "Jun 28, 2019")] <- 1
df_dates$SAR_SeaWatch123[in_interval(d, "Dec 31, 2019", "Feb 27, 2020")] <- 1
df_dates$SAR_SeaWatch123[in_interval(d, "Jun 07, 2020", "Jun 16, 2020")] <- 1
df_dates$SAR_SeaWatch123[in_interval(d, "Feb 25, 2021", "Mar 26, 2021")] <- 1

# SeaWatch 4
df_dates$SAR_SeaWatch4[in_interval(d, "Aug 15, 2020", "Sep 19, 2020")] <- 1
df_dates$SAR_SeaWatch4[in_interval(d, "Mar 03, 2021", "Oct 31, 2021")] <- 1

# Mission Lifeline
df_dates$SAR_Lifeline[in_interval(d, "Sep 14, 2017", "Jun 26, 2018")] <- 1
df_dates$SAR_Lifeline[in_interval(d, "Aug 26, 2019", "Sep 02, 2019")] <- 1

# Sea-Eye vessels
df_dates$SAR_SeaEye_TheSeaEye[in_interval(d, "Apr 19, 2016", "Jun 21, 2018")] <- 1
df_dates$SAR_SeaEye_TheSeefuchs[in_interval(d, "May 18, 2017", "Jun 06, 2018")] <- 1
df_dates$SAR_SeaEye_AlanKurdi[in_interval(d, "Dec 21, 2018", "May 04, 2020")] <- 1
df_dates$SAR_SeaEye_AlanKurdi[in_interval(d, "Sep 12, 2020", "Sep 25, 2020")] <- 1

# Open Arms vessels
df_dates$SAR_OpenArms_Astral[in_interval(d, "Jun 01, 2016", "Oct 31, 2021")] <- 1
df_dates$SAR_OpenArms_GolfoAzzuro[in_interval(d, "Dec 01, 2016", "Dec 01, 2017")] <- 1
df_dates$SAR_OpenArms_OpenArms[in_interval(d, "Dec 01, 2017", "Mar 16, 2018")] <- 1
df_dates$SAR_OpenArms_OpenArms[in_interval(d, "Apr 17, 2018", "Jan 13, 2019")] <- 1
df_dates$SAR_OpenArms_OpenArms[in_interval(d, "Apr 24, 2019", "Aug 20, 2019")] <- 1
df_dates$SAR_OpenArms_OpenArms[in_interval(d, "Aug 30, 2019", "Apr 17, 2021")] <- 1

# Mare Liberum
df_dates$SAR_MareLiberum[in_interval(d, "Aug 26, 2018", "Oct 31, 2021")] <- 1

# Mediterranea - multiple active periods
df_dates$SAR_Mediterranea[in_interval(d, "Oct 03, 2018", "Mar 18, 2019")] <- 1
df_dates$SAR_Mediterranea[in_interval(d, "Mar 28, 2019", "May 09, 2019")] <- 1
df_dates$SAR_Mediterranea[in_interval(d, "Jul 02, 2019", "Jul 06, 2019")] <- 1
df_dates$SAR_Mediterranea[in_interval(d, "Aug 23, 2019", "Sep 01, 2019")] <- 1
df_dates$SAR_Mediterranea[in_interval(d, "Feb 05, 2020", "Mar 18, 2020")] <- 1
df_dates$SAR_Mediterranea[in_interval(d, "Jun 10, 2020", "Sep 25, 2020")] <- 1

# SMH
df_dates$SAR_SMH[in_interval(d, "Oct 01, 2018", "Jan 17, 2019")] <- 1
df_dates$SAR_SMH[in_interval(d, "Apr 18, 2019", "May 07, 2020")] <- 1
df_dates$SAR_SMH[in_interval(d, "Dec 09, 2020", "Oct 31, 2021")] <- 1

# Louise Michel
df_dates$SAR_LouiseMichel[in_interval(d, "Aug 22, 2020", "Oct 22, 2020")] <- 1

# Refugee Rescue
df_dates$SAR_RefugeeRescue[in_interval(d, "Jan 15, 2016", "Aug 14, 2020")] <- 1

# MOAS
df_dates$SAR_MOAS[in_interval(d, "Aug 26, 2014", "Sep 06, 2017")] <- 1

# Jugend Rettet
df_dates$SAR_JugendRettet[in_interval(d, "Jul 24, 2016", "Aug 01, 2017")] <- 1

# MSF and SOS Mediterranee
df_dates$SAR_MSFandSOS[in_interval(d, "May 09, 2015", "Nov 19, 2018")] <- 1
df_dates$SAR_MSFandSOS[in_interval(d, "Jul 21, 2019", "Oct 31, 2021")] <- 1

# Save the Children
df_dates$SAR_SavetheChildren[in_interval(d, "Sep 08, 2016", "Oct 23, 2017")] <- 1

# MSF vessels
df_dates$SAR_MSF_BourbonArgos[in_interval(d, "May 09, 2015", "Aug 16, 2015")] <- 1
df_dates$SAR_MSF_BourbonArgos[in_interval(d, "Oct 03, 2015", "Jan 14, 2016")] <- 1
df_dates$SAR_MSF_BourbonArgos[in_interval(d, "May 06, 2016", "Nov 20, 2016")] <- 1
df_dates$SAR_MSF_Dignity1[in_interval(d, "Jun 13, 2015", "Nov 04, 2015")] <- 1
df_dates$SAR_MSF_Dignity1[in_interval(d, "Apr 22, 2016", "Oct 04, 2016")] <- 1
df_dates$SAR_MSF_VosPrudence[in_interval(d, "Mar 20, 2017", "Oct 05, 2017")] <- 1

# Resqship and Lifeboat
df_dates$SAR_Resqship[in_interval(d, "Apr 01, 2019", "Oct 28, 2019")] <- 1
df_dates$SAR_Lifeboat[in_interval(d, "Jul 01, 2016", "Sep 22, 2017")] <- 1

# Sum NGO SAR columns
df_dates$ngo_count <- rowSums(df_dates[, ngo_names], na.rm = TRUE)

fig_s1 <- ggplot(df_dates, aes(x = date, y = ngo_count)) +
  geom_step(linewidth = 0.4) +
  geom_vline(xintercept = as.numeric(ymd("2014-10-01")), col = "darkred", lty = 2) +
  geom_vline(xintercept = as.numeric(ymd("2017-02-01")), col = "darkred", lty = 2) +
  annotate("text", x = ymd("2014-10-01"), y = 12, label = "End of Mare Nostrum",
           angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = ymd("2017-02-01"), y = 12, label = "EU - Libyan Coast Guard",
           angle = 90, vjust = -0.5, size = 3) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(breaks = seq(0, 12, by = 2.5), limits = c(0, 12.5)) +
  labs(
    title = "Number of NGO-led search-and-rescue operations, 2009-2021",
    x = "Date", y = "Number"
  ) +
  theme_minimal()

ggsave("replicated_results/figures/supplementary/png/FigureS1_NGO_operations.png",
       fig_s1, width = 8, height = 5, dpi = 300)
ggsave("replicated_results/figures/supplementary/pdf/FigureS1_NGO_operations.pdf",
       fig_s1, width = 8, height = 5)

cat("[", format(Sys.time(), "%H:%M:%S"), "] Figure S1 saved.\n\n")

# ==============================================================================
# FIGURE S2: Decomposition of log of attempted crossings
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure S2...\n")

y <- ts(log(df$crossings_CMR), start = c(2009, 1), frequency = 12)
decomp <- decompose(y, type = "additive")

png("replicated_results/figures/supplementary/png/FigureS2_decomposition.png",
    width = 8, height = 10, units = "in", res = 300)
plot(decomp)
dev.off()

pdf("replicated_results/figures/supplementary/pdf/FigureS2_decomposition.pdf",
    width = 8, height = 10)
plot(decomp)
dev.off()

cat("[", format(Sys.time(), "%H:%M:%S"), "] Figure S2 saved.\n\n")

# ==============================================================================
# FIGURE S3: Death ratio over time
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure S3...\n")

df_ratio <- df %>%
  filter(!is.na(crossings_CMR) & crossings_CMR > 0) %>%
  mutate(
    arrivals_pushbacks = arrivals_CMR + LCG_pushbacks_count + TCG_pushbacks_count,
    death_ratio = dead_and_missing_Central_Mediterranean / (arrivals_pushbacks + dead_and_missing_Central_Mediterranean)
  ) %>%
  filter(!is.na(death_ratio) & !is.infinite(death_ratio))

fig_s3 <- ggplot(df_ratio, aes(x = date, y = death_ratio)) +
  geom_line(linewidth = 0.5) +
  scale_x_date(date_breaks = "1 year", date_labels = "Jan-%y",
               limits = c(ymd("2009-01-01"), ymd("2022-01-01"))) +
  scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.3, by = 0.1)) +
  labs(
    title = "Estimated number of deaths and missing migrants over the sum of arrivals and pushbacks in the CMR",
    x = "Date", y = "Ratio"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave("replicated_results/figures/supplementary/png/FigureS3_death_ratio.png",
       fig_s3, width = 10, height = 5, dpi = 300)
ggsave("replicated_results/figures/supplementary/pdf/FigureS3_death_ratio.pdf",
       fig_s3, width = 10, height = 5)

cat("[", format(Sys.time(), "%H:%M:%S"), "] Figure S3 saved.\n\n")

# ==============================================================================
# FIGURE S4: Structural breakpoints (deaths)
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure S4...\n")

y_deaths <- ts(df$dead_and_missing_Central_Mediterranean, start = c(2009, 1), frequency = 12)
y_deaths_log <- log(y_deaths + 0.1)

bp_deaths <- breakpoints(y_deaths_log ~ 1, h = 12)
ci_deaths <- confint(bp_deaths)

png("replicated_results/figures/supplementary/png/FigureS4_breakpoints_deaths.png",
    width = 8, height = 6, units = "in", res = 300)
plot(y_deaths_log, main = "Structural breakpoints in montly number of dead and missing in CMR (log)",
     ylab = "Log of number of attempted crossings", xlab = "Date")
lines(bp_deaths, col = "darkred", lty = 2)
lines(ci_deaths)
dev.off()

pdf("replicated_results/figures/supplementary/pdf/FigureS4_breakpoints_deaths.pdf", width = 8, height = 6)
plot(y_deaths_log, main = "Structural breakpoints in montly number of dead and missing in CMR (log)",
     ylab = "Log of number of attempted crossings", xlab = "Date")
lines(bp_deaths, col = "darkred", lty = 2)
lines(ci_deaths)
dev.off()

cat("[", format(Sys.time(), "%H:%M:%S"), "] Figure S4 saved.\n\n")

# ==============================================================================
# FIGURE S5: Structural breakpoints (crossings)
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure S5...\n")

y_crossings <- ts(df$crossings_CMR, start = c(2009, 1), frequency = 12)
y_crossings_log <- log(y_crossings)

bp_crossings <- breakpoints(y_crossings_log ~ 1, h = 12)
ci_crossings <- confint(bp_crossings)

png("replicated_results/figures/supplementary/png/FigureS5_breakpoints_crossings.png",
    width = 8, height = 6, units = "in", res = 300)
plot(y_crossings_log, main = "Structural breakpoints in montly number of attempted crossings (log)",
     ylab = "Log of number of attempted crossings", xlab = "Date")
lines(bp_crossings, col = "darkred", lty = 2)
lines(ci_crossings)
dev.off()

pdf("replicated_results/figures/supplementary/pdf/FigureS5_breakpoints_crossings.pdf", width = 8, height = 6)
plot(y_crossings_log, main = "Structural breakpoints in montly number of attempted crossings (log)",
     ylab = "Log of number of attempted crossings", xlab = "Date")
lines(bp_crossings, col = "darkred", lty = 2)
lines(ci_crossings)
dev.off()

cat("[", format(Sys.time(), "%H:%M:%S"), "] Figure S5 saved.\n\n")

# ==============================================================================
# FIGURES S6-S10: ORIGINAL VALIDATION WORKFLOW WITH ALL 8 BSTS SPECIFICATIONS
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Fitting default validation models (niter = 2500)...\n")
models_exist <- TRUE

if (!models_exist) {
  cat("    Pre-fitted models not found. Run 'reproduce_figures_full_model.R' first.\n")
  cat("    Skipping Figures S6-S10.\n\n")
} else {
  pre.period_mare_nostrum <- ymd(min(df_min_A$date), "2013-09-01")
  post.period_mare_nostrum <- ymd("2013-10-01", max(df_min_A$date))
  pre.period_sar_ngos <- ymd(min(df_min_A$date), "2014-10-01")
  post.period_sar_ngos <- ymd("2014-11-01", max(df_min_A$date))
  pre.period_sarlibya <- ymd(min(df_min_A$date), "2017-01-01")
  post.period_sar_libya <- ymd("2017-02-01", max(df_min_A$date))

  set.seed(270488)

  impact_marenostrum <- CausalImpact(
    df_min_A,
    pre.period = pre.period_mare_nostrum,
    post.period = post.period_mare_nostrum,
    alpha = 0.05,
    model.args = list(nseasons = NULL, dynamic.regression = FALSE, standardize.data = TRUE, max.flips = 100, niter = 2500)
  )

  impact_sarngos <- CausalImpact(
    df_min_A,
    pre.period = pre.period_sar_ngos,
    post.period = post.period_sar_ngos,
    alpha = 0.05,
    model.args = list(nseasons = NULL, dynamic.regression = FALSE, standardize.data = TRUE, max.flips = 100, niter = 2500)
  )

  impact_sarlibya <- CausalImpact(
    df_min_A,
    pre.period = pre.period_sarlibya,
    post.period = post.period_sar_libya,
    alpha = 0.05,
    model.args = list(nseasons = NULL, dynamic.regression = FALSE, standardize.data = TRUE, max.flips = 100, niter = 2500)
  )

  saveRDS(impact_marenostrum, "replicated_results/models/model_marenostrum_validation_default.RDS")
  saveRDS(impact_sarngos, "replicated_results/models/model_sarngos_validation_default.RDS")
  saveRDS(impact_sarlibya, "replicated_results/models/model_sarlibya_validation_default.RDS")

  # ============================================================================
  # DEFINE STATE SPECIFICATIONS FOR 8 MODEL COMPARISON
  # ============================================================================

  create_state_specs <- function(target) {
    sdy <- sd(target, na.rm = TRUE)
    ss <- list()
    sd.prior <- SdPrior(sigma.guess = 0.01 * sdy, upper.limit = sdy)

    state.spec_locallineartrend <- AddLocalLinearTrend(ss, target, sdy = sdy)
    state.spec_semilocallineartrend <- AddSemilocalLinearTrend(ss, target, sdy = sdy)
    state.spec_locallevel <- AddLocalLevel(ss, target, sigma.prior = sd.prior, sdy = sdy)

    state.spec_combined1 <- AddLocalLinearTrend(ss, target, sdy = sdy)
    state.spec_combined1 <- AddAutoAr(state.spec_combined1, target, sdy = sdy, lags = 3)

    state.spec_combined2 <- AddSemilocalLinearTrend(ss, target, sdy = sdy)
    state.spec_combined2 <- AddAutoAr(state.spec_combined2, target, sdy = sdy, lags = 3)

    state.spec_combined3 <- AddLocalLevel(ss, target, sigma.prior = sd.prior, sdy = sdy)
    state.spec_combined3 <- AddAutoAr(state.spec_combined3, target, sdy = sdy, lags = 3)

    state.spec_combined4 <- AddLocalLevel(ss, target, sigma.prior = sd.prior, sdy = sdy)
    state.spec_combined4 <- AddSemilocalLinearTrend(state.spec_combined4, target, sdy = sdy)
    state.spec_combined4 <- AddAutoAr(state.spec_combined4, target, sdy = sdy, lags = 3)

    list(
      LLT = state.spec_locallineartrend,
      SLLT = state.spec_semilocallineartrend,
      LL = state.spec_locallevel,
      `LLT+AR` = state.spec_combined1,
      `SLLT+AR` = state.spec_combined2,
      `LL+AR` = state.spec_combined3,
      `LL+SLLT+AR` = state.spec_combined4
    )
  }

  # ============================================================================
  # RUN ALL 8 MODEL SPECIFICATIONS FOR EACH INTERVENTION
  # ============================================================================

  run_all_models <- function(df_min_A, pre_end_idx, post_start_idx, intervention_name, default_impact, covar_drop_cols = c(1, 2, 3)) {
    cat("\n[", format(Sys.time(), "%H:%M:%S"), "] Processing:", intervention_name, "\n")

    target_y <- df_min_A$crossings_CMR[order(df_min_A$date)]
    post.period <- c(post_start_idx, nrow(df_min_A))
    post.period.response <- target_y[post.period[1]:post.period[2]]

    target_y_na <- target_y
    target_y_na[post.period[1]:post.period[2]] <- NA

    X_vars <- data.frame(lapply(df_min_A[, -covar_drop_cols], function(x) scale(x)))
    na_cols <- which(colSums(is.na(X_vars)) > 0)
    if (length(na_cols) > 0) X_vars <- X_vars[, -na_cols]
    covars_df <- as.matrix(X_vars)

    state_specs <- create_state_specs(target_y_na[1:pre_end_idx])

    impacts <- list()
    bsts_models <- list()
    model_names <- names(state_specs)

    for (i in seq_along(state_specs)) {
      model_name <- model_names[i]
      cat("    Running model:", model_name, "...")

      tryCatch({
        bsts_model <- bsts(
          target_y_na ~ covars_df,
          state.specification = state_specs[[i]],
          niter = 2500,
          max.flips = 100,
          seed = 270488,
          ping = 0
        )

        impact <- CausalImpact(
          bsts.model = bsts_model,
          post.period.response = post.period.response
        )

        impacts[[model_name]] <- impact
        bsts_models[[model_name]] <- bsts_model
        cat(" done\n")
      }, error = function(e) {
        cat(" ERROR:", conditionMessage(e), "\n")
      })
    }

    impacts[["Default"]] <- default_impact
    bsts_models[["Default"]] <- default_impact$model$bsts.model

    cat("[", format(Sys.time(), "%H:%M:%S"), "]", intervention_name, "complete.\n")

    list(impacts = impacts, bsts_models = bsts_models, dates = df_min_A$date)
  }

  # Run for all three interventions
  cat("\n========== RUNNING 8 MODEL SPECIFICATIONS ==========\n")

  results_marenostrum <- run_all_models(df_min_A, 26, 27, "Mare Nostrum", impact_marenostrum, covar_drop_cols = c(1, 2, 3))
  results_sarngos <- run_all_models(df_min_A, 39, 40, "NGO SAR", impact_sarngos, covar_drop_cols = c(1, 2))
  results_sarlibya <- run_all_models(df_min_A, 66, 67, "EU-Libya Cooperation", impact_sarlibya, covar_drop_cols = c(1, 2, 3))

  # ============================================================================
  # EXTRACT EFFECT SIZES
  # ============================================================================

  cat("\n[", format(Sys.time(), "%H:%M:%S"), "] Extracting effect sizes...\n")

  extract_effects <- function(impacts) {
    effects_list <- list()
    for (model_name in names(impacts)) {
      if (!is.null(impacts[[model_name]]$summary)) {
        effects_list[[model_name]] <- data.frame(
          model = model_name,
          AbsEffect = impacts[[model_name]]$summary$AbsEffect[1],
          AbsEffect.lower = impacts[[model_name]]$summary$AbsEffect.lower[1],
          AbsEffect.upper = impacts[[model_name]]$summary$AbsEffect.upper[1]
        )
      }
    }
    do.call(rbind, effects_list)
  }

  effs_marenostrum <- extract_effects(results_marenostrum$impacts)
  effs_sarngos <- extract_effects(results_sarngos$impacts)
  effs_sarlibya <- extract_effects(results_sarlibya$impacts)

  saveRDS(effs_marenostrum, "replicated_results/validation/effs_marenostrum.RDS")
  saveRDS(effs_sarngos, "replicated_results/validation/effs_sarngos.RDS")
  saveRDS(effs_sarlibya, "replicated_results/validation/effs_sarlibya.RDS")

  # ============================================================================
  # FIGURE S6: Cumulative Absolute Error Comparison (8 models)
  # ============================================================================

  cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure S6...\n")

  compare_marenostrum <- CompareBstsModels(results_marenostrum$bsts_models, main = "Mare Nostrum")
  compare_sarngos <- CompareBstsModels(results_sarngos$bsts_models, main = "NGO SAR")
  compare_sarlibya <- CompareBstsModels(results_sarlibya$bsts_models, main = "EU-Libya")

  dates_vec <- results_marenostrum$dates
  model_names_actual_order <- c("LLT", "SLLT", "LL", "LLT+AR", "SLLT+AR", "LL+AR", "LL+SLLT+AR", "Default")

  compare_to_long <- function(compare_mat, dates, model_names) {
    df <- as.data.frame(t(compare_mat))
    colnames(df) <- model_names[1:ncol(df)]
    df$date <- dates[1:nrow(df)]
    df %>%
      pivot_longer(cols = -date, names_to = "Models", values_to = "cum_error") %>%
      mutate(Models = factor(Models, levels = model_names))
  }

  compare_long_marenostrum <- compare_to_long(compare_marenostrum, dates_vec, model_names_actual_order)
  compare_long_sarngos <- compare_to_long(compare_sarngos, dates_vec, model_names_actual_order)
  compare_long_sarlibya <- compare_to_long(compare_sarlibya, dates_vec, model_names_actual_order)

  saveRDS(compare_long_marenostrum, "replicated_results/validation/compare_mod_long_marenostrum.RDS")
  saveRDS(compare_long_sarngos, "replicated_results/validation/compare_mod_long_sarngos.RDS")
  saveRDS(compare_long_sarlibya, "replicated_results/validation/compare_mod_long_sarlibya.RDS")

  model_colors <- c(
    "LLT" = "#E41A1C", "SLLT" = "#4DAF4A", "LL" = "#377EB8",
    "LLT+AR" = "#984EA3", "SLLT+AR" = "#FF7F00", "LL+AR" = "#A65628",
    "LL+SLLT+AR" = "#F781BF", "Default" = "#999999"
  )

  fig_s6a <- ggplot(compare_long_marenostrum, aes(x = date, y = cum_error, color = Models)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = model_colors) +
    theme_classic() +
    ylab("Cumulative Absolute Error") + xlab("Date") +
    ggtitle("A. Mare Nostrum") +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 1))

  fig_s6b <- ggplot(compare_long_sarngos, aes(x = date, y = cum_error, color = Models)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = model_colors) +
    theme_classic() +
    ylab("Cumulative Absolute Error") + xlab("Date") +
    ggtitle("B. NGOs search-and-rescue") +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 1))

  fig_s6c <- ggplot(compare_long_sarlibya, aes(x = date, y = cum_error, color = Models)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = model_colors) +
    theme_classic() +
    ylab("Cumulative Absolute Error") + xlab("Date") +
    ggtitle("C. EU and Libya cooperation") +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 1))

  fig_s6 <- grid.arrange(fig_s6a, fig_s6b, fig_s6c, ncol = 3,
                         top = "Model comparison: Cumulative one-step prediction errors")

  ggsave("replicated_results/figures/supplementary/png/FigureS6_model_validation.png",
         fig_s6, width = 15, height = 6, dpi = 300)
  ggsave("replicated_results/figures/supplementary/pdf/FigureS6_model_validation.pdf",
         fig_s6, width = 15, height = 6)

  cat("[", format(Sys.time(), "%H:%M:%S"), "] Figure S6 saved.\n")

  # ============================================================================
  # FIGURE S7: Effect Size Comparison (8 models)
  # ============================================================================

  cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure S7...\n")

  fig_s7a <- ggplot(effs_marenostrum, aes(x = model, y = AbsEffect,
                                           ymin = AbsEffect.lower, ymax = AbsEffect.upper)) +
    geom_point(size = 2) +
    geom_errorbar(width = 0.3) +
    geom_hline(yintercept = 0, lty = 2, color = "gray50") +
    ggtitle("A. Mare Nostrum") +
    ylab("Average AbsEffect") + xlab("Model") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  fig_s7b <- ggplot(effs_sarngos, aes(x = model, y = AbsEffect,
                                       ymin = AbsEffect.lower, ymax = AbsEffect.upper)) +
    geom_point(size = 2) +
    geom_errorbar(width = 0.3) +
    geom_hline(yintercept = 0, lty = 2, color = "gray50") +
    ggtitle("B. NGO's search-and-rescue") +
    ylab("Average AbsEffect") + xlab("Model") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  fig_s7c <- ggplot(effs_sarlibya, aes(x = model, y = AbsEffect,
                                        ymin = AbsEffect.lower, ymax = AbsEffect.upper)) +
    geom_point(size = 2) +
    geom_errorbar(width = 0.3) +
    geom_hline(yintercept = 0, lty = 2, color = "gray50") +
    ggtitle("C. EU and Libya cooperation") +
    ylab("Average AbsEffect") + xlab("Model") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  fig_s7 <- grid.arrange(fig_s7a, fig_s7b, fig_s7c, ncol = 3,
                         top = "Effect size comparison across different model specifications")

  ggsave("replicated_results/figures/supplementary/png/FigureS7_effect_sizes.png",
         fig_s7, width = 14, height = 5, dpi = 300)
  ggsave("replicated_results/figures/supplementary/pdf/FigureS7_effect_sizes.pdf",
         fig_s7, width = 14, height = 5)

  cat("[", format(Sys.time(), "%H:%M:%S"), "] Figure S7 saved.\n")

  # ============================================================================
  # FIGURE S8: ACF using bsts::AcfDist
  # ============================================================================

  cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure S8...\n")

  png("replicated_results/figures/supplementary/png/FigureS8_acf.png",
      width = 15, height = 5, units = "in", res = 300)
  par(mfrow = c(1, 3))
  AcfDist(impact_marenostrum$model$posterior.samples, main = "A. Mare Nostrum")
  AcfDist(impact_sarngos$model$posterior.samples, main = "B. NGOs")
  AcfDist(impact_sarlibya$model$posterior.samples, main = "C. EU and Libya cooperation")
  dev.off()

  pdf("replicated_results/figures/supplementary/pdf/FigureS8_acf.pdf", width = 15, height = 5)
  par(mfrow = c(1, 3))
  AcfDist(impact_marenostrum$model$posterior.samples, main = "A. Mare Nostrum")
  AcfDist(impact_sarngos$model$posterior.samples, main = "B. NGOs")
  AcfDist(impact_sarlibya$model$posterior.samples, main = "C. EU and Libya cooperation")
  dev.off()

  cat("[", format(Sys.time(), "%H:%M:%S"), "] Figure S8 saved.\n")

  # ============================================================================
  # FIGURE S9: QQ plots using bsts::qqdist
  # ============================================================================

  cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure S9...\n")

  png("replicated_results/figures/supplementary/png/FigureS9_qqplot.png",
      width = 15, height = 5, units = "in", res = 300)
  par(mfrow = c(1, 3))
  qqdist(impact_marenostrum$model$posterior.samples, main = "A. Mare Nostrum")
  qqdist(impact_sarngos$model$posterior.samples, main = "B. NGOs")
  qqdist(impact_sarlibya$model$posterior.samples, main = "C. EU and Libya cooperation")
  dev.off()

  pdf("replicated_results/figures/supplementary/pdf/FigureS9_qqplot.pdf", width = 15, height = 5)
  par(mfrow = c(1, 3))
  qqdist(impact_marenostrum$model$posterior.samples, main = "A. Mare Nostrum")
  qqdist(impact_sarngos$model$posterior.samples, main = "B. NGOs")
  qqdist(impact_sarlibya$model$posterior.samples, main = "C. EU and Libya cooperation")
  dev.off()

  cat("[", format(Sys.time(), "%H:%M:%S"), "] Figure S9 saved.\n")

  # ============================================================================
  # FIGURE S10: Coefficient inclusion probabilities
  # Shows top 30 covariates by inclusion probability for each intervention
  # ============================================================================

  cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure S10...\n")

  extract_coef_plot <- function(impact_obj, title, top_n = 30) {
    if (is.null(impact_obj$model$bsts.model)) return(NULL)

    bsts_model <- impact_obj$model$bsts.model
    coef_means <- colMeans(bsts_model$coefficients != 0)
    coef_values <- colMeans(bsts_model$coefficients)

    coef_df <- data.frame(
      variable = names(coef_means),
      inclusion_prob = coef_means,
      coef_value = coef_values
    ) %>%
      arrange(desc(inclusion_prob)) %>%
      head(top_n) %>%
      mutate(direction = ifelse(coef_value > 0, "positive", "negative"))

    if (nrow(coef_df) == 0) return(NULL)

    # Clean up variable names for display
    coef_df <- coef_df %>%
      mutate(variable = gsub("_", " ", variable),
             variable = gsub("lag 0[1-6]", "", variable),
             variable = str_trunc(variable, 50))

    ggplot(coef_df, aes(x = reorder(variable, inclusion_prob), y = inclusion_prob, fill = direction)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c("positive" = "steelblue", "negative" = "coral")) +
      labs(title = title, x = "", y = "Inclusion Probability") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 7),
            legend.position = "bottom",
            legend.title = element_blank())
  }

  fig_s10a <- extract_coef_plot(impact_marenostrum, "A. Mare Nostrum")
  fig_s10b <- extract_coef_plot(impact_sarngos, "B. NGOs search-and-rescue")
  fig_s10c <- extract_coef_plot(impact_sarlibya, "C. EU and Libya cooperation")

  if (!is.null(fig_s10a) && !is.null(fig_s10b) && !is.null(fig_s10c)) {
    fig_s10 <- grid.arrange(fig_s10a, fig_s10b, fig_s10c, ncol = 3,
                            top = "Covariate inclusion probabilities in BSTS models")

    ggsave("replicated_results/figures/supplementary/png/FigureS10_coefficients.png",
           fig_s10, width = 18, height = 10, dpi = 300)
    ggsave("replicated_results/figures/supplementary/pdf/FigureS10_coefficients.pdf",
           fig_s10, width = 18, height = 10)

    cat("[", format(Sys.time(), "%H:%M:%S"), "] Figure S10 saved.\n")
  } else {
    cat("[", format(Sys.time(), "%H:%M:%S"), "] Could not create Figure S10.\n")
  }
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=======================================================================\n")
cat("SUPPLEMENTARY FIGURES COMPLETE!\n")
cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=======================================================================\n\n")

cat("Generated figures:\n")
supp_files <- list.files("replicated_results/figures/supplementary/png", pattern = "\\.png$")
for (f in sort(supp_files)) {
  cat("  -", f, "\n")
}

cat("\nNOTE: Figure S11 (sensitivity analysis) requires re-running models\n")
cat("      without labor market indicators and is not included.\n")
