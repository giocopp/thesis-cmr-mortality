################################################################################
# Mortality Counterfactual Model
# ==============================
# Extension of Rodriguez Sanchez et al. (2023)
#
# Research question: Did the MoU / EU-Libya cooperation (Feb 2017) make the
# Central Mediterranean route more dangerous?
#
# Method: Same as original — BSTS + CausalImpact with spike-and-slab priors.
# Outcome: Mortality rate (deaths per 100 crossing attempts), log-transformed.
# Covariates: Exogenous push-pull factors + NEW sea condition variables (ERA5).
# Intervention: MoU / EU-Libya cooperation (Feb 2017 onwards).
#
# Additionally tests Mare Nostrum and NGO SAR periods, asking whether they
# REDUCED mortality (the reverse of the original paper's pull-factor question).
#
# TWO SPECIFICATIONS:
#   A. CURATED (primary): ~200-800 theoretically motivated predictors with
#      max.flips=-1 (full MCMC exploration). Focuses on sea conditions, weather,
#      key origin/transit country conflicts, destination unemployment, and
#      key exchange rates.
#   B. FULL (robustness): ~4,800 predictors (same filtering as original paper)
#      with max.flips=100. Tests whether results hold under kitchen-sink approach.
#
# PREREQUISITES:
#   - df_extended.RDS (built by build_mortality_dataset.R)
#   - R packages: CausalImpact, bsts, tidyverse, lubridate, scales, gridExtra
#
# USAGE:
#   source("Extension-1-BSTS-mortality/code/run_mortality_model.R")
################################################################################

cat("=", rep("=", 70), "\n", sep = "")
cat("MORTALITY COUNTERFACTUAL MODEL\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 71), "\n\n", sep = "")

# ==============================================================================
# 0. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(gridExtra)
  library(CausalImpact)
  library(bsts)
})

# --- Paths ---
THESIS_ROOT <- tryCatch({
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, "..", ".."))
}, error = function(e) {
  getwd()
})

DATA_DIR <- file.path(THESIS_ROOT, "Extension-1-BSTS-mortality", "data")
RESULTS_DIR <- file.path(THESIS_ROOT, "Extension-1-BSTS-mortality", "results")

dir.create(file.path(RESULTS_DIR, "models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "figures", "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "figures", "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "data"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Loading extended dataset...\n")

df <- readRDS(file.path(DATA_DIR, "df_extended.RDS"))

# Ensure outcome variables exist
df <- df %>%
  mutate(
    LCG_pushbacks_count = as.numeric(ifelse(is.na(LCG_pushbacks_count), 0, LCG_pushbacks_count)),
    TCG_pushbacks_count = as.numeric(ifelse(is.na(TCG_pushbacks_count), 0, TCG_pushbacks_count)),
    dead_and_missing_Central_Mediterranean = as.numeric(
      ifelse(is.na(dead_and_missing_Central_Mediterranean), 0,
             dead_and_missing_Central_Mediterranean)
    ),
    crossings_CMR = arrivals_CMR + LCG_pushbacks_count + TCG_pushbacks_count +
      dead_and_missing_Central_Mediterranean,
    mortality_rate_100 = (dead_and_missing_Central_Mediterranean / crossings_CMR) * 100
  )

cat("  Full dataset:", nrow(df), "rows x", ncol(df), "columns\n\n")


# ==============================================================================
# 2. VARIABLE FILTERING — CURATED SPECIFICATION (PRIMARY)
# ==============================================================================
#
# The original paper uses ~4,800 predictors with spike-and-slab variable
# selection (max.flips=100). This works well for crossings because many
# predictors (airport flows, seasonality, commodity prices) genuinely predict
# migration flows.
#
# For MORTALITY, these macro push-pull predictors have very weak signal —
# mortality depends on per-crossing risk factors (sea conditions, boat quality,
# rescue response) rather than migration decisions. Diagnostic analysis of the
# full-predictor model shows:
#   - Only 2 of 4,819 predictors selected (spurious: Algeria violence, soy prices)
#   - All 63 ERA5 sea condition variables have 0.000 inclusion probability
#   - Models run in ~25 sec (vs hours for crossings) = effectively empty model
#
# SOLUTION: Curate a smaller, theoretically motivated predictor set and use
# max.flips=-1 for exhaustive MCMC exploration at each iteration.
#
# Curated variables:
#   1. ERA5 sea conditions (wave height, wind, SST, cloud, fog) — base + lags 01-06
#   2. Weather in Italy/Malta/departure coast (temperature, precipitation, storms)
#   3. Conflicts in key Central Med origin/transit countries:
#      Libya (transit), Eritrea, Somalia, Nigeria, Sudan/South Sudan, Syria,
#      Tunisia, Egypt — UCDP events + disaster counts
#   4. Exchange rates for key origin currencies (EGP, NGN, TND, LYD, ETB, SOS,
#      SYP, DZD, MAD, GHS, GMD, KES)
#   5. Unemployment in key EU destination countries (Italy, Greece, Malta,
#      Germany, France, Spain, EU aggregate)
#   6. Key commodity prices (oil, energy index, food index)
#   7. Google Trends (employment/job searches in origin countries)
#   8. Time components (month, quarter, semester)
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Building CURATED predictor set...\n")

# First apply basic filtering: time window + remove high-order lags + remove
# endogenous variables (same removals as the full specification)
df_base <- df %>%
  filter(date >= "2011-02-01" & date < "2021-10-01") %>%
  dplyr::select(
    -c(
      # High-order lags (keep only 01-06)
      contains("lag_24"), contains("lag_23"), contains("lag_22"),
      contains("lag_21"), contains("lag_20"), contains("lag_19"),
      contains("lag_18"), contains("lag_17"), contains("lag_16"),
      contains("lag_15"), contains("lag_14"), contains("lag_13"),
      contains("lag_12"), contains("lag_11"), contains("lag_10"),
      contains("lag_09"), contains("lag_08"), contains("lag_07"),

      # Problematic variables (same as original)
      starts_with("airflow_Palestinian.Territories"),
      starts_with("asylum"),

      # Outcome-related: crossing volumes (ENDOGENOUS to MoU)
      any_of(c("arrivals_BSR", "arrivals_CMR", "arrivals_CRAG", "arrivals_EBR",
               "arrivals_EMR", "arrivals_OR", "arrivals_WAR", "arrivals_WBR",
               "arrivals_WMR", "crossings_CMR")),

      # Outcome-related: deaths on other routes
      any_of(c("dead_and_missing_Eastern_Mediterranean",
               "dead_and_missing_Western_Mediterranean")),

      # Endogenous: pushbacks (created by MoU)
      any_of(c("LCG_pushbacks_count", "TCG_pushbacks_count")),

      # Endogenous: geographic dispersion of deaths
      starts_with("sd_lat_"), starts_with("sd_lon_"),
      starts_with("frac_index_"),

      # Endogenous: SAR vessel and EU operation indicators
      starts_with("SAR_"), starts_with("FRONTEX_"), starts_with("EUNAVFOR_"),
      any_of(c("MARE_NOSTRUM", "COASTGUARD_LIBYA", "Extension_SAR_LIBYA")),

      # Not a predictor
      any_of(c("year"))
    )
  )

# Now select ONLY curated predictors from the base-filtered data
# Key Central Med origin/transit country patterns for UCDP + disasters
curated_countries <- c(
  "Libya", "Eritrea", "Somalia", "Nigeria", "Tunisia", "Egypt",
  "Sudan",           # matches both Sudan and South.Sudan / South Sudan
  "Syrian",          # matches Syrian.Arab.Republic
  "Syria"            # matches syria_trend_google
)
country_regex <- paste0("(", paste(curated_countries, collapse = "|"), ")")

# Key exchange rate currencies
curated_fx <- c("EGP", "NGN", "TND", "LYD", "ETB", "SOS", "SYP", "DZD",
                "MAD", "GHS", "GMD", "KES")
fx_regex <- paste0("^(", paste(curated_fx, collapse = "|"), ")_to_EURO")

# Key destination unemployment
curated_unemp <- c("ITALY", "GREECE", "MALTA", "GERMANY", "FRANCE", "SPAIN",
                    "eu_27", "euro_area_all")
unemp_regex <- paste0("^unem_(", paste(curated_unemp, collapse = "|"), ")")

# Key commodity prices
curated_commodities <- c("POILBRE", "PFOOD", "PNRG", "PFANDB")
commodity_regex <- paste0("^(", paste(curated_commodities, collapse = "|"), ")")

df_curated <- df_base %>%
  dplyr::select(
    date,
    dead_and_missing_Central_Mediterranean,
    mortality_rate_100,
    # 1. ERA5 sea conditions
    matches("^(wave_|wind_speed|sst_|cloud_cover|dewpoint_depression)"),
    # 2. Weather (Italy, Malta, central med, departure coast)
    matches("^(temperature_|precipitation_|daysstorm_)"),
    # 3. Key origin/transit country conflicts + disasters
    matches(country_regex),
    # 4. Key exchange rates
    matches(fx_regex),
    # 5. Key destination unemployment
    matches(unemp_regex),
    # 6. Key commodity prices
    matches(commodity_regex),
    # 7. Google Trends
    matches("^googlesearch_|syria_trend")
  )

cat("  Curated base dataset:", nrow(df_curated), "rows x", ncol(df_curated), "columns\n")


# ==============================================================================
# 2a. PREPARE CURATED MORTALITY RATE OUTCOME
# ==============================================================================

EPSILON <- 0.01

df_curated_mort <- data.frame(
  date = df_curated$date,
  mortality_rate = log(df_curated$mortality_rate_100 + EPSILON),
  dplyr::select(df_curated, -c(date, dead_and_missing_Central_Mediterranean,
                                 mortality_rate_100))
)

df_curated_mort <- df_curated_mort %>%
  mutate(
    month = month(date),
    semester = semester(date),
    quarter = quarter(date)
  )

df_curated_mort_clean <- df_curated_mort %>% na.omit()

cat("  Curated mortality rate: ", nrow(df_curated_mort_clean), " rows x ",
    ncol(df_curated_mort_clean), " columns\n", sep = "")
cat("  Curated predictors: ", ncol(df_curated_mort_clean) - 2,
    " (excluding date and outcome)\n\n", sep = "")


# ==============================================================================
# 2b. PREPARE CURATED DEATH COUNT OUTCOME (robustness)
# ==============================================================================

df_curated_deaths <- data.frame(
  date = df_curated$date,
  deaths_cmr = log(df_curated$dead_and_missing_Central_Mediterranean + 1),
  dplyr::select(df_curated, -c(date, dead_and_missing_Central_Mediterranean,
                                 mortality_rate_100))
)

df_curated_deaths <- df_curated_deaths %>%
  mutate(
    month = month(date),
    semester = semester(date),
    quarter = quarter(date)
  )

df_curated_deaths_clean <- df_curated_deaths %>% na.omit()


# ==============================================================================
# 2c. PREPARE FULL SPECIFICATION (for robustness comparison)
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Building FULL predictor set (robustness)...\n")

df_full_mort <- data.frame(
  date = df_base$date,
  mortality_rate = log(df_base$mortality_rate_100 + EPSILON),
  dplyr::select(df_base, -c(date, dead_and_missing_Central_Mediterranean,
                              mortality_rate_100))
)

df_full_mort <- df_full_mort %>%
  mutate(
    month = month(date),
    semester = semester(date),
    quarter = quarter(date)
  )

df_full_mort_clean <- df_full_mort %>% na.omit()

df_full_deaths <- data.frame(
  date = df_base$date,
  deaths_cmr = log(df_base$dead_and_missing_Central_Mediterranean + 1),
  dplyr::select(df_base, -c(date, dead_and_missing_Central_Mediterranean,
                              mortality_rate_100))
)

df_full_deaths <- df_full_deaths %>%
  mutate(
    month = month(date),
    semester = semester(date),
    quarter = quarter(date)
  )

df_full_deaths_clean <- df_full_deaths %>% na.omit()

cat("  Full mortality rate: ", nrow(df_full_mort_clean), " rows x ",
    ncol(df_full_mort_clean), " columns\n", sep = "")
cat("  Full predictors: ", ncol(df_full_mort_clean) - 2,
    " (excluding date and outcome)\n\n", sep = "")


# ==============================================================================
# 3. FIGURE 1: DESCRIPTIVE TIME SERIES — MORTALITY FOCUS
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Figure 1 (mortality descriptive)...\n")

intervention_dates <- list(
  mare_nostrum_start = ymd("2013-10-01"),
  mare_nostrum_end   = ymd("2014-10-01"),
  ngo_sar_start      = ymd("2014-11-01"),
  eu_libya_start     = ymd("2017-02-01"),
  lcg_sar_zone       = ymd("2017-08-01")
)

theme_paper <- function() {
  theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      axis.title = element_text(size = 14, color = "black"),
      axis.text  = element_text(size = 12, color = "black"),
      axis.line  = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

x_breaks <- ymd(c("2010-01-01", "2015-01-01", "2020-01-01"))
sar_start <- intervention_dates$mare_nostrum_start
sar_end   <- intervention_dates$eu_libya_start

# Panel A: Dead and missing over time
fig1a <- ggplot(df, aes(x = date, y = dead_and_missing_Central_Mediterranean)) +
  annotate("rect",
    xmin = intervention_dates$mare_nostrum_start,
    xmax = intervention_dates$mare_nostrum_end,
    ymin = -Inf, ymax = Inf, fill = "#B8D4E8", alpha = 0.35) +
  annotate("rect",
    xmin = intervention_dates$ngo_sar_start,
    xmax = intervention_dates$eu_libya_start,
    ymin = -Inf, ymax = Inf, fill = "#FFCCCC", alpha = 0.30) +
  annotate("rect",
    xmin = intervention_dates$eu_libya_start,
    xmax = intervention_dates$lcg_sar_zone,
    ymin = -Inf, ymax = Inf, fill = "#FFE4B5", alpha = 0.30) +
  annotate("rect",
    xmin = intervention_dates$lcg_sar_zone,
    xmax = max(df$date, na.rm = TRUE),
    ymin = -Inf, ymax = Inf, fill = "#F5DEB3", alpha = 0.30) +
  annotate("text", x = ymd("2014-04-01"), y = 1350, label = "EU\nMare Nostrum",
           size = 3.2, color = "blue") +
  annotate("text", x = ymd("2016-01-01"), y = 1350, label = "NGOs\nSAR",
           size = 3.2, color = "red") +
  annotate("text", x = ymd("2017-05-01"), y = 1350, label = "EU-LCG\nDeal",
           size = 3.2, color = "darkorange") +
  annotate("text", x = ymd("2019-06-01"), y = 1350, label = "Expansion\nLCG SAR-Zone",
           size = 3.2, color = "darkgreen") +
  geom_line(linewidth = 0.6, color = "black") +
  scale_x_date(breaks = x_breaks, date_labels = "%Y",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = comma) +
  labs(title = "A. Dead and missing in the Central Mediterranean",
       x = "Date", y = "Dead and missing (monthly count)") +
  theme_paper()

# Panel B: Mortality rate per 100 crossings
fig1b <- ggplot(df, aes(x = date, y = mortality_rate_100)) +
  annotate("rect",
    xmin = sar_start, xmax = sar_end,
    ymin = -Inf, ymax = Inf, fill = "#E6E6FA", alpha = 0.25) +
  annotate("text", x = ymd("2015-06-01"), y = 30,
           label = "Search-and-rescue", size = 4.5, color = "purple") +
  geom_line(linewidth = 0.6, color = "black") +
  scale_x_date(breaks = x_breaks, date_labels = "%Y",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(title = "B. Mortality rate per 100 attempted crossings",
       x = "Date", y = "Deaths per 100 crossings") +
  theme_paper()

figure1 <- grid.arrange(
  fig1a, fig1b, ncol = 1, heights = c(1, 1),
  bottom = grid::textGrob("Note: own calculations.",
                           x = 0.5, hjust = 0.5,
                           gp = grid::gpar(fontsize = 9))
)

ggsave(file.path(RESULTS_DIR, "figures", "png", "Figure1_mortality.png"),
       figure1, width = 8.5, height = 8.5, dpi = 300)
ggsave(file.path(RESULTS_DIR, "figures", "pdf", "Figure1_mortality.pdf"),
       figure1, width = 8.5, height = 8.5)

cat("  Figure 1 saved.\n\n")


# ==============================================================================
# 4. DEFINE INTERVENTION PERIODS (identical to original paper)
# ==============================================================================

pre.period_mare_nostrum  <- c(min(df_curated_mort_clean$date), ymd("2013-09-01"))
post.period_mare_nostrum <- c(ymd("2013-10-01"), max(df_curated_mort_clean$date))

pre.period_sar_ngos  <- c(min(df_curated_mort_clean$date), ymd("2014-10-01"))
post.period_sar_ngos <- c(ymd("2014-11-01"), max(df_curated_mort_clean$date))

pre.period_sarlibya   <- c(min(df_curated_mort_clean$date), ymd("2017-01-01"))
post.period_sar_libya <- c(ymd("2017-02-01"), max(df_curated_mort_clean$date))

# Seed (identical to original paper)
set.seed(270488)

# Model arguments — CURATED specification
# max.flips = -1 allows the MCMC to consider flipping every variable indicator
# at each iteration, giving much better mixing with the smaller predictor set
model_args_curated <- list(
  dynamic.regression = FALSE,
  standardize.data = TRUE,
  max.flips = -1,
  niter = 10000
)

# Model arguments — FULL specification (same as original paper)
model_args_full <- list(
  dynamic.regression = FALSE,
  standardize.data = TRUE,
  max.flips = 100,
  niter = 10000
)


# ==============================================================================
# 5. RUN CURATED MODELS — MORTALITY RATE (PRIMARY)
# ==============================================================================

cat(rep("=", 71), "\n", sep = "")
cat("CURATED SPECIFICATION — MORTALITY RATE\n")
cat(rep("=", 71), "\n\n", sep = "")

# --- Model A: Mare Nostrum — did it reduce mortality? ---
cat("[", format(Sys.time(), "%H:%M:%S"), "] Curated Model A: Mare Nostrum (mortality rate)...\n")
t1 <- Sys.time()

impact_cur_A <- CausalImpact(
  df_curated_mort_clean,
  pre.period = pre.period_mare_nostrum,
  post.period = post.period_mare_nostrum,
  alpha = 0.05,
  model.args = model_args_curated
)

t2 <- Sys.time()
cat("[", format(Sys.time(), "%H:%M:%S"), "] Curated Model A done in",
    round(difftime(t2, t1, units = "mins"), 1), "minutes.\n\n")
saveRDS(impact_cur_A,
        file.path(RESULTS_DIR, "models", "model_curated_mortality_marenostrum.RDS"))


# --- Model B: NGO SAR — did it reduce mortality? ---
cat("[", format(Sys.time(), "%H:%M:%S"), "] Curated Model B: NGO SAR (mortality rate)...\n")
t1 <- Sys.time()

impact_cur_B <- CausalImpact(
  df_curated_mort_clean,
  pre.period = pre.period_sar_ngos,
  post.period = post.period_sar_ngos,
  alpha = 0.05,
  model.args = model_args_curated
)

t2 <- Sys.time()
cat("[", format(Sys.time(), "%H:%M:%S"), "] Curated Model B done in",
    round(difftime(t2, t1, units = "mins"), 1), "minutes.\n\n")
saveRDS(impact_cur_B,
        file.path(RESULTS_DIR, "models", "model_curated_mortality_sarngos.RDS"))


# --- Model C: EU-Libya / MoU — did it increase mortality? ---
cat("[", format(Sys.time(), "%H:%M:%S"), "] Curated Model C: EU-Libya / MoU (mortality rate)...\n")
t1 <- Sys.time()

impact_cur_C <- CausalImpact(
  df_curated_mort_clean,
  pre.period = pre.period_sarlibya,
  post.period = post.period_sar_libya,
  alpha = 0.05,
  model.args = model_args_curated
)

t2 <- Sys.time()
cat("[", format(Sys.time(), "%H:%M:%S"), "] Curated Model C done in",
    round(difftime(t2, t1, units = "mins"), 1), "minutes.\n\n")
saveRDS(impact_cur_C,
        file.path(RESULTS_DIR, "models", "model_curated_mortality_sarlibya.RDS"))


# ==============================================================================
# 6. RUN CURATED MODELS — DEATH COUNT (robustness outcome)
# ==============================================================================

cat(rep("=", 71), "\n", sep = "")
cat("CURATED SPECIFICATION — DEATH COUNT (robustness)\n")
cat(rep("=", 71), "\n\n", sep = "")

pre.period_sarlibya_d   <- c(min(df_curated_deaths_clean$date), ymd("2017-01-01"))
post.period_sar_libya_d <- c(ymd("2017-02-01"), max(df_curated_deaths_clean$date))

cat("[", format(Sys.time(), "%H:%M:%S"), "] Curated Model C (death count)...\n")
t1 <- Sys.time()

impact_cur_C_deaths <- CausalImpact(
  df_curated_deaths_clean,
  pre.period = pre.period_sarlibya_d,
  post.period = post.period_sar_libya_d,
  alpha = 0.05,
  model.args = model_args_curated
)

t2 <- Sys.time()
cat("[", format(Sys.time(), "%H:%M:%S"), "] Curated deaths model done in",
    round(difftime(t2, t1, units = "mins"), 1), "minutes.\n\n")
saveRDS(impact_cur_C_deaths,
        file.path(RESULTS_DIR, "models", "model_curated_deaths_sarlibya.RDS"))


# ==============================================================================
# 7. RUN FULL SPECIFICATION — MODEL C ONLY (robustness predictor set)
# ==============================================================================

cat(rep("=", 71), "\n", sep = "")
cat("FULL SPECIFICATION — MODEL C ONLY (robustness)\n")
cat(rep("=", 71), "\n\n", sep = "")

cat("[", format(Sys.time(), "%H:%M:%S"), "] Full Model C: EU-Libya / MoU (mortality rate)...\n")
t1 <- Sys.time()

impact_full_C <- CausalImpact(
  df_full_mort_clean,
  pre.period = pre.period_sarlibya,
  post.period = post.period_sar_libya,
  alpha = 0.05,
  model.args = model_args_full
)

t2 <- Sys.time()
cat("[", format(Sys.time(), "%H:%M:%S"), "] Full Model C done in",
    round(difftime(t2, t1, units = "mins"), 1), "minutes.\n\n")
saveRDS(impact_full_C,
        file.path(RESULTS_DIR, "models", "model_full_mortality_sarlibya.RDS"))

cat("[", format(Sys.time(), "%H:%M:%S"), "] Full Model C (death count)...\n")
t1 <- Sys.time()

impact_full_C_deaths <- CausalImpact(
  df_full_deaths_clean,
  pre.period = pre.period_sarlibya_d,
  post.period = post.period_sar_libya_d,
  alpha = 0.05,
  model.args = model_args_full
)

t2 <- Sys.time()
cat("[", format(Sys.time(), "%H:%M:%S"), "] Full deaths model done in",
    round(difftime(t2, t1, units = "mins"), 1), "minutes.\n\n")
saveRDS(impact_full_C_deaths,
        file.path(RESULTS_DIR, "models", "model_full_deaths_sarlibya.RDS"))


cat("[", format(Sys.time(), "%H:%M:%S"), "] All models fitted.\n\n")


# ==============================================================================
# 8. EXTRACT RESULTS
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Extracting results...\n")

extract_results <- function(impact_obj, dates) {
  data.frame(
    date = dates,
    original = as.numeric(impact_obj$series$response),
    prediction = as.numeric(impact_obj$series$point.pred),
    prediction_lower = as.numeric(impact_obj$series$point.pred.lower),
    prediction_upper = as.numeric(impact_obj$series$point.pred.upper),
    pointwise_effect = as.numeric(impact_obj$series$point.effect),
    pointwise_effect_lower = as.numeric(impact_obj$series$point.effect.lower),
    pointwise_effect_upper = as.numeric(impact_obj$series$point.effect.upper),
    cumulative_effect = as.numeric(impact_obj$series$cum.effect),
    cumulative_effect_lower = as.numeric(impact_obj$series$cum.effect.lower),
    cumulative_effect_upper = as.numeric(impact_obj$series$cum.effect.upper)
  )
}

# Curated results
results_cur_A <- extract_results(impact_cur_A, df_curated_mort_clean$date)
results_cur_B <- extract_results(impact_cur_B, df_curated_mort_clean$date)
results_cur_C <- extract_results(impact_cur_C, df_curated_mort_clean$date)
results_cur_C_deaths <- extract_results(impact_cur_C_deaths, df_curated_deaths_clean$date)

# Full results
results_full_C <- extract_results(impact_full_C, df_full_mort_clean$date)
results_full_C_deaths <- extract_results(impact_full_C_deaths, df_full_deaths_clean$date)

saveRDS(
  list(
    curated_A = results_cur_A,
    curated_B = results_cur_B,
    curated_C = results_cur_C,
    curated_C_deaths = results_cur_C_deaths,
    full_C = results_full_C,
    full_C_deaths = results_full_C_deaths
  ),
  file.path(RESULTS_DIR, "data", "mortality_model_results.RDS")
)


# ==============================================================================
# 9. CREATE FIGURES — CURATED SPECIFICATION (PRIMARY)
# ==============================================================================

cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating figures...\n")

# Styling (matching original paper)
fill_pink <- "#FDE6E6"
fill_grey <- "#D9D9D9"
fill_ci   <- "#B7C9E2"
col_cf    <- "#6A6FB0"
col_smooth <- "grey40"
col_obs   <- "black"

theme_panel <- function() {
  theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      axis.title = element_text(size = 14, color = "black"),
      axis.text  = element_text(size = 12, color = "black"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.line  = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

moving_average <- function(x, k = 6) {
  n <- length(x)
  half <- floor(k / 2)
  vapply(seq_len(n), function(i) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    mean(x[lo:hi], na.rm = TRUE)
  }, numeric(1))
}

results_cur_A <- results_cur_A %>%
  mutate(date = as.Date(date), original_smooth = moving_average(original, k = 6))
results_cur_B <- results_cur_B %>%
  mutate(date = as.Date(date), original_smooth = moving_average(original, k = 6))
results_cur_C <- results_cur_C %>%
  mutate(date = as.Date(date), original_smooth = moving_average(original, k = 6))

date_min <- min(results_cur_A$date, na.rm = TRUE)
date_max <- max(results_cur_A$date, na.rm = TRUE)
x_breaks_2y <- seq(as.Date("2012-01-01"), date_max, by = "2 years")

mn_start  <- as.Date("2013-10-01")
mn_end    <- as.Date("2014-10-01")
ngo_start <- as.Date("2014-11-01")
ngo_end   <- as.Date("2017-02-01")
eu_start  <- as.Date("2017-02-01")

# ---- Figure 2: Counterfactual — Mortality Rate (Curated) ----

fig2a <- ggplot(results_cur_A, aes(x = date)) +
  annotate("rect", xmin = mn_start, xmax = mn_end,
           ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60) +
  annotate("rect", xmin = mn_end, xmax = date_max,
           ymin = -Inf, ymax = Inf, fill = fill_grey, alpha = 0.60) +
  geom_ribbon(aes(ymin = prediction_lower, ymax = prediction_upper),
              fill = fill_ci, alpha = 0.35) +
  geom_line(aes(y = original_smooth), color = col_smooth, linewidth = 0.8) +
  geom_line(aes(y = original), color = col_obs, linewidth = 0.8) +
  geom_line(aes(y = prediction), linetype = 2, color = col_cf, linewidth = 0.8) +
  geom_vline(xintercept = mn_start, color = "darkred", linetype = 2, linewidth = 0.8) +
  scale_x_date(breaks = x_breaks_2y, date_labels = "%y",
               limits = c(date_min, date_max), expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "A. Mare Nostrum (state-led SAR and\nanti-smuggler operations)",
       x = "", y = "") +
  theme_panel()

fig2b <- ggplot(results_cur_B, aes(x = date)) +
  annotate("rect", xmin = ngo_start, xmax = ngo_end,
           ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60) +
  annotate("rect", xmin = ngo_end, xmax = date_max,
           ymin = -Inf, ymax = Inf, fill = fill_grey, alpha = 0.60) +
  geom_ribbon(aes(ymin = prediction_lower, ymax = prediction_upper),
              fill = fill_ci, alpha = 0.35) +
  geom_line(aes(y = original_smooth), color = col_smooth, linewidth = 0.8) +
  geom_line(aes(y = original), color = col_obs, linewidth = 0.8) +
  geom_line(aes(y = prediction), linetype = 2, color = col_cf, linewidth = 0.8) +
  geom_vline(xintercept = ngo_start, color = "darkred", linetype = 2, linewidth = 0.8) +
  scale_x_date(breaks = x_breaks_2y, date_labels = "%y",
               limits = c(date_min, date_max), expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "B. NGOs (private-led SAR\nby various actors)",
       x = "", y = "Log mortality rate (deaths per 100 crossings)") +
  theme_panel()

fig2c <- ggplot(results_cur_C, aes(x = date)) +
  annotate("rect", xmin = eu_start, xmax = date_max,
           ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60) +
  geom_ribbon(aes(ymin = prediction_lower, ymax = prediction_upper),
              fill = fill_ci, alpha = 0.35) +
  geom_line(aes(y = original_smooth), color = col_smooth, linewidth = 0.8) +
  geom_line(aes(y = original), color = col_obs, linewidth = 0.8) +
  geom_line(aes(y = prediction), linetype = 2, color = col_cf, linewidth = 0.8) +
  geom_vline(xintercept = eu_start, color = "darkred", linetype = 2, linewidth = 0.8) +
  scale_x_date(breaks = x_breaks_2y, date_labels = "%y",
               limits = c(date_min, date_max), expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "C. EU-Libya cooperation / MoU\n(pushbacks and LCG SAR-zone extension)",
       x = "Date", y = "") +
  theme_panel()

note_fig2 <- "Note: own calculations. Curated specification. Outcome = log(mortality rate per 100 crossings + 0.01)."

figure2 <- grid.arrange(
  fig2a, fig2b, fig2c, ncol = 1,
  bottom = grid::textGrob(note_fig2, x = 0.01, hjust = 0,
                           gp = grid::gpar(fontsize = 9))
)

ggsave(file.path(RESULTS_DIR, "figures", "png", "Figure2_mortality_counterfactual.png"),
       figure2, width = 8, height = 10, dpi = 300)
ggsave(file.path(RESULTS_DIR, "figures", "pdf", "Figure2_mortality_counterfactual.pdf"),
       figure2, width = 8, height = 10)


# ---- Figure 3: Pointwise Effects — Mortality Rate (Curated) ----

fig3a <- ggplot(results_cur_A, aes(x = date)) +
  annotate("rect", xmin = mn_start, xmax = mn_end,
           ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60) +
  annotate("rect", xmin = mn_end, xmax = date_max,
           ymin = -Inf, ymax = Inf, fill = fill_grey, alpha = 0.60) +
  geom_ribbon(aes(ymin = pointwise_effect_lower, ymax = pointwise_effect_upper),
              fill = fill_ci, alpha = 0.35) +
  geom_line(aes(y = pointwise_effect), linetype = 2, color = col_cf, linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 0.8) +
  geom_vline(xintercept = mn_start, color = "darkred", linetype = 2, linewidth = 0.8) +
  scale_x_date(breaks = x_breaks_2y, date_labels = "%y",
               limits = c(date_min, date_max), expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "A. Mare Nostrum", x = "", y = "") +
  theme_panel()

fig3b <- ggplot(results_cur_B, aes(x = date)) +
  annotate("rect", xmin = ngo_start, xmax = ngo_end,
           ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60) +
  annotate("rect", xmin = ngo_end, xmax = date_max,
           ymin = -Inf, ymax = Inf, fill = fill_grey, alpha = 0.60) +
  geom_ribbon(aes(ymin = pointwise_effect_lower, ymax = pointwise_effect_upper),
              fill = fill_ci, alpha = 0.35) +
  geom_line(aes(y = pointwise_effect), linetype = 2, color = col_cf, linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 0.8) +
  geom_vline(xintercept = ngo_start, color = "darkred", linetype = 2, linewidth = 0.8) +
  scale_x_date(breaks = x_breaks_2y, date_labels = "%y",
               limits = c(date_min, date_max), expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "B. NGO search-and-rescue",
       x = "", y = "Diff. between observed and predicted") +
  theme_panel()

fig3c <- ggplot(results_cur_C, aes(x = date)) +
  annotate("rect", xmin = eu_start, xmax = date_max,
           ymin = -Inf, ymax = Inf, fill = fill_pink, alpha = 0.60) +
  geom_ribbon(aes(ymin = pointwise_effect_lower, ymax = pointwise_effect_upper),
              fill = fill_ci, alpha = 0.35) +
  geom_line(aes(y = pointwise_effect), linetype = 2, color = col_cf, linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 0.8) +
  geom_vline(xintercept = eu_start, color = "darkred", linetype = 2, linewidth = 0.8) +
  scale_x_date(breaks = x_breaks_2y, date_labels = "%y",
               limits = c(date_min, date_max), expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "C. EU-Libya cooperation / MoU", x = "Date", y = "") +
  theme_panel()

note_fig3 <- "Note: own calculations. Positive values = observed mortality HIGHER than predicted (route more dangerous)."

figure3 <- grid.arrange(
  fig3a, fig3b, fig3c, ncol = 1,
  bottom = grid::textGrob(note_fig3, x = 0.01, hjust = 0,
                           gp = grid::gpar(fontsize = 9))
)

ggsave(file.path(RESULTS_DIR, "figures", "png", "Figure3_mortality_pointwise.png"),
       figure3, width = 8, height = 10, dpi = 300)
ggsave(file.path(RESULTS_DIR, "figures", "pdf", "Figure3_mortality_pointwise.pdf"),
       figure3, width = 8, height = 10)

cat("  Figures 2 and 3 saved.\n\n")


# ==============================================================================
# 10. PRINT RESULTS SUMMARY
# ==============================================================================

cat(rep("=", 71), "\n", sep = "")
cat("MORTALITY MODEL RESULTS SUMMARY\n")
cat(rep("=", 71), "\n\n", sep = "")

print_model_summary <- function(label, impact_obj) {
  cat(label, "\n")
  cat("   Absolute Effect:", round(impact_obj$summary$AbsEffect[2], 4), "\n")
  cat("   Relative Effect:", round(impact_obj$summary$RelEffect[2] * 100, 2), "%\n")
  cat("   95% CI: [", round(impact_obj$summary$AbsEffect.lower[2], 4), ", ",
      round(impact_obj$summary$AbsEffect.upper[2], 4), "]\n", sep = "")
  cat("   p-value:", round(impact_obj$summary$p[2], 4), "\n")
  cat("   Significant:", ifelse(impact_obj$summary$p[2] < 0.05, "YES", "NO"), "\n")

  if (impact_obj$summary$AbsEffect[2] > 0) {
    cat("   Direction: Mortality HIGHER than predicted (route MORE dangerous)\n")
  } else {
    cat("   Direction: Mortality LOWER than predicted (route LESS dangerous)\n")
  }
  cat("\n")
}

cat("--- CURATED SPECIFICATION: Mortality Rate ---\n\n")
print_model_summary("A. MARE NOSTRUM (Oct 2013 - Oct 2014)", impact_cur_A)
print_model_summary("B. NGO SEARCH-AND-RESCUE (Nov 2014 onwards)", impact_cur_B)
print_model_summary("C. EU-LIBYA / MoU (Feb 2017 onwards)", impact_cur_C)

cat("--- CURATED SPECIFICATION: Death Count (robustness) ---\n\n")
print_model_summary("C. EU-LIBYA / MoU (Feb 2017 onwards) — death count", impact_cur_C_deaths)

cat("--- FULL SPECIFICATION: Model C only (robustness) ---\n\n")
print_model_summary("C. EU-LIBYA / MoU (Feb 2017 onwards) — mortality rate, full predictors", impact_full_C)
print_model_summary("C. EU-LIBYA / MoU (Feb 2017 onwards) — death count, full predictors", impact_full_C_deaths)

# Report inclusion probabilities for curated Model C
cat(rep("=", 71), "\n", sep = "")
cat("TOP PREDICTORS — CURATED MODEL C (mortality rate)\n")
cat(rep("=", 71), "\n\n", sep = "")

inc_probs <- colMeans(impact_cur_C$model$bsts.model$coefficients != 0)
top20 <- sort(inc_probs, decreasing = TRUE)[1:20]
for (i in seq_along(top20)) {
  cat(sprintf("  %2d. %-55s %.4f\n", i, names(top20)[i], top20[i]))
}

cat("\n")
cat(rep("=", 71), "\n", sep = "")
cat("INTERPRETATION GUIDE:\n")
cat("  - For Models A and B: a significant NEGATIVE effect would mean SAR\n")
cat("    reduced mortality (made the route safer).\n")
cat("  - For Model C: a significant POSITIVE effect would mean the MoU\n")
cat("    increased mortality (made the route more dangerous).\n")
cat(rep("=", 71), "\n\n", sep = "")

cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 71), "\n", sep = "")
