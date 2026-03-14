library(lubridate)
library(CausalImpact)
library(janitor)
library(countrycode)
library(tidyverse)
library(scales)
if (requireNamespace("imputeTS", quietly = TRUE)) {
  library(imputeTS)
} else {
  message("Package imputeTS not installed; continuing without it.")
}

# allow running from project root or from Replication/
if (file.exists("Replication/replicate_main.R") && file.exists("Original code and data/df.RDS")) {
  setwd("Replication")
}

# output folders for adapted run
dir.create("replicated_results/models", recursive = TRUE, showWarnings = FALSE)
dir.create("replicated_results/data", recursive = TRUE, showWarnings = FALSE)
dir.create("replicated_results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("replicated_results/figures/pdf", recursive = TRUE, showWarnings = FALSE)
dir.create("replicated_results/figures/png", recursive = TRUE, showWarnings = FALSE)

#read file
df <- readRDS(file="../Original code and data/df.RDS")
if (file.exists("../Original code and data/df_y_rec.RDS")) {
  df_y_rec <- readRDS(file="../Original code and data/df_y_rec.RDS")
  df <- left_join(df, df_y_rec, by="date")
}
if (!"y_rec" %in% names(df)) {
  df$y_rec <- 0
}

#Identifying relevant/important predictors to include in the model
df <- df %>% 
  mutate(LCG_pushbacks_count = as.numeric(ifelse(is.na(LCG_pushbacks_count), 0, LCG_pushbacks_count)),
         TCG_pushbacks_count = as.numeric(ifelse(is.na(TCG_pushbacks_count), 0, TCG_pushbacks_count)),
         dead_and_missing_Central_Mediterranean = as.numeric(ifelse(is.na(dead_and_missing_Central_Mediterranean), 0, dead_and_missing_Central_Mediterranean)),
         #crossings_CMR_r = arrivals_CMR + LCG_pushbacks_count + TCG_pushbacks_count + y_rec,
         crossings_CMR = arrivals_CMR + LCG_pushbacks_count + TCG_pushbacks_count + dead_and_missing_Central_Mediterranean,
         mortality_rate = (dead_and_missing_Central_Mediterranean/crossings_CMR)*1000)

#modeling
#Time series decomposition
#attempted crossings
y <- ts(df$crossings_CMR, start = c(2009,1), frequency=12)

stl_decomp1 <- stl(y, t.window = 24, s.window=12)
plot(stl_decomp1)
stl_decomp2 <- decompose(y, type="additive")
plot(stl_decomp2)

#install.packages("forecast")
library(forecast)
y_seasadj <- seasadj(stl_decomp2)
autoplot(cbind(y, y_seasadj))
y_star1 <- y - stl_decomp1$time.series[,1] - stl_decomp1$time.series[,2]
autoplot(cbind(y, y_star1))
y_star2 <- y - stl_decomp2$seasonal - stl_decomp2$trend
autoplot(cbind(y, y_star2))

#Structural changes in the series by means of chow test

#install.packages("strucchange")
library(strucchange)
bp_ts <- breakpoints(log(y) ~ 1, h=12)
summary(bp_ts)
ci_ts <- confint(bp_ts)
plot(log(y), 
     main = "Structural breakpoints in montly number of attempted crossings (log)",
     ylab = "Log of number of attempted crossings",
     xlab = "Date")
lines(bp_ts, col="darkred")
try(lines(ci_ts), silent = TRUE)

# no breakpoints are found for the mortality rate
library(strucchange)
bp_ts <- breakpoints(df$mortality_rate ~ 1, h=12)
summary(bp_ts)
ci_ts <- confint(bp_ts)
plot(df$y_rec,type="l")
lines(df$dead_and_missing_Central_Mediterranean, col="red")
try(lines(ci_ts), silent = TRUE)

#deaths
y_d1 <- ts(df$dead_and_missing_Central_Mediterranean, start = c(2009,1), frequency=12)
y_d2 <- ts(df$y_rec, start = c(2009,1), frequency=12)

stl_decomp1_d <- stl(y_d1, t.window = 24, s.window=12)
plot(stl_decomp1_d)
stl_decomp2_d <- decompose(y_d1, type="additive")
plot(stl_decomp2_d)

## seasonality decomposition of death rate
y_seasadj_d <- seasadj(stl_decomp2_d)
autoplot(cbind(y_d1, y_seasadj_d))
y_star1_d <- y_d1 - stl_decomp1_d$time.series[,1] - stl_decomp1_d$time.series[,2]
autoplot(cbind(y_d1, y_star1_d))
y_star2_d <- y_d1 - stl_decomp2_d$seasonal - stl_decomp2_d$trend
autoplot(cbind(y_d1, y_star2_d))

#Structural changes in the series by means of chow test on the death rate
bp_ts_d1 <- breakpoints(log(y_d1 + 1) ~ 1, h=12)
summary(bp_ts_d1)
ci_ts_d1 <- confint(bp_ts_d1)
plot(log(y_d1), main="Original data")
lines(bp_ts_d1)
lines(ci_ts_d1)

bp_ts_d2 <- breakpoints(log(y_d2 + 1) ~ 1, h=12)
summary(bp_ts_d2)
if (any(y_d2 != 0, na.rm = TRUE)) {
  ci_ts_d2 <- confint(bp_ts_d2)
  plot(log(y_d2), main="Reconstructed")
  lines(bp_ts_d2)
  try(lines(ci_ts_d2), silent = TRUE)
} else {
  message("Skipping reconstructed-deaths breakpoint plot: y_rec is all zeros in this data bundle.")
}

# no breakpoints are found for the mortality rate
library(strucchange)
bp_ts <- breakpoints(df$mortality_rate ~ 1, h=12)
summary(bp_ts)
ci_ts <- confint(bp_ts)
plot(df$mortality_rate)
lines(bp_ts)
try(lines(ci_ts), silent = TRUE)

# Estimating a first model

#dropping some of the variables before model selection
df_reduced <- df %>% 
  filter(date>="2011-02-01" & date<"2021-10-01") %>% 
  dplyr::select(-c(contains("lag_24",ignore.case = TRUE),
                   contains("lag_23",ignore.case = TRUE),
                   contains("lag_22",ignore.case = TRUE),
                   contains("lag_21",ignore.case = TRUE),
                   contains("lag_20",ignore.case = TRUE),
                   contains("lag_19",ignore.case = TRUE),
                   contains("lag_18",ignore.case = TRUE),
                   contains("lag_17",ignore.case = TRUE),
                   contains("lag_16",ignore.case = TRUE),
                   contains("lag_15",ignore.case = TRUE),
                   contains("lag_14",ignore.case = TRUE),
                   contains("lag_13",ignore.case = TRUE),
                   contains("lag_12",ignore.case = TRUE),
                   contains("lag_11",ignore.case = TRUE),
                   contains("lag_10",ignore.case = TRUE),
                   contains("lag_09",ignore.case = TRUE),
                   contains("lag_08",ignore.case = TRUE),
                   contains("lag_07",ignore.case = TRUE),
                   starts_with("airflow_Palestinian.Territories"),
                   starts_with("asylum"),
                   "arrivals_BSR","arrivals_CMR","arrivals_CRAG","arrivals_EBR","arrivals_EMR","arrivals_OR","arrivals_WAR","arrivals_WBR","arrivals_WMR",
                   "dead_and_missing_Eastern_Mediterranean","dead_and_missing_Central_Mediterranean","dead_and_missing_Western_Mediterranean",
                   "sd_lat__Eastern_Mediterranean","sd_lat__Central_Mediterranean","sd_lat__Western_Mediterranean","sd_lon__Eastern_Mediterranean","sd_lon__Central_Mediterranean","sd_lon__Western_Mediterranean",
                   "frac_index_2_to_10_deads_Eastern_Mediterranean","frac_index_2_to_10_deads_Central_Mediterranean","frac_index_2_to_10_deads_Western_Mediterranean",
                   "frac_index_less_than_1_dead_Eastern_Mediterranean","frac_index_less_than_1_dead_Central_Mediterranean","frac_index_less_than_1_dead_Western_Mediterranean",
                   "frac_index_more_than_10_deads_Eastern_Mediterranean","frac_index_more_than_10_deads_Central_Mediterranean","frac_index_more_than_10_deads_Western_Mediterranean",
                   "LCG_pushbacks_count","TCG_pushbacks_count",
                   "y_rec","mortality_rate"))

df_reduced_A <- data.frame(date=df_reduced$date,
                           crossings_CMR = log(df_reduced$crossings_CMR),
                           dplyr::select(df_reduced, -c(date,crossings_CMR)))

df_reduced_A <- df_reduced_A %>% 
  mutate(month = month(date),
         semester = semester(date),
         quarter = quarter(date))

df_min_A <- df_reduced_A %>% na.omit()

#Spike and slap prior to select variables
model <- logit.spike((crossings_CMR) ~ .,
                     data = dplyr::select(df_min_A, -c(date,month,semester,quarter)),
                     niter = 1000,
                     nthreads = 7,
                     seed = 270488)
smry.model <- summary(model)
plot(model, inc = 0.0001)
d_coefs <- data.frame(smry.model$coefficients)
d_coefs$names <- row.names(d_coefs)
d_coefs_best <- d_coefs %>% dplyr::filter(inc.prob >= 0.0001) %>% 
  dplyr::filter(names != "(Intercept)")
invisible(row.names(d_coefs %>% dplyr::filter(inc.prob >= 0.0001)))

df_min_A_ <- dplyr::select(df_min_A, c(date,crossings_CMR,
                                       starts_with("temperature"), starts_with("precipitation"), starts_with("daysstorm"),
                                       starts_with("ucdp_deaths_Syria"),
                                       airflow_LTU_lag_2,
                                       airflow_LUX_lag_2,
                                       PWOOLF_lag_05,
                                       airflow_NLD_lag_1,
                                       LBP_to_EURO_price_avg_lag_05,
                                       KES_to_EURO_price_avg_lag_04,
                                       airflow_QAT_lag_3,
                                       unem_GREECE_lag_04,
                                       num_expvio_Somalia,
                                       num_riots_Malawi_lag_06,
                                       PFSHMEAL,
                                       airflow_EST_lag_2,
                                       num_expvio_Cameroon_lag_04,
                                       disas_count_Botswana_lag_05,
                                       PWOOLF_lag_01,
                                       airflow_ITA_lag_2,
                                       airflow_MLI_lag_3,
                                       PZINC_lag_04,
                                       KHR_to_EURO_price_avg_lag_05,
                                       airflow_FIN_lag_4,
                                       num_battles_Togo_lag_04,
                                       airflow_GIN_lag_2,
                                       KES_to_EURO_price_avg_lag_02,
                                       NAD_to_EURO_price_avg_lag_01,
                                       PBARL_lag_03,
                                       airflow_ZMB_lag_6,
                                       airflow_GMB_lag_6,
                                       airflow_GAB_lag_6,
                                       airflow_ZWE_lag_2,
                                       airflow_SLE_lag_2,
                                       PCOFFOTM_lag_04,
                                       num_expvio_Egypt_lag_05,
                                       num_expvio_Burundi_lag_02,
                                       num_protest_Central.African.Republic_lag_01,
                                       TOP_to_EURO_price_avg_lag_06,
                                       BWP_to_EURO_price_avg_lag_06,
                                       TZS_to_EURO_price_avg_lag_03,
                                       airflow_BEN_lag_6,
                                       airflow_CMR_lag_5,
                                       airflow_GEO_lag_1))


'
starts_with("temperature"), starts_with("precipitation"), starts_with("daysstorm"),
                                       unem_ITALY_lag_06,
                                       PALUM_lag_05,
                                       airflow_EGY_lag_4,
                                       PWOOLC_lag_05,
                                       PALUM_lag_06,
                                       PNFUEL,
                                       PPOTASH_lag_02,
                                       PTOMATO_lag_05,
                                       PWOOLC_lag_04,
                                       IQD_to_EURO_price_avg_lag_03,
                                       airflow_TUN_lag_5,
                                       num_riots_Mauritania_lag_05,
                                       PEXGMETA_lag_05,
                                       num_expvio_Guinea.Bissau_lag_04,
                                       JOD_to_EURO_price_avg_lag_04,
                                       daysstorm_italy_lag_03,
                                       PFANDB,
                                       PNFUEL_lag_04,
                                       PAPPLE_lag_05,
                                       num_riots_Cameroon_lag_05,
                                       PWHEAMT_lag_04,
                                       airflow_SSD_lag_5,
                                       airflow_MRT_lag_4,
                                       disas_count_Viet.Nam_lag_05,
                                       PNGASEU_lag_03
'

df_reduced <- df_reduced %>% 
  mutate(crossings_CMR = log(df_reduced$crossings_CMR))

df_reduced_cc <- na.omit(df_reduced)
df_reduced_cc <- data.frame(date=df_reduced_cc$date,
                           crossings_CMR = df_reduced_cc$crossings_CMR,
                           dplyr::select(df_reduced_cc, -c(date,crossings_CMR))) %>% 
  mutate(month = month(date),
         semester = semester(date),
         quarter = quarter(date))

#Mare Nostrum Oct 18, 2013 - Oct 31, 2014
pre.period_mare_nostrum <- ymd(min(df_min_A$date), "2013-09-01")
post.period_mare_nostrum <- ymd("2013-10-01", max(df_min_A$date))
#SAR by NGOs
pre.period_sar_ngos <- ymd(min(df_min_A$date), "2014-10-01")
post.period_sar_ngos <- ymd("2014-11-01", max(df_min_A$date))
#EU Libya cooperation
pre.period_sarlibya <- ymd(min(df_min_A$date), "2017-01-01")
post.period_sar_libya <- ymd("2017-02-01", max(df_min_A$date))

set.seed(270488)

#Models
impact_marenostrum <- CausalImpact(df_min_A,
                                   pre.period = pre.period_mare_nostrum,
                                   post.period = post.period_mare_nostrum,
                                   alpha = 0.05,
                                   model.args = list(dynamic.regression=F, standardize.data=T, max.flips=100, niter=10000))
plot(impact_marenostrum, "original")
plot(impact_marenostrum$model$bsts.model, "coefficients", inc = 0.01)

impact_marenostrum$summary$AbsEffect
impact_marenostrum$summary$RelEffect

impact_sarngos <- CausalImpact(df_min_A,
                               pre.period = pre.period_sar_ngos,
                               post.period = post.period_sar_ngos,
                               alpha = 0.05,
                               model.args = list(dynamic.regression=F, standardize.data=T, max.flips=100, niter=10000))
plot(impact_sarngos, "original")
plot(impact_sarngos$model$bsts.model, "coefficients", inc=0.01)

impact_sarlibya <- CausalImpact(df_min_A,
                                pre.period = pre.period_sarlibya, post.period = post.period_sar_libya,
                                alpha = 0.05,
                                model.args = list(dynamic.regression=F, standardize.data=T, max.flips=100, niter=10000))
plot(impact_sarlibya)
plot(impact_sarlibya, "original")
plot(impact_sarlibya$model$bsts.model, "coefficients", inc=0.01)

saveRDS(impact_marenostrum, file="replicated_results/models/impact_marenostrum.RDS")
saveRDS(impact_sarngos, file="replicated_results/models/impact_sarngos.RDS")
saveRDS(impact_sarlibya, file="replicated_results/models/impact_sarlibya.RDS")

#Storing model results
model_df_results_A <- data.frame(Original = impact_marenostrum$series$response,
                                 Prediction = impact_marenostrum$series$point.pred,
                                 Prediction_lower = impact_marenostrum$series$point.pred.lower,
                                 Prediction_upper = impact_marenostrum$series$point.pred.upper,
                                 Pointwise_effect = impact_marenostrum$series$point.effect,      
                                 Pointwise_effect_lower = impact_marenostrum$series$point.effect.lower,
                                 Pointwise_effect_upper = impact_marenostrum$series$point.effect.upper,
                                 Cumulative_effect = impact_marenostrum$series$cum.effect,
                                 Cumulative_effect_lower = impact_marenostrum$series$cum.effect.lower,
                                 Cumulative_effect_upper = impact_marenostrum$series$cum.effect.upper)
model_df_results_A$date <- ymd(row.names(model_df_results_A))
model_df_results_A <- model_df_results_A %>% rename_with(tolower)

model_df_results_C <- data.frame(Original = impact_sarlibya$series$response,
                                 Prediction = impact_sarlibya$series$point.pred,
                                 Prediction_lower = impact_sarlibya$series$point.pred.lower,
                                 Prediction_upper = impact_sarlibya$series$point.pred.upper,
                                 Pointwise_effect = impact_sarlibya$series$point.effect,      
                                 Pointwise_effect_lower = impact_sarlibya$series$point.effect.lower,
                                 Pointwise_effect_upper = impact_sarlibya$series$point.effect.upper,
                                 Cumulative_effect = impact_sarlibya$series$cum.effect,
                                 Cumulative_effect_lower = impact_sarlibya$series$cum.effect.lower,
                                 Cumulative_effect_upper = impact_sarlibya$series$cum.effect.upper)
model_df_results_C$date <- ymd(row.names(model_df_results_C))
model_df_results_C <- model_df_results_C %>% rename_with(tolower)


model_df_results_B <- data.frame(Original = impact_sarngos$series$response,
                                 Prediction = impact_sarngos$series$point.pred,
                                 Prediction_lower = impact_sarngos$series$point.pred.lower,
                                 Prediction_upper = impact_sarngos$series$point.pred.upper,
                                 Pointwise_effect = impact_sarngos$series$point.effect,      
                                 Pointwise_effect_lower = impact_sarngos$series$point.effect.lower,
                                 Pointwise_effect_upper = impact_sarngos$series$point.effect.upper,
                                 Cumulative_effect = impact_sarngos$series$cum.effect,
                                 Cumulative_effect_lower = impact_sarngos$series$cum.effect.lower,
                                 Cumulative_effect_upper = impact_sarngos$series$cum.effect.upper)
model_df_results_B$date <- ymd(row.names(model_df_results_B))
model_df_results_B <- model_df_results_B %>% rename_with(tolower)

saveRDS(model_df_results_A, file="replicated_results/data/model_df_results_A_all.RDS")
saveRDS(model_df_results_B, file="replicated_results/data/model_df_results_B_all.RDS")
saveRDS(model_df_results_C, file="replicated_results/data/model_df_results_C_all.RDS")

#####################################################################
#####################################################################
##### DEATHS

y_d <- ts(df$dead_and_missing_Central_Mediterranean, start = c(2009,1), frequency=12)
library(strucchange)
bp_ts_d <- breakpoints(log(y_d+0.1) ~ 1, h=12)
summary(bp_ts_d)
ci_ts_d <- confint(bp_ts_d)
plot(log(y_d), 
     main = "Structural breakpoints in montly number of dead and missing in CMR",
     ylab = "Number of dead and missing",
     xlab = "Date")
lines(bp_ts_d, col="darkred")
lines(ci_ts_d)


df_reduced_d <- df %>% 
  filter(date>="2011-02-01" & date<"2021-10-01") %>% 
  dplyr::select(-c(contains("lag_24",ignore.case = TRUE),
                   contains("lag_23",ignore.case = TRUE),
                   contains("lag_22",ignore.case = TRUE),
                   contains("lag_21",ignore.case = TRUE),
                   contains("lag_20",ignore.case = TRUE),
                   contains("lag_19",ignore.case = TRUE),
                   contains("lag_18",ignore.case = TRUE),
                   contains("lag_17",ignore.case = TRUE),
                   contains("lag_16",ignore.case = TRUE),
                   contains("lag_15",ignore.case = TRUE),
                   contains("lag_14",ignore.case = TRUE),
                   contains("lag_13",ignore.case = TRUE),
                   contains("lag_12",ignore.case = TRUE),
                   contains("lag_11",ignore.case = TRUE),
                   contains("lag_10",ignore.case = TRUE),
                   contains("lag_09",ignore.case = TRUE),
                   contains("lag_08",ignore.case = TRUE),
                   contains("lag_07",ignore.case = TRUE),
                   starts_with("airflow_Palestinian.Territories"),
                   starts_with("asylum"),
                   "mortality_rate","y_rec",
                   "arrivals_BSR","arrivals_CRAG","arrivals_EBR","arrivals_EMR","arrivals_OR","arrivals_WAR","arrivals_WBR","arrivals_WMR",
                   "dead_and_missing_Eastern_Mediterranean","dead_and_missing_Western_Mediterranean",
                   "sd_lat__Eastern_Mediterranean","sd_lat__Central_Mediterranean","sd_lat__Western_Mediterranean","sd_lon__Eastern_Mediterranean","sd_lon__Central_Mediterranean","sd_lon__Western_Mediterranean",
                   "frac_index_2_to_10_deads_Eastern_Mediterranean","frac_index_2_to_10_deads_Central_Mediterranean","frac_index_2_to_10_deads_Western_Mediterranean",
                   "frac_index_less_than_1_dead_Eastern_Mediterranean","frac_index_less_than_1_dead_Central_Mediterranean","frac_index_less_than_1_dead_Western_Mediterranean",
                   "frac_index_more_than_10_deads_Eastern_Mediterranean","frac_index_more_than_10_deads_Central_Mediterranean","frac_index_more_than_10_deads_Western_Mediterranean"))

df_reduced_B <- data.frame(date=df_reduced_d$date,
                           dead_and_missing_Central_Mediterranean = log(df_reduced_d$dead_and_missing_Central_Mediterranean + 1),
                           dplyr::select(df_reduced_d, -c(date,dead_and_missing_Central_Mediterranean)))

df_reduced_B <- df_reduced_B %>% 
  mutate(month = month(date),
         semester = semester(date),
         quarter = quarter(date))

df_min_B <- df_reduced_B %>% na.omit()

# Deaths
model_d <- logit.spike((dead_and_missing_Central_Mediterranean) ~ .,
                     data = dplyr::select(df_min_B, -c(date,month,semester,quarter)),
                     niter = 1000,
                     nthreads = 7,
                     seed = 270488)
smry.model_d <- summary(model_d)
plot(model_d, inc = 0.1)
d_coefs <- data.frame(smry.model_d$coefficients)
d_coefs$names <- row.names(d_coefs)
d_coefs_best <- d_coefs %>% dplyr::filter(inc.prob >= 0.0001) %>% 
  dplyr::filter(names != "(Intercept)")
invisible(row.names(d_coefs %>% dplyr::filter(inc.prob >= 0.0001)))


impact_marenostrum_d <- CausalImpact(df_min_B,
                                   pre.period = pre.period_mare_nostrum,
                                   post.period = post.period_mare_nostrum,
                                   alpha = 0.05,
                                   model.args = list(dynamic.regression=F, standardize.data=T, max.flips=100, niter=10000))
plot(impact_marenostrum_d)
plot(impact_marenostrum_d$model$bsts.model, "coefficients", inc=0.1)

impact_sarngos_d <- CausalImpact(df_min_B,
                               pre.period = pre.period_sar_ngos,
                               post.period = post.period_sar_ngos,
                               alpha = 0.05,
                               model.args = list(dynamic.regression=F, standardize.data=T, max.flips=100, niter=10000))
plot(impact_sarngos_d)
plot(impact_sarngos$model$bsts.model, "coefficients", inc=0.001)

impact_sarlibya_d <- CausalImpact(df_min_B,
                                pre.period = pre.period_sarlibya, post.period = post.period_sar_libya,
                                alpha = 0.05,
                                model.args = list(dynamic.regression=F, standardize.data=T, max.flips=100, niter=10000))
plot(impact_sarngos_d)

#df_min_B <- df_reduced_B %>% na.omit()

#Spike and slap prior to select variables
model <- logit.spike((dead_and_missing_Central_Mediterranean) ~ .,
                     data = dplyr::select(df_reduced_B, -c(date,month,semester,quarter)),
                     niter = 1000,
                     nthreads = 7,
                     seed = 270488)

###############################################################################
# PLOTS - Figure 1, Figure 2, and Figure 3
###############################################################################

library(gridExtra)
library(grid)

###############################################################################
# FIGURE 1: Time-series of crossing attempts, mortality rate, and
#            main intervention periods, 2009-2021
###############################################################################

# Mortality rate per 100 attempted crossings (as in paper)
df <- df %>% mutate(mortality_rate_100 = (dead_and_missing_Central_Mediterranean / crossings_CMR) * 100)

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
  annotate("label", x = mn_start + lubridate::days(180), y = 29000,
           label = "EU\nMare Nostrum", color = "darkgreen", size = 2.8,
           fontface = "bold", fill = NA) +
  annotate("label", x = ngo_start + lubridate::days(400), y = 29000,
           label = "NGOs\nsearch-and-rescue", color = "darkorange3", size = 2.8,
           fontface = "bold", fill = NA) +
  annotate("label", x = eulcg_date + lubridate::days(40), y = 29000,
           label = "EU-LCG\nDeal", color = "red3", size = 2.8,
           fontface = "bold", fill = NA) +
  annotate("label", x = exp_start + lubridate::days(550), y = 29000,
           label = "Expansion\nLCG SAR-Zone", color = "brown", size = 2.8,
           fontface = "bold", fill = NA) +
  geom_line(linewidth = 0.4) +
  scale_x_date(breaks = seq(as.Date("2010-01-01"), as.Date("2020-01-01"), by = "5 years"),
               date_labels = "%Y") +
  scale_y_continuous(labels = scales::comma, breaks = seq(0, 30000, 10000)) +
  coord_cartesian(ylim = c(0, 31000)) +
  labs(title = "A. Attempted Crossings (arrivals, deaths, and pushbacks),\n    intervention periods, and structural breakpoints",
       x = "", y = "Attempted crossings") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        plot.margin = margin(5, 10, 5, 5))

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
  labs(title = "B. Mortality rate per 100 attempted crossings",
       x = "Date",
       y = expression(paste("Rate per 100 crossings", D[t]/Y[t], "*100"))) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        plot.margin = margin(5, 10, 5, 5))

figure1 <- arrangeGrob(
  fig1A, fig1B, ncol = 1, heights = c(1, 1),
  bottom = textGrob("Note: own calculations.", x = 0.02, hjust = 0,
                    gp = gpar(fontsize = 9, fontface = "italic"))
)

ggsave("replicated_results/figures/pdf/Figure1.pdf", figure1, width = 10, height = 10)
ggsave("replicated_results/figures/png/Figure1.png", figure1, width = 10, height = 10, dpi = 300)

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
    scale_x_date(breaks = x_breaks, date_labels = "%Y",
                 limits = c(date_min, date_max),
                 expand = expansion(mult = c(0.01, 0.01))) +
    labs(title = title,
         x = if (show_xlabel) "Date" else "",
         y = ylab) +
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
    scale_x_date(breaks = x_breaks, date_labels = "%Y",
                 limits = c(date_min, date_max),
                 expand = expansion(mult = c(0.01, 0.01))) +
    labs(title = title,
         x = if (show_xlabel) "Date" else "",
         y = ylab) +
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

#original CausalImpact plots
summary(impact_sarlibya)
plot(impact_sarlibya)
plot(impact_sarlibya$model$bsts.model, "coefficients")
summary(impact_sarngos)
plot(impact_sarngos)
plot(impact_sarngos$model$bsts.model, "coefficients")
summary(impact_marenostrum)
plot(impact_marenostrum)
plot(impact_marenostrum$model$bsts.model, "coefficients")


# Custom model
df_min <- df %>% filter(date >= "2011-02-01" & date < "2021-10-01") %>% na.omit()
required_custom <- c(
  "arrivals_CMR",
  "ALG_passengers_lag12", "EGY_passengers_lag12", "LYB_passengers_lag12", "MOR_passengers_lag12", "TUN_passengers_lag12",
  "unem_euro_area_all_lag1", "unem_euro_area_all_lag3", "unem_euro_area_all_lag6", "unem_euro_area_all_lag12",
  "num_battles_Libya", "num_battles_Libya_lag1", "num_battles_Libya_lag3", "num_battles_Libya_lag6", "num_battles_Libya_lag12",
  "num_battles_Somalia", "num_battles_Somalia_lag1", "num_battles_Somalia_lag3", "num_battles_Somalia_lag6", "num_battles_Somalia_lag12",
  "syria_trend_google", "syria_trend_google_lag1", "syria_trend_google_lag3", "syria_trend_google_lag6", "syria_trend_google_lag12",
  "oilprice_open_avg", "oilprice_open_avg_lag1", "oilprice_open_avg_lag3", "oilprice_open_avg_lag6", "oilprice_open_avg_lag12",
  "LDY_to_EURO_price_avg", "LDY_to_EURO_price_avg_lag1", "LDY_to_EURO_price_avg_lag3", "LDY_to_EURO_price_avg_lag6", "LDY_to_EURO_price_avg_lag12",
  "temperature_malta", "precipitation_malta", "daysstorm_malta",
  "temperature_italy", "precipitation_italy", "daysstorm_italy"
)
if (all(required_custom %in% names(df_min))) {
  attach(df_min)
  post.period2 <- c(29, 90)
  y <- df_min$arrivals_CMR
  post.period.response <- y[post.period2[1] : post.period2[2]]
  y[post.period2[1] : post.period2[2]] <- NA

  ss <- AddLocalLevel(list(), y)
  ss <- AddLocalLinearTrend(ss, y)
  ss <- AddSeasonal(list(), y, nseasons = 12, season.duration = 1)
  ss <- AddSeasonal(ss, y, nseasons = 4, season.duration = 3)
  ss <- AddTrig(ss, y, period = 3, frequencies = 1)
  ss <- AddAutoAr(ss, y, lags = 1)
  ss <- AddSemilocalLinearTrend(ss, y)

  #ss <- AddAr(ss, lags = 3, sigma.prior = SdPrior(3.0, 1.0))

  bsts.model <- bsts(y ~ 
                       ALG_passengers_lag12 + EGY_passengers_lag12 + LYB_passengers_lag12 + MOR_passengers_lag12 + TUN_passengers_lag12 + 
                       unem_euro_area_all_lag1 + unem_euro_area_all_lag3 + unem_euro_area_all_lag6 + unem_euro_area_all_lag12 + 
                       num_battles_Libya + num_battles_Libya_lag1 + num_battles_Libya_lag3 + num_battles_Libya_lag6 + num_battles_Libya_lag12 + 
                       num_battles_Somalia + num_battles_Somalia_lag1 + num_battles_Somalia_lag3 + num_battles_Somalia_lag6 + num_battles_Somalia_lag12 + 
                       syria_trend_google + syria_trend_google_lag1 + syria_trend_google_lag3 + syria_trend_google_lag6 + syria_trend_google_lag12 + 
                       oilprice_open_avg + oilprice_open_avg_lag1 + oilprice_open_avg_lag3 + oilprice_open_avg_lag6 + oilprice_open_avg_lag12 + 
                       LDY_to_EURO_price_avg + LDY_to_EURO_price_avg_lag1 + LDY_to_EURO_price_avg_lag3 + LDY_to_EURO_price_avg_lag6 + LDY_to_EURO_price_avg_lag12 + 
                       temperature_malta +  precipitation_malta + daysstorm_malta +  
                       temperature_italy + precipitation_italy + daysstorm_italy,  
                     state.specification = ss, 
                     family = "gaussian",
                     niter = 10000,
                     seed=270488)
  bsts.model1 <- bsts(y ~ temperature_malta +  precipitation_malta + daysstorm_malta +  
                        temperature_italy + precipitation_italy + daysstorm_italy,  
                      state.specification = ss, 
                      family = "gaussian",
                      niter = 1000)
  bsts.model2 <- bsts(y ~ oilprice_open_avg + oilprice_open_avg_lag1 + oilprice_open_avg_lag3 + oilprice_open_avg_lag6 + oilprice_open_avg_lag12 + 
                        LDY_to_EURO_price_avg + LDY_to_EURO_price_avg_lag1 + LDY_to_EURO_price_avg_lag3 + LDY_to_EURO_price_avg_lag6 + LDY_to_EURO_price_avg_lag12 + 
                        unem_euro_area_all_lag1 + unem_euro_area_all_lag3 + unem_euro_area_all_lag6 + unem_euro_area_all_lag12 + 
                        temperature_malta +  precipitation_malta + daysstorm_malta +  
                        temperature_italy + precipitation_italy + daysstorm_italy,  
                      state.specification = ss, 
                      family = "gaussian",
                      niter = 1000)
  bsts.model3 <- bsts(y ~
                        num_battles_Libya + num_battles_Libya_lag1 + num_battles_Libya_lag3 + num_battles_Libya_lag6 + num_battles_Libya_lag12 + 
                        num_battles_Somalia + num_battles_Somalia_lag1 + num_battles_Somalia_lag3 + num_battles_Somalia_lag6 + num_battles_Somalia_lag12 + 
                        syria_trend_google + syria_trend_google_lag1 + syria_trend_google_lag3 + syria_trend_google_lag6 + syria_trend_google_lag12 + 
                        oilprice_open_avg + oilprice_open_avg_lag1 + oilprice_open_avg_lag3 + oilprice_open_avg_lag6 + oilprice_open_avg_lag12 + 
                        unem_euro_area_all_lag1 + unem_euro_area_all_lag3 + unem_euro_area_all_lag6 + unem_euro_area_all_lag12 + 
                        LDY_to_EURO_price_avg + LDY_to_EURO_price_avg_lag1 + LDY_to_EURO_price_avg_lag3 + LDY_to_EURO_price_avg_lag6 + LDY_to_EURO_price_avg_lag12 + 
                        temperature_malta +  precipitation_malta + daysstorm_malta +  
                        temperature_italy + precipitation_italy + daysstorm_italy,  
                      state.specification = ss, 
                      family = "gaussian",
                      niter = 1000)
  plot(bsts.model, "components")
  plot(bsts.model, "coef")
  plot(bsts.model)

  CompareBstsModels(list("Model 1" = bsts.model1,
                         "Model 2" = bsts.model2,
                         "Model 3" = bsts.model3,
                         "Model full" = bsts.model),
                    colors = c("green", "red", "blue", "black"))

  impact2 <- CausalImpact(bsts.model = bsts.model,
                          post.period.response = post.period.response)
  plot(impact2)
  summary(impact2)

  errors <- bsts.prediction.errors(bsts.model, burn = 1000)
  PlotDynamicDistribution(errors$in.sample)

  pred.bsts <- predict(bsts.model, newdata=df_min[29:90,])
  plot(pred.bsts$median)
  try(detach(df_min), silent = TRUE)
} else {
  message("Skipping custom model block: required legacy predictors are not available in current df.RDS")
}
