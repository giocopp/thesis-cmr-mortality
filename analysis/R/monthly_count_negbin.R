# monthly_count_negbin.R
# ======================
# Monthly NegBin on death counts — the clean monthly complement to Model A.
#
# Why this model instead of the rate regression:
# - The rate regression (log(deaths/crossings + c)) has two problems:
#   (1) Kronmal/collider bias from conditioning on crossings
#       (dag_identification.qmd Section 5)
#   (2) Arbitrary constant c in log(0 + c) for zero-death months
#       (coefficient changes 2.4x between c=0.01 and c=1)
# - The NegBin on counts avoids both: no crossings denominator,
#   handles zeros natively through the count distribution
# - Captures the TOTAL effect (danger + volume) at monthly level
#   (estimands.qmd Section 3.4)
#
# Month-of-year FE are always included (seasonality).

library(MASS)
library(dplyr)
library(lubridate)

BASE_DIR <- file.path("replication", "rodriguez-sanchez", "extension")

# ============================================================
# 1. Load and prepare data
# ============================================================
df <- readRDS(file.path(BASE_DIR, "data", "df_extended.RDS"))

PRE_START <- as.Date("2011-02-01")
END_DATE  <- as.Date("2021-09-01")
MOU_DATE  <- as.Date("2017-02-01")

df <- df %>%
  filter(date >= PRE_START & date <= END_DATE) %>%
  mutate(
    deaths    = dead_and_missing_Central_Mediterranean,
    crossings = crossings_CMR,
    post_mou  = as.integer(date >= MOU_DATE),
    month_oy  = factor(month(date)),
    year_fac  = factor(year(date)),
    swh       = wave_height_central_med,
    wind      = wind_speed_central_med,
    swh_z     = scale(swh)[, 1],
    wind_z    = scale(wind)[, 1]
  )

cat("============================================================\n")
cat("MONTHLY COUNT MODEL: NEGATIVE BINOMIAL\n")
cat("============================================================\n\n")

cat("Analysis sample:", nrow(df), "months\n")
cat("  Pre-MoU:", sum(df$post_mou == 0), "months\n")
cat("  Post-MoU:", sum(df$post_mou == 1), "months\n")
cat("  Deaths: mean =", round(mean(df$deaths), 1),
    ", median =", median(df$deaths),
    ", max =", max(df$deaths),
    ", zeros =", sum(df$deaths == 0), "\n")
cat("  Variance/mean ratio:", round(var(df$deaths) / mean(df$deaths), 1),
    "(>>1 confirms overdispersion)\n\n")

# ============================================================
# 2. Helper: extract and display results
# ============================================================
nb_report <- function(model, label = "") {
  s <- summary(model)
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("Theta: %.3f (SE = %.3f)  |  AIC: %.1f  |  Log-lik: %.1f\n",
      model$theta, model$SE.theta, AIC(model), as.numeric(logLik(model))))

  coefs <- s$coefficients
  skip <- grepl("Intercept|month_oy|year_fac", rownames(coefs))
  coefs_show <- coefs[!skip, , drop = FALSE]

  cat(sprintf("\n%-35s %10s %10s %8s\n",
      "Coefficient", "Beta", "IRR", "p"))
  cat(paste(rep("-", 68), collapse = ""), "\n")
  for (i in seq_len(nrow(coefs_show))) {
    cat(sprintf("%-35s %+10.4f %10.4f %8.4f%s\n",
        rownames(coefs_show)[i],
        coefs_show[i, 1],
        exp(coefs_show[i, 1]),
        coefs_show[i, 4],
        ifelse(coefs_show[i, 4] < 0.01, " ***",
               ifelse(coefs_show[i, 4] < 0.05, " **",
                      ifelse(coefs_show[i, 4] < 0.1, " *", "")))))
  }
  invisible(model)
}


# ============================================================
# 3. NegBin specifications
# ============================================================

# --- Model 1: SWH x Post + month-of-year FE ---
m1 <- glm.nb(deaths ~ swh * post_mou + month_oy, data = df)
nb_report(m1, "Model 1: SWH x Post + month-of-year FE")

# --- Model 2: Wind x Post + month-of-year FE ---
m2 <- glm.nb(deaths ~ wind * post_mou + month_oy, data = df)
nb_report(m2, "Model 2: Wind x Post + month-of-year FE")

# --- Model 3: Joint SWH + Wind x Post ---
m3 <- glm.nb(deaths ~ swh * post_mou + wind * post_mou + month_oy, data = df)
nb_report(m3, "Model 3: SWH + Wind x Post (joint)")

# --- Model 4: Year FE + month-of-year FE (more demanding) ---
# Year FE absorb annual trends; Post main effect is mostly absorbed
# (only 2017 has both pre and post months).
# The interaction SWH x Post is still identified from within-year
# SWH variation across months.
m4 <- tryCatch(
  glm.nb(deaths ~ swh * post_mou + year_fac + month_oy, data = df),
  error = function(e) {
    cat("  Model 4 failed:", e$message, "\n")
    NULL
  }
)
if (!is.null(m4)) {
  nb_report(m4, "Model 4: SWH x Post + year FE + month-of-year FE")
}


# ============================================================
# 4. Overdispersion: NegBin vs Poisson
# ============================================================
cat("\n\n============================================================\n")
cat("OVERDISPERSION TEST: NegBin vs Poisson\n")
cat("============================================================\n\n")

p1 <- glm(deaths ~ swh * post_mou + month_oy, data = df, family = poisson)

ll_nb  <- as.numeric(logLik(m1))
ll_poi <- as.numeric(logLik(p1))
lr_stat <- 2 * (ll_nb - ll_poi)
# Boundary test: theta -> infinity under H0 (Poisson), one-sided
lr_p <- pchisq(lr_stat, df = 1, lower.tail = FALSE) / 2

cat(sprintf("LR test (H0: Poisson adequate): stat = %.3f, p = %.6f\n",
    lr_stat, lr_p))
cat(sprintf("NegBin theta: %.3f (lower = more overdispersion)\n", m1$theta))
cat(sprintf("AIC: Poisson = %.1f, NegBin = %.1f\n", AIC(p1), AIC(m1)))

# Poisson QMLE coefficients for comparison
cat("\nPoisson QMLE (SWH x Post):\n")
s_p <- summary(p1)$coefficients
int_p <- grep("swh:post_mou|post_mou:swh", rownames(s_p))
if (length(int_p) > 0) {
  cat(sprintf("  Beta = %+.4f, IRR = %.4f, p = %.4f\n",
      s_p[int_p, 1], exp(s_p[int_p, 1]), s_p[int_p, 4]))
  cat("  (Poisson SEs; use sandwich for QMLE robustness)\n")
}


# ============================================================
# 5. Crossings diagnostic: crossings ~ weather x post
# ============================================================
# Crossings as OUTCOME (not conditioning variable) is DAG-valid:
# {FE} is sufficient, no collider (dag_identification.qmd Section 12).
# Tests whether the MoU changed how weather affects crossing volume.

cat("\n\n============================================================\n")
cat("CROSSINGS DIAGNOSTIC (volume channel)\n")
cat("crossings ~ weather x post + month-of-year FE\n")
cat("============================================================\n")

mc1 <- glm.nb(crossings ~ swh * post_mou + month_oy, data = df)
nb_report(mc1, "Crossings ~ SWH x Post")

mc2 <- glm.nb(crossings ~ wind * post_mou + month_oy, data = df)
nb_report(mc2, "Crossings ~ Wind x Post")

cat("\n  Interpretation: if weather x post is null for crossings,\n")
cat("  then the death count interaction is primarily about DANGER,\n")
cat("  not volume (weather deterrence didn't change much post-MoU).\n")


# ============================================================
# 6. Comparison with rate regression
# ============================================================
cat("\n\n============================================================\n")
cat("COMPARISON: NegBin COUNTS vs LOG-RATE\n")
cat("============================================================\n\n")

df$log_rate_001 <- log(df$mortality_rate_100 + 0.01)
df$log_rate_1   <- log(df$mortality_rate_100 + 1)

m_rate1 <- lm(log_rate_001 ~ swh * post_mou + month_oy, data = df)
m_rate2 <- lm(log_rate_1   ~ swh * post_mou + month_oy, data = df)

get_int <- function(model) {
  s <- summary(model)$coefficients
  idx <- grep("swh:post_mou|post_mou:swh", rownames(s))
  if (length(idx) > 0) c(beta = s[idx, 1], p = s[idx, 4])
  else c(beta = NA, p = NA)
}

s_nb   <- summary(m1)$coefficients
idx_nb <- grep("swh:post_mou|post_mou:swh", rownames(s_nb))

cat(sprintf("%-35s %10s %8s\n", "Model", "Beta(SWH×Post)", "p"))
cat(paste(rep("-", 55), collapse = ""), "\n")
cat(sprintf("%-35s %+10.4f %8.4f\n", "NegBin (death counts)",
    s_nb[idx_nb, 1], s_nb[idx_nb, 4]))

r1 <- get_int(m_rate1)
cat(sprintf("%-35s %+10.4f %8.4f\n", "OLS log(rate + 0.01)",
    r1["beta"], r1["p"]))

r2 <- get_int(m_rate2)
cat(sprintf("%-35s %+10.4f %8.4f\n", "OLS log(rate + 1)",
    r2["beta"], r2["p"]))

cat("\n  The NegBin avoids:\n")
cat("  - Kronmal/collider bias (no crossings denominator)\n")
cat("  - Arbitrary log(Y+c) constant (coefficient changes 2.4x)\n")
cat("  It captures the total effect (danger + volume).\n")


# ============================================================
# 7. Summary table
# ============================================================
cat("\n\n============================================================\n")
cat("SUMMARY: ALL NegBin INTERACTION COEFFICIENTS\n")
cat("============================================================\n\n")

extract_int <- function(model, pattern) {
  s <- summary(model)$coefficients
  idx <- grep(pattern, rownames(s))
  if (length(idx) > 0) {
    list(beta = s[idx, 1], irr = exp(s[idx, 1]), p = s[idx, 4])
  } else {
    list(beta = NA, irr = NA, p = NA)
  }
}

cat(sprintf("%-40s %10s %8s %8s\n",
    "Specification", "Beta", "IRR", "p"))
cat(paste(rep("-", 70), collapse = ""), "\n")

specs <- list(
  list("SWH x Post (month-of-year FE)",        m1, "swh:post|post:swh"),
  list("Wind x Post (month-of-year FE)",        m2, "wind:post|post:wind"),
  list("SWH x Post (joint, month-of-year FE)",  m3, "swh:post|post:swh"),
  list("Wind x Post (joint, month-of-year FE)", m3, "wind:post|post:wind")
)
if (!is.null(m4)) {
  specs <- c(specs, list(
    list("SWH x Post (year + month-of-year FE)", m4, "swh:post|post:swh")
  ))
}

for (sp in specs) {
  r <- extract_int(sp[[2]], sp[[3]])
  cat(sprintf("%-40s %+10.4f %8.4f %8.4f%s\n",
      sp[[1]], r$beta, r$irr, r$p,
      ifelse(r$p < 0.01, " ***",
             ifelse(r$p < 0.05, " **",
                    ifelse(r$p < 0.1, " *", "")))))
}

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
