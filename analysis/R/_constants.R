# ── Project constants ─────────────────────────────────────────────────────────
# Sourced by _helpers.R; every analytical script picks these up via the
# single source(_helpers.R) call at the top of each file.

# ── Dates ─────────────────────────────────────────────────────────────────────
# MOU_DATE drives the analytical `post_mou` indicator built in
# 03_build/03_daily_panel.R. Set to the Italy-Libya MoU signing date.
MOU_DATE       <- as.Date("2017-02-02")
MOU_SIGN_DATE  <- as.Date("2017-02-02")

# Long-span analysis: 3-period political regime design aligned with the
# policy timeline (04_descriptive/05_policy_timeline.R).
#   1. SAR + border control     IOM_START   .. PERIOD_END_1
#   2. MoU + NGO containment    PERIOD_END_1+1 .. PERIOD_END_2
#   3. Lamorgese partial rollback PERIOD_END_2+1 .. PERIOD_END_3
DATE_START     <- as.Date("2013-01-01")
IOM_START      <- as.Date("2014-01-01")
PERIOD_END_1   <- as.Date("2017-02-01")  # last day before Italy-Libya MoU (2017-02-02)
PERIOD_END_2   <- as.Date("2020-10-20")  # last day before Lamorgese reforms (2020-10-21)
PERIOD_END_3   <- as.Date("2023-01-01")  # last day before Piantedosi decree (2023-01-02)

# Sub-period boundaries used for the period-1 robustness split
# (Mare Nostrum 2013-10-18 to 2014-10-31; Triton starts 2014-11-01).
MARE_NOSTRUM_END <- as.Date("2014-10-31")
TRITON_START     <- as.Date("2014-11-01")

PANEL_START    <- as.Date("2014-01-01")

# ── Geography ─────────────────────────────────────────────────────────────────
CMR_DEPARTURES         <- c("Libya", "Tunisia", "Algeria")
CMR_INCIDENT_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia")

# ── Estimation ────────────────────────────────────────────────────────────────
NW_LAG         <- 14   # Newey-West bandwidth for fixest::NW()
CROSSING_LAG   <- 14   # exogenous volume control: mean over [t-14, t-8]

# ── Plotting ──────────────────────────────────────────────────────────────────
LOESS_SPAN     <- 0.07
