# ── Project constants ─────────────────────────────────────────────────────────
# Sourced by _helpers.R; every analytical script picks these up via the
# single source(_helpers.R) call at the top of each file.

# ── Dates ─────────────────────────────────────────────────────────────────────
# MOU_DATE drives the analytical `post_mou` indicator built in
# 03_build/03_daily_panel.R. The signing date (Feb 2017) and the analytical
# cutoff (Jul 2017, allowing for operational implementation lag) are
# deliberately distinct; descriptive plot annotations use the signing date.
MOU_DATE       <- as.Date("2017-07-01")
MOU_SIGN_DATE  <- as.Date("2017-02-02")

# Long-span analysis (31): political-regime period boundaries.
DATE_START     <- as.Date("2013-01-01")
IOM_START      <- as.Date("2014-01-01")
PERIOD_END_1   <- as.Date("2017-01-31")  # end of Post-Arab Spring (Mare Nostrum + NGO SAR)
PERIOD_END_2   <- as.Date("2019-12-31")  # end of MoU + Salvini
PERIOD_END_3   <- as.Date("2022-10-21")  # end of partial rollback (Meloni sworn 10-22)

PANEL_START    <- as.Date("2014-01-01")

# ── Geography ─────────────────────────────────────────────────────────────────
CMR_DEPARTURES         <- c("Libya", "Tunisia", "Algeria")
CMR_INCIDENT_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia")

# ── Estimation ────────────────────────────────────────────────────────────────
NW_LAG         <- 14   # Newey-West bandwidth for fixest::NW()
CROSSING_LAG   <- 14   # exogenous volume control: mean over [t-14, t-8]

# ── Plotting ──────────────────────────────────────────────────────────────────
LOESS_SPAN     <- 0.07
