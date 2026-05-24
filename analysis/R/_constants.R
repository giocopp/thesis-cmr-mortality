# ── Project constants ────────────────────────────────────────────────────────

# ── Dates ────────────────────────────────────────────────────────────────────
MOU_DATE       <- as.Date("2017-02-02")
MOU_SIGN_DATE  <- as.Date("2017-02-02")

DATE_START     <- as.Date("2013-01-01")
IOM_START      <- as.Date("2014-01-01")
PERIOD_END_1   <- as.Date("2017-02-01")
PERIOD_END_2   <- as.Date("2020-10-20")
PERIOD_END_3   <- as.Date("2023-01-01")

MARE_NOSTRUM_END <- as.Date("2014-10-31")
TRITON_START     <- as.Date("2014-11-01")

PANEL_START    <- as.Date("2014-01-01")

# ── Geography ────────────────────────────────────────────────────────────────
CMR_DEPARTURES         <- c("Libya", "Tunisia", "Algeria")
CMR_INCIDENT_COUNTRIES <- c("Algeria", "Italy", "Libya", "Malta", "Tunisia")

# ── Estimation ───────────────────────────────────────────────────────────────
NW_LAG         <- 14
CROSSING_LAG   <- 14

# ── Plotting ─────────────────────────────────────────────────────────────────
LOESS_SPAN     <- 0.07
