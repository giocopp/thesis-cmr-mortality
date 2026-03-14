"""
Compute monthly average moon illumination fraction at Tripoli (32.9N, 13.1E).

Moon illumination affects nighttime visibility for both migrants navigating
and potential rescuers detecting boats. This is a genuinely exogenous predictor
of per-crossing mortality (detection/rescue probability mechanism).

Output: CSV with columns [date, moon_illumination_frac]
"""

import ephem
from datetime import datetime, timedelta
import csv
import os

# Tripoli, Libya — main departure corridor
OBSERVER_LAT = "32.9"
OBSERVER_LON = "13.1"

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "moon_illumination.csv")

# Date range matching the original dataset
START_YEAR, START_MONTH = 2009, 1
END_YEAR, END_MONTH = 2021, 10


def monthly_moon_illumination(year, month):
    """Compute average nightly moon illumination fraction for a given month."""
    observer = ephem.Observer()
    observer.lat = OBSERVER_LAT
    observer.lon = OBSERVER_LON
    observer.elevation = 0

    illuminations = []
    d = datetime(year, month, 1)

    # Determine last day of month
    if month == 12:
        next_month = datetime(year + 1, 1, 1)
    else:
        next_month = datetime(year, month + 1, 1)

    while d < next_month:
        # Compute at midnight local time (roughly 22:00-23:00 UTC for Tripoli)
        observer.date = ephem.Date(d)
        moon = ephem.Moon(observer)
        # moon.phase is 0-100, convert to 0-1 fraction
        illuminations.append(moon.phase / 100.0)
        d += timedelta(days=1)

    return sum(illuminations) / len(illuminations)


def main():
    print("Computing monthly moon illumination at Tripoli...")
    results = []

    year, month = START_YEAR, START_MONTH
    while (year, month) <= (END_YEAR, END_MONTH):
        illum = monthly_moon_illumination(year, month)
        date_str = f"{year}-{month:02d}-01"
        results.append((date_str, round(illum, 6)))
        print(f"  {date_str}: {illum:.4f}")

        # Next month
        if month == 12:
            year += 1
            month = 1
        else:
            month += 1

    # Write CSV
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["date", "moon_illumination_frac"])
        writer.writerows(results)

    print(f"\nSaved {len(results)} months to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
