#!/usr/bin/env python3
"""
download_era5_storm.py
======================
Download additional storm-intensity variables from ERA5 for event-level analysis.

Variables:
  1. i10fg — instantaneous 10m wind gust (atm grid, 0.25 deg)
  2. hmax  — maximum individual wave height (wave grid, 0.5 deg)

Same bounding box and period as download_era5_daily.py.
Output goes to the same directory with distinct filenames.

Prerequisites: same as download_era5_daily.py (cdsapi + CDS credentials).
"""

import cdsapi
from pathlib import Path

# ============================================================
# Configuration (matches download_era5_daily.py)
# ============================================================

AREA = [45, 4, 29, 26]  # [N, W, S, E]
YEARS = list(range(2014, 2026))
MONTHS = [f"{m:02d}" for m in range(1, 13)]
DAYS = [f"{d:02d}" for d in range(1, 32)]

BASE_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = BASE_DIR / "data" / "raw" / "era5"

DATASET = "reanalysis-era5-single-levels"


# ============================================================
# Download
# ============================================================

def download_year(client, year, output_dir):
    """Download one year of storm-intensity ERA5 data."""

    base_request = {
        "product_type": ["reanalysis"],
        "year": [str(year)],
        "month": MONTHS,
        "day": DAYS,
        "time": ["12:00"],
        "area": AREA,
        "data_format": "netcdf",
        "download_format": "unarchived",
    }

    files = {
        "atm_gust": (
            ["instantaneous_10m_wind_gust"],
            output_dir / f"era5_daily_cmr_atm_gust_{year}.nc",
        ),
        "wave_hmax": (
            ["maximum_individual_wave_height"],
            output_dir / f"era5_daily_cmr_wave_hmax_{year}.nc",
        ),
    }

    for label, (variables, output_file) in files.items():
        if output_file.exists():
            print(f"  {year}/{label}: already exists, skipping")
            continue

        print(f"  {year}/{label}: requesting...")
        request = {**base_request, "variable": variables}
        client.retrieve(DATASET, request, str(output_file))
        print(f"  {year}/{label}: saved to {output_file.name}")


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 60)
    print("ERA5 Storm Variables Download")
    print("  i10fg (wind gust, 0.25 deg)")
    print("  hmax  (max wave height, 0.5 deg)")
    print("=" * 60)

    print(f"\nArea: N={AREA[0]}, W={AREA[1]}, S={AREA[2]}, E={AREA[3]}")
    print(f"Years: {YEARS[0]}-{YEARS[-1]}")
    print(f"Output: {OUTPUT_DIR}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    client = cdsapi.Client()

    failed = []
    for year in YEARS:
        try:
            download_year(client, year, OUTPUT_DIR)
        except Exception as e:
            print(f"  {year}: FAILED - {e}")
            failed.append(year)

    print("\n" + "=" * 60)
    if failed:
        print(f"DONE with {len(failed)} failures: {failed}")
        print("Re-run the script to retry failed years.")
    else:
        print("DONE - all years downloaded successfully")
    print("=" * 60)

    for f in sorted(OUTPUT_DIR.glob("era5_daily_cmr_*_*.nc")):
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
