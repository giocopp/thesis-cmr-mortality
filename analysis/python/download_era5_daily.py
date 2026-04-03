#!/usr/bin/env python3
"""
download_era5_daily.py
======================
Download daily ERA5 reanalysis data for the Central Mediterranean Route.

Two separate requests per year (different native grids):
  1. Atmospheric (0.25 deg): u10, v10, SST, total precipitation
  2. Wave (0.5 deg): SWH, mean wave period

Bounding box: lat [29, 45], lon [4, 26]
  (CMR incidents after dropping Atlantic outliers, with 1-degree margin)

Period: 2014-01-01 to 2025-12-31

Output: data/raw/era5/era5_daily_cmr_atm_YYYY.nc
        data/raw/era5/era5_daily_cmr_wave_YYYY.nc

Prerequisites:
  pip install cdsapi
  Create ~/.cdsapirc with your CDS API credentials:
    url: https://cds.climate.copernicus.eu/api
    key: YOUR-PERSONAL-ACCESS-TOKEN
  Accept the ERA5 licence at:
    https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download
"""

import cdsapi
from pathlib import Path

# ============================================================
# Configuration
# ============================================================

AREA = [45, 4, 29, 26]  # [N, W, S, E]
YEARS = list(range(2008, 2026))
MONTHS = [f"{m:02d}" for m in range(1, 13)]
DAYS = [f"{d:02d}" for d in range(1, 32)]

BASE_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = BASE_DIR / "data" / "raw" / "era5"

DATASET = "reanalysis-era5-single-levels"

# Atmospheric variables (0.25 deg native grid)
ATM_VARIABLES = [
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "sea_surface_temperature",
    "total_precipitation",
]

# Wave variables (0.5 deg native grid) — MUST be separate request
WAVE_VARIABLES = [
    "significant_height_of_combined_wind_waves_and_swell",
    "mean_wave_period",
]


# ============================================================
# Download function
# ============================================================

def download_year(client, year, output_dir):
    """Download one year of daily ERA5 data (atm + wave separately)."""

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
        "atm": (ATM_VARIABLES, output_dir / f"era5_daily_cmr_atm_{year}.nc"),
        "wave": (WAVE_VARIABLES, output_dir / f"era5_daily_cmr_wave_{year}.nc"),
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
    print("ERA5 Daily Download for CMR Event-Level Analysis")
    print("=" * 60)

    print(f"\nAtmospheric variables ({len(ATM_VARIABLES)}):")
    for v in ATM_VARIABLES:
        print(f"  - {v}")
    print(f"\nWave variables ({len(WAVE_VARIABLES)}):")
    for v in WAVE_VARIABLES:
        print(f"  - {v}")

    print(f"\nArea: N={AREA[0]}, W={AREA[1]}, S={AREA[2]}, E={AREA[3]}")
    print(f"Years: {YEARS[0]}-{YEARS[-1]}")
    print(f"Time: 12:00 UTC (daily snapshot)")
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

    for f in sorted(OUTPUT_DIR.glob("era5_daily_cmr_*.nc")):
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
