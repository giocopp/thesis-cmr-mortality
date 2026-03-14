"""
Download ERA5 daily significant wave height for the Central Mediterranean.

Unlike the monthly means already in the dataset, daily data allows us to
compute extreme event statistics:
  - Monthly MAX wave height (captures lethal storm events)
  - Monthly SD of wave height (captures variability/instability)
  - Monthly count of days with SWH > 2m (dangerous conditions frequency)

These tail-event statistics are more relevant to mortality than monthly means,
which average away the extreme events that actually kill people.

Dataset: ERA5 single levels (reanalysis), daily
Variable: significant_height_of_combined_wind_waves_and_swell (swh)
Region: Central Mediterranean (31N-38N, 10E-20E)
Period: Jan 2009 - Oct 2021

Output: NetCDF file in Extension-2-new-data/data/
"""

import cdsapi
import os

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "era5_daily_waves_central_med.nc")

c = cdsapi.Client()

# Download daily mean wave data
# ERA5 wave variables are only available at certain times; we request 00:00
# and 12:00 to get a reasonable daily average
years = [str(y) for y in range(2009, 2022)]
months = [f"{m:02d}" for m in range(1, 13)]
days = [f"{d:02d}" for d in range(1, 32)]

print("Downloading ERA5 daily significant wave height...")
print(f"  Region: 31N-38N, 10E-20E")
print(f"  Period: 2009-01 to 2021-12")
print(f"  Output: {OUTPUT_FILE}")
print()

c.retrieve(
    "reanalysis-era5-single-levels",
    {
        "product_type": "reanalysis",
        "variable": "significant_height_of_combined_wind_waves_and_swell",
        "year": years,
        "month": months,
        "day": days,
        "time": ["00:00", "12:00"],
        "area": [38, 10, 31, 20],  # N, W, S, E
        "format": "netcdf",
    },
    OUTPUT_FILE,
)

print(f"\nSaved to {OUTPUT_FILE}")
