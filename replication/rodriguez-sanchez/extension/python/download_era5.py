"""
Download ERA5 Monthly Mean Data for Central Mediterranean Mortality Analysis
=============================================================================

Extension of Rodriguez Sanchez et al. (2023) — Mortality counterfactual model.

This script downloads ERA5 reanalysis monthly means for:
  1. Atmospheric variables over the Central Med crossing zone (0.25 deg grid)
  2. Ocean wave variables over the Central Med crossing zone (0.5 deg grid)
  3. Atmospheric variables over the North African departure coast (0.25 deg grid)

PREREQUISITES:
  1. Register (free) at https://cds.climate.copernicus.eu
  2. Accept the ERA5 Terms of Use on the dataset page:
     https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means
  3. Get your Personal Access Token from https://cds.climate.copernicus.eu/profile
  4. Create ~/.cdsapirc with:
       url: https://cds.climate.copernicus.eu/api
       key: YOUR-PERSONAL-ACCESS-TOKEN
  5. Install the CDS API client:
       pip install "cdsapi>=0.7.7"

USAGE:
  python download_era5.py

  Downloads are saved to ../data/era5/
  Each request takes approximately 10-30 minutes depending on CDS queue.
"""

import cdsapi
import os

# Output directory
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "era5")
os.makedirs(OUT_DIR, exist_ok=True)

# Study period: Jan 2009 – Oct 2021 (matching Rodriguez Sanchez et al.)
YEARS = [str(y) for y in range(2009, 2022)]
MONTHS = [f"{m:02d}" for m in range(1, 13)]

# Dataset name
DATASET = "reanalysis-era5-single-levels-monthly-means"

# Geographic bounding boxes [North, West, South, East]
CENTRAL_MED = [38, 10, 31, 20]       # Central Med crossing zone
NORTH_AFRICA_COAST = [34, 8, 30, 25]  # Libya/Tunisia departure coast

client = cdsapi.Client()


# =============================================================================
# REQUEST 1: Atmospheric variables — Central Mediterranean (0.25 deg grid)
# =============================================================================
print("=" * 70)
print("REQUEST 1/3: Atmospheric variables — Central Mediterranean")
print("=" * 70)

request_atmos = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "10m_u_component_of_wind",        # u10: east-west wind component
        "10m_v_component_of_wind",        # v10: north-south wind component
        "sea_surface_temperature",        # sst: critical for hypothermia survival
        "2m_temperature",                 # t2m: air temperature over sea
        "total_cloud_cover",              # tcc: visibility proxy
        "low_cloud_cover",                # lcc: fog proxy
        "2m_dewpoint_temperature",        # d2m: for fog index (t2m - d2m near 0 = fog)
    ],
    "year": YEARS,
    "month": MONTHS,
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": CENTRAL_MED,
}

outfile_1 = os.path.join(OUT_DIR, "era5_central_med_atmos_monthly.nc")
print(f"Downloading to: {outfile_1}")
client.retrieve(DATASET, request_atmos, outfile_1)
print(f"Done: {outfile_1}\n")


# =============================================================================
# REQUEST 2: Wave variables — Central Mediterranean (0.5 deg grid)
# MUST be separate from atmospheric — different spatial grid
# =============================================================================
print("=" * 70)
print("REQUEST 2/3: Wave variables — Central Mediterranean")
print("=" * 70)

request_waves = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "significant_height_of_combined_wind_waves_and_swell",  # swh: key capsizing predictor
        "mean_wave_period",                                      # mwp: short-period waves are more dangerous
        "mean_wave_direction",                                   # mwd: head-seas compound danger
    ],
    "year": YEARS,
    "month": MONTHS,
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": CENTRAL_MED,
}

outfile_2 = os.path.join(OUT_DIR, "era5_central_med_waves_monthly.nc")
print(f"Downloading to: {outfile_2}")
client.retrieve(DATASET, request_waves, outfile_2)
print(f"Done: {outfile_2}\n")


# =============================================================================
# REQUEST 3: Atmospheric variables — North African departure coast (0.25 deg)
# Libya/Tunisia coastal strip: weather conditions at departure
# =============================================================================
print("=" * 70)
print("REQUEST 3/3: Atmospheric variables — North Africa departure coast")
print("=" * 70)

request_coast = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "10m_u_component_of_wind",        # departure coast wind (u-component)
        "10m_v_component_of_wind",        # departure coast wind (v-component)
        "2m_temperature",                 # departure coast temperature
        "total_precipitation",            # departure coast precipitation
        "total_cloud_cover",              # departure coast visibility
    ],
    "year": YEARS,
    "month": MONTHS,
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": NORTH_AFRICA_COAST,
}

outfile_3 = os.path.join(OUT_DIR, "era5_north_africa_coast_monthly.nc")
print(f"Downloading to: {outfile_3}")
client.retrieve(DATASET, request_coast, outfile_3)
print(f"Done: {outfile_3}\n")


# =============================================================================
print("=" * 70)
print("ALL DOWNLOADS COMPLETE")
print(f"Files saved in: {OUT_DIR}")
print("=" * 70)
