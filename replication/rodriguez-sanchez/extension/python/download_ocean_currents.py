"""
Download monthly mean ocean surface currents from Copernicus Marine.

Dataset: MEDSEA_MULTIYEAR_PHY_006_004 (Mediterranean Sea Physics Reanalysis)
Variables: uo (eastward velocity), vo (northward velocity) at surface level
Region: Central Mediterranean (31N-38N, 10E-20E)
Period: Jan 2009 - Oct 2021

Ocean surface currents directly affect crossing duration, fuel consumption,
and drift trajectory. This is physically distinct from wave height and wind
speed, and perfectly exogenous to migration policy.

Output: NetCDF file in Extension-2-new-data/data/
"""

import copernicusmarine
import os

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "medsea_currents_central_med.nc")

print("Downloading Mediterranean Sea surface currents...")
print("  Dataset: MEDSEA_MULTIYEAR_PHY_006_004")
print("  Variables: uo (eastward), vo (northward)")
print("  Region: 31N-38N, 10E-20E (Central Mediterranean)")
print("  Depth: surface (0-1.5m)")
print("  Period: 2009-01 to 2021-10")
print()

copernicusmarine.subset(
    dataset_id="cmems_mod_med_phy-cur_my_4.2km_P1M-m",
    variables=["uo", "vo"],
    minimum_longitude=10,
    maximum_longitude=20,
    minimum_latitude=31,
    maximum_latitude=38,
    start_datetime="2009-01-01T00:00:00",
    end_datetime="2021-10-31T23:59:59",
    minimum_depth=0,
    maximum_depth=1.5018,
    output_filename="medsea_currents_central_med.nc",
    output_directory=OUTPUT_DIR,
    overwrite=True,
)

print(f"\nSaved to {OUTPUT_FILE}")
