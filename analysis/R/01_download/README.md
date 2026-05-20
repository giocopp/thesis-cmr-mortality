# Raw-data acquisition

The cleaning pipeline (`02_clean/`) reads from `data/raw/`. Sources:

| Source            | How to obtain                                                   | Snapshot used         |
|-------------------|-----------------------------------------------------------------|-----------------------|
| ERA5 (SWH)        | `01_era5_daily.py` and `02_era5_storm.py` (CDS API key required). | 2008-01-01 to 2025-?? |
| Frontex Themis    | FOI request (incident-level and monthly Themis tables).          | as supplied           |
| IOM MMP           | https://missingmigrants.iom.int — annual CSV download.           | 2014-2025             |
| UNITED            | https://united.unitedagainstrefugeedeaths.eu — incident export.  | 2013-2026             |
| UNHCR             | https://data.unhcr.org — monthly Italy arrivals.                 | 2014-                 |
| ACLED             | https://acleddata.com — Libya/Tunisia conflict events.           | 2010-                 |
| IMO SAR zones     | https://gisis.imo.org — GML polygons.                            | n/a                   |

Place files under `data/raw/{source}/` following the directory layout
expected by each cleaning script (see the header of each script in
`02_clean/`).
