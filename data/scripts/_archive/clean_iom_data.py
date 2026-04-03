#!/usr/bin/env python3
"""
IOM Data Cleaning Script
========================

Processes raw IOM Excel data files into clean, validated CSV files.

Input files (in data/raw/iom/):
    - ALL MED DATA 2010-2025_12.08.2025.xlsx  (Mediterranean crossings monthly TS)
    - IOM_MMP data_2014-2025_12.08.2025.xlsx  (Missing Migrants Project incidents)

Output files (in data/processed/):
    - med_crossings_monthlyTS.csv
    - iom_mmp_incidents_2014_2025_reg.csv
    - iom_mmp_2014_2025_all_types.csv
    - iom_mmp_index.csv

Fixes applied:
    - 2025 MMP column typo: 'No. dead/ missing' → 'No. dead/missing'
    - 2021 MMP column swap: 'Names' → renamed, 'No. minors' set to NaN
    - Junk rows (summary/blank) filtered from MMP data
    - Dates standardized to YYYY-MM-DD (imprecise dates preserved in
      'incident_date_raw' column)
    - Known coordinate error corrected: 2022.MMP0765 (Khums, Libya)
    - Crossings paths use hardcoded row boundaries (verified against Excel)

Known source-data issues NOT fixed here (require IOM confirmation):
    - 2021 route-level deaths (EMR/CMR/WMR/WAAR) in crossings file are
      byte-for-byte identical to 2014 — likely a data entry error in the
      original spreadsheet
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
import warnings

warnings.filterwarnings('ignore')


# =============================================================================
# Known coordinate corrections
# =============================================================================
# Format: Main ID → (correct_lat, correct_lon, source)
COORDINATE_FIXES = {
    '2022.MMP0765': (32.6486, 14.2714,
                     'Khums, Libya — original had Qatar coordinates'),
}


# =============================================================================
# Mediterranean crossings (monthly time series)
# =============================================================================

def clean_mediterranean_crossings(excel_path, output_dir):
    """
    Clean Mediterranean crossings data from the IOM Excel file.

    Uses hardcoded row boundaries (verified correct against the Excel
    structure as of 2026-03-14).
    """
    print("=" * 70)
    print("CLEANING MEDITERRANEAN CROSSINGS DATA")
    print("=" * 70)

    print(f"\nReading: {excel_path}")
    df_raw = pd.read_excel(
        excel_path,
        sheet_name='Crossings 2014-24 ALL ROUTE',
        header=None
    )
    print(f"  Raw shape: {df_raw.shape}")

    # Row boundaries per year (0-indexed rows in the header=None dataframe).
    # Each year block is 12 data rows + 1 TOTAL row + 1 separator.
    year_boundaries = {
        2014: (3, 14),   2015: (16, 27),  2016: (29, 40),
        2017: (42, 53),  2018: (55, 66),  2019: (68, 79),
        2020: (81, 92),  2021: (94, 105), 2022: (107, 118),
        2023: (120, 131), 2024: (133, 144), 2025: (146, 157),
    }

    month_map = {
        'January': 1, 'February': 2, 'March': 3, 'April': 4,
        'May': 5, 'June': 6, 'July': 7, 'August': 8,
        'September': 9, 'October': 10, 'November': 11, 'December': 12,
    }

    col_mapping = {
        2: 'waar',
        3: 'wmr_sea_arrivals',
        4: 'wmr_land_arrivals',
        5: 'total_sea_arrivals_waar_wmr',
        6: 'sea_arrivals_in_italy',
        7: 'sea_arrivals_in_malta',
        8: 'land_arrivals_in_greece',
        9: 'sea_arrivals_in_greece',
        10: 'sea_arrivals_in_cyprus',
        11: 'land_arrivals_in_cyprus',
        12: 'total_sea_arrivals_in_europe',
        13: 'total_arrivals_in_europe',
        14: 'interceptions_by_turkish_coast_guard',
        15: 'interceptions_by_libyan_coast_guard',
        16: 'interceptions_by_tunisian_coast_guard',
        17: 'interceptions_by_algerian_coast_guard',
        18: 'total_interceptions',
        19: 'emr',
        20: 'cmr',
        21: 'wmr',
        22: 'waar_1',
        23: 'total_deaths_on_maritime_routes_to_europe',
        24: 'total_attempted_crossings',
        25: 'mortality_rate_maritime_routes_to_europe',
        26: 'mortality_rate_1_in',
        27: 'emr_rate_of_death',
        28: 'cmr_rate_of_death',
        29: 'wmr_proportion_of_deaths_vs_arrivals',
        30: 'waar_proportion_of_deaths_vs_arrivals',
    }

    records = []
    for year, (start_row, end_row) in year_boundaries.items():
        year_data = df_raw.iloc[start_row:end_row + 1, :]

        for _, row in year_data.iterrows():
            month_name = row[1]
            if pd.isna(month_name) or 'TOTAL' in str(month_name):
                continue

            stripped = month_name.strip() if isinstance(month_name, str) else month_name
            month_num = month_map.get(stripped)
            if month_num is None:
                continue

            record = {
                'date': f"{year}-{month_num:02d}-01",
                'year': int(year),
                'month': int(month_num),
                'month_name': stripped,
            }
            for col_idx, col_name in col_mapping.items():
                record[col_name] = row[col_idx] if col_idx < len(row) else np.nan

            records.append(record)

    df_clean = pd.DataFrame(records)
    df_clean = df_clean.sort_values(['year', 'month']).reset_index(drop=True)

    print(f"\n  Rows: {len(df_clean)}")
    print(f"  Years: {df_clean['year'].min()}-{df_clean['year'].max()}")
    for y in sorted(df_clean['year'].unique()):
        n = (df_clean['year'] == y).sum()
        print(f"    {y}: {n} months")

    output_path = Path(output_dir) / 'med_crossings_monthlyTS.csv'
    df_clean.to_csv(output_path, index=False)
    print(f"\n  Saved: {output_path}")

    return df_clean


# =============================================================================
# IOM Missing Migrants Project (incident-level)
# =============================================================================

def _normalize_columns(df):
    """Fix known column name inconsistencies across year sheets."""
    renames = {}
    for col in df.columns:
        # 2025 typo: 'No. dead/ missing' → 'No. dead/missing'
        if col.strip() == 'No. dead/ missing':
            renames[col] = 'No. dead/missing'
    if renames:
        df = df.rename(columns=renames)
    return df


def _parse_date(raw_date):
    """
    Parse a single date value into (date_clean, date_raw, date_precision).

    Returns:
        date_clean: YYYY-MM-DD string or NaT
        date_raw:   original string (kept for imprecise dates)
        date_precision: 'day', 'month', 'imprecise', or 'unknown'
    """
    if pd.isna(raw_date):
        return pd.NaT, '', 'unknown'

    # Already a Timestamp from pandas Excel reading
    if isinstance(raw_date, pd.Timestamp):
        return raw_date.strftime('%Y-%m-%d'), str(raw_date.date()), 'day'

    s = str(raw_date).strip()

    if s.lower() in ('unknown', ''):
        return pd.NaT, s, 'unknown'

    # Timestamp string from concat: "YYYY-MM-DD HH:MM:SS"
    # Extract the date part directly (already in YYYY-MM-DD order).
    m = re.match(r'^(\d{4}-\d{2}-\d{2})\s+\d{2}:\d{2}:\d{2}', s)
    if m:
        return m.group(1), m.group(1), 'day'

    # Plain YYYY-MM-DD (no time component)
    m = re.match(r'^(\d{4})-(\d{2})-(\d{2})$', s)
    if m:
        return s, s, 'day'

    # DD.MM.YYYY format
    m = re.match(r'^(\d{1,2})\.(\d{1,2})\.(\d{4})$', s)
    if m:
        day, month, year = int(m.group(1)), int(m.group(2)), int(m.group(3))
        try:
            dt = pd.Timestamp(year=year, month=month, day=day)
            return dt.strftime('%Y-%m-%d'), s, 'day'
        except ValueError:
            return pd.NaT, s, 'imprecise'

    # Date ranges like '8-9.03.2023' → take first date
    m = re.match(r'^(\d{1,2})-\d{1,2}\.(\d{1,2})\.(\d{4})$', s)
    if m:
        day, month, year = int(m.group(1)), int(m.group(2)), int(m.group(3))
        try:
            dt = pd.Timestamp(year=year, month=month, day=day)
            return dt.strftime('%Y-%m-%d'), s, 'imprecise'
        except ValueError:
            return pd.NaT, s, 'imprecise'

    # YYYY-MM-?? or YYYY-MM-? (month-only precision)
    m = re.match(r'^(\d{4})-(\d{2})-\?+$', s)
    if m:
        year, month = int(m.group(1)), int(m.group(2))
        try:
            dt = pd.Timestamp(year=year, month=month, day=1)
            return dt.strftime('%Y-%m-%d'), s, 'month'
        except ValueError:
            return pd.NaT, s, 'imprecise'

    # ??/MM/YYYY or ??.MM.YYYY (day unknown)
    m = re.match(r'^\?\?[/.](\d{1,2})[/.](\d{4})$', s)
    if m:
        month, year = int(m.group(1)), int(m.group(2))
        try:
            dt = pd.Timestamp(year=year, month=month, day=1)
            return dt.strftime('%Y-%m-%d'), s, 'month'
        except ValueError:
            return pd.NaT, s, 'imprecise'

    # Fallback: try pandas parsing WITHOUT dayfirst to avoid swaps
    try:
        dt = pd.to_datetime(s)
        return dt.strftime('%Y-%m-%d'), s, 'day'
    except (ValueError, TypeError):
        return pd.NaT, s, 'imprecise'


def _standardize_dates(df):
    """Add standardized date columns to the MMP dataframe."""
    results = df['Incident date'].apply(_parse_date)
    df['incident_date_clean'] = results.apply(lambda x: x[0])
    df['incident_date_raw'] = results.apply(lambda x: x[1])
    df['incident_date_precision'] = results.apply(lambda x: x[2])
    return df


def _fix_coordinates(df):
    """Apply known coordinate corrections."""
    n_fixed = 0
    for main_id, (lat, lon, reason) in COORDINATE_FIXES.items():
        mask = df['Main ID'] == main_id
        if mask.any():
            df.loc[mask, 'Latitude'] = lat
            df.loc[mask, 'Longitude'] = lon
            n_fixed += mask.sum()
            print(f"    Fixed coordinates for {main_id}: {reason}")
    return df, n_fixed


def clean_iom_mmp_data(excel_path, output_dir):
    """
    Clean IOM Missing Migrants Project data.

    Reads each year sheet (2014-2025), normalizes columns, concatenates,
    filters junk rows, standardizes dates, and fixes known coordinate errors.
    """
    print("\n" + "=" * 70)
    print("CLEANING IOM MISSING MIGRANTS PROJECT DATA")
    print("=" * 70)

    excel_file = pd.ExcelFile(excel_path)
    print(f"\nReading: {excel_path}")
    print(f"  Sheets: {excel_file.sheet_names}")

    # --- Index sheet ---
    df_index = pd.read_excel(excel_path, sheet_name='Index')
    index_path = Path(output_dir) / 'iom_mmp_index.csv'
    df_index.to_csv(index_path, index=False)
    print(f"\n  Saved index: {index_path} ({len(df_index)} rows)")

    # --- Year sheets ---
    year_sheets = [s for s in excel_file.sheet_names if s.isdigit()]
    all_incidents = []

    for sheet in sorted(year_sheets, key=int):
        df_year = pd.read_excel(excel_path, sheet_name=sheet)
        df_year = _normalize_columns(df_year)
        n_before = len(df_year)

        # Filter junk rows: keep only rows with a non-null Main ID
        df_year = df_year[df_year['Main ID'].notna()].copy()

        # Also filter summary/total rows (Main ID containing 'TOTAL' etc.)
        df_year = df_year[
            ~df_year['Main ID'].astype(str).str.contains(
                'TOTAL|Total|total', na=False
            )
        ].copy()

        n_after = len(df_year)
        n_dropped = n_before - n_after
        suffix = f" (dropped {n_dropped} junk rows)" if n_dropped > 0 else ""
        print(f"  {sheet}: {n_after} rows{suffix}")

        all_incidents.append(df_year)

    df_all = pd.concat(all_incidents, ignore_index=True)
    print(f"\n  Combined: {len(df_all)} rows")

    # --- Standardize dates ---
    print("\n  Standardizing dates...")
    df_all = _standardize_dates(df_all)
    precision_counts = df_all['incident_date_precision'].value_counts()
    for prec, count in precision_counts.items():
        print(f"    {prec}: {count}")

    # --- Fix coordinates ---
    print("\n  Checking coordinates...")
    df_all, n_coord_fixed = _fix_coordinates(df_all)
    if n_coord_fixed == 0:
        print("    No coordinate fixes needed")

    # --- Save all incidents ---
    all_path = Path(output_dir) / 'iom_mmp_2014_2025_all_types.csv'
    df_all.to_csv(all_path, index=False)
    print(f"\n  Saved all types: {all_path} ({len(df_all)} rows)")

    # --- Filter to regular incidents ---
    df_regular = df_all[df_all['Incident Type'] == 'Incident'].copy()
    reg_path = Path(output_dir) / 'iom_mmp_incidents_2014_2025_reg.csv'
    df_regular.to_csv(reg_path, index=False)
    print(f"  Saved regular incidents: {reg_path} ({len(df_regular)} rows)")

    return {
        'index': df_index,
        'all_incidents': df_all,
        'regular_incidents': df_regular,
    }


# =============================================================================
# Validation
# =============================================================================

def validate_crossings(df):
    """Validate the crossings time series."""
    print("\n" + "=" * 70)
    print("VALIDATION: CROSSINGS")
    print("=" * 70)

    issues = []
    expected_years = list(range(2014, 2026))

    actual_years = sorted(df['year'].unique())
    if actual_years != expected_years:
        issues.append(f"Year mismatch: expected {expected_years}, got {actual_years}")
    else:
        print(f"  Years: {actual_years[0]}-{actual_years[-1]}")

    for year in expected_years:
        count = (df['year'] == year).sum()
        if count != 12:
            issues.append(f"Year {year}: expected 12 months, got {count}")

    dups = df.duplicated(subset=['year', 'month']).sum()
    if dups > 0:
        issues.append(f"Found {dups} duplicate year-month rows")

    if issues:
        print("\n  FAILED:")
        for issue in issues:
            print(f"    - {issue}")
        return False

    print(f"  {len(df)} rows, 12 months/year, no duplicates")
    print("  PASSED")
    return True


def validate_mmp(df):
    """Validate the MMP incident data."""
    print("\n" + "=" * 70)
    print("VALIDATION: MMP INCIDENTS")
    print("=" * 70)

    issues = []

    # No null Main IDs
    null_ids = df['Main ID'].isna().sum()
    if null_ids > 0:
        issues.append(f"{null_ids} rows with null Main ID")

    # No duplicate Main IDs
    dup_ids = df['Main ID'].duplicated().sum()
    if dup_ids > 0:
        issues.append(f"{dup_ids} duplicate Main IDs")

    # All should be Incident type
    non_incident = (df['Incident Type'] != 'Incident').sum()
    if non_incident > 0:
        issues.append(f"{non_incident} rows with Incident Type != 'Incident'")

    # Key columns present
    required = ['Incident date', 'Latitude', 'Longitude', 'No. dead',
                'No. dead/missing', 'Region of Incident', 'Route']
    for col in required:
        if col not in df.columns:
            issues.append(f"Missing column: {col}")

    # No. dead/missing should not be all-null for any year
    if 'Incident year' in df.columns and 'No. dead/missing' in df.columns:
        for year in sorted(df['Incident year'].dropna().unique()):
            yr_data = df[df['Incident year'] == year]['No. dead/missing']
            if yr_data.isna().all():
                issues.append(f"Year {int(year)}: No. dead/missing is all NaN")

    # Date standardization worked
    if 'incident_date_clean' in df.columns:
        null_dates = df['incident_date_clean'].isna().sum()
        pct = null_dates / len(df) * 100
        print(f"  Unparseable dates: {null_dates} ({pct:.1f}%)")

    if issues:
        print("\n  FAILED:")
        for issue in issues:
            print(f"    - {issue}")
        return False

    # Summary stats
    n_med = (df['Region of Incident'] == 'Mediterranean').sum()
    n_cmr = (df['Route'] == 'Central Mediterranean').sum()
    n_geocoded = df['Latitude'].notna().sum()
    print(f"  Total incidents: {len(df)}")
    print(f"  Mediterranean: {n_med}")
    print(f"  CMR: {n_cmr}")
    print(f"  Geocoded: {n_geocoded} ({n_geocoded/len(df)*100:.1f}%)")
    print("  PASSED")
    return True


# =============================================================================
# Main
# =============================================================================

def main():
    print("\n" + "=" * 70)
    print("IOM DATA CLEANING SCRIPT")
    print("=" * 70)

    base_dir = Path(__file__).parent
    raw_dir = base_dir / 'raw' / 'iom'
    output_dir = base_dir / 'processed'
    output_dir.mkdir(exist_ok=True)

    med_crossings_excel = raw_dir / 'ALL MED DATA 2010-2025_12.08.2025.xlsx'
    iom_mmp_excel = raw_dir / 'IOM_MMP data_2014-2025_12.08.2025.xlsx'

    # --- Crossings ---
    if med_crossings_excel.exists():
        df_crossings = clean_mediterranean_crossings(
            med_crossings_excel, output_dir
        )
        validate_crossings(df_crossings)
    else:
        print(f"\n  File not found: {med_crossings_excel}")

    # --- MMP ---
    if iom_mmp_excel.exists():
        mmp = clean_iom_mmp_data(iom_mmp_excel, output_dir)
        if mmp['regular_incidents'] is not None:
            validate_mmp(mmp['regular_incidents'])
    else:
        print(f"\n  File not found: {iom_mmp_excel}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)
    print(f"\nOutput: {output_dir}")
    for f in sorted(output_dir.glob("*.csv")):
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name} ({size_mb:.2f} MB)")


if __name__ == "__main__":
    main()
