# frontex_incidents_coords.RDS

## What this dataset is

This is the Frontex pad-194 Themis incident-level dataset (14,267 rows), enriched with detection and interception coordinates where possible. Each row is one boat event. The coordinates come from a separate Frontex extract (the "Triton 2014-2017" file) that covers a subset of the same incidents.

The base dataset (`frontex_incidents.RDS`, produced by `04_clean_frontex.R`) contains all original variables. This file adds five columns: four coordinate columns and a match-quality flag.

## How the coordinate matching works

The pad-194 file records every boat individually but has no coordinates. The Triton file covers November 2014 to September 2017 and has coordinates, but it often **aggregates** multiple boats into a single row. For example, if 5 boats were detected on the same day with the same transport type and interceptor, Triton stores them as one row with `n_incidents = 5` and a single coordinate.

The matching procedure (in `04b_merge_triton_coords.R`) works as follows:

1. **Group both datasets** by date, transport type, detected-by, intercepted-by, and SAR flag.
2. **Verify the match**: a group is accepted only if the number of pad rows equals Triton's `n_incidents`, and the sums of migrants and deaths agree exactly.
3. **Assign coordinates**: the (verified) Triton coordinate is assigned to every pad row in the matched group.

Rows outside the Triton date window (before Nov 2014 or after Sep 2017) have no coordinates and `match_quality = NA`.

## Variables

### From the original Frontex pad-194 data

| Variable | Description |
|---|---|
| `incident_id` | Unique incident identifier (Frontex IncidentNumber) |
| `date` | Detection date |
| `country_of_departure` | Country of departure (e.g., Libya, Tunisia) |
| `transport_type` | Vessel type as recorded by Frontex |
| `boat_category` | Simplified vessel category: Inflatable, Wooden, Metal, or Other |
| `num_persons` | Total persons on board |
| `num_deaths` | Number of deaths |
| `num_migrants` | Total irregular migrants |
| `sar_flag` | TRUE if search-and-rescue was involved, FALSE if not, NA if unknown |
| `in_op_area` | TRUE if the incident was inside the Frontex operation area |
| `operation_name` | Name of the Frontex operation (e.g., Triton, Themis) |
| `detected_by` | Who detected the vessel |
| `intercepted_by` | Who intercepted the vessel |
| `ngo_involved` | TRUE if an NGO vessel was involved in detection or interception |
| `event_type` | Detailed event classification (e.g., "SAR: NGO", "Not SAR: Coast Guard") |
| `event_type_agg` | Aggregated event type: "SAR" or "Not SAR" |

### Added by the Triton coordinate merge

| Variable | Description |
|---|---|
| `det_lat` | Detection latitude (WGS84). NA if no coordinate available. |
| `det_lon` | Detection longitude (WGS84). NA if no coordinate available. |
| `int_lat` | Interception latitude (WGS84). NA if no coordinate available. |
| `int_lon` | Interception longitude (WGS84). NA if no coordinate available. |
| `match_quality` | How the coordinate was obtained (see below). NA if not in Triton window. |
| `tri_n_inc` | Number of individual incidents that the matched Triton row aggregated. NA if unmatched. When > 1, the coordinate is shared across that many pad rows. |

## match_quality values

| Value | N rows | Meaning |
|---|---|---|
| `1-1_match` | 1,253 | One Triton row matched one pad row. The coordinate is specific to this event. `tri_n_inc = 1`. |
| `avg_match` | 1,402 | One Triton row matched multiple pad rows. The coordinate is the average stored in that Triton row, shared by all pad rows in the group. `tri_n_inc > 1`. |
| `amb_avg_match` | 243 | Multiple Triton rows matched the same group. The coordinate is a weighted centroid of those Triton rows. Only accepted if all Triton rows are within 28 km of each other (one ERA5 grid cell). |
| `no_coord` | 532 | The pad rows matched a Triton row, but Triton itself had no coordinates for that event. Match verified by metadata, but no lat/lon available. |
| `NA` | 10,837 | Row is outside the Triton date window or had no verifiable match. No coordinate information. |

## Summary

- 14,267 total rows (one per boat event)
- 2,898 rows (20.3%) have coordinates
- 11,369 rows (79.7%) have no coordinates (outside Triton window or unmatched)
