# Beam Data Analysis — WCTE T09 Beamline

ROOT-based C++ analysis for beamline hit reconstruction: time-of-flight, hodoscope
hit maps, trigger rate studies, and 2-D hit occupancy.

---

## Files

| File | Description |
|------|-------------|
| `basic_analysis.cpp` | Original analysis script — preserved as-is |
| `basic_analysis_v2.cpp` | Cleaned-up version: bug fixes, restructured, documented |

Always prefer `basic_analysis_v2.cpp` for new analysis work.  
The original is kept for reference and cross-checking.

---

## `basic_analysis.cpp`

The original production script used during the March–May 2025 beam tests.
Kept unmodified.  See `basic_analysis_v2.cpp` for the corrected and documented
version.

⚠️  **Known issues in the original** (fixed in v2):
- `tree->GetEntry(i)` was called **after** accessing branch data for TOF reconstruction
  at the first iteration, causing use of stale / uninitialized data.
- Channel constants and physics constants re-declared inside the event loop.
- Hardcoded absolute file path (`/Users/mgola/Documents/...`).
- Hardcoded momentum label (`"Beam Momentum: 500 MeV/c"`).
- Duplicate `SaveAs` call for `c_hodo_hits_2D`.

---

## `basic_analysis_v2.cpp`

Cleaned-up, fully documented version.

### What it does

1. **TOF (time-of-flight) mass identification** — measures time between T0 and T5
   (TOF) scintillators, converts to mass using the known flight path length; fills a
   1-D histogram to identify the beam species.
2. **Hodoscope hit-rate map** — counts T0-normalised hits per hodoscope channel
   (HD0–HD14) and displays as a 1-D bar chart and a 2-D occupancy map.
3. **Beam trigger plots** — histogram of total hits per trigger type.

### Key fixes over v1

| Fix | Description |
|-----|-------------|
| `GetEntry` order | `tree->GetEntry(i)` is now the **first** call in the event loop |
| Constants outside loop | `t0_channels`, `tof_channels`, `L_flight`, `c_light` are `const` and computed once |
| Relative file path | Input file path uses `beamline_run<run>_tuple_calib.root` (relative) |
| Dynamic run label | Momentum label reads the run string: `Form("Run %s", run_str.c_str())` |
| No duplicate SaveAs | Removed second `c_hodo_hits_2D->SaveAs(...)` call |
| Removed dead code | Removed large commented-out lead-glass and old Gaussian fit blocks |
| Section headers | Added comment headers for each analysis section |

### Usage

Compile and run with ROOT:

```bash
root -l -b -q 'basic_analysis_v2.cpp("1827")'
```

The argument is the run number string.  The macro expects the input file in the
working directory named `beamline_run<run>_tuple_calib.root`.

### Channel IDs

| Variable | Channels |
|----------|----------|
| `t0_channels` | {0, 1, 2, 3} |
| `tof_channels` | {48, 49, 50, …} (T5) |
| `hodo_channels` | HD0–HD6 = {32–38}, HD7 = 16, HD8–HD14 = {17–23} |

### Output

- `hodo_hit_rate_<run>.png` — 1-D per-channel hit rate
- `hodo_hits_2D_<run>.png` — 2-D occupancy map
- `tof_mass_<run>.png` — TOF mass histogram
- `trigger_rates_<run>.png` — trigger hit histogram

---

## Dependencies

- ROOT ≥ 6.26 with C++17 support
- Input: `beamline_run<run>_tuple_calib.root` (produced by the offline VME matcher)
