# README — `basic_analysis_v2.cpp`

## Purpose

Cleaned-up, fully documented version of `basic_analysis.cpp`.  Performs the same
analyses — TOF mass identification, hodoscope hit maps, trigger rate histograms —
with all known bugs fixed and dead code removed.

## Analyses performed

1. **TOF mass identification** — time between T0 and T5 scintillators, converted to
   mass using the measured flight path.  Used to confirm beam composition.
2. **Hodoscope hit-rate map** — T0-normalised hit rate per HD channel (HD0–HD14),
   displayed as a 1-D bar chart and a 2-D occupancy map.
3. **Trigger rate histogram** — total hits per beam trigger type.

## Usage

Compile and run with ROOT:

```bash
root -l -b -q 'basic_analysis_v2.cpp("1827")'
```

The argument is the run number string.  The macro reads:
```
beamline_run<run>_tuple_calib.root   (relative path, from working directory)
```

## Channel layout

| Constant | Value | Description |
|----------|-------|-------------|
| `t0_channels` | {0, 1, 2, 3} | T0 scintillator left/right × top/bottom |
| `tof_channels` | {48, 49, …} | T5 TOF scintillators |
| HD0–HD6 | {32–38} | Front hodoscope column |
| HD7 | 16 | Middle element |
| HD8–HD14 | {17–23} | Back hodoscope column |

## Physics constants

| Constant | Value | Description |
|----------|-------|-------------|
| `L_flight` | 8.63 m | T0 → T5 flight path length |
| `c_light` | 299.792 mm/ns | Speed of light |

## Outputs

| File | Description |
|------|-------------|
| `hodo_hit_rate_<run>.png` | 1-D hit rate per HD channel |
| `hodo_hits_2D_<run>.png` | 2-D occupancy map |
| `tof_mass_<run>.png` | TOF mass hypothesis histogram |
| `trigger_rates_<run>.png` | Trigger rate histogram |

## Fixes over v1

| Fix | Description |
|-----|-------------|
| `GetEntry` order | `tree->GetEntry(i)` is now **first** in the event loop |
| `const` constants | `t0_channels`, `tof_channels`, `L_flight`, `c_light` declared once outside loop |
| Relative path | Input file uses `beamline_run<run>_tuple_calib.root` (no hardcoded absolute path) |
| Dynamic label | Plot title reads `Form("Run %s", run_str.c_str())` |
| No duplicate SaveAs | Only one `c_hodo_hits_2D->SaveAs(...)` call |
| No dead code | Removed commented-out lead-glass calibration and old Gaussian fit blocks |
| Section comments | Clear section headers mark each analysis block |

## Dependencies

- ROOT ≥ 6.26 with C++17 support

## See also

- [`basic_analysis.cpp`](basic_analysis.cpp) — original (preserved, unmodified)
- [`README_basic_analysis.md`](README_basic_analysis.md) — original documentation
