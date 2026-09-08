# README — `basic_analysis.cpp`

## Purpose

Original ROOT/C++ analysis macro for the WCTE T09 beamline data.  Performs
time-of-flight (TOF) mass identification, hodoscope hit-rate maps, and beam
trigger rate histograms.

> ⚠️ **This is the preserved original.**  
> For corrected and documented analysis, use `basic_analysis_v2.cpp`.

## Known issues (fixed in v2)

| Issue | Description |
|-------|-------------|
| `GetEntry` bug | `tree->GetEntry(i)` called **after** reading branch data in the first loop iteration — uses stale / uninitialized values |
| Constants in loop | `t0_channels`, `tof_channels`, physics constants redeclared every iteration |
| Hardcoded file path | Absolute path `/Users/mgola/Documents/...` embedded in the ROOT file name |
| Hardcoded label | `"Beam Momentum: 500 MeV/c"` hard-coded regardless of run |
| Duplicate SaveAs | `c_hodo_hits_2D->SaveAs(...)` called twice |
| Dead code | Large commented-out blocks (lead-glass calibration, old Gaussian fit) |

## What it produces

- TOF mass histogram
- Hodoscope 1-D hit-rate bar chart
- Hodoscope 2-D occupancy map
- Beam trigger rate histogram

## Usage

```bash
root -l -b -q 'basic_analysis.cpp("1827")'
```

Input file must be present in the working directory:
```
beamline_run1827_tuple_calib.root
```

## Dependencies

- ROOT ≥ 6.26 with C++17

## See also

- [`basic_analysis_v2.cpp`](basic_analysis_v2.cpp) — corrected version
- [`README_basic_analysis_v2.md`](README_basic_analysis_v2.md) — v2 documentation
