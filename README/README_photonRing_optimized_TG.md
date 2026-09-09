# README — `photonRing_optimized_v2.py`

## Purpose

Event-display and Cherenkov ring reconstruction for WCTE.  Produces unrolled
cylindrical projections of PMT hits for individual beam events, and optionally
animates them as a GIF.  Intended for **visual inspection only** — not part of
the main analysis chain.

## What it does

1. Reads a single WCTE ROOT file (pre-`DataLoader` format, older runs).
2. Loads the WCTE detector geometry from `wcte_bldg157.geo` (AFS).
3. Loads per-channel timing offsets from a JSON calibration file.
4. For each event:
   - Applies a timing cut around the beam peak.
   - Renders an unrolled (φ, z) map of hit PMTs, coloured by charge.
   - Overlays the reconstructed Cherenkov ring.
5. Saves individual PNG frames and assembles them into a GIF.

## Inputs

| Parameter | Description |
|-----------|-------------|
| `root_file` | Path to WCTE ROOT file (single run) |
| `geo_file` | WCTE geometry file (AFS: `TaggedPhoton/TaggedPhoton.Geometry`) |
| `timing_offsets_json` | Per-channel timing offset JSON |
| `event_range` | Tuple `(start, stop)` — events to render |

## Outputs

```
event_display_<run>_ev<N>.png    # one per event
ring_animation_<run>.gif         # animated GIF
```

## Dependencies

```
uproot >= 5
awkward >= 2
numpy
matplotlib
imageio
Pillow
TaggedPhoton.Geometry  (CERN AFS package)
```

## Notes

- Requires CERN AFS access for the geometry file.
- Not parallelised — large event ranges may take several minutes.
- For the main analysis (charge extraction), use `wcte_charge_analysis_v2.py` instead.
