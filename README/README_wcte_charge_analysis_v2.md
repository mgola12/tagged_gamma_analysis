# README — `wcte_charge_analysis_v2.py`

## Purpose

Reads WCTE production ROOT files and extracts the total integrated charge for
every beam event, tagged by which hodoscope element (or element pair) recorded
an in-time positron hit.  This is **Step 1** of the tagged-photon analysis chain.

## Algorithm

1. Open the ROOT file via `uproot`.  Apply mPMT data-quality cuts through the
   local `DataLoader` (window mask, readout mask, good-PMT filter).
2. For each beam trigger event that passes T0 ≥ 1, T2, and NOT HC2:
   - Fit a Gaussian to the hit-time distribution to find the signal peak μ.
   - Sum all PMT charges in the window `[μ − 20, μ + 20]` ns.
3. Record the charge sum under every single-HD and HD-pair key that fired an
   in-time hodoscope hit for that event.
4. Write per-event charge dictionaries to JSON.

## Timing windows for hodoscope channels

| Group | Channels | Window relative to T0 |
|-------|----------|----------------------|
| A | HD0–HD3, HD8–HD11 | −50 to +100 ns |
| B | HD4–HD7, HD12–HD14 | +80 to +200 ns |

## Inputs

| Argument | Description |
|----------|-------------|
| `--run_patterns` | Glob(s) or explicit EOS path(s) to the merged production ROOT file |
| `--run_label` | Short tag for output filenames (e.g. `R1827`) |
| `--n_events` | Maximum number of events to process (default: all) |

## Outputs

```
per_event_total_charge_by_hd_<tag>.json        # {HD_key: {event_id: charge}}
per_event_total_charge_by_hd_pairs_<tag>.json  # {pair_key: {event_id: charge}}
total_charge_distribution_<tag>.png
charge_distribution_<tag>_<HD>.png             # one per single-HD channel
charge_distribution_<tag>_<HD>_<HD>.png        # one per HD pair
```

## Example

```bash
python3 wcte_charge_analysis_v2.py \
    --n_events 2000000 \
    --run_patterns "/eos/user/m/mogola/WCTE/1827/WCTE_merged_production_R1827.root" \
    --run_label R1827
```

## Dependencies

```
uproot >= 5
awkward >= 2
numpy
scipy
matplotlib
analysis_tools  (local package — DataLoader)
```

## Notes

- Requires EOS access (CERN AFS/EOS credentials or local mirror).
- `analysis_tools.DataLoader` applies per-run quality masks stored in the package;
  ensure the package is up-to-date with the run being processed.
- Run-to-momentum mapping:

  | Momentum | Runs |
  |----------|------|
  | 968 MeV/c | 1827, 1829, 1831 |
  | 774 MeV/c | 1812, 1814, 1820, 1825 |
  | 629 MeV/c | 1804, 1806, 1808, 1810 |
  | 483 MeV/c | 1786, 1788, 1790, 1792, 1794, 1796 |
