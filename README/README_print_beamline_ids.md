# README — `print_beamline_ids.py`

## Purpose

Data-quality utility.  Scans the first N events of a beamline ROOT file and
prints every unique `beamline_id` / `beamline_id_name` pair found in the data.
Use this to verify the channel map of a new run before running the full analysis.

## Expected output

```
beamline_id   beamline_id_name
-----------   ----------------
0             T0_left_top
1             T0_left_bot
2             T0_right_top
3             T0_right_bot
8             T2
11            HC2
16            HD7
17            HD8
...
```

## CLI

```bash
python3 print_beamline_ids.py \
    --run_file /eos/.../WCTE_offline_R1827S0_VME_matched.root \
    --n_events 5000
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--run_file` | *(required)* | Path to offline VME-matched ROOT file |
| `--n_events` | 5000 | Number of events to scan |

## Known channel IDs

| ID | Detector | Role |
|----|----------|------|
| 0–3 | T0 (L/R × top/bot) | Beam trigger |
| 8 | T2 | Downstream trigger |
| 11 | HC2 | Halo veto |
| 16 | HD7 | Hodoscope (middle) |
| 17–23 | HD8–HD14 | Hodoscope (back) |
| 32–38 | HD0–HD6 | Hodoscope (front) |

## Dependencies

```
uproot >= 5
awkward >= 2
```
