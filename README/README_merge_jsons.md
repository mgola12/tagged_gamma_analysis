# README — `merge_jsons.py`

## Purpose

Alternative merger for per-run JSON charge files.  Uses an integer offset
(`run_index × 10_000_000`) rather than a run-number prefix to de-duplicate event
IDs.  This is **Step 2b** (alternative to `merge_json_by_momentum.py`).

## De-duplication strategy

```python
# Run index 0 → events 0, 1, 2, …
# Run index 1 → events 10_000_000, 10_000_001, …
# Run index 2 → events 20_000_000, …
offset = run_index * 10_000_000
merged[str(int(event_id) + offset)] = charge
```

## Configuration

Edit at the top of the script:

```python
INPUT_FILES = [
    "per_event_total_charge_by_hd_R1827.json",
    "per_event_total_charge_by_hd_R1829.json",
    ...
]
OUTPUT_PREFIX = "merged_hd"   # → merged_hd_<label>.json
LABEL = "968MeVc"
```

## Outputs

```
merged_hd_<label>.json
merged_pairs_<label>.json
```

## Usage

```bash
python3 merge_jsons.py
```

## Notes

- **Prefer `merge_json_by_momentum.py`** — its run-prefixed keys make it easy to
  trace any event back to its source run.
- This script is kept for backward compatibility with earlier analysis versions.
- With ≤6 runs per momentum the offset scheme is safe (no event numbers exceed
  10 million per run in practice).
