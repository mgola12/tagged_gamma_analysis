# README — `merge_json_by_momentum.py`

## Purpose

Merges all per-run JSON charge files for a given beam momentum into a single
combined JSON.  Event keys are prefixed with the run number so entries from
different runs never collide.  This is **Step 2a** of the analysis chain.

## Key design

Event IDs are made globally unique by prepending the run number:

```
"R1827_12345"   ← event 12345 from run 1827
"R1829_12345"   ← same event number but a different run
```

This allows full traceability back to the source run.

## Configuration

Edit these variables at the top of the script before running:

```python
DATA_DIR = "/path/to/per_run_jsons/"   # directory containing per-run JSONs
MOMENTUM_RUNS = {
    "968MeVc": ["R1827", "R1829", "R1831"],
    "774MeVc": ["R1812", "R1820", "R1825"],
    "629MeVc": ["R1804", "R1806", "R1808", "R1810"],
    "483MeVc": ["R1786", "R1788", "R1790", "R1792", "R1794", "R1796"],
}
OUTPUT_DIR = "Check_mPMT_qualityCuts/"
```

> ⚠️ Run 1814 (774 MeV/c) is absent from EOS — omit it from `MOMENTUM_RUNS`.

## Inputs

Per-run JSON files produced by `wcte_charge_analysis_v2.py`:
```
per_event_total_charge_by_hd_<run_label>.json
per_event_total_charge_by_hd_pairs_<run_label>.json
```

## Outputs

```
Check_mPMT_qualityCuts/merged_hd_<momentum>MeVc.json
Check_mPMT_qualityCuts/merged_hd_pairs_<momentum>MeVc.json
```

## Usage

```bash
python3 merge_json_by_momentum.py
```

## Notes

- Preferred over `merge_jsons.py` because run-prefixed keys preserve run traceability.
- Downstream scripts (`overlay_normalize_final_v2.py`, `expectVsMeasured_*.py`) all
  expect the `merged_hd_<momentum>MeVc.json` format produced by this script.
