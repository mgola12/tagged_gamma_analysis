"""
merge_json_by_momentum.py

Merges per-run JSON charge files into one file per beam momentum.
Event keys are prefixed with the run number (e.g. "R1827_12345") so
entries from different runs never collide.

Output files (written to the same directory as the input files):
  merged_hd_{momentum}MeVc.json        — single-hit
  merged_hd_pairs_{momentum}MeVc.json  — matched pairs
"""

import json
import os

DATA_DIR = "/Users/mgola/Downloads/Tagged_gamma_16072026/Check_mPMT_qualityCuts"

MOMENTUM_RUNS = {
    "968": [1827, 1829, 1831],
    "774": [1812, 1814, 1820, 1825],
    "629": [1804, 1806, 1808, 1810],
    "483": [1786, 1788, 1790, 1792, 1794, 1796],
}

SUFFIX = "_start0_n2000000"

for momentum, runs in MOMENTUM_RUNS.items():
    merged_single = {}
    merged_pairs  = {}

    for run in runs:
        tag = f"R{run}{SUFFIX}"

        # ── Single-hit ────────────────────────────────────────────────────────
        single_path = os.path.join(DATA_DIR, f"per_event_total_charge_by_hd_{tag}.json")
        if os.path.exists(single_path):
            with open(single_path) as f:
                data = json.load(f)
            for hd_key, event_dict in data.items():
                merged_single.setdefault(hd_key, {})
                for evt_id, charge in event_dict.items():
                    merged_single[hd_key][f"R{run}_{evt_id}"] = charge
            print(f"  [{momentum} MeV/c] loaded single: {single_path}")
        else:
            print(f"  WARNING: missing {single_path}")

        # ── Pairs ─────────────────────────────────────────────────────────────
        pairs_path = os.path.join(DATA_DIR, f"per_event_total_charge_by_hd_pairs_{tag}.json")
        if os.path.exists(pairs_path):
            with open(pairs_path) as f:
                data = json.load(f)
            for pair_key, event_dict in data.items():
                merged_pairs.setdefault(pair_key, {})
                for evt_id, charge in event_dict.items():
                    merged_pairs[pair_key][f"R{run}_{evt_id}"] = charge
            print(f"  [{momentum} MeV/c] loaded pairs : {pairs_path}")
        else:
            print(f"  WARNING: missing {pairs_path}")

    # ── Write merged files ────────────────────────────────────────────────────
    out_single = os.path.join(DATA_DIR, f"merged_hd_{momentum}MeVc.json")
    with open(out_single, "w") as f:
        json.dump(merged_single, f, indent=2)
    total_single = sum(len(v) for v in merged_single.values())
    print(f"  → Saved: {out_single}  ({total_single} total events across all HD channels)")

    out_pairs = os.path.join(DATA_DIR, f"merged_hd_pairs_{momentum}MeVc.json")
    with open(out_pairs, "w") as f:
        json.dump(merged_pairs, f, indent=2)
    total_pairs = sum(len(v) for v in merged_pairs.values())
    print(f"  → Saved: {out_pairs}  ({total_pairs} total events across all pairs)\n")

print("Done.")
