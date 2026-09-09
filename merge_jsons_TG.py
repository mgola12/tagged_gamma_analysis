"""
merge_jsons.py

Merges per-run JSON files into one JSON per beam momentum,
for both single-HD and HD-pairs.

Event IDs from different runs are offset by run_index * 10_000_000
to avoid collisions.
"""

import json
import os

FOLDER = "/Users/mgola/Downloads/Tagged_gamma_16072026/Check_mPMT_qualityCuts"

MOMENTUM_RUNS = {
    "968MeVc": [1827, 1829, 1831],
    "774MeVc": [1812, 1814, 1820, 1825],
    "629MeVc": [1804, 1806, 1808, 1810],
    "483MeVc": [1786, 1788, 1790, 1792, 1794, 1796],
}

def merge_files(run_numbers, prefix):
    merged = {}
    for run_idx, run in enumerate(run_numbers):
        fname = os.path.join(FOLDER, f"{prefix}_R{run}_start0_n2000000.json")
        if not os.path.exists(fname):
            print(f"  ⚠ Missing: {fname}")
            continue
        with open(fname) as f:
            data = json.load(f)
        offset = run_idx * 10_000_000
        n_events = 0
        for hd_key, event_dict in data.items():
            if hd_key not in merged:
                merged[hd_key] = {}
            for eid, charge in event_dict.items():
                merged[hd_key][str(int(eid) + offset)] = charge
                n_events += 1
        print(f"  R{run}: {len(data)} HD keys, {n_events} total entries added (offset={offset})")
    return merged

for label, runs in MOMENTUM_RUNS.items():
    print(f"\n{'='*55}")
    print(f"Momentum: {label}")

    # Single-hit
    print("  Single-HD:")
    merged_hd = merge_files(runs, "per_event_total_charge_by_hd")
    out_hd = os.path.join(FOLDER, f"merged_hd_{label}.json")
    with open(out_hd, "w") as f:
        json.dump(merged_hd, f, indent=2)
    print(f"  ✅ Saved: {out_hd}")
    for key, v in sorted(merged_hd.items()):
        print(f"     {key:6s}: {len(v):>7} events")

    # Pairs
    print("  Pairs:")
    merged_pairs = merge_files(runs, "per_event_total_charge_by_hd_pairs")
    out_pairs = os.path.join(FOLDER, f"merged_pairs_{label}.json")
    with open(out_pairs, "w") as f:
        json.dump(merged_pairs, f, indent=2)
    print(f"  ✅ Saved: {out_pairs}")
    for key, v in sorted(merged_pairs.items()):
        print(f"     {key:12s}: {len(v):>7} events")

print(f"\n{'='*55}")
print("Done.")
