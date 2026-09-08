"""
wcte_charge_analysis_v2.py

Produces per-run:
  - Single-hit charge distributions per hodoscope element (JSON + PNG)
  - Matched-pair charge distributions HD0+HD8 .. HD6+HD14 (JSON + PNG)
  - Total charge distribution (PNG)

Changes vs v1:
  1. Data-quality cuts via DataLoader:
       - loader.apply_mPMT_data_quality_cuts()   → window_data_quality_mask + hit_pmt_readout_mask
     DataLoader.iterate() is used so the cut is applied batch-by-batch automatically.

  2. Offline beam-event selection (mirrors the hardware trigger):
       - T0 fired  : at least one of IDs {0, 1, 2, 3} present in beamline_pmt_tdc_ids
       - T2 fired  : ID 8 present in beamline_pmt_tdc_ids
       - HC2 veto  : ID 11 must NOT be present in beamline_pmt_tdc_ids

Usage:
    python3 wcte_charge_analysis_v2.py --n_events 2000000 --start_event 0 \\
        --run_patterns "/eos/.../production_v1_0/1827/WCTE_merged_production_R1827.root" \\
        --run_label "R1827"
"""

import uproot
import numpy as np
import glob
import awkward as ak
import json
import argparse
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
from analysis_tools import DataLoader

# ── Beamline trigger / veto IDs ───────────────────────────────────────────────
# From compute_tagging_efficiency.py (confirmed against beamline mapping)
T0_IDS = {0, 1, 2, 3}   # T0-0L, T0-1L, T0-0R, T0-1R  (require ≥1)
T2_ID  = 8               # T2 scintillator
HC2_ID = 11              # HC2 halo counter — veto

# ── Hodoscope channel IDs ─────────────────────────────────────────────────────
# Front: HD0=32, HD1=33, HD2=34, HD3=35, HD4=36, HD5=37, HD6=38, HD7=16
# Back:  HD8=17, HD9=18, HD10=19, HD11=20, HD12=21, HD13=22, HD14=23
ALL_HD_IDS = {32, 33, 34, 35, 36, 37, 38, 16, 17, 18, 19, 20, 21, 22, 23}

hd_id_map = {
    32: "HD0",  33: "HD1",  34: "HD2",  35: "HD3",  36: "HD4",
    37: "HD5",  38: "HD6",  16: "HD7",  17: "HD8",  18: "HD9",
    19: "HD10", 20: "HD11", 21: "HD12", 22: "HD13", 23: "HD14"
}

# Hodoscope pairs: front HD_i (32+i) + back HD_{i+8} (17+i), i = 0..6
# HD7 (ID 16) has no back partner — single-hit only
PAIR_INDICES = list(range(7))

def pair_name(i):
    return f"HD{i}+HD{i+8}"

def gauss_function(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# ── Argument parser ────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--n_events",     type=int, default=100000)
parser.add_argument("--start_event",  type=int, default=0)
parser.add_argument("--run_patterns", type=str,
    default="/eos/experiment/wcte/data/2025_commissioning/processed_offline_data/"
            "production_v1_0/1827/WCTE_merged_production_R1827.root")
parser.add_argument("--run_label",    type=str, default="Run")
args      = parser.parse_args()
run_label = args.run_label
chunk_tag = f"{run_label}_start{args.start_event}_n{args.n_events}"

# ── Collect run files ──────────────────────────────────────────────────────────
run_file_list = []
for pattern in args.run_patterns.split(","):
    matched = glob.glob(pattern.strip())
    if not matched:
        print(f"No files matched pattern: {pattern.strip()}")
    run_file_list.extend(matched)

if not run_file_list:
    print("No run files found. Exiting.")
    sys.exit(1)

print(f"Run files ({len(run_file_list)}):")
for f in run_file_list:
    print(f"  {f}")

# ── Branches (DataLoader will append mask branches automatically) ─────────────
branches_needed = [
    "hit_mpmt_slot_ids", "hit_pmt_position_ids",
    "hit_pmt_charges", "hit_pmt_calibrated_times",
    "beamline_pmt_tdc_ids", "beamline_pmt_tdc_times",
]

# ── Accumulators ───────────────────────────────────────────────────────────────
hd_charge_map_by_event   = {}
pair_charge_map_by_event = {}
charges_to_plot          = []

global_event_id = 0

# ── Event loop — one file at a time ───────────────────────────────────────────
for file_path in run_file_list:
    print(f"\n{'='*60}")
    print(f"Processing: {file_path}")

    # ── DataLoader: enable all three layers of quality cuts ───────────────────
    loader = DataLoader(file_path, branches_to_load=list(branches_needed))
    loader.apply_mPMT_data_quality_cuts()   # window_data_quality_mask + hit_pmt_readout_mask

    # ── Good-PMT channel set for this run ─────────────────────────────────────
    good_slots, good_positions = loader.get_good_wcte_pmts()
    good_channel_set = set(int(s) * 100 + int(p)
                           for s, p in zip(good_slots, good_positions))
    print(f"  Good PMT channels: {len(good_channel_set)}")

    # ── Iterate in batches (cuts applied inside DataLoader.iterate) ───────────
    for batch in loader.iterate(
            verbose=True,
            entry_start=args.start_event,
            entry_stop=args.start_event + args.n_events):

        for i_local in range(len(batch)):
            event = batch[i_local]

            # ── WCTE hit arrays ────────────────────────────────────────────────
            if ak.count(event['hit_pmt_charges']) == 0:
                continue

            channel_ids_flat = np.asarray(
                ak.to_numpy(event['hit_mpmt_slot_ids'] * 100
                            + event['hit_pmt_position_ids']))
            digi_hit_charge  = np.asarray(ak.to_numpy(event['hit_pmt_charges']))
            digi_hit_time    = np.asarray(ak.to_numpy(event['hit_pmt_calibrated_times']))

            # ── Good-PMT mask ──────────────────────────────────────────────────
            good_mask = np.array([cid in good_channel_set for cid in channel_ids_flat])
            if not np.any(good_mask):
                continue
            channel_ids_flat = channel_ids_flat[good_mask]
            digi_hit_charge  = digi_hit_charge[good_mask]
            hit_time         = digi_hit_time[good_mask]

            # ── Beamline TDC map ───────────────────────────────────────────────
            tdc_ids   = ak.to_list(event["beamline_pmt_tdc_ids"])
            tdc_times = ak.to_list(event["beamline_pmt_tdc_times"])
            id_time_map  = {cid: t for cid, t in zip(tdc_ids, tdc_times) if t is not None}
            tdc_ids_set  = set(id_time_map.keys())

            # ── Offline beam-event selection: T0 ∧ T2 ∧ ¬HC2 ─────────────────
            if not (T0_IDS & tdc_ids_set):        # require ≥1 T0 channel
                continue
            if T2_ID not in tdc_ids_set:           # require T2
                continue
            if HC2_ID in tdc_ids_set:              # veto HC2
                continue

            # ── Gaussian fit on timing distribution ────────────────────────────
            bins      = np.linspace(1675, 1760, 80)
            y, edges  = np.histogram(hit_time, bins=bins)
            x_centers = 0.5 * (edges[:-1] + edges[1:])
            try:
                popt, _ = curve_fit(gauss_function, x_centers, y, p0=[100, 1710, 3])
                A, mu, sigma = popt
            except RuntimeError as e:
                print(f"  Event {global_event_id}: Gaussian fit failed — {e}")
                global_event_id += 1
                continue

            # ── Charge sum in timing window (mu ± 20 ns), good PMTs only ──────
            intime_mask = (hit_time >= mu - 20) & (hit_time <= mu + 20)
            charge_sum  = float(np.sum(digi_hit_charge[intime_mask]))
            charges_to_plot.append(charge_sum)

            # ── Single-hit: any HD channel present in TDC map ─────────────────
            for cid in ALL_HD_IDS:
                if cid in tdc_ids_set:
                    hd_key = hd_id_map[cid]
                    hd_charge_map_by_event.setdefault(hd_key, {})[str(global_event_id)] = charge_sum

            # ── Pairs: front HD_i AND back HD_{i+8} both present ──────────────
            for i in PAIR_INDICES:
                if (32 + i) in tdc_ids_set and (17 + i) in tdc_ids_set:
                    key = pair_name(i)
                    pair_charge_map_by_event.setdefault(key, {})[str(global_event_id)] = charge_sum

            if global_event_id % 10000 == 0:
                print(f"  Processed event {global_event_id}")
            global_event_id += 1

# ── Save JSON ──────────────────────────────────────────────────────────────────
single_json_path = f"per_event_total_charge_by_hd_{chunk_tag}.json"
with open(single_json_path, "w") as f:
    json.dump(hd_charge_map_by_event, f, indent=2)
print(f"\nSaved: {single_json_path}")

pair_json_path = f"per_event_total_charge_by_hd_pairs_{chunk_tag}.json"
with open(pair_json_path, "w") as f:
    json.dump(pair_charge_map_by_event, f, indent=2)
print(f"Saved: {pair_json_path}")

# ── Summary ────────────────────────────────────────────────────────────────────
print(f"\n{'='*55}")
print("Single-hit counts (after all cuts + T0∧T2∧¬HC2 selection):")
for hd_key in sorted(hd_charge_map_by_event):
    print(f"  {hd_key:6s} : {len(hd_charge_map_by_event[hd_key])} events")
print("\nMatched pair counts:")
for i in PAIR_INDICES:
    key = pair_name(i)
    print(f"  {key:12s} : {len(pair_charge_map_by_event.get(key, {}))}")
print(f"{'='*55}")

# ── Total charge distribution ──────────────────────────────────────────────────
if charges_to_plot:
    y_c, edges_c = np.histogram(charges_to_plot, bins=100)
    x_c = 0.5 * (edges_c[:-1] + edges_c[1:])
    plt.figure()
    plt.hist(charges_to_plot, bins=100, histtype='step', linewidth=1.5, label="Data")
    try:
        popt_c, _ = curve_fit(gauss_function, x_c, y_c,
                               p0=[np.max(y_c), x_c[np.argmax(y_c)], np.std(charges_to_plot)/2])
        x_fit = np.linspace(min(x_c), max(x_c), 1000)
        plt.plot(x_fit, gauss_function(x_fit, *popt_c), color="orange", label="Gaussian fit")
        plt.text(0.95, 0.95, f"mu={popt_c[1]:.0f}\nsigma={popt_c[2]:.0f}",
                 transform=plt.gca().transAxes, va='top', ha='right')
    except RuntimeError:
        pass
    plt.xlabel("Total Charge (all PMTs)")
    plt.ylabel("Event Count")
    plt.title(f"Total PMT Charge Distribution — {run_label}")
    plt.legend()
    plt.savefig(f"total_charge_distribution_{chunk_tag}.png", dpi=150)
    plt.close()

# ── Per-HD plots ───────────────────────────────────────────────────────────────
for hd_key, event_charge_dict in hd_charge_map_by_event.items():
    charges = list(event_charge_dict.values())
    if not charges:
        continue
    plt.figure()
    plt.hist(charges, bins=100, histtype='step', linewidth=1.5)
    plt.xlabel("Total PMT Charge (per event)")
    plt.ylabel("Counts")
    plt.title(f"Total Charge Distribution — {hd_key} — {run_label}")
    plt.xlim(0, 1.0e6)
    plt.savefig(f"charge_distribution_{chunk_tag}_{hd_key}.png", dpi=150)
    plt.close()
    print(f"Saved: charge_distribution_{chunk_tag}_{hd_key}.png")

# ── Per-pair plots ─────────────────────────────────────────────────────────────
for i in PAIR_INDICES:
    key     = pair_name(i)
    charges = list(pair_charge_map_by_event.get(key, {}).values())
    if not charges:
        continue
    plt.figure()
    plt.hist(charges, bins=100, histtype='step', linewidth=1.5, color='darkorange')
    plt.xlabel("Total PMT Charge (per event)")
    plt.ylabel("Counts")
    plt.title(f"Total Charge Distribution — {key} — {run_label}")
    plt.xlim(0, 1.0e6)
    plt.savefig(f"charge_distribution_{chunk_tag}_{key.replace('+', '_')}.png", dpi=150)
    plt.close()
    print(f"Saved: charge_distribution_{chunk_tag}_{key.replace('+', '_')}.png")
