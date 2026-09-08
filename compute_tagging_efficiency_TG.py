"""
compute_tagging_efficiency.py

Computes the tagged photon beam efficiency in three steps:

  N_beam     : events where T0 (any of IDs 0,1,2,3) AND T2 (ID 8) fire
               AND HC2 (ID 11) does NOT fire in beamline_pmt_tdc_ids,
               AND the WCTE detector has valid hits with non-zero charge

  N_tagged   : N_beam events where a valid in-time hodoscope pair fires:
               front HD_i (beamline_id 32+i) AND back HD_{i+8} (beamline_id 17+i)
               both appear in-time relative to T0, for i in 0..6

  N_detected : N_tagged events where the WCTE charge exceeds the
               per-channel threshold (μ − 3σ) calibrated at 774 MeV/c

Three efficiencies per pair:
  eps_hodo  = N_tagged   / N_beam      (geometric acceptance of hodoscope)
  eps_det   = N_detected / N_tagged    (photon detection eff. given a tag)
  eps_total = N_detected / N_beam      (overall tagging efficiency)

Beamline ID table (from beamline_pmt_tdc_ids):
  T0-0L=0, T0-1L=1, T0-0R=2, T0-1R=3
  T2=8,  HC-2=11
  HD7=16
  HD8=17,  HD9=18,  HD10=19, HD11=20, HD12=21, HD13=22, HD14=23
  HD0=32,  HD1=33,  HD2=34,  HD3=35,  HD4=36,  HD5=37,  HD6=38
  TOF-0..15 = IDs 48-63  (NOT hodoscope)

Hodoscope timing groups (t_HD - t_T0):
  Group A (HD0-HD3, HD8-HD11): window [-50, +100] ns
  Group B (HD4-HD7, HD12-HD14): window [+80, +200] ns

Per-channel detection threshold at 774 MeV/c (mu - 3*sigma from Gaussian fits):
  HD2=30841, HD3=149661, HD4=223430, HD5=271783, HD6=309697
  HD10=80176, HD11=173245, HD12=236336, HD13=280368, HD14=314689
  HD0, HD1, HD7, HD8, HD9: no signal at 774 MeV/c

Usage:
  python3 compute_tagging_efficiency.py \\
      --run_file /eos/experiment/wcte/data/.../WCTE_offline_R1827S0_VME_matched.root \\
      --n_events 200000 --run_label R1827_774MeVc
"""

import uproot
import numpy as np
import awkward as ak
import argparse
import json
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ── Beamline channel IDs ──────────────────────────────────────────────────────
T0_IDS = {0, 1, 2, 3}   # T0-0L, T0-1L, T0-0R, T0-1R
T2_ID  = 8
HC2_ID = 11

# Hodoscope pairs: pair index i (0..6)
#   front: HD_i      → beamline_id = 32 + i   (HD0=32 .. HD6=38)
#   back:  HD_{i+8}  → beamline_id = 17 + i   (HD8=17 .. HD14=23)
PAIR_INDICES = list(range(7))

def pair_name(i):
    return f"HD{i}+HD{i+8}"

# Single HD channels (no veto applied — veto already handled in hardware trigger)
SINGLE_HD_MAP = {
    32: "HD0",  33: "HD1",  34: "HD2",  35: "HD3",
    36: "HD4",  37: "HD5",  38: "HD6",  16: "HD7",
    17: "HD8",  18: "HD9",  19: "HD10", 20: "HD11",
    21: "HD12", 22: "HD13", 23: "HD14",
}
ALL_SINGLE_HD_KEYS = [
    "HD0","HD1","HD2","HD3","HD4","HD5","HD6","HD7",
    "HD8","HD9","HD10","HD11","HD12","HD13","HD14",
]

# ── Hodoscope in-time window (T0-subtracted) ──────────────────────────────────
# Group A: HD0-HD3 (IDs 32-35), HD8-HD11 (IDs 17-20)  → peak ~25-37 ns after T0
# Group B: HD4-HD6 (IDs 36-38), HD7 (ID 16), HD12-HD14 (IDs 21-23) → peak ~126-143 ns
GROUP_A_IDS = {32, 33, 34, 35, 17, 18, 19, 20}
GROUP_B_IDS = {36, 37, 38, 16, 21, 22, 23}

INTIME_A_LO, INTIME_A_HI = -50.0, 100.0   # ns, T0-relative, Group A
INTIME_B_LO, INTIME_B_HI =  80.0, 200.0   # ns, T0-relative, Group B

def is_intime(cid, dt):
    """Return True if dt = t_channel - t_T0 is within the in-time window for cid."""
    if cid in GROUP_A_IDS:
        return INTIME_A_LO < dt < INTIME_A_HI
    if cid in GROUP_B_IDS:
        return INTIME_B_LO < dt < INTIME_B_HI
    return False

# ── Per-channel detection thresholds at 774 MeV/c (mu - 3*sigma) ─────────────
# Derived from Gaussian fits to WCTE total charge distributions.
# HD0, HD1, HD7, HD8, HD9 have no signal at 774 MeV/c — set to None.
THRESH_774 = {
    "HD0":  None,
    "HD1":  None,
    "HD2":  30_841,
    "HD3":  149_661,
    "HD4":  223_430,
    "HD5":  271_783,
    "HD6":  309_697,
    "HD7":  None,
    "HD8":  None,
    "HD9":  None,
    "HD10": 80_176,
    "HD11": 173_245,
    "HD12": 236_336,
    "HD13": 280_368,
    "HD14": 314_689,
}

# For pairs: average of front and back channel thresholds.
# Pair i has no threshold if either front or back is None.
PAIR_THRESH_774 = {}
for _i in PAIR_INDICES:
    _fn, _bn = f"HD{_i}", f"HD{_i+8}"
    _ft, _bt = THRESH_774[_fn], THRESH_774[_bn]
    if _ft is not None and _bt is not None:
        PAIR_THRESH_774[pair_name(_i)] = (_ft + _bt) / 2.0

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# ── Arguments ─────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--run_file", type=str,
    default="/eos/experiment/wcte/data/2025_commissioning/offline_data_vme_match/WCTE_offline_R1827S0_VME_matched.root")
parser.add_argument("--n_events",  type=int, default=200000)
parser.add_argument("--run_label", type=str, default="R1827_774MeVc")
args = parser.parse_args()

print(f"\nRun: {args.run_label}")
print(f"File: {args.run_file}")
print(f"Events: {args.n_events}")
print(f"Group A in-time: [{INTIME_A_LO}, {INTIME_A_HI}] ns  (T0-relative)")
print(f"Group B in-time: [{INTIME_B_LO}, {INTIME_B_HI}] ns  (T0-relative)")
print(f"Pair thresholds at 774 MeV/c (avg front+back, mu-3sigma):")
for _i in PAIR_INDICES:
    _pn = pair_name(_i)
    _th = PAIR_THRESH_774.get(_pn)
    print(f"  {_pn:12s}: {int(_th) if _th is not None else 'N/A'}")

# ── Timing calibration ────────────────────────────────────────────────────────
with open("SR20250401183945_ER20250401183945_ST20250401183945_ET20250401183945.json") as f:
    timing_data = json.load(f)
timing_offset_map = {entry["channel_id"]: entry["timing_offset"]
                     for entry in timing_data["data"]}

# ── Branches needed ───────────────────────────────────────────────────────────
branches = [
    "hit_mpmt_slot_ids", "hit_mpmt_card_ids",
    "hit_pmt_position_ids", "hit_pmt_charges", "hit_pmt_times",
    "beamline_pmt_tdc_ids", "beamline_pmt_tdc_times",
]

CHUNK_SIZE = 10000

print(f"\nReading up to {args.n_events} events in chunks of {CHUNK_SIZE} from:\n  {args.run_file}")
with uproot.open(args.run_file) as f:
    tree = f["WCTEReadoutWindows"]
    n    = min(args.n_events, tree.num_entries)
print(f"Total events to process: {n}")

# ── Counters ──────────────────────────────────────────────────────────────────
N_beam = 0

N_tagged      = {pair_name(i): 0 for i in PAIR_INDICES}
N_detected    = {pair_name(i): 0 for i in PAIR_INDICES}
N_tagged_any  = 0
N_detected_any = 0

N_tagged_single    = {hd: 0 for hd in ALL_SINGLE_HD_KEYS}
N_detected_single  = {hd: 0 for hd in ALL_SINGLE_HD_KEYS}
N_tagged_single_any   = 0
N_detected_single_any = 0

# ── Event loop ────────────────────────────────────────────────────────────────
print("\nRunning event loop ...")
with uproot.open(args.run_file) as f:
    tree = f["WCTEReadoutWindows"]
    for chunk_start in range(0, n, CHUNK_SIZE):
        chunk_stop = min(chunk_start + CHUNK_SIZE, n)
        if chunk_start % 50000 == 0:
            print(f"  {chunk_start}/{n}")
        chunk = tree.arrays(branches, entry_start=chunk_start,
                            entry_stop=chunk_stop, library="ak")
        for event in chunk:

            # ── Beam event: T0 AND T2 fired, HC2 did not fire ─────────────
            tdc_ids_all = set(ak.to_list(event["beamline_pmt_tdc_ids"]))

            if not (T0_IDS & tdc_ids_all):
                continue
            if T2_ID not in tdc_ids_all:
                continue
            if HC2_ID in tdc_ids_all:
                continue

            # ── Valid WCTE hits required ───────────────────────────────────
            ids_sel = event["hit_mpmt_card_ids"] < 120
            if ak.sum(ids_sel) == 0:
                continue

            # ── WCTE charge in timing window ───────────────────────────────
            channel_ids_flat = np.hstack(
                (event["hit_mpmt_slot_ids"] * 100 + event["hit_pmt_position_ids"])[ids_sel])
            digi_charge    = np.hstack(event["hit_pmt_charges"][ids_sel])
            digi_time      = np.hstack(event["hit_pmt_times"][ids_sel])
            offsets        = np.array([timing_offset_map.get(int(c), 0.0)
                                       for c in channel_ids_flat])
            corrected_time = digi_time - offsets

            bins_t = np.linspace(1675, 1760, 80)
            y_t, edges_t = np.histogram(corrected_time, bins=bins_t)
            x_t = 0.5 * (edges_t[:-1] + edges_t[1:])
            try:
                popt, _ = curve_fit(gauss, x_t, y_t, p0=[np.max(y_t), 1710, 3])
                mu_wcte = popt[1]
            except RuntimeError:
                continue

            if mu_wcte < 1675 or mu_wcte > 1760:
                continue

            t_mask = (corrected_time >= mu_wcte - 45) & (corrected_time <= mu_wcte + 45)
            charge = float(np.sum(digi_charge[t_mask]))

            if charge == 0.0:
                continue

            N_beam += 1

            # ── Build full beamline TDC map ────────────────────────────────
            tdc_ids_list   = ak.to_list(event["beamline_pmt_tdc_ids"])
            tdc_times_list = ak.to_list(event["beamline_pmt_tdc_times"])
            id_time_map    = {cid: t for cid, t in zip(tdc_ids_list, tdc_times_list)
                              if t is not None}

            # T0 reference time (earliest T0 hit)
            t0_hits = [id_time_map[cid] for cid in T0_IDS if cid in id_time_map]
            t0_ref  = min(t0_hits)

            # ── Single-channel efficiency ──────────────────────────────────
            any_single_fired    = False
            any_single_detected = False
            for cid, hd in SINGLE_HD_MAP.items():
                if cid not in id_time_map:
                    continue
                dt = id_time_map[cid] - t0_ref
                if not is_intime(cid, dt):
                    continue
                N_tagged_single[hd] += 1
                any_single_fired = True
                thresh = THRESH_774[hd]
                if thresh is not None and charge > thresh:
                    N_detected_single[hd] += 1
                    any_single_detected = True

            if any_single_fired:
                N_tagged_single_any += 1
                if any_single_detected:
                    N_detected_single_any += 1

            # ── Pair efficiency ────────────────────────────────────────────
            any_pair_fired    = False
            any_pair_detected = False
            for i in PAIR_INDICES:
                front_id = 32 + i
                back_id  = 17 + i

                front_ok = (front_id in id_time_map and
                            is_intime(front_id, id_time_map[front_id] - t0_ref))
                back_ok  = (back_id in id_time_map and
                            is_intime(back_id, id_time_map[back_id] - t0_ref))

                if not (front_ok and back_ok):
                    continue

                pname = pair_name(i)
                N_tagged[pname] += 1
                any_pair_fired = True

                thresh = PAIR_THRESH_774.get(pname)
                if thresh is not None and charge > thresh:
                    N_detected[pname] += 1
                    any_pair_detected = True

            if any_pair_fired:
                N_tagged_any += 1
                if any_pair_detected:
                    N_detected_any += 1

# ── Results ───────────────────────────────────────────────────────────────────
print(f"\n{'='*80}")
print(f"Run: {args.run_label}   |   N_beam = {N_beam}")
print(f"{'='*80}\n")

def pct(num, den):
    if den == 0:
        return 0.0, 0.0
    p   = num / den
    err = np.sqrt(p * (1 - p) / den)
    return 100 * p, 100 * err

# ── Per-pair table ────────────────────────────────────────────────────────────
print(f"PAIR EFFICIENCY (threshold = mu-3sigma at 774 MeV/c)")
print(f"{'Pair':12s}  {'Threshold':>10}  {'N_tagged':>9}  {'eps_hodo (%)':>14}  "
      f"{'N_detected':>11}  {'eps_det (%)':>12}  {'eps_total (%)':>14}")
print("-" * 100)

for i in PAIR_INDICES:
    pname  = pair_name(i)
    thresh = PAIR_THRESH_774.get(pname)
    thresh_s = f"{int(thresh):,}" if thresh is not None else "N/A"
    nt     = N_tagged[pname]
    nd     = N_detected[pname]
    eh, eh_e = pct(nt, N_beam)
    ed, ed_e = pct(nd, nt)
    et, et_e = pct(nd, N_beam)
    print(f"{pname:12s}  {thresh_s:>10}  {nt:>9}  "
          f"{eh:>10.2f}±{eh_e:.2f}  "
          f"{nd:>11}  "
          f"{ed:>8.2f}±{ed_e:.2f}  "
          f"{et:>10.2f}±{et_e:.2f}")

print("-" * 100)
eh_t, eh_t_e = pct(N_tagged_any, N_beam)
ed_t, ed_t_e = pct(N_detected_any, N_tagged_any)
et_t, et_t_e = pct(N_detected_any, N_beam)
print(f"{'ALL PAIRS':12s}  {'':>10}  {N_tagged_any:>9}  "
      f"{eh_t:>10.2f}±{eh_t_e:.2f}  "
      f"{N_detected_any:>11}  "
      f"{ed_t:>8.2f}±{ed_t_e:.2f}  "
      f"{et_t:>10.2f}±{et_t_e:.2f}")

# ── Single-channel table ──────────────────────────────────────────────────────
print(f"\nSINGLE HD CHANNEL EFFICIENCY (threshold = mu-3sigma at 774 MeV/c)")
print(f"{'HD':8s}  {'Threshold':>10}  {'N_tagged':>9}  {'eps_hodo (%)':>14}  "
      f"{'N_detected':>11}  {'eps_det (%)':>12}  {'eps_total (%)':>14}")
print("-" * 100)

for hd in ALL_SINGLE_HD_KEYS:
    thresh   = THRESH_774[hd]
    thresh_s = f"{int(thresh):,}" if thresh is not None else "N/A"
    nt       = N_tagged_single[hd]
    nd       = N_detected_single[hd]
    eh, eh_e = pct(nt, N_beam)
    ed, ed_e = pct(nd, nt)
    et, et_e = pct(nd, N_beam)
    print(f"{hd:8s}  {thresh_s:>10}  {nt:>9}  "
          f"{eh:>10.2f}±{eh_e:.2f}  "
          f"{nd:>11}  "
          f"{ed:>8.2f}±{ed_e:.2f}  "
          f"{et:>10.2f}±{et_e:.2f}")

print("-" * 100)
eh_s, eh_s_e = pct(N_tagged_single_any, N_beam)
ed_s, ed_s_e = pct(N_detected_single_any, N_tagged_single_any)
et_s, et_s_e = pct(N_detected_single_any, N_beam)
print(f"{'ALL SINGLE':8s}  {'':>10}  {N_tagged_single_any:>9}  "
      f"{eh_s:>10.2f}±{eh_s_e:.2f}  "
      f"{N_detected_single_any:>11}  "
      f"{ed_s:>8.2f}±{ed_s_e:.2f}  "
      f"{et_s:>10.2f}±{et_s_e:.2f}")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*80}")
print(f"Summary:")
print(f"  N_beam                         = {N_beam}")
print(f"  N_tagged   (any pair)          = {N_tagged_any}   "
      f"eps_hodo = {eh_t:.2f} ± {eh_t_e:.2f} %")
print(f"  N_detected (any pair)          = {N_detected_any}   "
      f"eps_det  = {ed_t:.2f} ± {ed_t_e:.2f} %")
print(f"  eps_total  (pairs)             = {et_t:.2f} ± {et_t_e:.2f} %")
print(f"  N_tagged   (any single HD)     = {N_tagged_single_any}   "
      f"eps_hodo = {eh_s:.2f} ± {eh_s_e:.2f} %")
print(f"  N_detected (any single HD)     = {N_detected_single_any}   "
      f"eps_det  = {ed_s:.2f} ± {ed_s_e:.2f} %")
print(f"  eps_total  (single HD)         = {et_s:.2f} ± {et_s_e:.2f} %")
print(f"{'='*80}")

# ── Plot: eps_det per pair ────────────────────────────────────────────────────
pair_names_with_thresh = [pair_name(i) for i in PAIR_INDICES if pair_name(i) in PAIR_THRESH_774]
eps_det_vals = [pct(N_detected[pn], N_tagged[pn])[0] for pn in pair_names_with_thresh]
eps_det_errs = [pct(N_detected[pn], N_tagged[pn])[1] for pn in pair_names_with_thresh]

if pair_names_with_thresh:
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(pair_names_with_thresh))
    ax.bar(x, eps_det_vals, yerr=eps_det_errs, capsize=5, color="steelblue",
           alpha=0.8, ecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(pair_names_with_thresh, rotation=30, ha="right")
    ax.set_ylabel("Detection efficiency (%)\n[N_detected / N_tagged]", fontsize=12)
    ax.set_title(f"Photon detection efficiency per pair — {args.run_label}\n"
                 f"(threshold = μ−3σ at 774 MeV/c)", fontsize=11)
    ax.set_ylim(0, 110)
    ax.grid(True, axis="y", linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(f"tagging_efficiency_{args.run_label}.png", dpi=150)
    plt.close(fig)
    print(f"\nSaved: tagging_efficiency_{args.run_label}.png")
