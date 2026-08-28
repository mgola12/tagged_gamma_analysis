import uproot
import numpy as np
import glob
import awkward as ak
import json
import argparse
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit

import sys
from pathlib import Path

root = Path(__file__).resolve().parent.parent

paths = [
    root / "LicketyFit-main",
    root / "LicketyFit-main" / "LicketyFit",
    root / "Geometry", 
     root / "Geometry" / "Geometry"
]

for p in paths:
    if p.exists() and str(p) not in sys.path:
        sys.path.insert(0, str(p))

from EventDisplay import EventDisplay
import imageio.v2 as imageio
import os
from PIL import Image
import matplotlib.colors as colors
sys.path.insert(0, "/home/photogrammetry/rituparna/Geometry-main/Geometry")
from Geometry.WCD import WCD
from Geometry.Device import Device

# --- Utility Functions ---
def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, p0, p):
    v2 = (p - np.array(p0)).flatten()
    v2_u = unit_vector(v2)
    v1_u = unit_vector(v1)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def gauss_function(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# --- Constants ---
vg = 2.20027795333758801e8 * 100 / 1.e10  # mm/ns
origin = np.array([0.0, -200.0, -1520.0 + 188.0])  # mm
#HD_IDS = {33, 17, 18}  # HD1, HD8, HD9
#VETO_HD_IDS = {32, 17}
VETO_HD_IDS = {888}
ALL_HD_IDS = {32, 33, 34, 35, 36, 37, 38, 16, 17, 18, 19, 20, 21, 22, 23}
TOF_IDS = set(range(48, 64))  # TOF channels

hd_id_map = {
    32: "HD0", 33: "HD1", 34: "HD2", 35: "HD3", 36: "HD4",
    37: "HD5", 38: "HD6", 16: "HD7", 17: "HD8", 18: "HD9",
    19: "HD10", 20: "HD11", 21: "HD12", 22: "HD13", 23: "HD14"
}

# --- Argument Parser ---
parser = argparse.ArgumentParser()
parser.add_argument("--n_events", type=int, default=-1, help="Number of events to process")
parser.add_argument("--event_display", action="store_true", help="Enable event display")
parser.add_argument("--start_event", type=int, default=0, help="Starting event index")
args = parser.parse_args()
start_event = args.start_event
chunk_tag = f"start{start_event}"
# --- Geometry ---
my_hall = Device.open_file('/data/rituparna/Geometry/examples/wcte_bldg157.geo')
my_wcte = my_hall.wcds[0]
'''
'''

geometry_cache = {}

for slot_id, mpmt in enumerate(my_wcte.mpmts):
    for pos_id, pmt in enumerate(mpmt.pmts):
        try:
            placement = pmt.get_placement(place_info="design")
        except Exception as e:
            print(f"⚠ PMT {slot_id},{pos_id} get_placement failed: {type(e).__name__}: {e}")
            continue

        loc = (
            placement.get("location")
            or placement.get("loc")
            or placement.get("position")
            or placement.get("pos")
        )

        if loc is None:
            print(f"⚠ PMT {slot_id},{pos_id} placement keys: {list(placement.keys())}")
            continue

        loc = np.asarray(loc, dtype=float)
        theta = angle_between(np.array([0, 0, 1]), origin, loc)
        phi = 0.0
        geometry_cache[(slot_id, pos_id)] = (loc, theta, phi)


# --- Timing offsets ---
with open("SR20250401183945_ER20250401183945_ST20250401183945_ET20250401183945.json", "r") as f:
    timing_data = json.load(f)
timing_offset_map = {entry["channel_id"]: entry["timing_offset"] for entry in timing_data["data"]}

# --- File loading ---
file_pattern = "/data/2025_commissioning/WCTE_offline_R1827S0_VME_matched.root"
file_list = glob.glob(file_pattern)
file_path = file_list[0]

branches_needed = [
    "hit_mpmt_slot_ids",
    "hit_mpmt_card_ids",
    "hit_pmt_position_ids",
    "hit_pmt_charges",
    "hit_pmt_times",
    "beamline_pmt_tdc_ids",
    "beamline_pmt_tdc_times"
]
KNOWN_HD_IDS = set(hd_id_map.keys())
charges_to_plot = []
hd_charge_map_by_event = {}
event_charges_map = {}
all_event_charges = {}
per_event_total_charge_map = {}

# Total charge accumulator for HD channels (excluding veto)
hd_charge_distribution = {hd_id_map[cid]: [] for cid in ALL_HD_IDS - VETO_HD_IDS}
hd_total_charge_map = {hd_id_map[cid]: 0.0 for cid in ALL_HD_IDS - VETO_HD_IDS}
hd_charge_map_by_event = {}
hd_energy_map_by_event = {}
# --- Event Loop ---p
with uproot.open(file_path) as root_file:
    tree = root_file["WCTEReadoutWindows"]
    total_entries = tree.num_entries
    print(f"Total entries in tree: {total_entries}")
    start_event = args.start_event
    if args.n_events > 0:
        end_event = min(start_event + args.n_events, total_entries)
    else:
        end_event = total_entries
    print(f"Processing events from {start_event} to {end_event - 1}")
    #events = tree.arrays(branches_needed, entry_start=event_id, entry_stop=event_id + 1, library="ak")
    events = tree.arrays(branches_needed, entry_start=start_event, entry_stop=end_event, library="ak")
    for i_event, event in enumerate(events):
        event_id = start_event + i_event
    #for event_id in range(n_to_run):
        if event_id % 1000 == 0:
            print(f"\n▶ Processing event {event_id}")

        if len(events) == 0:
            print(f"⚠ No event found at index {event_id}")
            continue
        #event = events[0]
        event = events[i_event]

        channel_ids = (event['hit_mpmt_slot_ids'] * 100) + event['hit_pmt_position_ids']
        ids_selection = (event['hit_mpmt_card_ids'] < 120)

        if ak.sum(ids_selection) == 0:
            print(f"⚠ Event {event_id}: No valid hits after card ID selection. Skipping.")
            continue
        channel_ids_flat = np.hstack(channel_ids[ids_selection])
        digi_hit_pmt = np.hstack(((event['hit_mpmt_slot_ids']*19)+event['hit_pmt_position_ids'])[ids_selection])
        digi_hit_charge = np.hstack(event['hit_pmt_charges'][ids_selection])
        digi_hit_time = np.hstack(event['hit_pmt_times'][ids_selection])

        hit_pmt_position_list = []
        hit_pmt_angle_list = []
        hit_pmt_phi_angle_list = []
        hit_pmt_charge_list = []
        hit_pmt_time_list = []
        for i, pmt in enumerate(digi_hit_pmt):
            slot_ids = pmt // 19
            position_ids = pmt % 19
            pmt_position, theta, phi = geometry_cache.get((slot_ids, position_ids), (None, None, None))
            if pmt_position is None:
                continue
            theta = angle_between(np.array([0., 0., 1.]), origin, pmt_position)
            pmt_tof = np.linalg.norm(pmt_position - origin) / vg

            sin_theta = np.sin(theta)
            if sin_theta != 0:
                cos_phi = np.cos(angle_between(np.array([1., 0., 0.]), origin, pmt_position)) / sin_theta
                cos_phi = np.clip(cos_phi, -1.0, 1.0)
                phi = np.arccos(cos_phi) if pmt_position[1] >= origin[1] else 2 * np.pi - np.arccos(cos_phi)
            else:
                phi = 0.0

            hit_pmt_position_list.append(pmt_position)
            hit_pmt_angle_list.append(theta)
            hit_pmt_phi_angle_list.append(phi)
            hit_pmt_charge_list.append(digi_hit_charge[i])
            hit_pmt_time_list.append(digi_hit_time[i] - pmt_tof)

        timing_offsets = np.array([timing_offset_map.get(cid, 0.0) for cid in channel_ids_flat])
        corrected_digi_hit_time = digi_hit_time - timing_offsets

        bins = np.linspace(1675, 1760, 80)
        y, edges = np.histogram(corrected_digi_hit_time, bins=bins)
        x_centers = 0.5 * (edges[:-1] + edges[1:])
        try:
            popt, _ = curve_fit(gauss_function, x_centers, y, p0=[100, 1710, 3])
            A, mu, sigma = popt
            params_text = f"A = {A:.2f}\nμ = {mu:.2f}\nσ = {sigma:.2f}"
        except RuntimeError as e:
            print(f"⚠ Event {event_id}: Gaussian fit failed — {e}")
            continue


        tdc_ids = ak.to_list(event["beamline_pmt_tdc_ids"])
        tdc_times = ak.to_list(event["beamline_pmt_tdc_times"])

        
        id_time_map = {
            cid: time for cid, time in zip(tdc_ids, tdc_times)
            if time is not None
        }

        #print(f"Event {event_id}: id_time_map keys = {sorted(id_time_map.keys())}")
        
        '''    
        '''
        veto_hit = VETO_HD_IDS.intersection(id_time_map)
        tof_hit = TOF_IDS.intersection(id_time_map)

        
        #if len(veto_hit) > 0 or len(tof_hit) > 0:
        #    print(f"⛔ Skipping: Veto HD = {veto_hit}, TOF hit = {tof_hit}")
        #    continue
        '''
        '''
        #Selected time window before summing charge
        time_min = mu - 2 * sigma
        time_max = mu + 2 * sigma
        valid_mask = (
            (corrected_digi_hit_time >= time_min) &
            (corrected_digi_hit_time <= time_max)
        )


        charge_sum = float(np.sum(digi_hit_charge[valid_mask]))
        
        #charge_sum = float(np.sum(event["hit_pmt_charges"])) # simple summing the charge
        '''
        For Subtracting the dark rate per event each PMT
        Q_DARK_1PE = 140  # Typical 1 PE charge from your Q1PE range (130–150)
        dark_charge_subtracted = []
        
        for q in digi_hit_charge:
            if q < 200:  # Likely a dark hit (unphysical low charge)
                dark_charge_subtracted.append(0.0)
            else:
                dark_charge_subtracted.append(q - Q_DARK_1PE)
        
        charge_sum = float(np.sum(dark_charge_subtracted))
        '''
        charges_to_plot.append(charge_sum)
        per_event_total_charge_map[str(event_id)] = charge_sum
        
        for cid in ALL_HD_IDS - VETO_HD_IDS:
            if cid in id_time_map:
                hd_key = hd_id_map.get(cid, f"HD_{cid}")
                if hd_key not in hd_charge_map_by_event:
                    hd_charge_map_by_event[hd_key] = {}
                hd_charge_map_by_event[hd_key][str(event_id)] = charge_sum
                if hd_key not in hd_energy_map_by_event:
                    hd_energy_map_by_event[hd_key] = {}
                hd_energy_map_by_event[hd_key][str(event_id)] = 1.06395e-3 * charge_sum - 30.13
                
        
        '''
        '''
        #print(f"✅ Saved: event_display_entry{event_id}.png & _Corrected.png")

'''
'''

with open(f"per_event_total_charge_by_hd_all_Events_Run1827.json", "w") as f:
    json.dump(hd_charge_map_by_event, f, indent=2)

with open(f"per_event_total_energy_by_hd_all_Events_Run1827.json", "w") as f:
    json.dump(hd_energy_map_by_event, f, indent=2)


#print(f"✅ Saved: per_event_total_charge_by_hd_all_Events_Run1808_{chunk_tag}.json")

#Charge sum only
plt.figure()
bins = 100
y, edges = np.histogram(charges_to_plot, bins=bins)
x_centers = 0.5 * (edges[:-1] + edges[1:])

A_init = np.max(y)
mu_init = x_centers[np.argmax(y)]
sigma_init = np.std(x_centers)

if not np.isfinite(sigma_init) or sigma_init <= 0:
    sigma_init = 1.0

try:
    popt, _ = curve_fit(
        gauss_function,
        x_centers,
        y,
        p0=[A_init, mu_init, sigma_init],
        bounds=([0, -np.inf, 1e-6], [np.inf, np.inf, np.inf]),
        maxfev=10000
    )
    A_fit, mu_fit, sigma_fit = popt
    plt.hist(charges_to_plot, bins=100, histtype='step', linewidth=1.5, label="Data")
    x_fit = np.linspace(min(x_centers), max(x_centers), 1000)
    y_fit = gauss_function(x_fit, *popt)
    plt.plot(x_fit, y_fit, label="Gaussian fit", color="orange")

    fit_text = f"$\\mu = {mu_fit:.2f}$\n$\\sigma = {sigma_fit:.2f}$"
    plt.text(0.95, 0.95, fit_text, transform=plt.gca().transAxes,
             verticalalignment='top', horizontalalignment='right')

except RuntimeError as e:
    print(f"⚠ Gaussian fit failed: {e}")
    plt.hist(charges_to_plot, bins=100, histtype='step', linewidth=1.5, label="Data")



plt.xlabel("Total Charge (all PMTs)")
plt.ylabel("Event Count")
plt.title("Total PMT Charge Distribution (All Events)")
plt.legend()
        # plt.savefig("total_charge_distribution_all_events.png", dpi=150)
plt.close()
#print("✅ Saved: total_charge_distribution_all_events_{chunk_tag}.png (with Gaussian fit)")


    
# Plot histograms: one per HD
for hd_key, event_charge_dict in hd_charge_map_by_event.items():
    charges = list(event_charge_dict.values())
    if len(charges) == 0:
        continue
    plt.figure()
    plt.hist(charges, bins=100, histtype='step', linewidth=1.5)
    plt.xlabel("Total PMT Charge (per event)")
    plt.ylabel("Counts")
    plt.title(f"Total Charge Distribution — {hd_key}")
    plt.xlim(0, 1.0e6)  
    plt.savefig(f"charge_distribution_Run_1827_{hd_key}.pdf", dpi=150)
    plt.close()
    
for hd_key, event_energy_dict in hd_energy_map_by_event.items():
    energies = list(event_energy_dict.values())
    if len(energies) == 0:
        continue
    plt.figure()
    plt.hist(energies, bins=100, histtype='step', linewidth=1.5)
    plt.xlabel("Photon Energy (MeV)")
    plt.ylabel("Counts")
    plt.title(f"Photon Energy Distribution — {hd_key}")
    plt.xlim(0, 1500)  
    plt.savefig(f"energy_distribution_Run_1827_{hd_key}.pdf", dpi=150)
    plt.close()
