# Tagged Gamma Analysis — WCTE T09 Beamline

Analysis scripts for the tagged photon beam facility at the CERN T09 beamline.  
Dataset: March–May 2025, four beam momenta (483 / 629 / 774 / 968 MeV/c).

---

## Repository layout

```
Tagged_gamma_16072026/
├── wcte_charge_analysis_v2.py       # Step 1 — per-run charge extraction
├── merge_json_by_momentum.py        # Step 2a — merge JSONs (run-prefixed keys)
├── merge_jsons.py                   # Step 2b — merge JSONs (integer-offset keys)
├── overlay_normalize_final_v2.py    # Step 3a — overlay single-HD charge distributions
├── overlay_pairs_normalize.py       # Step 3b — overlay HD-pair charge distributions
├── photonRing_optimized_v2.py       # Event display / Cherenkov ring reconstruction
├── compute_tagging_efficiency.py    # Tagging efficiency per pair and single channel
├── print_beamline_ids.py            # Data quality: list beamline IDs and names
└── analysis_tools/                  # Shared DataLoader (quality cuts, good-PMT list)
```

---

## Physics overview

A positron beam from T09 passes through trigger scintillators (T0, T2) and emits a
bremsstrahlung photon before reaching the Halbach permanent magnet.  The magnet bends
the residual positron into the 15-element hodoscope; the photon continues to WCTE.

A hit on hodoscope element **HD_i** tags the photon energy as:
```
E_γ = p_beam − ⟨p_positron⟩_i
```
where `⟨p_positron⟩_i` is the mean positron momentum for that hodoscope element.

**Beamline channel IDs** (beamline_pmt_tdc_ids):

| ID | Detector |
|----|----------|
| 0–3 | T0 scintillator (left/right × top/bottom) |
| 8 | T2 trigger scintillator |
| 11 | HC2 halo counter (veto) |
| 16 | HD7 |
| 17–23 | HD8–HD14 (back column) |
| 32–38 | HD0–HD6 (front column) |
| 48–63 | TOF / T5 channels (not hodoscope) |

**Offline beam-event selection:** T0 ≥ 1 hit **AND** T2 hit **AND** HC2 = 0.

---

## Step-by-step workflow

### 1 · Extract per-event charges — `wcte_charge_analysis_v2.py`

Reads WCTE production ROOT files from EOS.  For each beam event it:
- Applies mPMT data-quality cuts via `DataLoader` (window mask + readout mask).
- Restricts to good-PMT channels (from `DataLoader.get_good_wcte_pmts()`).
- Fits a Gaussian to the event timing distribution to find the signal peak μ.
- Sums charge in the window `[μ − 20, μ + 20]` ns.
- Records the per-event charge sum keyed by single-HD hit or HD-pair coincidence.

**Outputs** (written to the working directory):
```
per_event_total_charge_by_hd_<tag>.json       # single-HD → {event_id: charge}
per_event_total_charge_by_hd_pairs_<tag>.json # HD-pair   → {event_id: charge}
total_charge_distribution_<tag>.png
charge_distribution_<tag>_<HD>.png            # one per HD channel
charge_distribution_<tag>_<HD>_<HD>.png       # one per pair
```

**Example:**
```bash
python3 wcte_charge_analysis_v2.py \
    --n_events 2000000 \
    --run_patterns "/eos/.../1827/WCTE_merged_production_R1827.root" \
    --run_label R1827
```

**Run → momentum mapping:**

| Momentum | Runs |
|----------|------|
| 968 MeV/c | 1827, 1829, 1831 |
| 774 MeV/c | 1812, 1814, 1820, 1825 (run 1814 missing from EOS) |
| 629 MeV/c | 1804, 1806, 1808, 1810 |
| 483 MeV/c | 1786, 1788, 1790, 1792, 1794, 1796 |

---

### 2a · Merge JSON files by momentum — `merge_json_by_momentum.py`

Merges all per-run JSON files for a given momentum into a single file.  
Event keys are prefixed with the run number (`R1827_12345`) so entries from
different runs never collide.

**Outputs** (in `Check_mPMT_qualityCuts/`):
```
merged_hd_<momentum>MeVc.json
merged_hd_pairs_<momentum>MeVc.json
```

Configure `DATA_DIR` and `MOMENTUM_RUNS` at the top of the script; then run:
```bash
python3 merge_json_by_momentum.py
```

---

### 2b · Alternative merge — `merge_jsons.py`

Same purpose as 2a but uses an integer offset (`run_index × 10_000_000`) rather
than a run-number prefix to de-duplicate event IDs.  Outputs:
```
merged_hd_<label>.json
merged_pairs_<label>.json
```

Both merging scripts produce equivalent data; `merge_json_by_momentum.py` is
preferred because the prefix makes it easier to trace events back to their run.

---

### 3a · Overlay single-HD charge distributions — `overlay_normalize_final_v2.py`

Loads the merged JSON files for one or more beam momenta, plots unit-area normalised
charge histograms for each hodoscope element, and overlays a Gaussian fit per momentum.

**Default:** runs at 968 MeV/c only.  Use `--momenta` to add more:
```bash
# All four momenta
python3 overlay_normalize_final_v2.py \
    --momenta "1000 MeV" "800 MeV" "650 MeV" "500 MeV"

# Custom output directory
python3 overlay_normalize_final_v2.py \
    --momenta "1000 MeV" "800 MeV" \
    --output-dir my_overlay_plots
```

**Skip channels** (beam does not hit these elements at a given momentum) are
configured per momentum inside `BEAM_CONFIG`; they are omitted from that momentum's
histogram rather than plotted as empty.

**Outputs:** one PNG per HD channel in the output directory.

---

### 3b · Overlay HD-pair charge distributions — `overlay_pairs_normalize.py`

Same as 3a but for matched-pair coincidences (HD_i + HD_{i+8}).  
The Gaussian fit is restricted to `charge > 50 000 a.u.` to avoid the zero-spike
at the left edge pulling the fit.

```bash
python3 overlay_pairs_normalize.py          # all four momenta
python3 overlay_pairs_normalize.py --momenta "1000 MeV" "800 MeV"
```

**Pair momentum ranges** (from hodoscope calculator):

| Pair | p-range (MeV/c) |
|------|-----------------|
| HD1+HD9  | 694–727 |
| HD2+HD10 | 560–585 |
| HD3+HD11 | 469–490 |
| HD4+HD12 | 405–422 |
| HD5+HD13 | 356–371 |
| HD6+HD14 | 318–331 |

---

### Auxiliary: Event display — `photonRing_optimized_v2.py`

Reads a single ROOT file (early-format, pre-`DataLoader`) and renders unrolled
cylindrical projections of WCTE hits.  Requires CERN AFS access for the geometry
file (`wcte_bldg157.geo`) and a per-channel timing offset JSON.  Intended for
visual inspection only; not part of the main analysis chain.

---

### Auxiliary: Tagging efficiency — `compute_tagging_efficiency.py`

Computes three efficiency figures per hodoscope pair and single channel:

| Efficiency | Definition |
|-----------|------------|
| ε_hodo | N_tagged / N_beam (geometric acceptance) |
| ε_det  | N_detected / N_tagged (photon detection given tag) |
| ε_total | N_detected / N_beam |

Detection uses the per-channel WCTE charge threshold μ − 3σ derived at 774 MeV/c.
In-time hodoscope hits are required relative to T0 (Group A: −50 to +100 ns;
Group B: +80 to +200 ns).

```bash
python3 compute_tagging_efficiency.py \
    --run_file /eos/.../WCTE_offline_R1827S0_VME_matched.root \
    --n_events 200000 --run_label R1827_774MeVc
```

---

### Auxiliary: Print beamline IDs — `print_beamline_ids.py`

Scans the first N events of a ROOT file and prints every unique `beamline_id` /
`beamline_id_name` pair.  Useful for verifying the channel map of a new run.

```bash
python3 print_beamline_ids.py \
    --run_file /eos/.../WCTE_offline_R1827S0_VME_matched.root \
    --n_events 5000
```

---

## Dependencies

```
uproot >= 5
awkward >= 2
numpy
scipy
matplotlib
analysis_tools   (local package in this directory)
```

For `photonRing_optimized_v2.py` additionally:
```
imageio
Pillow
TaggedPhoton.Geometry  (CERN AFS: /afs/cern.ch/user/m/mogola/condor_jobs/TaggedPhoton)
```
