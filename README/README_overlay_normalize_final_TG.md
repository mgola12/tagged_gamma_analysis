# README — `overlay_normalize_final_v2.py`

## Purpose

Loads merged JSON charge files, plots unit-area normalised charge histograms for
each hodoscope element (single-HD), and overlays a Gaussian fit per beam momentum.
This is **Step 3a** of the analysis chain and is the source of the Gaussian fit
parameters (μ, σ) used downstream by `expectVsMeasured_*.py`.

## Algorithm

1. Load `merged_hd_<momentum>MeVc.json` for each requested momentum.
2. For each HD channel and each momentum, histogram the per-event charge values.
3. Normalise the histogram to unit area.
4. Fit a Gaussian in the range `[charge_min, charge_max]` — configurable per
   channel to exclude the low-charge noise peak.
5. Overlay all momenta on one canvas per HD channel and annotate with fit results.

## Configuration

At the top of the script:

```python
DATA_DIR = "Check_mPMT_qualityCuts/"   # location of merged JSONs
BEAM_CONFIG = {
    "968 MeV": {"skip_channels": [], "color": "red"},
    "774 MeV": {"skip_channels": [0, 1, 14], "color": "orange"},
    ...
}
```

`skip_channels` lists HD indices that the beam cannot illuminate at that momentum
(geometric acceptance); they are omitted from that momentum's histogram.

## CLI

```bash
# Default: 968 MeV/c only
python3 overlay_normalize_final_v2.py

# All four momenta
python3 overlay_normalize_final_v2.py \
    --momenta "968 MeV" "774 MeV" "629 MeV" "483 MeV"

# Custom output directory
python3 overlay_normalize_final_v2.py \
    --momenta "968 MeV" "774 MeV" \
    --output-dir overlay_plots/
```

## Outputs

One PNG per hodoscope channel in the output directory:
```
charge_overlay_HD0.png
charge_overlay_HD1.png
...
charge_overlay_HD14.png
```

## Fit parameters extracted

The Gaussian μ (peak position in QDC a.u.) and σ are read off these plots and
hard-coded into `expectVsMeasured_final_v13_singleOnly.py` (arrays `Q_s`, `err_s`).

## Dependencies

```
numpy
scipy
matplotlib
json (stdlib)
```
