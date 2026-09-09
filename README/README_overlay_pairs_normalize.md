# README — `overlay_pairs_normalize.py`

## Purpose

Same as `overlay_normalize_final_v2.py` but for **matched-pair coincidences**
(HD_i + HD_{i+8}).  Produces per-pair overlaid charge histograms with Gaussian
fits.  This is **Step 3b** of the analysis chain.

## Pairs and momentum ranges

| Pair | Positron momentum range (MeV/c) | Tagged photon energy (MeV) |
|------|----------------------------------|---------------------------|
| HD0 + HD8  | 759–792 | 176–209 |
| HD1 + HD9  | 694–727 | 241–274 |
| HD2 + HD10 | 560–585 | 383–408 |
| HD3 + HD11 | 469–490 | 478–499 |
| HD4 + HD12 | 405–422 | 546–563 |
| HD5 + HD13 | 356–371 | 597–612 |
| HD6 + HD14 | 318–331 | 637–650 |

## Algorithm

Same normalise-and-fit procedure as `overlay_normalize_final_v2.py` but:
- Input JSON: `merged_hd_pairs_<momentum>MeVc.json`
- Fit range starts at `charge > 50 000 a.u.` to exclude the zero-charge spike
  caused by events where one HD element fired but no photon was detected.

## CLI

```bash
# All four momenta
python3 overlay_pairs_normalize.py

# Subset
python3 overlay_pairs_normalize.py --momenta "968 MeV" "774 MeV"
```

## Outputs

One PNG per hodoscope pair:
```
charge_overlay_pair_HD0_HD8.png
charge_overlay_pair_HD1_HD9.png
...
charge_overlay_pair_HD6_HD14.png
```

## Dependencies

```
numpy
scipy
matplotlib
json (stdlib)
```
