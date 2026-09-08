# T10 Background / Calibration Analysis

Scripts for tagged-photon flux measurement, WCTE energy calibration, and
expected-vs-measured energy comparison.  These codes consume merged JSON files
produced by the upstream `Tagged_gamma_16072026` analysis chain.

---

## Files

| Script | Purpose |
|--------|---------|
| `taggedPhotonFlux_v2.py` | Tagged γ per-spill bar chart for all four momenta |
| `expectVsMeasured_final_v13_singleOnly.py` | Expected vs measured energy — single HD, three calibrations |
| `expectVsMeasured_final_v14_highE_only.py` | Expected vs measured energy — single HD, high-E calibration only *(new)* |
| `new_QDC_to_eV_from_Claude_v2.py` | QDC → MeV/c calibration from positron-beam data |
| `new_QDC_to_eV_from_Claude_v3.py` | Same with major + minor grid lines added *(new)* |

---

## `taggedPhotonFlux_v2.py`

Produces the tagged-photon flux bar chart (Fig. 6 in the paper).

**What it does:**
- Divides total tagged counts per hodoscope channel by the number of beam spills
  to obtain **tagged γ per spill**.
- Colour-codes bars by beam momentum (blue=483, green=629, orange=774, red=968).
- Adds per-bar value labels and a momentum summary table inset.
- Y-axis is logarithmic.

**Key inputs (hard-coded at the top of the script):**

```python
SPILLS = {483: 3881, 629: 2851, 774: 1137, 968: 1982}
```

> ⚠️  The 774 MeV/c spill count (1137) is a lower bound — run 1814 is missing
> from EOS and is not included.

**Output:** `taggedPhotonFlux_v2.png`

```bash
python3 taggedPhotonFlux_v2.py
```

---

## `new_QDC_to_eV_from_Claude_v2.py` / `v3`

Fits the WCTE total charge (QDC, a.u.) vs positron beam momentum (MeV/c) with
two models:

| Model | Formula |
|-------|---------|
| Linear (free) | p = m × Q + c |
| Quadratic (origin) | p = m₁ × Q + m₂ × Q² |

Data points are from positron calibration runs covering 150–1000 MeV/c.
Errors on QDC peak positions from Gaussian fits are used as weights.

**High-E calibration slope** (used downstream):
```
m_highE = 1.01933e-3  MeV per a.u.   (free linear fit, high-energy runs)
```

**v3** adds major + minor grid lines to the plot for easier reading.

**Output:** `calibration_linear_vs_quadratic.png`

```bash
python3 new_QDC_to_eV_from_Claude_v3.py   # recommended (has grid)
```

---

## `expectVsMeasured_final_v13_singleOnly.py`

Compares **expected** photon energies (from hodoscope geometry) with **measured**
WCTE energies (Gaussian-fit peak of total charge, converted via calibration) for
all single-HD channels at all four momenta.

Produces three output plots, one per calibration function:

| Tag | Calibration |
|-----|-------------|
| `_linear` | E = 1.0632×10⁻³ × Q − 29.44 |
| `_highE` | E = 1.01933×10⁻³ × Q |
| `_quad` | E = 9.023×10⁻⁴ × Q + 1.535×10⁻¹⁰ × Q² |

For each calibration a weighted linear fit is shown both free and through the origin.
Statistical and systematic (2.25%) errors are combined in quadrature.

**Output:** `wcte_measured_vs_expected_v13_singleOnly_measured_{linear,highE,quad}.png`

```bash
python3 expectVsMeasured_final_v13_singleOnly.py
```

---

## `expectVsMeasured_final_v14_highE_only.py` *(new)*

Simplified version of v13 that produces **a single clean plot** using only the
high-E calibration (E = 1.01933×10⁻³ × Q), which gives the best linearity over
the full 110–709 MeV range.

Improvements over v13:
- One output file instead of three.
- Reference diagonal (y = x) added for visual comparison.
- Y-axis label includes the calibration formula for self-documentation.

**Output:** `wcte_measured_vs_expected_v14_highE_only.png`

```bash
python3 expectVsMeasured_final_v14_highE_only.py
```

---

## Data flow

```
Tagged_gamma_16072026/
└── Check_mPMT_qualityCuts/
    ├── merged_hd_968MeVc.json   ──┐
    ├── merged_hd_774MeVc.json     │  → expectVsMeasured_*.py
    ├── merged_hd_629MeVc.json     │     (Gaussian fits done upstream by
    └── merged_hd_483MeVc.json   ──┘      overlay_normalize_final_v2.py)

Per-channel Gaussian fit results (mu, sigma)
    → hard-coded Q_s / err_s arrays in expectVsMeasured_*.py

Positron calibration runs (ROOT)
    → Gaussian-fit peak positions hard-coded in new_QDC_to_eV_from_Claude_v*.py
```

---

## Dependencies

```
numpy
scipy
matplotlib
```
