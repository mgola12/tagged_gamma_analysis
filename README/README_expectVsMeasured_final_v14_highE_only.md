# README — `expectVsMeasured_final_v14_highE_only.py`

## Purpose

Simplified, single-output version of `expectVsMeasured_final_v13_singleOnly.py`.
Uses **only the high-energy calibration** (E = 1.01933×10⁻³ × Q), which gives
the best linearity across the full 110–709 MeV tagged-photon energy range.

## Why high-E only?

The high-E calibration is fitted from positron beam runs at ≥362 MeV/c, covering
the energy range of interest for the tagged photon paper.  The full-range linear
and quadratic calibrations include low-energy runs and can introduce systematic
biases at the high-energy end.

## Improvements over v13

| Feature | v13 | v14 |
|---------|-----|-----|
| Number of output plots | 3 | 1 |
| Reference diagonal (y = x) | ✗ | ✓ |
| Y-axis formula label | ✗ | ✓ |
| Grid lines | ✗ | ✓ |

## Calibration used

```
E_measured [MeV] = 1.01933 × 10⁻³ × Q [a.u.]
```

No intercept (origin-forced), fitted from positron calibration data at
362–968 MeV/c.

## Input data (hard-coded)

```python
# Per momentum, per single-HD: (Q_peak [a.u.], err [a.u.], E_expected [MeV])
Q_s   = {968: [...], 774: [...], 629: [...], 483: [...]}
err_s = {968: [...], 774: [...], 629: [...], 483: [...]}
exp_s = {968: [...], 774: [...], 629: [...], 483: [...]}
```

## Output

```
wcte_measured_vs_expected_v14_highE_only.png
```

## Usage

```bash
python3 expectVsMeasured_final_v14_highE_only.py
```

## Dependencies

```
numpy
scipy
matplotlib
```
