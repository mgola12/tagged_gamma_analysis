# README — `new_QDC_to_eV_from_Claude_v3.py`

## Purpose

Updated version of `new_QDC_to_eV_from_Claude_v2.py` with **major and minor grid
lines** added to the calibration plot.  All physics and fit logic is identical to v2.

## Changes vs v2

```python
# Added to the plot section:
ax.grid(True, which='major', linestyle='--', linewidth=0.7, alpha=0.6)
ax.grid(True, which='minor', linestyle=':',  linewidth=0.4, alpha=0.4)
ax.minorticks_on()
```

## Models and coefficients

See [`README_new_QDC_to_eV_from_Claude_v2.md`](README_new_QDC_to_eV_from_Claude_v2.md)
for the full description.  Coefficients are identical between v2 and v3.

## Output

```
calibration_linear_vs_quadratic.png
```

## Usage

```bash
python3 new_QDC_to_eV_from_Claude_v3.py
```

## Dependencies

```
numpy
scipy
matplotlib
```

> **Recommended version** — use v3 for all final paper figures.
