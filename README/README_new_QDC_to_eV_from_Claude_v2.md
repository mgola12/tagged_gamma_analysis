# README — `new_QDC_to_eV_from_Claude_v2.py`

## Purpose

Derives the QDC → MeV/c calibration by fitting two models to positron calibration
run data.  The resulting coefficients are used downstream by `expectVsMeasured_*.py`.

## Models fitted

| Model | Formula |
|-------|---------|
| Linear (free) | p = m × Q + c |
| Quadratic (origin) | p = m₁ × Q + m₂ × Q² |

Fits are weighted by the statistical errors on the QDC Gaussian peak positions.

## Calibration data points

12 positron beam settings covering 150–968 MeV/c.  The QDC peak positions are
obtained from Gaussian fits to the total WCTE charge distributions.

## Calibration coefficients (from v2)

```
Linear fit:
  m = 1.01933 × 10⁻³ MeV/a.u.   (slope)
  c = −29.44 MeV                 (intercept)

High-E linear (fitted to ≥362 MeV/c runs only, in v3):
  m_highE = 1.01933 × 10⁻³ MeV/a.u.   (≈ same; intercept ≈ 0)

Quadratic:
  m₁ = 9.023 × 10⁻⁴ MeV/a.u.
  m₂ = 1.535 × 10⁻¹⁰ MeV/a.u.²
```

## Output

```
calibration_linear_vs_quadratic.png
```

## Usage

```bash
python3 new_QDC_to_eV_from_Claude_v2.py
```

## Notes

- v3 (`new_QDC_to_eV_from_Claude_v3.py`) is identical but adds major + minor grid
  lines for easier reading.  Use v3 for final plots.

## Dependencies

```
numpy
scipy
matplotlib
```
