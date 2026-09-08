# README — `expectVsMeasured_final_v13_singleOnly.py`

## Purpose

Compares **expected** tagged photon energies (from hodoscope geometry) with
**measured** WCTE energies (from Gaussian fits to the total charge distribution)
for all single-HD channels at all four beam momenta.

Produces three separate output plots, one per calibration function.

## Calibrations compared

| Label | Formula | Notes |
|-------|---------|-------|
| `linear` | E = 1.0632×10⁻³ × Q − 29.44 MeV | Full range, free linear |
| `highE` | E = 1.01933×10⁻³ × Q | High-energy runs, origin-forced |
| `quad` | E = 9.023×10⁻⁴ × Q + 1.535×10⁻¹⁰ × Q² | Quadratic, origin-forced |

## Input data (hard-coded arrays)

```python
# Beam momenta
momenta = [968, 774, 629, 483]   # MeV/c

# Single-HD channels: Gaussian fit peak positions Q [a.u.] per (momentum, channel)
Q_s     = { ... }    # from overlay_normalize_final_v2.py fits

# Statistical uncertainties on Q_s
err_s   = { ... }

# Expected energies from hodoscope geometry
exp_s   = { ... }
```

## Fits performed

For each calibration, the script fits:
- **Free linear**: E_measured = a × E_expected + b
- **Origin-forced linear**: E_measured = a × E_expected

Both use χ²-weighted least squares.  Statistical and systematic (2.25%) errors
are combined in quadrature for the weights.

## Outputs

```
wcte_measured_vs_expected_v13_singleOnly_measured_linear.png
wcte_measured_vs_expected_v13_singleOnly_measured_highE.png
wcte_measured_vs_expected_v13_singleOnly_measured_quad.png
```

## Usage

```bash
python3 expectVsMeasured_final_v13_singleOnly.py
```

## Dependencies

```
numpy
scipy
matplotlib
```

## See also

- `expectVsMeasured_final_v14_highE_only.py` — single clean plot using highE only
- `new_QDC_to_eV_from_Claude_v3.py` — source of calibration coefficients
