# README — `taggedPhotonFlux_v2.py`

## Purpose

Produces the **tagged-photon flux bar chart** (Fig. 6 in the paper): tagged photon
counts per beam spill for each hodoscope channel, at all four beam momenta.

## What it does

1. Reads the total tagged-photon count per HD channel from hard-coded summary
   dictionaries (derived from the merged JSON files).
2. Divides by the number of beam spills to obtain **tagged γ per spill**.
3. Plots a grouped bar chart, colour-coded by momentum.
4. Adds a logarithmic y-axis, per-bar value labels, and a summary table inset.

## Key inputs (hard-coded)

```python
# Number of recorded beam spills per momentum
SPILLS = {483: 3881, 629: 2851, 774: 1137, 968: 1982}

# Tagged counts per HD channel per momentum (from merged JSONs)
COUNTS = { ... }
```

> ⚠️ The 774 MeV/c spill count (1137) excludes run 1814, which is missing from
> EOS.  The flux value for 774 MeV/c is therefore a **lower bound**.

## Output

```
taggedPhotonFlux_v2.png
```

## Usage

```bash
python3 taggedPhotonFlux_v2.py
```

## Dependencies

```
numpy
matplotlib
```
