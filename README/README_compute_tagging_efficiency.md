# README — `compute_tagging_efficiency.py`

## Purpose

Computes the tagged-photon **tagging efficiency** for each hodoscope element
(single-HD and HD-pair), broken down into three components:

| Efficiency | Symbol | Definition |
|-----------|--------|------------|
| Hodoscope efficiency | ε_hodo | N_tagged / N_beam |
| Detection efficiency | ε_det | N_detected / N_tagged |
| Total efficiency | ε_total | N_detected / N_beam |

ε_hodo measures geometric acceptance; ε_det measures the probability that a
tagged photon produces a detectable WCTE signal above threshold.

## Algorithm

1. Loop over beam events (T0 ≥ 1, T2, NOT HC2).
2. For each event, find all in-time hodoscope hits (Group A: −50 to +100 ns;
   Group B: +80 to +200 ns relative to T0).
3. Check whether the WCTE total charge exceeds the detection threshold
   `μ − 3σ` derived from 774 MeV/c Gaussian fits.
4. Accumulate counts per channel to compute the three efficiencies.

## Thresholds (774 MeV/c)

The `μ − 3σ` thresholds are hard-coded from Gaussian fits at 774 MeV/c.
Update these if re-running at a different reference momentum.

## CLI

```bash
python3 compute_tagging_efficiency.py \
    --run_file /eos/.../WCTE_offline_R1827S0_VME_matched.root \
    --n_events 200000 \
    --run_label R1827_774MeVc
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--run_file` | *(required)* | Offline ROOT file |
| `--n_events` | all | Events to process |
| `--run_label` | `""` | Tag for output filenames |

## Outputs

```
tagging_efficiency_<label>.png     # bar chart per channel
tagging_efficiency_<label>.csv     # table: channel, N_beam, N_tag, N_det, efficiencies
```

## Dependencies

```
uproot >= 5
awkward >= 2
numpy
matplotlib
analysis_tools  (DataLoader)
```
