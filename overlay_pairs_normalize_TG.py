#!/usr/bin/env python3
"""
Overlay HD-pair charge distributions with normalization and Gaussian fits.

Examples
--------
Run with all four beam momenta (default):
    python3 overlay_pairs_normalize.py

Run with a subset:
    python3 overlay_pairs_normalize.py --momenta "1000 MeV" "800 MeV"
"""

import argparse
import json
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit


# ---------- Gaussian ----------
def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


# ---------- Beam configuration ----------
BEAM_CONFIG = {
    "1000 MeV": {
        "file": (
            "/Users/mgola/Downloads/Tagged_gamma_16072026/Check_mPMT_qualityCuts/"
            "merged_pairs_968MeVc.json"
        ),
        "color": "tab:blue",
        "skip_pairs": set(),
    },
    "800 MeV": {
        "file": (
            "/Users/mgola/Downloads/Tagged_gamma_16072026/Check_mPMT_qualityCuts/"
            "merged_pairs_774MeVc.json"
        ),
        "color": "tab:orange",
        "skip_pairs": {"HD1+HD9"},
    },
    "650 MeV": {
        "file": (
            "/Users/mgola/Downloads/Tagged_gamma_16072026/Check_mPMT_qualityCuts/"
            "merged_pairs_629MeVc.json"
        ),
        "color": "tab:green",
        "skip_pairs": {"HD1+HD9", "HD2+HD10"},
    },
    "500 MeV": {
        "file": (
            "/Users/mgola/Downloads/Tagged_gamma_16072026/Check_mPMT_qualityCuts/"
            "merged_pairs_483MeVc.json"
        ),
        "color": "tab:red",
        "skip_pairs": {"HD1+HD9", "HD2+HD10", "HD3+HD11", "HD4+HD12"},
    },
}

# Corrected beam momentum display labels
DISPLAY_LABEL = {
    "1000 MeV": "968 MeV/c",
    "800 MeV":  "774 MeV/c",
    "650 MeV":  "629 MeV/c",
    "500 MeV":  "483 MeV/c",
}

# Momentum range for each pair (from hodoscope calculator)
PAIR_PRANGE = {
    "HD1+HD9":   (694, 727),
    "HD2+HD10":  (560, 585),
    "HD3+HD11":  (469, 490),
    "HD4+HD12":  (405, 422),
    "HD5+HD13":  (356, 371),
    "HD6+HD14":  (318, 331),
}

DEFAULT_MOMENTA = ["1000 MeV", "800 MeV", "650 MeV", "500 MeV"]

BIN_EDGES = np.linspace(0, 1.4e6, 500)


# Fit only bins above this charge to avoid the zero-spike pulling the Gaussian.
# The full distribution (including the spike and tail) is still displayed.
FIT_CHARGE_MIN = 50_000


def parse_args():
    parser = argparse.ArgumentParser(
        description="Overlay HD-pair charge distributions with normalization and Gaussian fits."
    )
    parser.add_argument(
        "--momenta",
        nargs="+",
        default=DEFAULT_MOMENTA,
        choices=list(BEAM_CONFIG.keys()),
        help=f"Beam momenta to include (default: all four).",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory for plots.",
    )
    return parser.parse_args()


def load_datasets(selected_labels):
    datasets = []
    for label in selected_labels:
        cfg = BEAM_CONFIG[label]
        with open(cfg["file"]) as f:
            data = json.load(f)
        datasets.append((label, data, cfg["color"], cfg["skip_pairs"]))
    return datasets


def fit_and_plot_pair(pair_key, datasets, output_dir):
    fig, ax = plt.subplots(figsize=(7, 5))

    legend_handles = []
    legend_labels  = []
    stats_lines    = []

    pmin, pmax = PAIR_PRANGE.get(pair_key, (0, 0))
    print(f"\n----- PAIR {pair_key}  ({pmin}–{pmax} MeV/c) -----")

    for label, data, color, skip_pairs in datasets:
        if pair_key in skip_pairs:
            continue

        charges = list(map(float, data.get(pair_key, {}).values()))
        if not charges:
            continue

        charges = np.array(charges)

        # ---------- Histogram ----------
        counts, edges = np.histogram(charges, bins=BIN_EDGES)
        N_total = np.sum(counts)
        if N_total == 0:
            continue

        # ---------- Normalization ----------
        counts_norm = counts / N_total
        x = 0.5 * (edges[:-1] + edges[1:])

        # ---------- Initial guesses (use only peak region to avoid zero-spike) ----------
        fit_mask = x >= FIT_CHARGE_MIN
        x_fit_region      = x[fit_mask]
        counts_fit_region = counts_norm[fit_mask]
        if counts_fit_region.sum() == 0:
            continue
        p0 = [
            np.max(counts_fit_region),
            float(np.median(charges[charges >= FIT_CHARGE_MIN])) if np.any(charges >= FIT_CHARGE_MIN) else 4e5,
            float(np.std(charges[charges  >= FIT_CHARGE_MIN]))   if np.any(charges >= FIT_CHARGE_MIN) else 5e4,
        ]

        # ---------- Gaussian fit (restricted to peak region) ----------
        fit_ok = False
        popt   = None
        try:
            popt, pcov = curve_fit(gauss, x_fit_region, counts_fit_region, p0=p0, maxfev=20000)
            A, mu, sigma = popt
            perr = np.sqrt(np.diag(pcov))
            mu_err, sigma_err = perr[1], perr[2]
            fit_ok = True
            print(
                f"  {label:>7}  mu = {mu:10.1f} +/- {mu_err:8.1f}   "
                f"sigma = {sigma:10.1f} +/- {sigma_err:8.1f}   N = {N_total}"
            )
        except RuntimeError:
            print(f"  {label:>7}  fit failed   N = {N_total}")
        except Exception as e:
            print(f"  {label:>7}  fit error: {e}")

        # ---------- Plot normalized histogram ----------
        n, b, patches = ax.hist(
            edges[:-1],
            bins=edges,
            weights=counts_norm,
            histtype="step",
            linewidth=1.5,
            color=color,
        )
        display_label = DISPLAY_LABEL.get(label, label)
        legend_handles.append(patches[0])
        legend_labels.append(display_label)

        # ---------- Plot fit ----------
        if fit_ok:
            yfit = gauss(x, *popt)
            ax.plot(x, yfit, "--", color=color)
            stats_lines.append(
                (f"{display_label}: μ={popt[1]:.1f}, σ={popt[2]:.1f}", color)
            )

    # ---------- Stats text ----------
    for i, (text_str, color) in enumerate(stats_lines):
        ax.text(
            0.99, 0.97 - 0.065 * i,
            text_str,
            transform=ax.transAxes,
            fontsize=15,
            color=color,
            ha="right",
            va="top",
            fontweight="bold",
        )

    # ---------- Legend ----------
    if legend_handles:
        header_handle = Line2D([0], [0], color="black", lw=1.5)
        handles = [header_handle] + legend_handles
        labels  = [f"PAIR {pair_key}"] + legend_labels

        leg = ax.legend(
            handles, labels,
            loc="lower right",
            frameon=True,
            prop={"size": 12, "weight": "bold"},
            handlelength=2.0,
            handletextpad=0.6,
            borderaxespad=0.8,
        )
        for txt in leg.get_texts():
            txt.set_ha("left")

    ax.set_xlabel("Total WCTE charge (a.u.)", fontsize=20, fontweight="bold")
    ax.set_ylabel("Normalized counts",        fontsize=20, fontweight="bold")
    ax.set_xlim(BIN_EDGES[0], BIN_EDGES[-1])
    fig.tight_layout()

    safe_key = pair_key.replace("+", "_")
    out_path = os.path.join(output_dir, f"overlay_charge_distribution_{safe_key}.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  ✅ Saved: {out_path}")


def main():
    args     = parse_args()
    selected = args.momenta

    output_dir = args.output_dir or (
        "overlaid_pairs_charge_plots_" + "_".join(m.replace(" ", "") for m in selected)
    )
    os.makedirs(output_dir, exist_ok=True)

    datasets = load_datasets(selected)

    # All pair keys present across selected datasets
    all_pair_keys = sorted(
        set().union(*[set(d) for _, d, _, _ in datasets]),
        key=lambda k: int(k.split("+")[0].replace("HD", ""))
    )

    print("\n=== Fitted parameters (mu, sigma) per HD pair ===")
    for pair_key in all_pair_keys:
        fit_and_plot_pair(pair_key, datasets, output_dir)

    print(f"\nOverlay plots saved to: {output_dir}/")
    print("Done.")


if __name__ == "__main__":
    main()
