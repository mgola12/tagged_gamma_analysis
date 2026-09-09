#!/usr/bin/env python3
"""
Overlay HD charge distributions with normalization and Gaussian fits.

By default this only processes the 1000 MeV beam file. To include more
beam momenta later, either pass --momenta on the command line or change
the DEFAULT_MOMENTA list below.

Examples
--------
Run with just 1000 MeV (default):
    python3 overlay_normalize_final_v2.py

Run with 1000 MeV and 800 MeV overlaid:
    python3 overlay_normalize_final_v2.py --momenta "1000 MeV" "800 MeV"

Run with all four:
    python3 overlay_normalize_final_v2.py --momenta "1000 MeV" "800 MeV" "650 MeV" "500 MeV"
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
# To add a new beam momentum, just add a new entry here.
BEAM_CONFIG = {
    "1000 MeV": {
        "file": (
            "/Users/mgola/Downloads/Tagged_gamma_16072026/Check_mPMT_qualityCuts/"
            "merged_hd_968MeVc.json"
        ),
        "color": "tab:blue",
        "skip_channels": set(),
    },
    "800 MeV": {
        "file": (
            "/Users/mgola/Downloads/Tagged_gamma_16072026/Check_mPMT_qualityCuts/"
            "merged_hd_774MeVc.json"
        ),
        "color": "tab:orange",
        "skip_channels": {"HD1", "HD9"},
    },
    "650 MeV": {
        "file": (
            "/Users/mgola/Downloads/Tagged_gamma_16072026/Check_mPMT_qualityCuts/"
            "merged_hd_629MeVc.json"
        ),
        "color": "tab:green",
        "skip_channels": {"HD1", "HD2", "HD9", "HD10"},
    },
    "500 MeV": {
        "file": (
            "/Users/mgola/Downloads/Tagged_gamma_16072026/Check_mPMT_qualityCuts/"
            "merged_hd_483MeVc.json"
        ),
        "color": "tab:red",
        "skip_channels": {"HD1", "HD2", "HD3", "HD4", "HD9", "HD10", "HD11", "HD12"},
    },
}

# Corrected beam momentum values shown on plots (legend + stats text only).
# Internal keys above (file paths, --momenta, skip_channels) are unaffected.
DISPLAY_LABEL = {
    "1000 MeV": "968 MeV/c",
    "800 MeV": "774 MeV/c",
    "650 MeV": "629 MeV/c",
    "500 MeV": "483 MeV/c",
}

# Change this if you'd rather hard-code the default selection instead of
# always passing --momenta on the command line.
DEFAULT_MOMENTA = ["1000 MeV"]

BIN_EDGES = np.linspace(0, 1.4e6, 500)



def parse_args():
    parser = argparse.ArgumentParser(
        description="Overlay HD charge distributions with normalization and Gaussian fits."
    )
    parser.add_argument(
        "--momenta",
        nargs="+",
        default=DEFAULT_MOMENTA,
        choices=list(BEAM_CONFIG.keys()),
        help=f"Beam momenta to include (default: {DEFAULT_MOMENTA}).",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Output directory for plots. If omitted, an automatic name is "
            "built from the selected momenta."
        ),
    )
    return parser.parse_args()


def load_datasets(selected_labels):
    """Load only the JSON files needed for the selected momenta."""
    datasets = []
    for label in selected_labels:
        cfg = BEAM_CONFIG[label]
        with open(cfg["file"]) as f:
            data = json.load(f)
        datasets.append((label, data, cfg["color"], cfg["skip_channels"]))
    return datasets


def fit_and_plot_hd(hd_key, datasets, output_dir):
    fig, ax = plt.subplots(figsize=(7, 5))

    legend_handles = []
    legend_labels = []
    stats_lines = []

    print(f"\n----- ELEMENT {hd_key} -----")

    for label, data, color, skip_channels in datasets:
        if hd_key in skip_channels:
            continue

        charges = list(map(float, data.get(hd_key, {}).values()))
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

        # ---------- Initial guesses ----------
        p0 = [
            np.max(counts_norm),
            float(np.mean(charges)),
            float(np.std(charges)) if np.std(charges) > 0 else 1.0,
        ]

        # ---------- Gaussian Fit ----------
        fit_ok = False
        popt = None
        try:
            popt, pcov = curve_fit(gauss, x, counts_norm, p0=p0, maxfev=20000)
            A, mu, sigma = popt
            perr = np.sqrt(np.diag(pcov))
            mu_err, sigma_err = perr[1], perr[2]
            fit_ok = True
            print(
                f"{label:>7}  mu = {mu:10.1f} +/- {mu_err:8.1f}   "
                f"sigma = {sigma:10.1f} +/- {sigma_err:8.1f}"
            )
        except RuntimeError:
            print(f"{label:>7}  fit failed")
        except Exception as e:
            print(f"{label:>7}  fit error: {e}")

        # ---------- Plot normalized histogram ----------
        n, b, patches = plt.hist(
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
            plt.plot(x, yfit, "--", color=color)

            stats_lines.append((f"{display_label}: μ={popt[1]:.1f}, σ={popt[2]:.1f}", color))

    # ---------- Stats text box, placed INSIDE the axes (top area) ----------
    # Right-aligned to match the right edge of the legend box below.
    for i, (text_str, color) in enumerate(stats_lines):
        ax.text(
            0.99,
            0.97 - 0.065 * i,
            text_str,
            transform=ax.transAxes,
            fontsize=15,
            color=color,
            ha="right",
            va="top",
            fontweight="bold",
        )

    # ---------- Legend, placed INSIDE the axes (lower right) ----------
    if legend_handles:
        header_handle = Line2D([0], [0], color="black", lw=1.5)
        handles = [header_handle] + legend_handles
        labels = [f"ELEMENT {hd_key}"] + legend_labels

        leg = ax.legend(
            handles,
            labels,
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
    ax.set_ylabel("Normalized counts", fontsize=20, fontweight="bold")
    ax.set_xlim(BIN_EDGES[0], BIN_EDGES[-1])
    fig.tight_layout()

    out_path = os.path.join(output_dir, f"overlay_charge_distribution_{hd_key}.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    selected = args.momenta

    output_dir = args.output_dir or (
        "overlaid_hd_charge_plots_with_fits_" + "_".join(m.replace(" ", "") for m in selected)
    )
    os.makedirs(output_dir, exist_ok=True)

    datasets = load_datasets(selected)

    # All HD keys present across just the selected datasets
    all_hd_keys = sorted(set().union(*[set(d) for _, d, _, _ in datasets]))

    print("\n=== Fitted parameters (mu, sigma) per HD ===")
    for hd_key in all_hd_keys:
        fit_and_plot_hd(hd_key, datasets, output_dir)

    print(f"\nOverlay plots saved to: {output_dir}")
    print("Done - normalized histograms and Gaussian fits generated.")


if __name__ == "__main__":
    main()
