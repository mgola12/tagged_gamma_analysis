"""
expectVsMeasured_final_v14_highE_only.py

Plots measured vs expected photon energy using the high-E calibration only:
    E_meas [MeV] = 1.01933e-3 × Q [a.u.]

Single hodoscope channels only (no pairs). A single linear fit is performed
through the origin (forced) and a free fit (slope + intercept) are both shown.

Calibration source: WCTE positron-beam calibration runs (150–1000 MeV/c),
  fit from new_QDC_to_eV_from_Claude_v2.py.

Outputs
-------
  wcte_measured_vs_expected_v14_highE_only.png

Data (mu from Gaussian fits to per-channel charge distributions)
----------------------------------------------------------------
  968 MeV/c → HD1,HD2,HD3,HD4,HD5,HD6,HD7,HD9,HD10,HD11,HD12,HD13,HD14
  774 MeV/c → HD2,HD3,HD4,HD5,HD6,HD7,HD10,HD11,HD12,HD13,HD14
  629 MeV/c → HD3,HD4,HD5,HD6,HD7,HD11,HD12,HD13,HD14
  483 MeV/c → HD5,HD6,HD7,HD13,HD14
"""

import numpy as np
import matplotlib.pyplot as plt

# ── High-E calibration: E [MeV] = m * Q [a.u.] ──────────────────────────────
# Slope derived from positron-beam calibration (intercept forced to zero)
CAL_SLOPE = 1.01933e-3   # MeV per a.u.

def cal_highE(Q):
    return CAL_SLOPE * Q

# ── Raw Gaussian-fit charge means (a.u.) — Single HD channels ────────────────
# Ordering follows the beamline hodoscope layout (front HD0–HD6, back HD8–HD14)
Q_s = {
    968: np.array([194403.5, 345129.2, 438007.2, 502805.5, 550406.6, 586695.5, 615371.6,
                   262638.2, 380635.0, 458709.0, 515420.1, 558904.7, 591555.1]),
    774: np.array([161225.8, 260361.6, 327921.4, 376591.0, 414649.3, 444823.7,
                   200146.5, 281855.3, 340763.9, 384903.6, 419710.1]),
    629: np.array([118523.5, 188992.0, 239067.0, 277623.4, 308417.0,
                   141729.9, 202508.2, 247751.3, 283217.1]),
    483: np.array([100098.8, 140723.1, 172736.1, 109444.1, 146348.3]),
}

# ── Expected photon energies [MeV] ───────────────────────────────────────────
exp_s = {
    968: np.array([238.5, 387, 485, 554, 605.5, 645.5, 677,
                   297, 421, 505.5, 567.5, 614, 651]),
    774: np.array([193, 291, 360, 411.5, 451.5, 483,
                   227, 311.5, 373.5, 420, 457]),
    629: np.array([146, 215, 266.5, 306.5, 338,
                   166.5, 228.5, 275, 312]),
    483: np.array([120.5, 160.5, 192, 129, 166]),
}

# ── Statistical errors on mu [a.u.], converted to MeV via CAL_SLOPE ─────────
err_s = {
    968: np.array([2502.3, 612.1, 646.4, 514.0, 587.2, 633.0, 717.0,
                   842.2, 687.3, 633.7, 558.9, 580.6, 632.7]) * CAL_SLOPE,
    774: np.array([1196.4, 479.7, 456.8, 407.1, 402.7, 442.8,
                   726.8, 473.0, 483.4, 381.1, 386.8])        * CAL_SLOPE,
    629: np.array([1007.5, 370.6, 313.4, 332.0, 304.9,
                   632.0, 352.9, 317.6, 326.8])                * CAL_SLOPE,
    483: np.array([495.7, 258.8, 233.4, 418.0, 250.1])        * CAL_SLOPE,
}

# ── Build calibrated arrays ───────────────────────────────────────────────────
MOMENTA   = [968, 774, 629, 483]
MOM_COLOR = {968: 'tab:blue', 774: 'tab:green', 629: 'tab:orange', 483: 'tab:red'}

expected_all  = np.concatenate([exp_s[bm]             for bm in MOMENTA])
measured_all  = np.concatenate([cal_highE(Q_s[bm])    for bm in MOMENTA])
err_stat_all  = np.concatenate([err_s[bm]             for bm in MOMENTA])
syst_all      = 0.0225 * measured_all               # 2.25% systematic
err_total_all = np.sqrt(err_stat_all**2 + syst_all**2)
labels_all    = np.concatenate([[bm] * len(exp_s[bm]) for bm in MOMENTA])

# ── Linear fits ──────────────────────────────────────────────────────────────
w = 1.0 / err_stat_all**2

# Free fit: y = slope * x + intercept
coeffs, cov    = np.polyfit(expected_all, measured_all, 1, w=1.0 / err_stat_all, cov=True)
slope, intercept = coeffs
slope_err      = np.sqrt(cov[0, 0])
intercept_err  = np.sqrt(cov[1, 1])

# Fit through origin: y = slope0 * x
slope0     = np.sum(w * expected_all * measured_all) / np.sum(w * expected_all**2)
slope0_err = 1.0 / np.sqrt(np.sum(w * expected_all**2))

print("=" * 65)
print("High-E calibration only  [E = 1.01933e-3 × Q]")
print(f"  Free   : slope = {slope:.4f} ± {slope_err:.4f}  "
      f"intercept = {intercept:.2f} ± {intercept_err:.2f}")
print(f"  Origin : slope = {slope0:.4f} ± {slope0_err:.4f}")
print("=" * 65)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 7))

bm_handles = []
for bm in MOMENTA:
    mask = labels_all == bm
    eb = ax.errorbar(
        expected_all[mask], measured_all[mask],
        yerr=err_total_all[mask],
        fmt='o', color=MOM_COLOR[bm],
        ecolor='gray', capsize=3, markersize=7,
        label=f'{bm} MeV/c',
    )
    bm_handles.append(eb)

x_line = np.linspace(min(expected_all) * 0.95, max(expected_all) * 1.02, 300)
sign   = '+' if intercept >= 0 else '−'

line_free, = ax.plot(
    x_line, slope * x_line + intercept, 'r-', linewidth=2,
    label=(f'Free fit: y = ({slope:.3f}±{slope_err:.3f})x '
           f'{sign} ({abs(intercept):.1f}±{intercept_err:.1f})'),
)
line_origin, = ax.plot(
    x_line, slope0 * x_line, 'b--', linewidth=2,
    label=f'Origin fit: y = ({slope0:.3f}±{slope0_err:.3f})x',
)

# Reference diagonal
ax.plot(x_line, x_line, 'k:', linewidth=1, alpha=0.4, label='y = x')

# Legends
leg1 = ax.legend(
    handles=bm_handles, fontsize=16, loc='upper left',
    title='Beam Momentum', title_fontsize=16,
    frameon=True, markerscale=1.4,
    prop={'size': 16, 'weight': 'bold'},
)
leg1.get_title().set_fontweight('bold')
ax.add_artist(leg1)

ax.legend(
    handles=[line_free, line_origin],
    fontsize=13, loc='lower right',
    frameon=True, prop={'size': 13, 'weight': 'bold'},
)

ax.set_xlabel("Expected Photon Energy (MeV)", fontsize=18, fontweight='bold')
ax.set_ylabel("WCTE Measured Energy (MeV)\n[High-E calibration: E = 1.01933×10⁻³ × Q]",
              fontsize=16, fontweight='bold')
ax.tick_params(axis='both', labelsize=14)
ax.grid(True, linestyle='--', alpha=0.35)

plt.tight_layout()
fname = "wcte_measured_vs_expected_v14_highE_only.png"
plt.savefig(fname, dpi=150)
plt.close()
print(f"Saved: {fname}")
