"""
new_QDC_to_eV_from_Claude_v3.py

Converts WCTE total integrated charge (QDC, in arbitrary units) to beam momentum
(MeV/c) using data from positron-beam calibration runs.

Two fit models are compared:
  • Linear (free):      p [MeV/c] = m × Q + c
  • Quadratic (origin): p [MeV/c] = m₁ × Q + m₂ × Q²

The resulting calibration coefficients are used as the charge-to-energy conversion
in expectVsMeasured_final_v14_highE_only.py (and v13).

Changes vs v2
-------------
  • Major and minor grid lines added to the plot for easier reading.
  • Minor tick marks enabled on both axes.

Output
------
  calibration_linear_vs_quadratic.png
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ── Data: positron beam calibration points ────────────────────────────────────
# beam_momenta [MeV/c] are the known beam settings.
# qdc_means    [a.u.]  are the Gaussian-fit peak positions of the WCTE total
#                       charge distribution for each setting.
# dc_mean_errors       are the statistical uncertainties on qdc_means.
beam_momenta = np.array([150, 160, 180, 190, 192.65, 202.65, 212.65,
                          362.13, 434.82, 628.82, 774.38, 968.41])

qdc_means = np.array([
    173052.2, 179251.2, 185473.1, 203765.4,
    209774.6, 224115.6, 229145.7, 365814.0,
    440663.0, 628137.0, 758746.0, 930547.0
])

dc_mean_errors = np.array([
    166.3, 166.6, 163.2, 189.0,
    174.7, 188.3, 188.7, 200.0,
    205.0, 250.0, 226.0, 230.0
])

# ── Fit models ────────────────────────────────────────────────────────────────
def linear_free(x, m, c):
    return m * x + c

def quadratic_origin(x, m1, m2):
    return m1 * x + m2 * x**2

# Free linear fit
popt_lin, pcov_lin   = curve_fit(linear_free, qdc_means, beam_momenta,
                                  sigma=dc_mean_errors, absolute_sigma=True)
m_lin, c_lin         = popt_lin
m_lin_err, c_lin_err = np.sqrt(np.diag(pcov_lin))

# Quadratic fit through origin
popt_quad, pcov_quad   = curve_fit(quadratic_origin, qdc_means, beam_momenta,
                                    sigma=dc_mean_errors, absolute_sigma=True,
                                    p0=[1e-3, 1e-10])
m1, m2               = popt_quad
m1_err, m2_err       = np.sqrt(np.diag(pcov_quad))

print("=" * 65)
print(f"Linear fit (free):  m = {m_lin:.5e} ± {m_lin_err:.3e}  "
      f"({m_lin_err / m_lin * 100:.2f}%)")
print(f"                    c = {c_lin:.3f} ± {c_lin_err:.3f}")
print(f"Quadratic fit: m1 = {m1:.5e} ± {m1_err:.3e}")
print(f"               m2 = {m2:.5e} ± {m2_err:.3e}")
print("=" * 65)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))

x_fit = np.linspace(0, max(qdc_means) * 1.05, 500)
cmap  = plt.cm.get_cmap("tab10", len(beam_momenta))

sign = '+' if c_lin >= 0 else '−'
ax.plot(x_fit, linear_free(x_fit, m_lin, c_lin), 'k--', linewidth=1.5,
        label=f"Linear: m = {m_lin:.4e}\n         c = {sign}{abs(c_lin):.2f}")
ax.plot(x_fit, quadratic_origin(x_fit, m1, m2), 'r-', linewidth=1.5,
        label=f"Quadratic: m₁ = {m1:.3e}\n           m₂ = {m2:.3e}")

for i, (q, p, e) in enumerate(zip(qdc_means, beam_momenta, dc_mean_errors)):
    ax.errorbar(q, p, xerr=e, fmt='o', color=cmap(i),
                ecolor='black', elinewidth=1, capsize=3, label=f"{p:.0f} MeV/c")

ax.set_xlim(left=0)
ax.set_ylim(bottom=0)

# ── Grid: major + minor ───────────────────────────────────────────────────────
ax.grid(True, which='major', linestyle='--', linewidth=0.7, alpha=0.6)
ax.grid(True, which='minor', linestyle=':',  linewidth=0.4, alpha=0.4)
ax.minorticks_on()

ax.set_xlabel("Total WCTE Charge (a.u.)", fontsize=13, fontweight='bold')
ax.set_ylabel("Beam Momentum (MeV/c)",    fontsize=13, fontweight='bold')

ax.legend(
    fontsize=10,
    ncol=2,
    loc='upper left',
    markerscale=1.3,
    labelspacing=0.6,
    handletextpad=0.6,
    columnspacing=1.0,
    borderpad=0.8,
)

plt.tight_layout()
plt.savefig("calibration_linear_vs_quadratic.png", dpi=150)
plt.show()
