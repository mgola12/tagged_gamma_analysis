import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# -----------------------------------------
# INPUT ARRAYS
# -----------------------------------------

# Total spills per beam momentum (from count_spills.py)
# 483: runs 1786+1788+1790+1792+1794+1796
# 629: runs 1804+1806+1808+1810
# 774: runs 1812+1820+1825  (run 1814 missing from EOS)
# 968: runs 1827+1829+1831
SPILLS = {483: 3881, 629: 2851, 774: 1137, 968: 1982}

beam_momenta = [
    483, 483, 483, 483, 483, 483, '', '',
    629, 629, 629, 629, 629, 629, 629, 629, 629, '', '',
    774, 774, 774, 774, 774, 774, 774, 774, 774, 774, 774, '', '',
    968, 968, 968, 968, 968, 968, 968, 968, 968, 968, 968, 968, 968, '', ''
]

hd_labels = [
    "HD5 (110 MeV)", "HD6 (156 MeV)", "HD7 (194 MeV)", "HD12 (52 MeV)", "HD13 (121 MeV)", "HD14 (162 MeV)", "", "",
    "HD3 (131 MeV)", "HD4 (211 MeV)", "HD5 (268 MeV)", "HD6 (312 MeV)", "HD7 (349 MeV)", "HD11 (159 MeV)", "HD12 (227 MeV)", "HD13 (278 MeV)", "HD14 (318 MeV)", "", "",
    "HD2 (185 MeV)", "HD3 (297 MeV)", "HD4 (373 MeV)", "HD5 (428 MeV)", "HD6 (471 MeV)", "HD7 (507 MeV)", "HD10 (230 MeV)", "HD11 (322 MeV)", "HD12 (388 MeV)", "HD13 (437 MeV)", "HD14 (477 MeV)", "", "",
    "HD1 (235 MeV)", "HD2 (404 MeV)", "HD3 (509 MeV)", "HD4 (580 MeV)", "HD5 (634 MeV)", "HD6 (675 MeV)", "HD7 (709 MeV)", "HD9 (311 MeV)", "HD10 (444 MeV)", "HD11 (532 MeV)", "HD12 (595 MeV)", "HD13 (643 MeV)", "HD14 (680 MeV)", "", ""
]

tagged_counts = [
    245643, 122818, 72910, 691877, 186850, 97443, np.nan, np.nan,
    285951, 107970, 60658, 38427, 26011, 193626, 89940, 49324, 31101, np.nan, np.nan,
    179651, 64235, 32948, 20913, 13900, 9865, 113104, 49839, 28335, 17446, 11459, np.nan, np.nan,
    516961, 137984, 67964, 39094, 26727, 19054, 13881, 253153, 102289, 55185, 34283, 22478, 15769, np.nan, np.nan
]

# -----------------------------------------
# COLOR ASSIGNMENT BASED ON MOMENTUM
# -----------------------------------------

colors = []
for bm in beam_momenta:
    if bm == 483:
        colors.append("tab:blue")
    elif bm == 629:
        colors.append("tab:green")
    elif bm == 774:
        colors.append("tab:orange")
    elif bm == 968:
        colors.append("tab:red")
    else:
        colors.append("darkgray")  # blanks ("")
        

# -----------------------------------------
# CONVERT TO PHOTONS PER SPILL
# -----------------------------------------

photons_per_spill = []
for bm, count in zip(beam_momenta, tagged_counts):
    if bm == '' or np.isnan(count):
        photons_per_spill.append(np.nan)
    else:
        photons_per_spill.append(count / SPILLS[bm])

# -----------------------------------------
# PLOT
# -----------------------------------------

x = np.arange(len(hd_labels))
fig, ax = plt.subplots(figsize=(16, 6))

bars = ax.bar(x, photons_per_spill, color=colors)

ax.set_ylabel("Tagged γ per spill", fontsize=20, fontweight='bold')
ax.set_title("")

ax.set_yscale('log')
ax.set_ylim(0.1, 1e4)

# -----------------------------------------
# VALUES ABOVE BARS = photons per spill
# -----------------------------------------

for i, val in enumerate(photons_per_spill):
    if not np.isnan(val):
        ax.text(i, val * 1.15, f"{val:.1f}",
                ha='center', va='bottom', fontsize=8, rotation=90)

# -----------------------------------------
# X-AXIS LABELS
# -----------------------------------------

ax.set_xticks(x)
ax.set_xticklabels(hd_labels, rotation=90, fontsize=9)

ax.set_xlabel("Tagged γ Energy for Each Hodoscope Channel",
              fontsize=20, fontweight='bold')

# -----------------------------------------
# REMOVE SECONDARY AXIS
# -----------------------------------------

ax.yaxis.grid(True, linestyle='--', alpha=0.6)

# -----------------------------------------
# LEGEND
# -----------------------------------------

legend_patches = [
    mpatches.Patch(color="tab:blue", label="483 MeV/c"),
    mpatches.Patch(color="tab:green", label="629 MeV/c"),
    mpatches.Patch(color="tab:orange", label="774 MeV/c"),
    mpatches.Patch(color="tab:red", label="968 MeV/c"),
]
legend = ax.legend(
    handles=legend_patches,
    title="Beam Momentum:",
    loc="upper right",
    fontsize=12,
    frameon=False
)
legend.get_title().set_fontweight('bold')
legend.get_title().set_fontsize(13)

# -----------------------------------------
# INSET SUMMARY TABLE
# -----------------------------------------
# Summed events per calibrated beam momentum and photon energy range.

table_data = [
    ["$p_{\\mathrm{beam}}$ (MeV/$c$)", "Total Events", r"$E_\gamma$ Range (MeV)"],
    ["483",   "7172235",  " 52 – 194"],
    ["629",   "4774947",  "131 – 349"],
    ["774",   "2499834",  "185 – 507"],
    ["968",   "3087574",  "235 – 709"],
]

row_colors = [
    ["#f0f0f0", "#f0f0f0", "#f0f0f0"],   # header row
    ["#cce5ff", "#cce5ff", "#cce5ff"],   # 483 MeV/c (blue tint)
    ["#d4edda", "#d4edda", "#d4edda"],   # 629 MeV/c (green tint)
    ["#ffecd2", "#ffecd2", "#ffecd2"],   # 774 MeV/c (orange tint)
    ["#f8d7da", "#f8d7da", "#f8d7da"],   # 968 MeV/c (red tint)
]

# Place the table in the upper-left area of the axes
tbl = ax.table(
    cellText  = table_data,
    cellColours = row_colors,
    cellLoc   = "center",
    loc       = "upper left",
    bbox      = [0.34, 0.73, 0.32, 0.25],   # [x0, y0, width, height] in axes coords — raised to clear orange bars
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(10)

# Bold header row
for col in range(3):
    tbl[(0, col)].set_text_props(fontweight='bold')

# Bold first column (momentum values)
for row in range(1, 5):
    tbl[(row, 0)].set_text_props(fontweight='bold')

plt.tight_layout()
plt.savefig("taggedPhotonFlux_v2.png", dpi=150, bbox_inches="tight")
plt.show()
