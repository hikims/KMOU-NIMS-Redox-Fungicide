import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# ==============================
# Path
# ==============================
BASE = "/Users/mkim/Desktop/Research/fungicide_jh"

DATA_PATH = os.path.join(
    BASE,
    "threshold_parameter",
    "summary",
    "threshold_estimates.txt"
)

OUT_DIR = os.path.join(BASE, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ==============================
# Publication style
# ==============================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 13,
    "axes.labelsize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 13,
    "axes.linewidth": 1.0,
    "lines.linewidth": 2.2,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "ps.fonttype": 42,
    "pdf.fonttype": 42,
})

# ==============================
# Load data
# columns:
# group factor dA_factor dF_factor alpha_factor A0_lower A0_upper threshold_mid
# ==============================
data = np.genfromtxt(
    DATA_PATH,
    comments="#",
    dtype=None,
    encoding=None
)

group = data["f0"]
factor = data["f1"]
A0_lower = data["f5"]
A0_upper = data["f6"]
threshold_mid = data["f7"]

# ==============================
# Plot settings
# ==============================
groups = [
    ("dA", r"$d_A$", "o"),
    ("dF", r"$d_F$", "s"),
    ("alpha", r"$\alpha$", "^"),
]

colors = {
    "dA": "#1f77b4",
    "dF": "#ff7f0e",
    "alpha": "#2ca02c",
}

# Maximum tested dose used in threshold search
A0_MAX = 3.0

# y-position for not-found crosses and annotation
NOT_FOUND_Y = 3.08
TEXT_Y = 3.35

# ==============================
# Figure
# ==============================
fig, ax = plt.subplots(figsize=(6.8, 4.8))
fig.patch.set_facecolor("white")

# ==============================
# Plot
# ==============================
for g, label, marker in groups:
    mask = group == g

    x = factor[mask]
    y = threshold_mid[mask]
    y_low = A0_lower[mask]
    y_high = A0_upper[mask]

    order = np.argsort(x)
    x = x[order]
    y = y[order]
    y_low = y_low[order]
    y_high = y_high[order]

    # Found threshold:
    # finite values below A0_MAX are plotted as circles/squares/triangles with lines.
    found = (~np.isnan(y)) & (y < A0_MAX)

    # Not found:
    # NaN values or values at/above the tested maximum are plotted as crosses.
    not_found = ~found
    
    found_factors = set(np.round(x[found], 6).tolist())
    not_found = not_found & np.array(
        [round(float(xx), 6) not in found_factors for xx in x],
        dtype=bool,
    )

    # ------------------------------
    # Plot found thresholds
    # ------------------------------
    if np.any(found):
        yerr = np.vstack([
            y[found] - y_low[found],
            y_high[found] - y[found]
        ])

        ax.errorbar(
            x[found],
            y[found],
            yerr=yerr,
            marker=marker,
            markersize=6.8,
            linewidth=2.2,
            elinewidth=1.3,
            capsize=4,
            capthick=1.2,
            color=colors[g],
            markerfacecolor=colors[g],
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=label,
            zorder=4,
        )

    # ------------------------------
    # Plot not-found cases
    # ------------------------------
    if np.any(not_found):
        ax.scatter(
            x[not_found],
            np.full(np.sum(not_found), NOT_FOUND_Y),
            marker="x",
            s=85,
            linewidths=2.0,
            color=colors[g],
            zorder=5,
        )

# ==============================
# Reference line and annotation
# ==============================
ax.axhline(
    A0_MAX,
    linestyle="--",
    linewidth=1.0,
    color="0.45",
    zorder=1,
)

ax.text(
    0.31,
    TEXT_Y,
    "not found within tested range",
    fontsize=12,
    color="0.20",
    ha="left",
    va="center",
)

# ==============================
# Axes
# ==============================
ax.set_xlabel("Multiplicative factor")
ax.set_ylabel(r"Estimated threshold $A_{0,\mathrm{emp}}^*$")

ax.set_xlim(0.20, 2.05)
ax.set_ylim(1.0, 3.5)

ax.set_xticks([0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00])
ax.set_yticks(np.arange(1.0, 3.51, 0.5))

ax.grid(
    True,
    which="major",
    linestyle="-",
    linewidth=0.5,
    color="0.88",
    zorder=0,
)

for spine in ax.spines.values():
    spine.set_linewidth(1.0)

# ==============================
# Legend
# ==============================
leg = ax.legend(
    loc="upper right",
    bbox_to_anchor=(0.98, 0.72),
    frameon=True,
    fancybox=False,
    facecolor="white",
    edgecolor="0.75",
    handlelength=2.0,
    borderpad=0.45,
)

leg.get_frame().set_linewidth(0.7)

# ==============================
# Layout
# ==============================
fig.subplots_adjust(
    left=0.15,
    right=0.98,
    bottom=0.16,
    top=0.96,
)

# ==============================
# Save
# ==============================
png_path = os.path.join(OUT_DIR, "fig5.png")
eps_path = os.path.join(OUT_DIR, "fig5.eps")

fig.savefig(
    png_path,
    format="png",
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)

fig.savefig(
    eps_path,
    format="eps",
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)

plt.close(fig)

print("Saved:")
print(" -", png_path)
print(" -", eps_path)