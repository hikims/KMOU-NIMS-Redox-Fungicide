import os
import warnings
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=UserWarning)

# ==============================
# Path
# ==============================
BASE = "/Users/mkim/Desktop/Research/fungicide_jh"

DATA_PATH = os.path.join(
    BASE,
    "threshold_parameter",
    "summary",
    "threshold_scan_results.txt"
)

OUT_DIR = os.path.join(BASE, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ==============================
# Config
# ------------------------------
# The paper plots the NORMALIZED aggregate burden
#     M_bar(T_obs) = M_tot(T_obs) / |Omega|,   |Omega| = 100 for Omega = [0,10]^2
# (paper notation: \overline{M}).
#
# Set OMEGA_AREA according to what column `total_F_t100` of the data file
# actually stores:
#   * column already holds the normalized average  -> OMEGA_AREA = 1.0
#   * column holds the raw integral  M_tot = \int F dx  -> OMEGA_AREA = 100.0
#
# Quick check: at the lowest dose (A0 = 1) the curves should land in the
# ~1e-2 .. 1e-1 range (close to the persistence plateau F_inf = 0.1).
# If they instead come out around ~1e0 .. 1e1, the column is the raw
# integral and OMEGA_AREA must be set to 100.0.
# ==============================
OMEGA_AREA = 1.0   # <-- verify against threshold_scan_results.txt

# ==============================
# Publication style
# ==============================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.labelsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "axes.linewidth": 1.0,
    "lines.linewidth": 2.0,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "ps.fonttype": 42,
    "pdf.fonttype": 42,
})

# ==============================
# Load data
# Actual column order in threshold_scan_results.txt:
#   0:tag  1:A0  2:dA_factor  3:dF_factor  4:alpha_factor
#   5:dA   6:dF  7:alpha
#   8:total_F_t20   9:area_t20
#   10:total_F_t50  11:area_t50
#   12:total_F_t100 13:area_t100
# The empirical threshold itself is defined through the infected-area
# criterion (Section: numerical method); these figures show the
# normalized aggregate burden M_bar(T_obs) as a complementary measure only.
# ==============================
data = np.loadtxt(DATA_PATH, comments="#", dtype=str)

print("Data shape:", data.shape)

tag = data[:, 0]
A0 = data[:, 1].astype(float)

dA_factor = data[:, 2].astype(float)
dF_factor = data[:, 3].astype(float)
alpha_factor = data[:, 4].astype(float)

# Raw column 12 = total_F_t100 (T_obs = 100).
total_F_raw = data[:, 12].astype(float)

# Normalized aggregate burden  M_bar(T_obs) = M_tot(T_obs) / |Omega|.
Mbar_tobs = total_F_raw / OMEGA_AREA

# ==============================
# Helper functions
# ==============================
def positive_only(y):
    """Replace non-positive values by NaN for log-scale plotting."""
    return np.where(y > 0, y, np.nan)


def unique_by_A0(x, y):
    """Sort by A0 and remove duplicated A0 values."""
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    x_unique, unique_idx = np.unique(x, return_index=True)
    y_unique = y[unique_idx]

    return x_unique, y_unique


def get_mask(group_name, fac):
    if group_name == "dA":
        return np.char.startswith(tag, "dA_") & np.isclose(dA_factor, fac)

    if group_name == "dF":
        return np.char.startswith(tag, "dF_") & np.isclose(dF_factor, fac)

    if group_name == "alpha":
        return np.char.startswith(tag, "alpha_") & np.isclose(alpha_factor, fac)

    raise ValueError(f"Unknown group_name: {group_name}")


def plot_group(group_name, factors_to_plot, parameter_label, output_base):
    """
    group_name:
        "alpha" or "dF"
    parameter_label:
        r"\alpha" or r"d_F"
    output_base:
        "fig7" or "fig8"
    """

    fig, ax = plt.subplots(figsize=(6.7, 4.6))
    fig.patch.set_facecolor("white")

    colors = {
        0.25: "#1f77b4",
        0.5:  "#ff7f0e",
        1.0:  "#2ca02c",
        1.5:  "#d62728",
        2.0:  "#9467bd",
    }

    markers = {
        0.25: "o",
        0.5:  "s",
        1.0:  "^",
        1.5:  "D",
        2.0:  "v",
    }

    for fac in factors_to_plot:
        mask = get_mask(group_name, fac)

        x = A0[mask]
        y = Mbar_tobs[mask]

        if len(x) == 0:
            print(f"No data for {group_name}, factor={fac}")
            continue

        x, y = unique_by_A0(x, y)
        y_plot = positive_only(y)

        label = rf"${parameter_label}\times {fac:g}$"

        ax.plot(
            x,
            y_plot,
            marker=markers.get(fac, "o"),
            markersize=4.8,
            linewidth=2.0,
            color=colors.get(fac, None),
            markerfacecolor=colors.get(fac, None),
            markeredgecolor="white",
            markeredgewidth=0.45,
            label=label,
            zorder=3,
        )

    # ==============================
    # Axes
    # ==============================
    ax.set_yscale("log")

    ax.set_xlabel(r"Initial fungicide dose $A_0$")
    ax.set_ylabel("Normalized aggregate burden\n"
              r"$\overline{M}(T_{\mathrm{obs}})$")

    # y-axis lower bound kept low so small values are not clipped
    ax.set_ylim(1e-13, 2e-1)

    # x-axis matches the scanned dose interval [1, 3]
    ax.set_xlim(0.9, 3.1)
    ax.set_xticks([1.0, 1.5, 2.0, 2.5, 3.0])

    # Cleaner grid
    ax.grid(
        True,
        which="major",
        linestyle="-",
        linewidth=0.5,
        color="0.88",
        zorder=0,
    )

    ax.grid(
        True,
        which="minor",
        linestyle=":",
        linewidth=0.35,
        color="0.92",
        zorder=0,
    )

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    # ==============================
    # Legend
    # ==============================
    leg = ax.legend(
        loc="lower right",
        frameon=True,
        fancybox=False,
        facecolor="white",
        edgecolor="0.75",
        handlelength=2.0,
        borderpad=0.45,
        labelspacing=0.35,
    )

    leg.get_frame().set_linewidth(0.7)

    # ==============================
    # Layout
    # ==============================
    fig.subplots_adjust(
        left=0.16,
        right=0.98,
        bottom=0.20,
        top=0.97,
    )

    # ==============================
    # Save
    # ==============================
    png_path = os.path.join(OUT_DIR, f"{output_base}.png")
    eps_path = os.path.join(OUT_DIR, f"{output_base}.eps")

    fig.savefig(
        png_path,
        format="png",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.08,
        facecolor="white",
        edgecolor="none",
    )

    fig.savefig(
        eps_path,
        format="eps",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.08,
        facecolor="white",
        edgecolor="none",
    )

    plt.close(fig)

    print("Saved:")
    print(" -", png_path)
    print(" -", eps_path)


# ==============================
# Generate Figures
# ==============================
factors_to_plot = [0.25, 0.5, 1.0, 1.5, 2.0]

# Figure 7: varying alpha
plot_group(
    group_name="alpha",
    factors_to_plot=factors_to_plot,
    parameter_label=r"\alpha",
    output_base="fig7"
)

# Figure 8: varying dF
plot_group(
    group_name="dF",
    factors_to_plot=factors_to_plot,
    parameter_label=r"d_F",
    output_base="fig8"
)