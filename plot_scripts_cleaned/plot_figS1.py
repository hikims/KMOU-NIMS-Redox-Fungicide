#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure S1: Estimated empirical fixed-time suppression threshold
A0_emp*(T_obs; I*) as a function of multiplicative factors applied to
dA, dF, and alpha.

Filled markers : threshold found within the tested dose range.
Crosses        : the criterion I(T_obs) < I* was not met up to A0_max = 3.

INPUT
    results/threshold_parameter/summary/threshold_estimates.txt
    Whitespace-separated columns:
        group  factor  dA_factor  dF_factor  alpha_factor
        A0_lower  A0_upper  threshold_mid

OUTPUT (written to results/figures/)
    figS1.png, figS1.eps

USAGE
    python scripts/plot_figS1.py
"""

import os
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# ============================================================
# Repository-relative paths
# ============================================================
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_DATA = (
    REPO_ROOT / "results" / "threshold_parameter" / "summary" / "threshold_estimates.txt"
)
DEFAULT_OUTDIR = REPO_ROOT / "results" / "figures"


# ============================================================
# Publication style
# ============================================================
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


# ============================================================
# Plot configuration
# ============================================================
GROUPS = [
    ("dA", r"$d_A$", "o"),
    ("dF", r"$d_F$", "s"),
    ("alpha", r"$\alpha$", "^"),
]

COLORS = {
    "dA": "#1f77b4",
    "dF": "#ff7f0e",
    "alpha": "#2ca02c",
}

# Maximum tested dose used in the threshold search
A0_MAX = 3.0

# Vertical position for "not found" crosses and the annotation text
NOT_FOUND_Y = 3.08
TEXT_Y = 3.35


# ============================================================
# Main
# ============================================================
def main():
    data_file = DEFAULT_DATA
    out_dir = DEFAULT_OUTDIR
    out_dir.mkdir(parents=True, exist_ok=True)

    if not data_file.exists():
        raise FileNotFoundError(f"Missing data file: {data_file}")

    # columns: group  factor  dA_factor  dF_factor  alpha_factor
    #          A0_lower  A0_upper  threshold_mid
    data = np.genfromtxt(data_file, comments="#", dtype=None, encoding=None)

    group = data["f0"]
    factor = data["f1"]
    A0_lower = data["f5"]
    A0_upper = data["f6"]
    threshold_mid = data["f7"]

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    fig.patch.set_facecolor("white")

    for g, label, marker in GROUPS:
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

        # "Found" = finite threshold strictly below A0_MAX
        found = (~np.isnan(y)) & (y < A0_MAX)

        # Exclude any factor that already has a "found" entry from "not found"
        found_factors = set(np.round(x[found], 6).tolist())
        not_found = (~found) & np.array(
            [round(float(xx), 6) not in found_factors for xx in x],
            dtype=bool,
        )

        if np.any(found):
            yerr = np.vstack([y[found] - y_low[found], y_high[found] - y[found]])
            ax.errorbar(
                x[found], y[found], yerr=yerr,
                marker=marker, markersize=6.8,
                linewidth=2.2, elinewidth=1.3,
                capsize=4, capthick=1.2,
                color=COLORS[g],
                markerfacecolor=COLORS[g],
                markeredgecolor="white",
                markeredgewidth=0.7,
                label=label,
                zorder=4,
            )

        if np.any(not_found):
            ax.scatter(
                x[not_found],
                np.full(np.sum(not_found), NOT_FOUND_Y),
                marker="x", s=85, linewidths=2.0,
                color=COLORS[g],
                zorder=5,
            )

    # Reference line and "not found" annotation
    ax.axhline(A0_MAX, linestyle="--", linewidth=1.0, color="0.45", zorder=1)
    ax.text(0.31, TEXT_Y, "not found within tested range",
            fontsize=12, color="0.20", ha="left", va="center")

    ax.set_xlabel("Multiplicative factor")
    ax.set_ylabel(r"Estimated threshold $A_{0,\mathrm{emp}}^*$")
    ax.set_xlim(0.20, 2.05)
    ax.set_ylim(1.0, 3.5)
    ax.set_xticks([0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00])
    ax.set_yticks(np.arange(1.0, 3.51, 0.5))
    ax.grid(True, which="major", linestyle="-",
            linewidth=0.5, color="0.88", zorder=0)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    leg = ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.98, 0.72),
        frameon=True, fancybox=False,
        facecolor="white", edgecolor="0.75",
        handlelength=2.0, borderpad=0.45,
    )
    leg.get_frame().set_linewidth(0.7)

    fig.subplots_adjust(left=0.15, right=0.98, bottom=0.16, top=0.96)

    png_path = out_dir / "figS1.png"
    eps_path = out_dir / "figS1.eps"
    fig.savefig(png_path, format="png", dpi=600,
                bbox_inches="tight", facecolor="white", edgecolor="none")
    fig.savefig(eps_path, format="eps", dpi=300,
                bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)

    print("Saved:")
    print(" -", png_path)
    print(" -", eps_path)


if __name__ == "__main__":
    main()
