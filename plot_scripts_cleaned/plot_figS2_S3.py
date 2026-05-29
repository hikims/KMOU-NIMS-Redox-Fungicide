#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures S2 and S3: Dose-response curves of the normalized aggregate burden
Mbar(T_obs) at T_obs = 100.

    Figure S2 : five values of the fungicide consumption rate alpha
    Figure S3 : five values of the fungal diffusivity dF

These dose-response curves are complementary aggregate measures: the
empirical suppression threshold itself is defined through the infected-area
criterion of Section 5.1 of the paper, not through Mbar.

INPUT
    results/threshold_parameter/summary/threshold_scan_results.txt
    Column order:
        0:tag  1:A0  2:dA_factor  3:dF_factor  4:alpha_factor
        5:dA   6:dF  7:alpha
        8:total_F_t20    9:area_t20
       10:total_F_t50   11:area_t50
       12:total_F_t100  13:area_t100

    Mbar(T_obs) = total_F(T_obs) / |Omega|
    Set OMEGA_AREA according to whether column 12 already stores Mbar (1.0)
    or the raw integral total_F (100.0 for Omega = [0,10]^2).

OUTPUT (written to results/figures/)
    figS2.png, figS2.eps        # varying alpha
    figS3.png, figS3.eps        # varying dF

USAGE
    python scripts/plot_figS2_S3.py
"""

import os
import warnings
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


warnings.filterwarnings("ignore", category=UserWarning)


# ============================================================
# Repository-relative paths
# ============================================================
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_DATA = (
    REPO_ROOT / "results" / "threshold_parameter" / "summary" / "threshold_scan_results.txt"
)
DEFAULT_OUTDIR = REPO_ROOT / "results" / "figures"


# ============================================================
# Configuration
# ============================================================
# Set OMEGA_AREA = 1.0 if column 12 already stores the normalized average,
# or 100.0 if it stores the raw integral total_F (Omega = [0,10]^2).
OMEGA_AREA = 1.0


# ============================================================
# Publication style
# ============================================================
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


# ============================================================
# Helpers
# ============================================================
def positive_only(y: np.ndarray) -> np.ndarray:
    """Replace non-positive values by NaN for log-scale plotting."""
    return np.where(y > 0, y, np.nan)


def unique_by_A0(x: np.ndarray, y: np.ndarray):
    """Sort by A0 and remove duplicated A0 values."""
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    x_unique, unique_idx = np.unique(x, return_index=True)
    y_unique = y[unique_idx]
    return x_unique, y_unique


def get_mask(tag, dA_factor, dF_factor, alpha_factor, group_name, fac):
    if group_name == "dA":
        return np.char.startswith(tag, "dA_") & np.isclose(dA_factor, fac)
    if group_name == "dF":
        return np.char.startswith(tag, "dF_") & np.isclose(dF_factor, fac)
    if group_name == "alpha":
        return np.char.startswith(tag, "alpha_") & np.isclose(alpha_factor, fac)
    raise ValueError(f"Unknown group_name: {group_name}")


def plot_group(
    group_name: str,
    factors_to_plot,
    parameter_label: str,
    output_base: str,
    tag, A0, dA_factor, dF_factor, alpha_factor, Mbar_tobs,
    out_dir: Path,
):
    """
    group_name      : "alpha" or "dF"
    parameter_label : raw LaTeX symbol, e.g. r"\\alpha" or r"d_F"
    output_base     : output file stem, e.g. "figS2" or "figS3"
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
        0.25: "o", 0.5: "s", 1.0: "^", 1.5: "D", 2.0: "v",
    }

    for fac in factors_to_plot:
        mask = get_mask(tag, dA_factor, dF_factor, alpha_factor, group_name, fac)
        x = A0[mask]
        y = Mbar_tobs[mask]
        if len(x) == 0:
            print(f"No data for {group_name}, factor={fac}")
            continue

        x, y = unique_by_A0(x, y)
        y_plot = positive_only(y)
        label = rf"${parameter_label}\times {fac:g}$"

        ax.plot(
            x, y_plot,
            marker=markers.get(fac, "o"),
            markersize=4.8,
            linewidth=2.0,
            color=colors.get(fac),
            markerfacecolor=colors.get(fac),
            markeredgecolor="white",
            markeredgewidth=0.45,
            label=label,
            zorder=3,
        )

    ax.set_yscale("log")
    ax.set_xlabel(r"Initial fungicide dose $A_0$")
    ax.set_ylabel("Normalized aggregate burden\n"
                  r"$\overline{M}(T_{\mathrm{obs}})$")
    ax.set_ylim(1e-13, 2e-1)
    ax.set_xlim(0.9, 3.1)
    ax.set_xticks([1.0, 1.5, 2.0, 2.5, 3.0])
    ax.grid(True, which="major", linestyle="-",
            linewidth=0.5, color="0.88", zorder=0)
    ax.grid(True, which="minor", linestyle=":",
            linewidth=0.35, color="0.92", zorder=0)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    leg = ax.legend(
        loc="lower right",
        frameon=True, fancybox=False,
        facecolor="white", edgecolor="0.75",
        handlelength=2.0, borderpad=0.45, labelspacing=0.35,
    )
    leg.get_frame().set_linewidth(0.7)

    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.20, top=0.97)

    png_path = out_dir / f"{output_base}.png"
    eps_path = out_dir / f"{output_base}.eps"
    fig.savefig(png_path, format="png", dpi=600,
                bbox_inches="tight", pad_inches=0.08,
                facecolor="white", edgecolor="none")
    fig.savefig(eps_path, format="eps", dpi=300,
                bbox_inches="tight", pad_inches=0.08,
                facecolor="white", edgecolor="none")
    plt.close(fig)

    print(f"Saved {output_base}:")
    print(" -", png_path)
    print(" -", eps_path)


# ============================================================
# Main
# ============================================================
def main():
    data_file = DEFAULT_DATA
    out_dir = DEFAULT_OUTDIR
    out_dir.mkdir(parents=True, exist_ok=True)

    if not data_file.exists():
        raise FileNotFoundError(f"Missing data file: {data_file}")

    data = np.loadtxt(data_file, comments="#", dtype=str)
    print("Data shape:", data.shape)

    tag = data[:, 0]
    A0 = data[:, 1].astype(float)
    dA_factor = data[:, 2].astype(float)
    dF_factor = data[:, 3].astype(float)
    alpha_factor = data[:, 4].astype(float)

    # Column 12 = total_F at T_obs = 100
    total_F_raw = data[:, 12].astype(float)
    Mbar_tobs = total_F_raw / OMEGA_AREA

    factors_to_plot = [0.25, 0.5, 1.0, 1.5, 2.0]

    # Figure S2: varying alpha
    plot_group(
        group_name="alpha",
        factors_to_plot=factors_to_plot,
        parameter_label=r"\alpha",
        output_base="figS2",
        tag=tag, A0=A0,
        dA_factor=dA_factor, dF_factor=dF_factor, alpha_factor=alpha_factor,
        Mbar_tobs=Mbar_tobs,
        out_dir=out_dir,
    )

    # Figure S3: varying dF
    plot_group(
        group_name="dF",
        factors_to_plot=factors_to_plot,
        parameter_label=r"d_F",
        output_base="figS3",
        tag=tag, A0=A0,
        dA_factor=dA_factor, dF_factor=dF_factor, alpha_factor=alpha_factor,
        Mbar_tobs=Mbar_tobs,
        out_dir=out_dir,
    )


if __name__ == "__main__":
    main()
