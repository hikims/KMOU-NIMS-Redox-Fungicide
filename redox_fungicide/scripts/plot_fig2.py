#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2: Dose-time map of finite-time practical suppression.

Two panels:
    (a) Normalized aggregate burden  log10 Mbar(t)
    (b) Infected-area fraction       A(t)  in %

The dashed white contour in panel (a) marks Mbar = 2e-2, the representative
finite-time suppression level discussed in Section 5.2 of the paper.

INPUT
    results/dose_time_map/B1_A0_Tobs_map.txt
    Expected columns:
        A0  time  total_A  total_F  M_inf  area_percent
        min_A  max_A  min_R  max_R  min_F  max_F
        log10_total_F  log10_M_inf

OUTPUT (written to results/figures/)
    fig2.pdf, fig2.eps, fig2.png            : two-panel figure
    fig2_Mbar_only.pdf/.eps/.png            : single-panel variant
    fig2_summary_table.csv                  : representative values for A0 = 0..5

DOMAIN
    Omega = [0,10] x [0,10],  |Omega| = 100
    Mbar(t) = total_F(t) / |Omega|

USAGE
    python scripts/plot_fig2.py
    python scripts/plot_fig2.py --data path/to/B1_A0_Tobs_map.txt
    python scripts/plot_fig2.py --outdir path/to/output --contour 2e-2
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ============================================================
# Repository-relative paths
# ============================================================
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_DATA = REPO_ROOT / "results" / "dose_time_map" / "B1_A0_Tobs_map.txt"
DEFAULT_OUTDIR = REPO_ROOT / "results" / "figures"


# ============================================================
# Constants
# ============================================================
DOMAIN_AREA = 100.0
MBAR_CONTOUR = 2.0e-2


# ============================================================
# Style
# ============================================================
def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 18,
        "axes.labelsize": 25,
        "axes.titlesize": 24,
        "axes.titleweight": "bold",
        "axes.linewidth": 1.7,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.fontsize": 17,
        "legend.frameon": True,
        "legend.fancybox": False,
        "legend.framealpha": 1.0,
        "legend.edgecolor": "0.80",
        "lines.linewidth": 3.0,
        "lines.markersize": 8.0,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "mathtext.fontset": "dejavusans",
    })


# ============================================================
# Data loading
# ============================================================
def read_b1_data(data_file: Path) -> pd.DataFrame:
    cols = [
        "A0", "time", "total_A", "total_F", "M_inf", "area_percent",
        "min_A", "max_A", "min_R", "max_R", "min_F", "max_F",
        "log10_total_F", "log10_M_inf",
    ]
    df = pd.read_csv(
        data_file,
        sep=r"\s+",
        comment="#",
        names=cols,
        engine="python",
    )
    if df.empty:
        raise ValueError(f"No numerical data were read from {data_file}")

    df["Mbar"] = df["total_F"] / DOMAIN_AREA
    df["log10_Mbar"] = np.log10(np.maximum(df["Mbar"].to_numpy(), 1e-300))

    df = (
        df.sort_values(["A0", "time"])
          .drop_duplicates(subset=["A0", "time"], keep="last")
          .reset_index(drop=True)
    )
    return df


def pivot_quantity(df: pd.DataFrame, quantity: str) -> pd.DataFrame:
    return (
        df.pivot_table(
            index="time", columns="A0", values=quantity, aggfunc="mean",
        )
        .sort_index()
        .sort_index(axis=1)
    )


# ============================================================
# Summary table
# ============================================================
def make_summary_table(df: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    target_doses = [0, 1, 2, 3, 4, 5]
    available_A0 = np.array(sorted(df["A0"].unique()), dtype=float)
    rows = []

    for target in target_doses:
        if available_A0.size == 0:
            continue
        A0_val = available_A0[np.argmin(np.abs(available_A0 - target))]
        sub = df[np.isclose(df["A0"], A0_val)].copy()
        if sub.empty:
            continue

        idx_min_mbar = sub["Mbar"].idxmin()
        row_min_mbar = sub.loc[idx_min_mbar]
        idx_min_minf = sub["M_inf"].idxmin()
        row_min_minf = sub.loc[idx_min_minf]
        idx_t100 = (sub["time"] - 100.0).abs().idxmin()
        row_t100 = sub.loc[idx_t100]

        rows.append({
            "A0": A0_val,
            "min_Mbar": row_min_mbar["Mbar"],
            "time_min_Mbar": row_min_mbar["time"],
            "min_Minf": row_min_minf["M_inf"],
            "time_min_Minf": row_min_minf["time"],
            "Mbar_at_100": row_t100["Mbar"],
            "Minf_at_100": row_t100["M_inf"],
            "area_percent_at_100": row_t100["area_percent"],
        })

    summary = pd.DataFrame(rows)
    out_file = out_dir / "fig2_summary_table.csv"
    summary.to_csv(out_file, index=False)
    print(f"Saved summary table: {out_file}")
    return summary


# ============================================================
# Plot helpers
# ============================================================
def add_suppression_contour(ax, A0_grid, T_grid, Z_mbar, contour_level):
    zmin = np.nanmin(Z_mbar)
    zmax = np.nanmax(Z_mbar)
    if not (zmin <= contour_level <= zmax):
        print(
            f"[warning] contour level {contour_level:g} is outside "
            f"data range [{zmin:g}, {zmax:g}]. No contour was drawn."
        )
        return
    ax.contour(
        A0_grid, T_grid, Z_mbar,
        levels=[contour_level],
        colors="white",
        linewidths=2.4,
        linestyles="--",
    )


def add_contour_legend(ax) -> None:
    handle = Line2D(
        [0], [0],
        color="black", linewidth=3.0, linestyle="--",
        label=r"$\overline{M}=2\times10^{-2}$",
    )
    leg = ax.legend(
        handles=[handle],
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        frameon=True, fancybox=False, framealpha=1.0,
        facecolor="white", edgecolor="0.75",
        borderpad=0.6, handlelength=2.6,
    )
    leg.get_frame().set_linewidth(1.3)


def style_axes(ax) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(1.7)
    ax.tick_params(axis="both", which="major", width=1.5, length=6)


def style_colorbar(cbar, label: str) -> None:
    cbar.set_label(label, fontsize=21, labelpad=12)
    cbar.ax.tick_params(labelsize=16, width=1.3, length=5)
    cbar.outline.set_linewidth(1.3)


def save_figure(fig, out_dir: Path, stem: str) -> None:
    for ext in ["pdf", "eps", "png"]:
        out_file = out_dir / f"{stem}.{ext}"
        fig.savefig(out_file, dpi=400, bbox_inches="tight")
        print(f"Saved figure: {out_file}")


# ============================================================
# Plotting
# ============================================================
def make_two_panel_figure(df, out_dir, contour_level=MBAR_CONTOUR):
    p_log_mbar = pivot_quantity(df, "log10_Mbar")
    p_mbar = pivot_quantity(df, "Mbar")
    p_area = pivot_quantity(df, "area_percent")

    A0s = p_log_mbar.columns.to_numpy(dtype=float)
    times = p_log_mbar.index.to_numpy(dtype=float)
    A0_grid, T_grid = np.meshgrid(A0s, times)

    fig = plt.figure(figsize=(17.2, 6.2))
    gs = fig.add_gridspec(
        nrows=1, ncols=5,
        width_ratios=[1.0, 0.045, 0.34, 1.0, 0.045],
        wspace=0.20,
    )

    ax0 = fig.add_subplot(gs[0, 0])
    cax0 = fig.add_subplot(gs[0, 1])
    spacer = fig.add_subplot(gs[0, 2])
    spacer.axis("off")
    ax1 = fig.add_subplot(gs[0, 3])
    cax1 = fig.add_subplot(gs[0, 4])

    # Panel (a): normalized aggregate burden
    im0 = ax0.pcolormesh(
        A0_grid, T_grid, p_log_mbar.to_numpy(),
        shading="auto", cmap="viridis",
    )
    cbar0 = fig.colorbar(im0, cax=cax0)
    style_colorbar(cbar0, r"$\log_{10}\overline{M}(t)$")
    add_suppression_contour(ax0, A0_grid, T_grid, p_mbar.to_numpy(), contour_level)
    add_contour_legend(ax0)
    ax0.set_xlabel(r"Initial fungicide dose $A_0$", fontsize=25)
    ax0.set_ylabel(r"Observation time $t$", fontsize=25)
    ax0.set_title(r"(a) Normalized aggregate burden", fontsize=24, fontweight="bold", pad=14)
    ax0.set_xlim(float(np.nanmin(A0s)), float(np.nanmax(A0s)))
    ax0.set_ylim(float(np.nanmin(times)), float(np.nanmax(times)))
    ax0.set_xticks(np.arange(0, 5.1, 1.0))
    ax0.set_yticks(np.arange(0, 101, 20))
    style_axes(ax0)

    # Panel (b): infected-area fraction
    im1 = ax1.pcolormesh(
        A0_grid, T_grid, p_area.to_numpy(),
        shading="auto", cmap="viridis",
    )
    cbar1 = fig.colorbar(im1, cax=cax1)
    style_colorbar(cbar1, r"Infected-area fraction $\mathcal{A}(t)$ (%)")
    ax1.set_xlabel(r"Initial fungicide dose $A_0$", fontsize=25)
    ax1.set_ylabel(r"Observation time $t$", fontsize=25)
    ax1.set_title(r"(b) Infected-area fraction", fontsize=24, fontweight="bold", pad=14)
    ax1.set_xlim(float(np.nanmin(A0s)), float(np.nanmax(A0s)))
    ax1.set_ylim(float(np.nanmin(times)), float(np.nanmax(times)))
    ax1.set_xticks(np.arange(0, 5.1, 1.0))
    ax1.set_yticks(np.arange(0, 101, 20))
    style_axes(ax1)

    fig.subplots_adjust(left=0.065, right=0.965, bottom=0.17, top=0.86)

    save_figure(fig, out_dir, "fig2")
    plt.close(fig)


def make_mbar_only_figure(df, out_dir, contour_level=MBAR_CONTOUR):
    p_log_mbar = pivot_quantity(df, "log10_Mbar")
    p_mbar = pivot_quantity(df, "Mbar")

    A0s = p_log_mbar.columns.to_numpy(dtype=float)
    times = p_log_mbar.index.to_numpy(dtype=float)
    A0_grid, T_grid = np.meshgrid(A0s, times)

    fig, ax = plt.subplots(figsize=(8.8, 6.2), constrained_layout=True)
    im = ax.pcolormesh(
        A0_grid, T_grid, p_log_mbar.to_numpy(),
        shading="auto", cmap="viridis",
    )
    cbar = fig.colorbar(im, ax=ax, pad=0.035)
    style_colorbar(cbar, r"$\log_{10}\overline{M}(t)$")
    add_suppression_contour(ax, A0_grid, T_grid, p_mbar.to_numpy(), contour_level)
    add_contour_legend(ax)
    ax.set_xlabel(r"Initial fungicide dose $A_0$", fontsize=25)
    ax.set_ylabel(r"Observation time $t$", fontsize=25)
    ax.set_title(r"Dose--time map of finite-time suppression",
                 fontsize=24, fontweight="bold", pad=14)
    ax.set_xlim(float(np.nanmin(A0s)), float(np.nanmax(A0s)))
    ax.set_ylim(float(np.nanmin(times)), float(np.nanmax(times)))
    ax.set_xticks(np.arange(0, 5.1, 1.0))
    ax.set_yticks(np.arange(0, 101, 20))
    style_axes(ax)

    save_figure(fig, out_dir, "fig2_Mbar_only")
    plt.close(fig)


# ============================================================
# Main
# ============================================================
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Figure 2: dose-time map.")
    parser.add_argument("--data", type=str, default=None,
                        help=f"Path to B1_A0_Tobs_map.txt (default: {DEFAULT_DATA})")
    parser.add_argument("--outdir", type=str, default=None,
                        help=f"Output directory (default: {DEFAULT_OUTDIR})")
    parser.add_argument("--contour", type=float, default=MBAR_CONTOUR,
                        help="Contour level for Mbar (default: 2e-2).")
    return parser.parse_args()


def main() -> None:
    set_publication_style()
    args = parse_args()

    data_file = Path(args.data).expanduser() if args.data else DEFAULT_DATA
    if not data_file.exists():
        raise FileNotFoundError(
            f"Data file not found: {data_file}\n"
            f"Provide it with --data path/to/B1_A0_Tobs_map.txt"
        )

    out_dir = Path(args.outdir).expanduser() if args.outdir else DEFAULT_OUTDIR
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Data file       : {data_file}")
    print(f"Output directory: {out_dir}")
    print(f"Contour level   : {args.contour:g}")

    df = read_b1_data(data_file)

    print("\nData summary")
    print(f"  A0 range     : {df['A0'].min():.3f} to {df['A0'].max():.3f}")
    print(f"  time range   : {df['time'].min():.3f} to {df['time'].max():.3f}")
    print(f"  Mbar range   : {df['Mbar'].min():.6e} to {df['Mbar'].max():.6e}")

    make_two_panel_figure(df, out_dir, contour_level=args.contour)
    make_mbar_only_figure(df, out_dir, contour_level=args.contour)
    summary = make_summary_table(df, out_dir)

    print("\nRepresentative summary")
    if not summary.empty:
        print(summary.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()
