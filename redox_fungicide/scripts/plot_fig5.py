#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 5: Spatial dynamics of the fungal density F(x,t).

A 3 x 4 grid of color maps:
    Rows: three fungicide doses (A0 = 0, 2.25, 5)
    Cols: four observation times (t = 0, 10, 50, 100)

The t = 0 column is constructed from the exact initial condition described
in Section 5.1 of the paper: a circular patch centered at (5,5) with
radius 2, where F = 0.5 inside and F = 0 outside. This avoids any noise
from the very first saved snapshot.

INPUT
    results/spatial_snapshots/A0_<A0>/field_<tag>.dat
    Accepted file name patterns:
        field_t000.dat / field_t010.dat / field_t050.dat / field_t100.dat
        field_0000.dat / field_0100.dat / field_0500.dat / field_1000.dat
    Each file has five whitespace-separated columns:
        x  y  A  F  R

OUTPUT (written to results/figures/)
    fig5.png, fig5.pdf

USAGE
    python scripts/plot_fig5.py
"""

import os
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


# ============================================================
# Repository-relative paths
# ============================================================
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_SNAP_DIR = REPO_ROOT / "results" / "spatial_snapshots"
DEFAULT_OUTDIR = REPO_ROOT / "results" / "figures"


# ============================================================
# Publication style
# ============================================================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 15,
    "axes.titlesize": 18,
    "axes.labelsize": 18,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "axes.linewidth": 0.8,
    "ps.fonttype": 42,
    "pdf.fonttype": 42,
})


# ============================================================
# Grid size (must match the simulation)
# ============================================================
Nx, Ny = 500, 500


# ============================================================
# Row / column configuration
# ============================================================
A0_rows = [
    ("No fungicide\n$A_0=0$", "A0_0.000"),
    ("Intermediate\n$A_0=2.25$", "A0_2.250"),
    ("High fungicide\n$A_0=5$", "A0_5.000"),
]

time_files = [
    (r"$t=0$",   ["field_t000.dat", "field_0000.dat"]),
    (r"$t=10$",  ["field_t010.dat", "field_0100.dat"]),
    (r"$t=50$",  ["field_t050.dat", "field_0500.dat"]),
    (r"$t=100$", ["field_t100.dat", "field_1000.dat"]),
]


# ============================================================
# Helpers
# ============================================================
def find_existing_file(folder: Path, candidates):
    """Return the first existing file in `folder` whose name matches any candidate."""
    for fname in candidates:
        path = folder / fname
        if path.exists():
            return path
    return None


def load_F(filepath: Path) -> np.ndarray:
    """Load the F field from a 5-column snapshot file (x y A F R)."""
    data = np.loadtxt(filepath)
    F = data[:, 3].reshape(Nx, Ny)
    return F


def make_physical_initial_patch() -> np.ndarray:
    """
    Build the exact t = 0 initial condition on the physical domain.
    The paper uses: center = (5, 5), radius = 2, F = 0.5 inside, F = 0 outside.
    """
    x = np.linspace(0.0, 10.0, Nx)
    y = np.linspace(0.0, 10.0, Ny)
    X, Y = np.meshgrid(x, y, indexing="ij")

    F0 = np.zeros((Nx, Ny))
    mask = (X - 5.0) ** 2 + (Y - 5.0) ** 2 <= 2.0 ** 2
    F0[mask] = 0.5
    return F0


# ============================================================
# Figure
# ============================================================
def main():
    snap_dir = DEFAULT_SNAP_DIR
    out_dir = DEFAULT_OUTDIR
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 4, figsize=(12.2, 8.4))
    fig.patch.set_facecolor("white")

    im = None
    missing_files = []

    F0_exact = make_physical_initial_patch()

    for row, (row_label, A0_folder) in enumerate(A0_rows):
        folder = snap_dir / A0_folder

        for col, (tlabel, candidates) in enumerate(time_files):
            ax = axes[row, col]

            if col == 0:
                # Use the exact analytical initial patch for t = 0
                F = F0_exact
            else:
                path = find_existing_file(folder, candidates)
                if path is None:
                    missing_files.append(
                        str(folder / " OR ".join(candidates))
                    )
                    F = None
                else:
                    F = load_F(path)

            if F is not None:
                img = ax.imshow(
                    F.T,
                    origin="lower",
                    cmap="viridis",
                    vmin=0.0,
                    vmax=0.5,
                    aspect="equal",
                    interpolation="nearest",
                    extent=[0, 10, 0, 10],
                )
                if im is None:
                    im = img
            else:
                ax.text(0.5, 0.5, "missing",
                        ha="center", va="center", fontsize=15,
                        transform=ax.transAxes)

            if row == 0:
                ax.set_title(tlabel, fontsize=18, pad=10)
            if col == 0:
                ax.set_ylabel(row_label, fontsize=18, rotation=90, labelpad=8)

            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_linewidth(0.8)

    # Layout
    fig.subplots_adjust(
        left=0.08, right=0.86,
        bottom=0.07, top=0.92,
        wspace=0.12, hspace=0.20,
    )

    # Shared colorbar
    cax = fig.add_axes([0.89, 0.20, 0.025, 0.60])
    if im is not None:
        cbar = fig.colorbar(im, cax=cax)
    else:
        norm = Normalize(vmin=0.0, vmax=0.5)
        sm = ScalarMappable(norm=norm, cmap="viridis")
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)

    cbar.set_label(r"Fungal density $F$", fontsize=18, labelpad=12)
    cbar.ax.tick_params(labelsize=14)

    # Save
    png_path = out_dir / "fig5.png"
    pdf_path = out_dir / "fig5.pdf"
    fig.savefig(png_path, format="png", dpi=600,
                bbox_inches="tight", facecolor="white", edgecolor="none")
    fig.savefig(pdf_path, format="pdf",
                bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)

    print("Saved:")
    print(" -", png_path)
    print(" -", pdf_path)

    if missing_files:
        print("\nMissing files:")
        for f in missing_files:
            print(" -", f)


if __name__ == "__main__":
    main()
