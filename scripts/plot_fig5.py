#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 5: Spatial dynamics of the fungal density F(x,t).

Notes:
  * t = 0 panels are loaded from actual simulation snapshots, not reconstructed.
  * Physical coordinates are the solver's node coordinates x = i*dx,
    with dx = Lx/(Nx-1) = 10/499 taken from src/parameters.h.
  * A thin white contour marks F = F_thr = 1e-3.
  * PDF/EPS outputs embed TrueType fonts.

INPUT (written by scripts/spatial_dynamics_exec.sh)
    results/fig5_rerun_A0_0_2p25_5/data/
          A0_0.000/
            field_0000.dat
            field_0100.dat
            field_0500.dat
            field_1000.dat
          A0_2.250/
            field_0000.dat
            field_0100.dat
            field_0500.dat
            field_1000.dat
          A0_5.000/
            field_0000.dat
            field_0100.dat
            field_0500.dat
            field_1000.dat

Each snapshot file should have columns:
    x  y  A  F  R

OUTPUT (written to results/figures/)
    fig5.pdf, fig5.eps, fig5.png

USAGE
    python scripts/plot_fig5.py
    python scripts/plot_fig5.py --snap-dir results/spatial_snapshots
    python scripts/plot_fig5.py --no-contour
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


# ============================================================
# Repository-relative paths
# ============================================================
# Every script in this repository resolves its inputs and outputs relative to
# the repository root, so `python scripts/<name>.py` behaves identically on
# any machine and from any working directory.  Set RESULTS_ROOT
# if the results tree is stored outside the repository, e.g.
#     export RESULTS_ROOT=/scratch/redox_fungicide/results
# This is the same variable the *_exec.sh drivers accept.
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
RESULTS_DIR = Path(
    os.environ.get("RESULTS_ROOT", REPO_ROOT / "results")
).expanduser().resolve()
DEFAULT_OUTDIR = RESULTS_DIR / "figures"

# Snapshot roots, nearest first.  `spatial_snapshots` is the layout produced
# by spatial_dynamics_exec.sh; the second entry is the dedicated Figure 5
# re-run used for the submitted version.
DEFAULT_SNAP_DIRS = [
    RESULTS_DIR / "fig5_rerun_A0_0_2p25_5" / "data",
    RESULTS_DIR / "spatial_snapshots",
]


# ============================================================
# Publication style
# ============================================================

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "dejavusans",
    "text.usetex": False,
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
    "svg.fonttype": "none",
})


# ============================================================
# Grid and figure parameters
# ============================================================

NX, NY = 500, 500
# Grid spacing, matching src/parameters.h:  dx = Lx / (Nx - 1) = 10/499.
LX = LY = 10.0
H = LX / (NX - 1)

F_THR = 1.0e-3
VMIN, VMAX = 0.0, 0.5

# initial.cpp evaluates the initial patch at the grid NODES, x = i*dx, so the
# plotted coordinates are node coordinates, not cell centres.
x_nodes = np.arange(NX) * H
y_nodes = np.arange(NY) * H
Xc, Yc = np.meshgrid(x_nodes, y_nodes, indexing="xy")

# imshow extent is measured to the outer edges of the first and last cells.
EXTENT = [-0.5 * H, LX + 0.5 * H, -0.5 * H, LY + 0.5 * H]


# ============================================================
# Rows and columns
# ============================================================

A0_ROWS = [
    ("No fungicide\n$A_0=0$", 0.0),
    ("Intermediate\n$A_0=2.25$", 2.25),
    ("High fungicide\n$A_0=5$", 5.0),
]

TIME_SPECS = [
    (r"$t=0$", 0),
    (r"$t=10$", 10),
    (r"$t=50$", 50),
    (r"$t=100$", 100),
]


# ============================================================
# Flexible file discovery
# ============================================================

def a0_folder_candidates(a0: float) -> list[str]:
    """Possible folder names for a given dose."""
    if np.isclose(a0, 0.0):
        return [
            "A0_0.000", "A0_0", "A0_0.0", "A0_0p000",
            "0", "0.000", "0p000",
        ]
    if np.isclose(a0, 2.25):
        return [
            "A0_2.250", "A0_2.25", "A0_2p250", "A0_2p25",
            "2.250", "2.25", "2p250", "2p25",
        ]
    if np.isclose(a0, 5.0):
        return [
            "A0_5.000", "A0_5", "A0_5.0", "A0_5p000",
            "5", "5.000", "5p000",
        ]

    return [
        f"A0_{a0:.3f}",
        f"A0_{a0:.2f}",
        f"A0_{a0:g}",
        f"{a0:.3f}",
        f"{a0:g}",
    ]


def time_file_candidates(t: int) -> list[str]:
    """Possible snapshot filenames for each displayed time."""
    if t == 0:
        return [
            "field_0000.dat",
            "field_t000.dat",
            "field_000.dat",
            "field_t0.dat",
            "field_0.dat",
            "F_0000.dat",
            "F_t000.dat",
        ]

    return [
        f"field_{10*t:04d}.dat",   # t=10 -> field_0100.dat
        f"field_t{t:03d}.dat",     # t=10 -> field_t010.dat
        f"field_{t:04d}.dat",
        f"field_{t:03d}.dat",
        f"field_{t}.dat",
        f"F_{10*t:04d}.dat",
        f"F_t{t:03d}.dat",
    ]


def find_snapshot_file(snap_dir: Path, a0: float, t: int) -> Optional[Path]:
    """Find the snapshot file for a given dose and time."""

    for folder_name in a0_folder_candidates(a0):
        folder = snap_dir / folder_name
        for filename in time_file_candidates(t):
            candidate = folder / filename
            if candidate.exists():
                return candidate

    # Fallback: recursive search inside snap_dir.
    # This helps if the folder names differ slightly.
    expected_names = set(time_file_candidates(t))
    dose_tokens = [name.lower() for name in a0_folder_candidates(a0)]

    if snap_dir.exists():
        matches = []
        for p in snap_dir.rglob("*.dat"):
            if p.name not in expected_names:
                continue
            full_text = "/".join(part.lower() for part in p.parts)
            if any(tok in full_text for tok in dose_tokens):
                matches.append(p)

        if matches:
            matches.sort(key=lambda p: (len(p.parts), str(p)))
            return matches[0]

    return None


def count_available_panels(snap_dir: Path) -> int:
    """Count how many of the 3 x 4 required panels are available."""
    if not snap_dir.exists():
        return 0

    count = 0
    for _, a0 in A0_ROWS:
        for _, t in TIME_SPECS:
            if find_snapshot_file(snap_dir, a0, t) is not None:
                count += 1
    return count


def choose_snap_dir(cli_snap_dir: Optional[str]) -> Path:
    """Choose the best snapshot directory."""

    if cli_snap_dir is not None:
        p = Path(cli_snap_dir).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(f"--snap-dir does not exist: {p}")
        return p

    scored = []
    for p in DEFAULT_SNAP_DIRS:
        if p.exists():
            scored.append((count_available_panels(p), p))

    print("Snapshot-directory candidates:")
    if scored:
        scored.sort(key=lambda item: item[0], reverse=True)
        for count, path in scored:
            print(f"  {count:2d}/12 panels: {path}")

        best_count, best_path = scored[0]
        if best_count > 0:
            return best_path

    print("  No valid panels found in default candidates.")
    return DEFAULT_SNAP_DIRS[0]


# ============================================================
# Loading and plotting helpers
# ============================================================

def load_F(filepath: Path) -> np.ndarray:
    """Load fungal density F from a snapshot file.

    Expected columns:
        x  y  A  F  R

    The fourth column, index 3, is F.
    """
    data = np.loadtxt(filepath)

    if data.ndim != 2:
        raise ValueError(f"Expected 2D array in {filepath}, got shape {data.shape}")

    if data.shape[1] < 4:
        raise ValueError(
            f"Expected at least 4 columns in {filepath}; "
            f"snapshot format should be x y A F R. Got shape {data.shape}"
        )

    n_rows = data.shape[0]
    expected = NX * NY

    if n_rows == expected:
        return data[:, 3].reshape(NX, NY)

    # Allow a different square grid, but warn through the returned shape.
    side = int(round(np.sqrt(n_rows)))
    if side * side == n_rows:
        print(f"[warning] {filepath} has {side} x {side} grid, not {NX} x {NY}.")
        return data[:, 3].reshape(side, side)

    raise ValueError(
        f"Cannot reshape F field from {filepath}. "
        f"Expected {expected} rows or a square number of rows; got {n_rows}."
    )


def add_detection_contour(ax, F: np.ndarray) -> None:
    """Add white contour at F = F_THR if the level exists in the panel."""
    if F.shape != (NX, NY):
        return

    fmin = np.nanmin(F)
    fmax = np.nanmax(F)

    if fmin <= F_THR <= fmax:
        ax.contour(
            Xc,
            Yc,
            F.T,
            levels=[F_THR],
            colors="0.92",
            linewidths=0.45,
        )


def save_figure(fig, out_dir: Path, stem: str) -> None:
    """Save PNG, PDF, and EPS."""
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = [
        (out_dir / f"{stem}.png", "png", 600),
        (out_dir / f"{stem}.pdf", "pdf", None),
        (out_dir / f"{stem}.eps", "eps", 300),
    ]

    print("Saved:")
    for path, fmt, dpi in outputs:
        kwargs = dict(
            format=fmt,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
        if dpi is not None:
            kwargs["dpi"] = dpi

        fig.savefig(path, **kwargs)
        print(" -", path)


# ============================================================
# Main plotting routine
# ============================================================

def make_figure(snap_dir: Path, out_dir: Path, add_contour: bool = True) -> None:
    available = count_available_panels(snap_dir)

    if available < 12:
        missing = []
        for _, a0 in A0_ROWS:
            for _, t in TIME_SPECS:
                p = find_snapshot_file(snap_dir, a0, t)
                if p is None:
                    missing.append((a0, t))

        msg = [
            f"Only {available}/12 required panels were found in:",
            f"  {snap_dir}",
            "",
            "Missing panels:",
        ]
        for a0, t in missing:
            msg.append(f"  A0={a0:g}, t={t}")

        msg.append("")
        msg.append("Run scripts/spatial_dynamics_exec.sh first, or pass --snap-dir")
        msg.append("pointing at a directory that contains the A0_0.000, A0_2.250")
        msg.append("and A0_5.000 folders.")

        raise FileNotFoundError("\n".join(msg))

    print(f"Using snapshot directory: {snap_dir}")
    print(f"Output directory: {out_dir}")

    fig, axes = plt.subplots(3, 4, figsize=(12.2, 8.4))
    fig.patch.set_facecolor("white")

    im = None
    loaded_files: list[Path] = []

    for row, (row_label, a0) in enumerate(A0_ROWS):
        for col, (time_label, t) in enumerate(TIME_SPECS):
            ax = axes[row, col]
            path = find_snapshot_file(snap_dir, a0, t)

            if path is None:
                # This should not happen because we checked above.
                ax.text(
                    0.5,
                    0.5,
                    "missing",
                    ha="center",
                    va="center",
                    fontsize=15,
                    transform=ax.transAxes,
                )
                continue

            F = load_F(path)
            loaded_files.append(path)

            img = ax.imshow(
                F.T,
                origin="lower",
                cmap="viridis",
                vmin=VMIN,
                vmax=VMAX,
                aspect="equal",
                interpolation="nearest",
                extent=EXTENT,
            )

            if add_contour:
                add_detection_contour(ax, F)

            if im is None:
                im = img

            if row == 0:
                ax.set_title(time_label, fontsize=18, pad=10)

            if col == 0:
                ax.set_ylabel(
                    row_label,
                    fontsize=18,
                    rotation=90,
                    labelpad=8,
                )

            ax.set_xticks([])
            ax.set_yticks([])

            for spine in ax.spines.values():
                spine.set_linewidth(0.8)

    # Layout: reserve right margin for shared colorbar.
    fig.subplots_adjust(
        left=0.08,
        right=0.86,
        bottom=0.08,
        top=0.92,
        wspace=0.12,
        hspace=0.20,
    )

    cax = fig.add_axes([0.89, 0.20, 0.025, 0.60])

    if im is not None:
        cbar = fig.colorbar(im, cax=cax)
    else:
        norm = Normalize(vmin=VMIN, vmax=VMAX)
        sm = ScalarMappable(norm=norm, cmap="viridis")
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)

    cbar.set_label(r"Fungal density $F$", fontsize=18, labelpad=12)
    cbar.ax.tick_params(labelsize=14)

    if add_contour:
        fig.text(
            0.50,
            0.030,
            r"White contour: $F=F_{\rm thr}=10^{-3}$",
            ha="center",
            va="bottom",
            fontsize=10.5,
            color="0.25",
        )

    save_figure(fig, out_dir, "fig5")
    plt.close(fig)

    print("")
    print("Loaded files:")
    for p in loaded_files:
        print(" -", p)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Figure 5 spatial snapshots."
    )

    parser.add_argument(
        "--snap-dir",
        type=str,
        default=None,
        help="Snapshot directory. Default: results/fig5_rerun_A0_0_2p25_5/data.",
    )

    parser.add_argument(
        "--outdir",
        type=str,
        default=str(DEFAULT_OUTDIR),
        help="Output directory. Default: results/figures.",
    )

    parser.add_argument(
        "--no-contour",
        action="store_true",
        help="Disable the white F=1e-3 contour.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    snap_dir = choose_snap_dir(args.snap_dir)
    out_dir = Path(args.outdir).expanduser().resolve()

    make_figure(
        snap_dir=snap_dir,
        out_dir=out_dir,
        add_contour=not args.no_contour,
    )


if __name__ == "__main__":
    main()