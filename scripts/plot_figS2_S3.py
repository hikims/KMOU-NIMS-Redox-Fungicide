#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures S2 and S3: dose-response curves of the normalized aggregate burden.

Figure S2:
    variation of the fungicide-consumption rate alpha.

Figure S3:
    variation of the fungal diffusivity dF.

The empirical suppression threshold itself is defined through the infected-area
criterion, not through the aggregate burden. These two figures are supporting
dose-response summaries for the aggregate burden at T_obs.

Input:
        threshold_parameter/summary/threshold_scan_results.txt

The practical column order is:
    tag  A0  dA_factor  dF_factor  alpha_factor  dA  dF  alpha
    total_F_t20  area_t20  total_F_t50  area_t50  total_F_t100  area_t100

The tag has the form:
    dA_0.25x_A0_1.00
    dF_1.50x_A0_2.00
    alpha_2.00x_A0_4.80

Outputs:
    figS2.png, figS2.pdf, figS2.eps
    figS3.png, figS3.pdf, figS3.eps
    figS2_data_used.csv
    figS3_data_used.csv

Usage:
    python plot_figS2_S3.py

Optional:
    python plot_figS2_S3.py --x-left 0.50 --x-right 5.80
    python plot_figS2_S3.py --x-left 0.50 --x-right 5.50
    python plot_figS2_S3.py --omega-area 100
"""

from __future__ import annotations

import argparse
import os
import re
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


warnings.filterwarnings("ignore", category=UserWarning)


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


def _first_existing(candidates):
    """First candidate that exists; otherwise the first, for clear errors."""
    candidates = [Path(p) for p in candidates]
    for path in candidates:
        if path.exists():
            return path
    return candidates[0]


# The parameter-dependence scan is written by scripts/threshold_exec.sh
# (driver main_supp_figS1_threshold).  The second candidate is the output of
# the legacy sensitivity driver, kept so an older results tree still works.
# --data overrides both.
DEFAULT_DATA = _first_existing([
    RESULTS_DIR / "threshold_parameter" / "summary" / "threshold_scan_results.txt",
    RESULTS_DIR / "sensitivity" / "summary" / "sensitivity_results.txt",
])


# ============================================================
# Configuration
# ============================================================

DEFAULT_TOBS = 100

# Keep this as 1.0 if total_F_t100 already stores the normalized aggregate burden.
# Use --omega-area 100 if total_F_t100 stores the raw integral over Omega=[0,10]^2.
DEFAULT_OMEGA_AREA = 1.0

# Reviewer-polished x-axis range:
#   left  = 0.50: starts at A0 = 0.5 as requested.
#   right = 5.80: creates blank space for the lower-right legend.
# Tick labels are shown only from A0 = 0.5 to A0 = 5.0.
DEFAULT_X_LEFT = 0.50
DEFAULT_X_RIGHT = 5.80
DEFAULT_X_TICK_MIN = 0.50
DEFAULT_X_TICK_MAX = 5.00

FACTORS_TO_PLOT = [0.25, 0.50, 1.00, 1.50, 2.00]

COLORS = {
    0.25: "#1f77b4",
    0.50: "#ff7f0e",
    1.00: "#2ca02c",
    1.50: "#d62728",
    2.00: "#9467bd",
}

MARKERS = {
    0.25: "o",
    0.50: "s",
    1.00: "^",
    1.50: "D",
    2.00: "v",
}

TAG_PATTERN = re.compile(
    r"^(dA|DA|dF|DF|alpha)_([0-9]+(?:\.[0-9]+)?)x_A0_([0-9]+(?:\.[0-9]+)?)$"
)


# ============================================================
# Publication style
# ============================================================

def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "text.usetex": False,

        "font.size": 12,
        "axes.labelsize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 13,
        "axes.linewidth": 1.0,
        "lines.linewidth": 2.4,

        "figure.dpi": 300,
        "savefig.dpi": 600,
        "ps.fonttype": 42,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })


# ============================================================
# Data parsing
# ============================================================

def normalize_group(group: str) -> str:
    if group in {"DA", "dA"}:
        return "dA"
    if group in {"DF", "dF"}:
        return "dF"
    if group == "alpha":
        return "alpha"
    return group


def parse_tag(tag: str) -> tuple[Optional[str], float, float]:
    """
    Parse group, multiplicative factor, and A0 from a tag.

    Example:
        dF_1.50x_A0_2.00 -> ("dF", 1.50, 2.00)
    """
    m = TAG_PATTERN.match(str(tag))

    if m is None:
        return None, np.nan, np.nan

    group = normalize_group(m.group(1))
    factor = float(m.group(2))
    a0_tag = float(m.group(3))

    return group, factor, a0_tag


def read_scan_file(data_file: Path) -> pd.DataFrame:
    """
    Read threshold_scan_results.txt robustly.

    The legacy header can be misleading, so the file is read line-by-line and
    assigned the practical column order used by the plotting scripts.
    """
    if not data_file.exists():
        raise FileNotFoundError(f"Missing data file: {data_file}")

    cols = [
        "tag",
        "A0",
        "dA_factor",
        "dF_factor",
        "alpha_factor",
        "dA",
        "dF",
        "alpha",
        "total_F_t20",
        "area_t20",
        "total_F_t50",
        "area_t50",
        "total_F_t100",
        "area_t100",
    ]

    rows = []
    skipped = 0

    with open(data_file, "r", encoding="utf-8", errors="ignore") as f:
        for line_number, line in enumerate(f, start=1):
            s = line.strip()

            if not s:
                continue
            if s.startswith("#"):
                continue

            parts = s.split()

            if len(parts) < len(cols):
                skipped += 1
                print(
                    f"[warning] skipping line {line_number}: "
                    f"expected at least {len(cols)} columns, got {len(parts)}"
                )
                continue

            if len(parts) > len(cols):
                parts = parts[:len(cols)]

            rows.append(parts)

    if not rows:
        raise ValueError(f"No numerical rows were read from {data_file}")

    if skipped:
        print(f"[warning] skipped {skipped} malformed rows.")

    df = pd.DataFrame(rows, columns=cols)

    parsed = [parse_tag(tag) for tag in df["tag"].astype(str)]
    parsed_df = pd.DataFrame(parsed, columns=["tag_group", "tag_factor", "A0_tag"])

    df = pd.concat([df, parsed_df], axis=1)

    for col in df.columns:
        if col not in {"tag", "tag_group"}:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["A0", "tag_group", "tag_factor"]).copy()

    mismatch = np.abs(df["A0"] - df["A0_tag"]) > 5.0e-4
    if mismatch.any():
        print(f"[warning] {int(mismatch.sum())} rows have A0 values inconsistent with the tag.")

    # Keep the last occurrence of duplicate deterministic scans.
    df["_row_order"] = np.arange(len(df))
    df = (
        df.sort_values("_row_order")
          .drop_duplicates(subset=["tag_group", "tag_factor", "A0"], keep="last")
          .drop(columns=["_row_order"])
          .sort_values(["tag_group", "tag_factor", "A0"])
          .reset_index(drop=True)
    )

    return df


# ============================================================
# Curve extraction
# ============================================================

def positive_only(y: np.ndarray) -> np.ndarray:
    """Replace non-positive values by NaN for log-scale plotting."""
    return np.where(y > 0.0, y, np.nan)


def factor_label(factor: float) -> str:
    if np.isclose(factor, round(factor)):
        return f"{int(round(factor))}"
    return f"{factor:g}"


def total_column_for_tobs(tobs: int) -> str:
    return f"total_F_t{int(tobs)}"


def extract_curve(
    df: pd.DataFrame,
    group_name: str,
    factor: float,
    total_col: str,
    omega_area: float,
) -> pd.DataFrame:
    """
    Extract one factor curve for group_name in {"alpha", "dF"}.

    Uses the tag-parsed group and factor. This avoids relying on an old or
    misleading header line in the scan file.
    """
    sub = df[
        (df["tag_group"] == group_name)
        & np.isclose(df["tag_factor"], factor)
    ].copy()

    if sub.empty:
        return sub

    if total_col not in sub.columns:
        raise ValueError(f"Column {total_col!r} is missing from the scan file.")

    sub = sub.dropna(subset=["A0", total_col]).copy()
    sub = (
        sub.sort_values("A0")
           .drop_duplicates(subset=["A0"], keep="last")
           .reset_index(drop=True)
    )

    sub["Mbar_Tobs"] = sub[total_col] / omega_area
    sub["Mbar_plot"] = positive_only(sub["Mbar_Tobs"].to_numpy(dtype=float))

    return sub


# ============================================================
# Plotting helpers
# ============================================================

def save_all(fig, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    for ext, dpi in [("png", 600), ("pdf", None), ("eps", 300)]:
        path = out_dir / f"{stem}.{ext}"

        kwargs = dict(
            format=ext,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )

        if dpi is not None:
            kwargs["dpi"] = dpi

        fig.savefig(path, **kwargs)
        print(f"Saved: {path}")


def choose_xlim(
    all_x: list[np.ndarray],
    x_left: Optional[float],
    x_right: Optional[float],
) -> tuple[float, float]:
    finite_x = []

    for x in all_x:
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        if x.size:
            finite_x.append(x)

    if not finite_x:
        return DEFAULT_X_LEFT, DEFAULT_X_RIGHT

    xcat = np.concatenate(finite_x)
    xmin = float(np.nanmin(xcat))
    xmax = float(np.nanmax(xcat))

    span = max(xmax - xmin, 1.0)

    if x_left is None:
        left = max(DEFAULT_X_LEFT, xmin - 0.03 * span)
    else:
        left = float(x_left)

    if x_right is None:
        right = xmax + 0.03 * span
    else:
        right = float(x_right)

    if right <= left:
        right = left + 1.0

    return left, right


def set_xticks_from_xlim(
    ax,
    left: float,
    right: float,
    tick_min: float = DEFAULT_X_TICK_MIN,
    tick_max: float = DEFAULT_X_TICK_MAX,
) -> None:
    """
    Set clean half-step ticks.

    The x-axis range may extend beyond A0=5 to create blank space for the
    legend, but tick labels are shown only from tick_min to tick_max.
    """
    tick_start = max(tick_min, np.ceil(left * 2.0) / 2.0)
    tick_end = min(tick_max, np.floor(right * 2.0) / 2.0)

    ticks = np.arange(tick_start, tick_end + 1.0e-9, 0.5)

    if ticks.size > 0:
        ax.set_xticks(ticks)


def plot_group(
    df: pd.DataFrame,
    group_name: str,
    parameter_label: str,
    output_base: str,
    total_col: str,
    omega_area: float,
    out_dir: Path,
    x_left: Optional[float],
    x_right: Optional[float],
    x_tick_min: float,
    x_tick_max: float,
    y_min: float,
    y_max: float,
    save_csv: bool = True,
) -> None:
    """
    Plot one supplement figure.

    group_name:
        "alpha" for Fig. S2
        "dF"    for Fig. S3

    parameter_label:
        r"\\alpha" or r"d_F"
    """
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    fig.patch.set_facecolor("white")

    used_frames = []
    all_x_for_xlim = []

    for factor in FACTORS_TO_PLOT:
        sub = extract_curve(
            df=df,
            group_name=group_name,
            factor=factor,
            total_col=total_col,
            omega_area=omega_area,
        )

        if sub.empty:
            print(f"[warning] no data for {group_name} factor {factor:g}")
            continue

        x = sub["A0"].to_numpy(dtype=float)
        y = sub["Mbar_plot"].to_numpy(dtype=float)

        all_x_for_xlim.append(x)

        label = rf"${parameter_label}\times {factor_label(factor)}$"

        ax.plot(
            x,
            y,
            marker=MARKERS[factor],
            markersize=6.5,
            linewidth=2.4,
            color=COLORS[factor],
            markerfacecolor=COLORS[factor],
            markeredgecolor="white",
            markeredgewidth=0.6,
            label=label,
            zorder=3,
        )

        sub = sub.copy()
        sub["figure"] = output_base
        sub["plotted_group"] = group_name
        sub["plotted_factor"] = factor
        used_frames.append(sub)

    if not used_frames:
        raise RuntimeError(f"No curves were plotted for {group_name}.")

    left, right = choose_xlim(all_x_for_xlim, x_left=x_left, x_right=x_right)

    ax.set_xlim(left, right)
    set_xticks_from_xlim(
        ax,
        left,
        right,
        tick_min=x_tick_min,
        tick_max=x_tick_max,
    )

    ax.set_yscale("log")
    ax.set_ylim(y_min, y_max)

    ax.set_xlabel(r"Initial fungicide dose $A_0$", fontsize=22)
    ax.set_ylabel(
        r"Normalized aggregate burden" + "\n" + r"$\overline{M}(T_{\rm obs})$",
        fontsize=22,
    )

    ax.grid(True, which="major", color="0.88", linestyle="-", linewidth=0.6)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    # With x_right = 5.80, the lower-right legend sits in the blank region
    # and covers fewer data curves.
    leg = ax.legend(
        loc="lower right",
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="0.75",
        handlelength=2.0,
        borderpad=0.55,
    )
    leg.get_frame().set_linewidth(0.7)

    fig.subplots_adjust(left=0.15, right=0.98, bottom=0.17, top=0.96)

    save_all(fig, out_dir, output_base)
    plt.close(fig)

    if save_csv:
        used = pd.concat(used_frames, ignore_index=True)
        csv_path = out_dir / f"{output_base}_data_used.csv"
        used.to_csv(csv_path, index=False)
        print(f"Saved: {csv_path}")


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Figures S2 and S3 from threshold_scan_results.txt."
    )

    parser.add_argument(
        "--data",
        type=str,
        default=str(DEFAULT_DATA),
        help="Path to threshold_scan_results.txt.",
    )

    parser.add_argument(
        "--outdir",
        type=str,
        default=str(DEFAULT_OUTDIR),
        help="Output directory.",
    )

    parser.add_argument(
        "--tobs",
        type=int,
        default=DEFAULT_TOBS,
        help="Observation time. Uses total_F_t<Tobs>. Default: 100.",
    )

    parser.add_argument(
        "--omega-area",
        type=float,
        default=DEFAULT_OMEGA_AREA,
        help=(
            "Normalization area. Use 1.0 if total_F_t100 already stores Mbar; "
            "use 100.0 if it stores the raw integral over Omega=[0,10]^2."
        ),
    )

    parser.add_argument(
        "--x-left",
        type=float,
        default=DEFAULT_X_LEFT,
        help=(
            "Left x-axis limit. Default: 0.50."
        ),
    )

    parser.add_argument(
        "--x-right",
        type=float,
        default=DEFAULT_X_RIGHT,
        help=(
            "Right x-axis limit. Default: 5.80, chosen to create blank space "
            "for the lower-right legend."
        ),
    )

    parser.add_argument(
        "--x-tick-min",
        type=float,
        default=DEFAULT_X_TICK_MIN,
        help="Smallest x tick to show. Default: 0.50.",
    )

    parser.add_argument(
        "--x-tick-max",
        type=float,
        default=DEFAULT_X_TICK_MAX,
        help="Largest x tick to show. Default: 5.00.",
    )

    parser.add_argument(
        "--y-min",
        type=float,
        default=1.0e-13,
        help="Lower y-axis limit for the log scale.",
    )

    parser.add_argument(
        "--y-max",
        type=float,
        default=1.0e-1,
        help="Upper y-axis limit for the log scale.",
    )

    parser.add_argument(
        "--no-csv",
        action="store_true",
        help="Do not save the data-used CSV files.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    set_publication_style()

    data_file = Path(args.data).expanduser().resolve()
    out_dir = Path(args.outdir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    total_col = total_column_for_tobs(args.tobs)

    print(f"Results root : {RESULTS_DIR}")
    print(f"Data file    : {data_file}")
    print(f"Output dir   : {out_dir}")
    print(f"T_obs        : {args.tobs}")
    print(f"Total column : {total_col}")
    print(f"Omega area   : {args.omega_area}")
    print(f"x range      : {args.x_left} to {args.x_right}")
    print(f"x ticks      : {args.x_tick_min} to {args.x_tick_max}")

    df = read_scan_file(data_file)

    if total_col not in df.columns:
        raise ValueError(
            f"Column {total_col!r} not found. "
            f"Available total_F columns: {[c for c in df.columns if c.startswith('total_F_t')]}"
        )

    print(f"Rows loaded  : {len(df)}")
    print("Groups       :", ", ".join(sorted(df["tag_group"].dropna().unique())))
    print("Factors      :", ", ".join(f"{x:g}" for x in sorted(df["tag_factor"].dropna().unique())))

    # Figure S2: alpha variation.
    plot_group(
        df=df,
        group_name="alpha",
        parameter_label=r"\alpha",
        output_base="figS2",
        total_col=total_col,
        omega_area=args.omega_area,
        out_dir=out_dir,
        x_left=args.x_left,
        x_right=args.x_right,
        x_tick_min=args.x_tick_min,
        x_tick_max=args.x_tick_max,
        y_min=args.y_min,
        y_max=args.y_max,
        save_csv=not args.no_csv,
    )

    # Figure S3: dF variation.
    plot_group(
        df=df,
        group_name="dF",
        parameter_label=r"d_F",
        output_base="figS3",
        total_col=total_col,
        omega_area=args.omega_area,
        out_dir=out_dir,
        x_left=args.x_left,
        x_right=args.x_right,
        x_tick_min=args.x_tick_min,
        x_tick_max=args.x_tick_max,
        y_min=args.y_min,
        y_max=args.y_max,
        save_csv=not args.no_csv,
    )


if __name__ == "__main__":
    main()