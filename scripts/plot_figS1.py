#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure S1: empirical fixed-time critical dose from threshold scan results.

This script reads the scan-level file

        threshold_parameter/summary/threshold_scan_results.txt

and reconstructs the first bracketing interval for the infected-area criterion
at T_obs = 100.

The scan file is expected to have the practical column order

    tag  A0  dA_factor  dF_factor  alpha_factor  dA  dF  alpha
    total_F_t20  area_t20  total_F_t50  area_t50  total_F_t100  area_t100

where the parameter group and multiplicative factor are parsed from the tag:

    dA_0.25x_A0_1.00      -> group=dA,    factor=0.25
    dF_1.50x_A0_2.00      -> group=dF,    factor=1.50
    alpha_2.00x_A0_4.80   -> group=alpha, factor=2.00

Outputs:
    figS1.png
    figS1.pdf
    figS1.eps
    figS1_thresholds_from_scan.csv
    figS1_caption_suggestion.txt

Usage:
    python plot_figS1.py

Optional:
    python plot_figS1.py --data /path/to/threshold_scan_results.txt
    python plot_figS1.py --outdir /path/to/figures
    python plot_figS1.py --criterion 1e-8 --tobs 100
"""

from __future__ import annotations

import argparse
import re
import os
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

GROUP_ORDER = ["dA", "dF", "alpha"]

GROUP_LABELS = {
    "dA": r"$d_A$",
    "dF": r"$d_F$",
    "alpha": r"$\alpha$",
}

GROUP_MARKERS = {
    "dA": "o",
    "dF": "s",
    "alpha": "^",
}

GROUP_COLORS = {
    "dA": "#1f77b4",
    "dF": "#ff7f0e",
    "alpha": "#2ca02c",
}

GROUP_X_OFFSET = {
    "dA": -0.018,
    "dF": 0.000,
    "alpha": 0.018,
}

TAG_PATTERN = re.compile(
    r"^(dA|dF|alpha)_([0-9]+(?:\.[0-9]+)?)x_A0_([0-9]+(?:\.[0-9]+)?)$"
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
        "axes.labelsize": 15,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 10.5,
        "axes.linewidth": 1.0,
        "lines.linewidth": 2.0,

        "figure.dpi": 300,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


# ============================================================
# Data loading and parsing
# ============================================================

def parse_tag(tag: str) -> tuple[str, float, float]:
    """
    Return:
        group, factor, A0_from_tag

    Example:
        dA_0.25x_A0_1.00 -> ("dA", 0.25, 1.00)
    """
    m = TAG_PATTERN.match(str(tag))
    if m is None:
        raise ValueError(f"Could not parse tag: {tag!r}")

    group = m.group(1)
    factor = float(m.group(2))
    a0_from_tag = float(m.group(3))

    return group, factor, a0_from_tag


def read_scan_file(data_file: Path) -> pd.DataFrame:
    """
    Read the scan-level data file and add group/factor columns from the tag.

    This function reads the file line-by-line rather than relying entirely on
    the legacy header, because some older scan files have a header that does
    not exactly match the practical column order.
    """
    if not data_file.exists():
        raise FileNotFoundError(f"Missing scan file: {data_file}")

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
                # Keep only the practical columns used for this plot.
                parts = parts[:len(cols)]

            rows.append(parts)

    if not rows:
        raise ValueError(f"No numerical rows were read from {data_file}")

    if skipped:
        print(f"[warning] skipped {skipped} malformed rows.")

    df = pd.DataFrame(rows, columns=cols)

    parsed = []
    bad_tags = []

    for tag in df["tag"].astype(str):
        try:
            parsed.append(parse_tag(tag))
        except ValueError:
            parsed.append((None, np.nan, np.nan))
            bad_tags.append(tag)

    if bad_tags:
        print("[warning] Unrecognized tags were skipped:")
        for tag in bad_tags[:10]:
            print("  ", tag)
        if len(bad_tags) > 10:
            print(f"  ... and {len(bad_tags) - 10} more")

    parsed_df = pd.DataFrame(parsed, columns=["group", "factor", "A0_tag"])
    df = pd.concat([df, parsed_df], axis=1)
    df = df.dropna(subset=["group", "factor", "A0"]).copy()

    # Convert numeric columns robustly.
    for col in df.columns:
        if col not in {"tag", "group"}:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Remove rows that lack A0 or the area column.
    df = df.dropna(subset=["A0"]).copy()

    # Check tag A0 against numeric A0, but do not stop on tiny formatting differences.
    mismatch = np.abs(df["A0"] - df["A0_tag"]) > 5.0e-4
    if mismatch.any():
        print(f"[warning] {int(mismatch.sum())} rows have A0 values inconsistent with the tag.")

    # Keep the last occurrence of duplicate scans.
    df["_row_order"] = np.arange(len(df))

    df = (
        df.sort_values("_row_order")
          .drop_duplicates(subset=["group", "factor", "A0"], keep="last")
          .drop(columns=["_row_order"])
          .sort_values(["group", "factor", "A0"])
          .reset_index(drop=True)
    )

    return df


# ============================================================
# Threshold reconstruction
# ============================================================

def area_column_for_tobs(tobs: int) -> str:
    col = f"area_t{int(tobs)}"
    return col


def compute_one_threshold(
    sub: pd.DataFrame,
    criterion: float,
    area_col: str,
) -> dict:
    """
    Compute the first bracketing interval for one group/factor.

    The suppression criterion is:
        area_col <= criterion

    Output status:
        found  : crossing found inside tested range
        below  : criterion is already satisfied at the smallest tested dose
        above  : criterion is not satisfied up to the largest tested dose
    """
    sub = sub.sort_values("A0").drop_duplicates(subset=["A0"], keep="last").copy()

    A0 = sub["A0"].to_numpy(dtype=float)
    area = sub[area_col].to_numpy(dtype=float)

    good = np.isfinite(A0) & np.isfinite(area)
    A0 = A0[good]
    area = area[good]

    if len(A0) == 0:
        return {
            "status": "missing",
            "A0_lower": np.nan,
            "A0_upper": np.nan,
            "threshold_mid": np.nan,
            "A0_min_tested": np.nan,
            "A0_max_tested": np.nan,
            "area_at_min": np.nan,
            "area_at_max": np.nan,
            "n_points": 0,
        }

    A0_min = float(A0[0])
    A0_max = float(A0[-1])
    area_min_dose = float(area[0])
    area_max_dose = float(area[-1])

    satisfied = area <= criterion

    if satisfied[0]:
        return {
            "status": "below",
            "A0_lower": np.nan,
            "A0_upper": A0_min,
            "threshold_mid": np.nan,
            "A0_min_tested": A0_min,
            "A0_max_tested": A0_max,
            "area_at_min": area_min_dose,
            "area_at_max": area_max_dose,
            "n_points": len(A0),
        }

    if not np.any(satisfied):
        return {
            "status": "above",
            "A0_lower": A0_max,
            "A0_upper": np.nan,
            "threshold_mid": np.nan,
            "A0_min_tested": A0_min,
            "A0_max_tested": A0_max,
            "area_at_min": area_min_dose,
            "area_at_max": area_max_dose,
            "n_points": len(A0),
        }

    # First index where the criterion is satisfied.
    idx = int(np.argmax(satisfied))

    if idx == 0:
        # Should already be handled above, but keep this for safety.
        return {
            "status": "below",
            "A0_lower": np.nan,
            "A0_upper": A0_min,
            "threshold_mid": np.nan,
            "A0_min_tested": A0_min,
            "A0_max_tested": A0_max,
            "area_at_min": area_min_dose,
            "area_at_max": area_max_dose,
            "n_points": len(A0),
        }

    lower = float(A0[idx - 1])
    upper = float(A0[idx])
    mid = 0.5 * (lower + upper)

    return {
        "status": "found",
        "A0_lower": lower,
        "A0_upper": upper,
        "threshold_mid": mid,
        "A0_min_tested": A0_min,
        "A0_max_tested": A0_max,
        "area_at_min": area_min_dose,
        "area_at_max": area_max_dose,
        "n_points": len(A0),
    }


def reconstruct_thresholds(
    df: pd.DataFrame,
    criterion: float,
    tobs: int,
) -> pd.DataFrame:
    area_col = area_column_for_tobs(tobs)

    if area_col not in df.columns:
        raise ValueError(
            f"Requested T_obs={tobs}, but column {area_col!r} is not present. "
            f"Available area columns are: {[c for c in df.columns if c.startswith('area_t')]}"
        )

    records = []

    for group in GROUP_ORDER:
        df_g = df[df["group"] == group]

        if df_g.empty:
            print(f"[warning] no rows for group {group!r}")
            continue

        factors = sorted(df_g["factor"].dropna().unique())

        for factor in factors:
            sub = df_g[np.isclose(df_g["factor"], factor)].copy()

            out = compute_one_threshold(
                sub=sub,
                criterion=criterion,
                area_col=area_col,
            )

            # Representative parameter values from first row.
            first = sub.iloc[0]

            record = {
                "group": group,
                "factor": float(factor),
                "dA_factor": float(first["dA_factor"]),
                "dF_factor": float(first["dF_factor"]),
                "alpha_factor": float(first["alpha_factor"]),
                "dA": float(first["dA"]),
                "dF": float(first["dF"]),
                "alpha": float(first["alpha"]),
                "criterion": criterion,
                "T_obs": tobs,
            }
            record.update(out)
            records.append(record)

    result = pd.DataFrame(records)

    result = (
        result.sort_values(["group", "factor"])
              .reset_index(drop=True)
    )

    return result


# ============================================================
# Plot helpers
# ============================================================

def make_group_legend_handles() -> list[Line2D]:
    handles = []

    for group in GROUP_ORDER:
        handles.append(
            Line2D(
                [0],
                [0],
                color=GROUP_COLORS[group],
                marker=GROUP_MARKERS[group],
                linestyle="-",
                linewidth=2.2,
                markersize=7.0,
                markerfacecolor=GROUP_COLORS[group],
                markeredgecolor="white",
                markeredgewidth=0.7,
                label=GROUP_LABELS[group],
            )
        )

    return handles


def make_status_legend_handles() -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker="o",
            color="0.25",
            linestyle="None",
            markersize=7.0,
            label="bracketed threshold",
        ),
        Line2D(
            [0],
            [0],
            marker="v",
            color="0.25",
            linestyle="None",
            markersize=8.0,
            label="below tested range",
        ),
        Line2D(
            [0],
            [0],
            marker="x",
            color="0.25",
            linestyle="None",
            markersize=8.0,
            markeredgewidth=1.8,
            label="above tested range",
        ),
    ]


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


def write_caption_suggestion(out_dir: Path) -> None:
    caption = r"""Suggested caption for Fig. S1:

Estimated infected-area critical dose
$\widehat A^{\mathcal I}_{0,c}$ at $T_{\rm obs}=100$
under multiplicative perturbations of $d_A$, $d_F$, and $\alpha$.
Markers show the midpoint of the first bracketing interval in the dose scan,
and vertical error bars denote the corresponding bracket.
The downward triangle indicates that the infected-area criterion was already
satisfied at the smallest tested dose. The criterion is
$\mathcal I(T_{\rm obs})\le 10^{-8}\%$, which corresponds to a
zero-detectable-cell criterion on the $500\times500$ grid.
The symbol $\times$ is reserved for cases above the tested range.
"""

    path = out_dir / "figS1_caption_suggestion.txt"
    path.write_text(caption, encoding="utf-8")
    print(f"Saved: {path}")


# ============================================================
# Plotting
# ============================================================

def plot_thresholds(
    thresholds: pd.DataFrame,
    out_dir: Path,
    y_max: Optional[float] = None,
    show_bottom_note: bool = True,
) -> None:
    if thresholds.empty:
        raise ValueError("No thresholds to plot.")

    found_df = thresholds[thresholds["status"] == "found"]

    if found_df.empty:
        raise ValueError("No bracketed thresholds were found.")

    global_min_tested = float(np.nanmin(thresholds["A0_min_tested"]))
    global_max_found = float(np.nanmax(found_df["A0_upper"]))

    if y_max is None:
        # Keep a little space above the highest found threshold.
        y_max = max(5.2, global_max_found + 0.35)

    below_y = max(0.05, global_min_tested - 0.10)
    above_y = min(y_max - 0.10, global_max_found + 0.25)

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    fig.patch.set_facecolor("white")

    # ------------------------------------------------------------
    # Plot each parameter family.
    # ------------------------------------------------------------
    for group in GROUP_ORDER:
        sub = thresholds[thresholds["group"] == group].copy()

        if sub.empty:
            continue

        sub = sub.sort_values("factor")

        color = GROUP_COLORS[group]
        marker = GROUP_MARKERS[group]
        offset = GROUP_X_OFFSET[group]

        # Bracketed thresholds.
        found = sub[sub["status"] == "found"].copy()

        if not found.empty:
            x = found["factor"].to_numpy(dtype=float) + offset
            y = found["threshold_mid"].to_numpy(dtype=float)
            y_lower = found["A0_lower"].to_numpy(dtype=float)
            y_upper = found["A0_upper"].to_numpy(dtype=float)

            yerr = np.vstack([
                np.maximum(0.0, y - y_lower),
                np.maximum(0.0, y_upper - y),
            ])

            ax.errorbar(
                x,
                y,
                yerr=yerr,
                marker=marker,
                markersize=7.2,
                linewidth=2.4,
                elinewidth=1.2,
                capsize=4,
                capthick=1.1,
                color=color,
                markerfacecolor=color,
                markeredgecolor="white",
                markeredgewidth=0.7,
                zorder=4,
            )

        # Below tested range.
        below = sub[sub["status"] == "below"].copy()

        if not below.empty:
            x = below["factor"].to_numpy(dtype=float) + offset
            ax.scatter(
                x,
                np.full_like(x, below_y, dtype=float),
                marker="v",
                s=85,
                color=color,
                edgecolors="white",
                linewidths=0.6,
                zorder=5,
            )

        # Above tested range.
        above = sub[sub["status"] == "above"].copy()

        if not above.empty:
            x = above["factor"].to_numpy(dtype=float) + offset
            ax.scatter(
                x,
                np.full_like(x, above_y, dtype=float),
                marker="x",
                s=90,
                linewidths=2.0,
                color=color,
                zorder=5,
            )

    # ------------------------------------------------------------
    # Reference line for the lowest tested dose.
    # ------------------------------------------------------------
    ax.axhline(
        global_min_tested,
        linestyle="--",
        linewidth=0.9,
        color="0.72",
        zorder=1,
    )

    ax.text(
        0.205,
        global_min_tested + 0.035,
        r"lowest tested $A_0$",
        ha="left",
        va="bottom",
        fontsize=9.5,
        color="0.35",
    )

    # ------------------------------------------------------------
    # Axes and legends.
    # ------------------------------------------------------------
    ax.set_xlabel("Multiplicative factor", fontsize=20)
    ax.set_ylabel(
        r"Critical dose $\widehat{A}^{\mathcal{I}}_{0,c}$",
        fontsize=20,
    )

    ax.set_xlim(0.20, 2.06)
    ax.set_ylim(0.0, y_max)

    ax.set_xticks([0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00])
    ax.set_yticks(np.arange(0.0, y_max + 1.0e-9, 1.0))

    ax.grid(True, which="major", color="0.88", linestyle="-", linewidth=0.6)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    # Legend 1: parameter groups.
    leg1 = ax.legend(
        handles=make_group_legend_handles(),
        loc="upper left",
        bbox_to_anchor=(0.025, 0.975),
        frameon=True,
        fancybox=False,
        facecolor="white",
        edgecolor="0.75",
        borderpad=0.45,
    )
    leg1.get_frame().set_linewidth(0.7)
    ax.add_artist(leg1)

    # Legend 2: status markers.
    leg2 = ax.legend(
        handles=make_status_legend_handles(),
        loc="lower right",
        bbox_to_anchor=(0.985, 0.03),
        frameon=True,
        fancybox=False,
        facecolor="white",
        edgecolor="0.75",
        borderpad=0.45,
    )
    leg2.get_frame().set_linewidth(0.7)

    # Short in-figure note only.  The long zero-cell explanation is written
    # to figS1_caption_suggestion.txt and should go in the manuscript caption.
    if show_bottom_note:
        fig.text(
            0.5,
            0.020,
            r"Error bars: first bracketing interval; $\blacktriangledown$ below tested range.",
            ha="center",
            va="bottom",
            fontsize=10.5,
            color="0.25",
        )
        fig.subplots_adjust(left=0.12, right=0.985, bottom=0.17, top=0.965)
    else:
        fig.subplots_adjust(left=0.12, right=0.985, bottom=0.12, top=0.965)

    save_all(fig, out_dir, "figS1")
    plt.close(fig)


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Figure S1 from threshold_scan_results.txt."
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
        help="Output directory for figS1.",
    )

    parser.add_argument(
        "--criterion",
        type=float,
        default=1.0e-8,
        help="Infected-area criterion in percent. Default: 1e-8.",
    )

    parser.add_argument(
        "--tobs",
        type=int,
        default=100,
        help="Observation time. Uses area_t<Tobs>. Default: 100.",
    )

    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="Optional y-axis maximum. Default: automatic, at least 5.2.",
    )

    parser.add_argument(
        "--no-bottom-note",
        action="store_true",
        help="Remove the short note below the figure.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    set_publication_style()

    data_file = Path(args.data).expanduser().resolve()
    out_dir = Path(args.outdir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading scan data: {data_file}")

    scan_df = read_scan_file(data_file)

    print(f"Rows loaded: {len(scan_df)}")
    print("Groups:", ", ".join(sorted(scan_df["group"].dropna().unique())))
    print("Factors:", ", ".join(f"{x:g}" for x in sorted(scan_df["factor"].dropna().unique())))

    thresholds = reconstruct_thresholds(
        scan_df,
        criterion=args.criterion,
        tobs=args.tobs,
    )

    out_csv = out_dir / "figS1_thresholds_from_scan.csv"
    thresholds.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")

    print("\nReconstructed thresholds:")
    cols_to_print = [
        "group",
        "factor",
        "status",
        "A0_lower",
        "A0_upper",
        "threshold_mid",
        "A0_min_tested",
        "A0_max_tested",
        "area_at_min",
        "area_at_max",
    ]
    print(thresholds[cols_to_print].to_string(index=False))

    plot_thresholds(
        thresholds=thresholds,
        out_dir=out_dir,
        y_max=args.ymax,
        show_bottom_note=not args.no_bottom_note,
    )

    write_caption_suggestion(out_dir)


if __name__ == "__main__":
    main()