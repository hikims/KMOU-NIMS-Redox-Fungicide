#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 6: Recovered-state front propagation.

Numerical verification of the spreading speed of Theorem 3.5,

    c_inf = 2 * sqrt( dF * (1 - B) / (1 + theta^m) ),   B = beta * theta^m.

Every other analytical result in the paper has a verification panel
(Proposition 3.3 -> Fig 4a, redox recovery -> Fig 4b, Theorem 3.4 -> Fig 4d).
This figure supplies the missing one for the spreading speed.

Panels:
    (a) Equivalent front radius r(t), the prediction r(0) + c_inf * t, and
        the Bramson-corrected linear fit.
    (b) Instantaneous speed dr/dt against c_inf, with the fit window shaded.
    (c) Front radius for three doses.  A finite dose translates the front
        trajectory in time but leaves its asymptotic slope unchanged, which
        is the content of Theorems 3.4 and 3.5 taken together.

------------------------------------------------------------------------
INPUT -- two layouts, the scalar one strongly preferred
------------------------------------------------------------------------
(1) SCALAR DIAGNOSTIC  (~160 lines per run)

    results/front_speed/front_beta_<beta:.3f>_A0_<A0:.3f>.dat

        # t     n_lev[0]  n_lev[1]  n_lev[2]
        0.0     31417     31417     31417
        1.0     ...

    n_lev[k] = number of grid cells with F >= LEVELS[k].  One column is
    enough; the extra columns enable the level-sensitivity check reported
    in the caption.

    C++ side, inside the existing output block:

        const double lev[3] = {0.02, 0.05, 0.08};
        long cnt[3] = {0, 0, 0};
        for (int i = 0; i < Nx*Ny; ++i)
            for (int k = 0; k < 3; ++k)
                if (F[i] >= lev[k]) ++cnt[k];
        fprintf(f_front, "%.6f %ld %ld %ld\\n", t, cnt[0], cnt[1], cnt[2]);

(2) SNAPSHOT FALLBACK

    results/spatial_snapshots/A0_<A0:.3f>/field_XXXX.dat,  t = XXXX / 10.
    Works, but needs the fields dumped densely in time (~1.6 GB), whereas
    the scalar diagnostic needs a few kilobytes.

------------------------------------------------------------------------
REQUIRED RUN SETTINGS -- these differ from the Figure 5 runs
------------------------------------------------------------------------
    beta     = 10          (B = 0.9)
    A0       = 0, 2.25, 5
    t_final  = 160         (Figure 5 stopped at 100 -- too short)
    dt_out   = 1.0         (Figure 5 saved 4 times -- cannot fit a slope)
    grid     = 500 x 500, h = 0.02

OUTPUT (written to results/figures/)
    fig6.png, fig6.eps
    fig6_front_speed.csv     : fitted speeds and deviations

USAGE
    python scripts/plot_fig6.py                 # beta = 10
    python scripts/plot_fig6.py --beta 8        # secondary check
    python scripts/plot_fig6.py --scales        # print Table 3 and exit
"""

import argparse
import csv
import re
import os
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


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

DEFAULT_FRONT_DIR = RESULTS_DIR / "front_speed"
DEFAULT_SNAP_DIR = RESULTS_DIR / "spatial_snapshots"


# ============================================================
# Model parameters (keep in sync with configs/parameters.h)
# ============================================================
THETA = 0.30
M_EXP = 2.0
D_F = 1.0e-3

# Discretisation
NX = NY = 500
LX = 10.0
# Grid spacing, matching src/parameters.h:  dx = Lx / (Nx - 1).
# The solver uses a node-based grid, so this is 10/499 = 0.0200401, not 0.02.
# Using 0.02 here understates the front radius, and hence the fitted speed,
# by 0.2 %.  The value is only a fallback: when the data file header carries
# dx and dy, radius_from_counts uses those instead.
H = LX / (NX - 1)
R_INIT = 2.0             # initial lesion radius
N_DIM = 2                # spatial dimension, for the curvature correction

# Super-level thresholds written by the C++ diagnostic, in column order.
LEVELS = [0.02, 0.05, 0.08]

A0_LIST = [0.0, 2.25, 5.0]

# Reference radius used to quantify the treatment-induced front delay
R_DELAY_REF = 3.0

# beta values tabulated in Table 3 of the manuscript
TABLE_BETAS = [5.0, 8.0, 10.0, 11.0]

# Saved formats.  main.tex includes fig6.eps, so EPS must stay; PDF is added
# because it is what most viewers handle best, and PNG for quick inspection.
OUTPUT_FORMATS = ("pdf", "eps", "png")


# ============================================================
# Publication style
# ============================================================
# Matched to plot_fig2.py so that the two multi-panel figures of the main
# text are typographically identical.  With FIG_WIDTH_IN below and inclusion
# at \textwidth in elsarticle 3p (~468 pt), the on-page sizes come out at
# roughly 9.4 pt for axis labels and 9.1 pt for panel titles, i.e. the same
# as the body text.
#
# font.sans-serif is pinned explicitly.  Without it matplotlib falls back to
# DejaVu Sans, which is why Figures 3 and 4 render in a different typeface
# from Figures 2 and 5 in the submitted version.
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
    "legend.fontsize": 15,
    "legend.frameon": True,
    "legend.fancybox": False,
    "legend.framealpha": 1.0,
    "legend.edgecolor": "0.80",
    "lines.linewidth": 3.0,
    "lines.markersize": 8.0,
    "axes.grid": True,
    "axes.axisbelow": False,
    "grid.color": "0.88",
    "grid.alpha": 1.0,
    "grid.linewidth": 0.8,
    "grid.linestyle": "-",
    "figure.dpi": 120,
    "savefig.dpi": 400,
    "mathtext.fontset": "dejavusans",
    # Embed fonts as Type 42 (TrueType) so the text layer survives
    # typesetting.  Verify with:  pdftotext fig6.pdf -
    "ps.fonttype": 42,
    "pdf.fonttype": 42,
})

# Width in inches.  Chosen to match plot_fig2.py, so both figures scale to
# \textwidth by the same factor and their text ends up the same size.
FIG_WIDTH_IN = 17.2
FIG_HEIGHT_IN = 5.6

COLOR_MEAS = "#c0392b"
COLOR_PRED = "k"
COLOR_FIT = "#1f6fbf"
DOSE_COLORS = ["#c0392b", "#e08214", "#1f6fbf"]


def blend_with_white(color, amount):
    """EPS-safe replacement for alpha blending (PostScript has no alpha)."""
    rgb = np.array(mcolors.to_rgb(color))
    return tuple((1.0 - amount) * np.ones(3) + amount * rgb)


SHADE_FIT = blend_with_white(COLOR_FIT, 0.10)


def set_legend(ax, loc="best", ncol=1):
    leg = ax.legend(
        loc=loc, ncol=ncol, frameon=True, framealpha=1.0,
        facecolor="white", edgecolor="0.80", fancybox=False,
        borderpad=0.4, handlelength=1.8, labelspacing=0.35,
        handletextpad=0.6,
    )
    leg.get_frame().set_linewidth(1.0)
    return leg


# ============================================================
# Analytic predictions (Table 3 of the manuscript)
# ============================================================
def predictions(beta: float) -> dict:
    B = beta * THETA ** M_EXP
    if B >= 1.0:
        raise ValueError(
            f"beta = {beta:g} gives B = {B:.3f} >= 1: the recovered state is "
            f"extinction and no front exists."
        )

    a_inf = (1.0 - B) / (1.0 + THETA ** M_EXP)
    c_inf = 2.0 * np.sqrt(D_F * a_inf)
    lam_star = np.sqrt(a_inf / D_F)
    width = np.sqrt(D_F / a_inf)

    return {
        "beta": beta,
        "B": B,
        "a_inf": a_inf,
        "c_inf": c_inf,
        "lam_star": lam_star,
        "bramson": 3.0 / (2.0 * lam_star),
        "width": width,
        "width_cells": width / H,
        "t_wall": (0.5 * LX - R_INIT) / c_inf,
        "F_lev": 0.5 * (1.0 - B),
    }


def print_scale_table() -> None:
    """Reproduce Table 3 of the manuscript."""
    print(f"{'beta':>6} {'B':>7} {'a_inf':>10} {'c_inf':>10} "
          f"{'width/h':>9} {'t_wall':>8}")
    for beta in TABLE_BETAS:
        p = predictions(beta)
        print(f"{p['beta']:6.0f} {p['B']:7.2f} {p['a_inf']:10.4f} "
              f"{p['c_inf']:10.5f} {p['width_cells']:9.1f} {p['t_wall']:8.0f}")
    print("\nOnly beta with a resolved front (width/h >~ 4) and t_wall inside "
          "the simulated window\nare usable. beta = 10 is the primary case; "
          "beta = 8 is a truncated secondary check.")


def fit_window(pred: dict):
    """Fit interval: inside the domain, past the initial transient."""
    t_hi = min(0.9 * pred["t_wall"], 140.0)
    t_lo = max(40.0, 0.25 * t_hi)
    return t_lo, t_hi


# ============================================================
# Loading
# ============================================================
def radius_from_counts(counts, cell_area=None) -> np.ndarray:
    """
    Equivalent radius of a super-level set occupying `counts` cells.

    `cell_area` should be dx*dy taken from the file header.  H is used only
    for snapshot-fallback data, which carries no header.
    """
    if cell_area is None:
        cell_area = H * H
    return np.sqrt(np.asarray(counts, dtype=float) * cell_area / np.pi)


def series_radius(rec: dict, col: int) -> np.ndarray:
    """Radius for one run: use the simulation's own column when present."""
    if rec.get("radii") is not None and rec["radii"].shape[1] > col:
        return rec["radii"][:, col]

    meta = rec.get("meta", {})
    dx = meta.get("dx")
    dy = meta.get("dy", dx)
    cell_area = (dx * dy) if dx else None
    c_idx = min(col, rec["counts"].shape[1] - 1)
    return radius_from_counts(rec["counts"][:, c_idx], cell_area)


def parse_header(path: Path) -> dict:
    """
    Read the '# key = value' header written by main6.cpp.

    The grid spacing and the threshold levels come from the simulation, not
    from constants duplicated here, so this script cannot disagree with the
    run that produced the data.
    """
    meta = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            body = line[1:].strip()
            if "=" not in body:
                continue
            key, _, val = body.partition("=")
            key = key.strip()
            val = val.split("#")[0].strip()
            try:
                nums = [float(v) for v in val.split()]
            except ValueError:
                meta[key] = val
                continue
            meta[key] = nums[0] if len(nums) == 1 else nums
    return meta


def load_front_scalar(front_dir: Path, beta: float, A0: float):
    path = front_dir / f"front_beta_{beta:.3f}_A0_{A0:.3f}.dat"
    if not path.exists():
        return None

    meta = parse_header(path)
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]

    levels = meta.get("levels")
    if levels is None:
        # Pre-header file: counts only, thresholds unknown.
        return {"t": data[:, 0], "counts": data[:, 1:],
                "radii": None, "meta": meta, "source": str(path)}

    if not isinstance(levels, list):
        levels = [levels]
    n_lev = len(levels)

    counts = data[:, 1:1 + n_lev]
    radii = (data[:, 1 + n_lev:1 + 2 * n_lev]
             if data.shape[1] >= 1 + 2 * n_lev else None)

    return {"t": data[:, 0], "counts": counts, "radii": radii,
            "meta": meta, "source": str(path)}


def load_front_fields(snap_dir: Path, A0: float, F_lev: float):
    folder = snap_dir / f"A0_{A0:.3f}"
    if not folder.is_dir():
        return None

    rx = re.compile(r"field_(\d+)\.dat$")
    records = []

    for path in sorted(folder.glob("field_*.dat")):
        m = rx.search(path.name)
        if m is None:
            continue
        t = int(m.group(1)) / 10.0
        F = np.loadtxt(path)[:, 3]
        records.append((t, float(np.count_nonzero(F >= F_lev))))

    if len(records) < 3:
        return None

    records.sort()
    return {
        "t": np.array([r[0] for r in records]),
        "counts": np.array([[r[1]] for r in records]),
        "radii": None,
        "meta": {},
        "source": str(folder / "field_*.dat"),
    }


def load_front(front_dir, snap_dir, beta, A0, F_lev):
    rec = load_front_scalar(front_dir, beta, A0)
    if rec is not None:
        return rec

    rec = load_front_fields(snap_dir, A0, F_lev)
    if rec is not None:
        print(f"  [fallback] A0 = {A0:g}: counted from field snapshots "
              f"({rec['t'].size} times). The scalar diagnostic is preferred.")
        return rec

    return None


def primary_column(pred: dict, meta: dict) -> int:
    """
    Column of the threshold closest to F_lev = (1 - B)/2.

    main6.cpp records `primary_level_index` explicitly; fall back to matching
    against its `levels` list, and to the module default only for files that
    predate the header.
    """
    if "primary_level_index" in meta:
        return int(meta["primary_level_index"])

    levels = meta.get("levels", LEVELS)
    if not isinstance(levels, list):
        levels = [levels]
    return int(np.argmin(np.abs(np.array(levels) - pred["F_lev"])))


# ============================================================
# Speed estimation
# ============================================================
def bramson_fit(t, r, pred, t_lo, t_hi):
    """
    For compactly supported data the KPP front obeys

        r(t) = c_inf * t - (3 / 2 lambda*) log t + O(1),

    so a raw finite-time slope underestimates c_inf.  Fitting
    r(t) + (3 / 2 lambda*) log t against t removes that bias.
    """
    mask = (t >= t_lo) & (t <= t_hi) & (r > 0.0) & (t > 0.0)
    if np.count_nonzero(mask) < 4:
        return None, mask

    y = r[mask] + pred["bramson"] * np.log(t[mask])
    c, b = np.polyfit(t[mask], y, 1)
    return {"c": float(c), "intercept": float(b)}, mask


def raw_slope(t, r, t_lo, t_hi):
    """Uncorrected slope, reported to show the size of the Bramson bias."""
    mask = (t >= t_lo) & (t <= t_hi) & (r > 0.0)
    if np.count_nonzero(mask) < 4:
        return None
    c, _ = np.polyfit(t[mask], r[mask], 1)
    return float(c)


def predicted_rdot(t, r, pred):
    """
    Finite-time front speed including the two known corrections:

        rdot(t) = c_inf - (3 / 2 lambda*) / t - (n-1) dF / r(t)

    The second term is the Bramson logarithmic delay, the third the radial
    curvature correction.

    Both are documented in Section 5.6.  Plotting this instead of the bare
    horizontal line c_inf is what turns the residual gap from an unexplained
    discrepancy into a quantitative match.
    """
    out = np.full_like(np.asarray(t, dtype=float), np.nan)
    ok = (np.asarray(t) > 0) & (np.asarray(r) > 0)
    tt = np.asarray(t, dtype=float)[ok]
    rr = np.asarray(r, dtype=float)[ok]
    out[ok] = pred["c_inf"] - pred["bramson"] / tt - (N_DIM - 1) * D_F / rr
    return out


def corrected_speed_estimate(t, r, pred, t_lo, t_hi):
    """
    Speed estimate with both finite-time corrections removed:

        c = < rdot(t) + (3 / 2 lambda*) / t + (n-1) dF / r(t) >

    averaged over the fit window.  Unlike the Bramson-only linear fit this
    also removes the curvature term, which for a radial front does not vanish
    and biases the linear fit low by about dF / <r>.
    """
    good = r > 0.0
    if np.count_nonzero(good) < 4:
        return None

    rdot = np.full_like(r, np.nan)
    rdot[good] = np.gradient(r[good], t[good])

    mask = (t >= t_lo) & (t <= t_hi) & good & (t > 0)
    if np.count_nonzero(mask) < 4:
        return None

    vals = (rdot[mask]
            + pred["bramson"] / t[mask]
            + (N_DIM - 1) * D_F / r[mask])
    return {"c": float(np.mean(vals)), "std": float(np.std(vals))}


def crossing_time(t, r, r_ref):
    """First time at which r(t) reaches r_ref, by linear interpolation."""
    above = np.where(r >= r_ref)[0]
    if above.size == 0 or above[0] == 0:
        return None
    i = above[0]
    r0, r1 = r[i - 1], r[i]
    if r1 == r0:
        return float(t[i])
    return float(t[i - 1] + (r_ref - r0) * (t[i] - t[i - 1]) / (r1 - r0))


def level_sensitivity(rec, pred, t_lo, t_hi):
    """Fitted speed for each available threshold, for the caption."""
    counts = rec["counts"]
    if counts.shape[1] < 2:
        return None

    meta = rec.get("meta", {})
    levels = meta.get("levels", LEVELS)
    if not isinstance(levels, list):
        levels = [levels]

    out = []
    for k in range(counts.shape[1]):
        lev = levels[k] if k < len(levels) else float("nan")
        fit_k, _ = bramson_fit(rec["t"], series_radius(rec, k),
                               pred, t_lo, t_hi)
        if fit_k is not None:
            out.append((lev, fit_k["c"]))
    return out


def report_grid(meta: dict, pred: dict) -> None:
    """Cross-check the file header against the constants in this script."""
    if not meta:
        print("  [note] no header in the data file; falling back to the "
              "constants defined in this script")
        return

    dx = meta.get("dx")
    if dx is not None:
        print(f"  grid (file)   = {int(meta.get('Nx', 0))} x "
              f"{int(meta.get('Ny', 0))}, dx = {dx:.6f}, "
              f"Lx = {meta.get('Lx', float('nan')):.4f}")
        if abs(dx - H) > 1e-9:
            print(f"  [note] the run used dx = {dx:.6f}, not H = {H:.6f} as set "
                  f"in this script. The file value is used.")

    c_file = meta.get("c_inf")
    if c_file is not None and abs(c_file - pred["c_inf"]) > 1e-9:
        print(f"  [warning] c_inf from the file ({c_file:.6f}) differs from "
              f"the value recomputed here ({pred['c_inf']:.6f}). Check that "
              f"theta, m and dF match configs/parameters.h.")


# ============================================================
# Panels
# ============================================================
def panel_a(ax, t, r, pred, fit, mask):
    ax.plot(t, r, color=COLOR_MEAS, lw=3.0,
            label=r"measured", zorder=4)

    t_line = np.linspace(t.min(), t.max(), 300)
    ax.plot(t_line, R_INIT + pred["c_inf"] * t_line,
            color=COLOR_PRED, ls="--", lw=2.0,
            label=r"$r(0)+c_\infty t$", zorder=5)

    if fit is not None:
        t_fit = t[mask]
        r_fit = (fit["c"] * t_fit + fit["intercept"]
                 - pred["bramson"] * np.log(t_fit))
        ax.plot(t_fit, r_fit, color=COLOR_FIT, ls=":", lw=2.6,
                label=rf"fit  $c={fit['c']:.5f}$", zorder=6)

    ax.axhline(0.5 * LX, color="0.45", ls="-.", lw=1.4,
               label=r"domain edge")

    ax.set_title(r"(a)  Front radius", fontweight="bold")
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"$r(t)$")
    ax.set_ylim(0.0, 0.5 * LX * 1.08)
    set_legend(ax, loc="upper left")


def panel_b(ax, t, r, pred, fit, t_lo, t_hi, t_start=5.0):
    good = r > 0.0
    rdot = np.full_like(r, np.nan)
    if np.count_nonzero(good) > 2:
        rdot[good] = np.gradient(r[good], t[good])

    # The first few time units are front formation, not propagation: the
    # patch relaxes from F = 0.5 to the plateau 1 - B before a KPP front
    # exists.  Starting the axis at t_start keeps that spike from squashing
    # the range that matters.
    shown = t >= t_start

    ax.plot(t[shown], rdot[shown], color=COLOR_MEAS, lw=2.6,
            label=r"measured", zorder=4)

    # Below the formation time the prediction is meaningless (r is still
    # collapsing, so dF/r blows up); clip it rather than let it run off the
    # bottom of the panel.
    pred_rdot = predicted_rdot(t, r, pred)
    visible = shown & np.isfinite(pred_rdot) & (pred_rdot > 0.35 * pred["c_inf"])
    ax.plot(t[visible], pred_rdot[visible],
            color="#2a8a3e", ls="-", lw=2.4,
            label=r"finite-$t$ prediction", zorder=5)

    ax.axhline(pred["c_inf"], color=COLOR_PRED, ls="--", lw=2.0,
               label=rf"$c_\infty={pred['c_inf']:.5f}$", zorder=5)

    ax.axvspan(t_lo, t_hi, color=SHADE_FIT, zorder=0)
    # Label the shaded band at the top.  At the bottom it collided with the
    # legend, which sits lower right because that is the only empty corner:
    # the curves rise from lower left and c_inf caps the top.
    ax.text(0.5 * (t_lo + t_hi), 0.965, "fit window",
            transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=15, color="0.30")

    ax.set_title(r"(b)  Instantaneous speed", fontweight="bold")
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"$\dot r(t)$")
    ax.set_xlim(0.0, t.max() * 1.02)

    finite = rdot[shown][np.isfinite(rdot[shown])]
    lo = min(finite.min(), 0.6 * pred["c_inf"]) if finite.size else 0.0
    ax.set_ylim(max(0.0, 0.9 * lo), 1.12 * pred["c_inf"])
    set_legend(ax, loc="lower right")


def choose_delay_reference(series):
    """
    Largest round radius every dose actually reaches.

    A fixed reference (r = 3) is not reached by the highest dose within the
    simulated window, so its delay would simply be missing from the figure.
    """
    reachable = [float(np.max(r)) for _, _, r in series if np.size(r)]
    if not reachable:
        return R_DELAY_REF

    cap = min(reachable) - 0.05
    for cand in (3.5, 3.0, 2.75, 2.5, 2.25, 2.0):
        if cand <= cap:
            return cand
    return max(0.5 * cap, 0.1)


def panel_c(ax, series, r_ref_level):
    delays = {}
    t_ref0 = None

    for k, (A0, t, r) in enumerate(series):
        color = DOSE_COLORS[k % len(DOSE_COLORS)]
        emerged = r > 0.0

        # While the super-level set is empty the lesion is suppressed below
        # the measurement threshold; r is reported as zero there.
        ax.plot(t[~emerged], np.zeros(np.count_nonzero(~emerged)),
                color=color, lw=1.6, ls=":", zorder=2)
        ax.plot(t[emerged], r[emerged], color=color, lw=3.0,
                label=rf"$A_0={A0:g}$", zorder=4)

        tc = crossing_time(t, r, r_ref_level)
        if tc is not None:
            if k == 0:
                t_ref0 = tc
            elif t_ref0 is not None:
                delays[A0] = tc - t_ref0

    ax.axhline(r_ref_level, color="0.45", ls="--",
               lw=1.6, label=rf"$r={r_ref_level:g}$ reference")

    ax.set_title(r"(c)  Dose-induced front delay", fontweight="bold")
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"$r(t)$")
    ax.set_ylim(0.0, 0.5 * LX * 1.08)
    set_legend(ax, loc="upper left")

    return delays


# ============================================================
# Main
# ============================================================
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Plot Figure 6: recovered-state front speed.")
    ap.add_argument("--beta", type=float, default=10.0,
                    help="beta of the run to analyse (default: 10)")
    ap.add_argument("--frontdir", type=str, default=None,
                    help=f"Scalar diagnostic directory (default: {DEFAULT_FRONT_DIR})")
    ap.add_argument("--snapdir", type=str, default=None,
                    help=f"Snapshot directory for the fallback "
                         f"(default: {DEFAULT_SNAP_DIR})")
    ap.add_argument("--outdir", type=str, default=None,
                    help=f"Output directory (default: {DEFAULT_OUTDIR})")
    ap.add_argument("--scales", action="store_true",
                    help="print the Table 3 front scales and exit")
    return ap.parse_args()


def main():
    args = parse_args()

    if args.scales:
        print_scale_table()
        return

    beta = args.beta
    pred = predictions(beta)
    t_lo, t_hi = fit_window(pred)

    front_dir = Path(args.frontdir).expanduser() if args.frontdir else DEFAULT_FRONT_DIR
    snap_dir = Path(args.snapdir).expanduser() if args.snapdir else DEFAULT_SNAP_DIR
    out_dir = Path(args.outdir).expanduser() if args.outdir else DEFAULT_OUTDIR
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("Figure 6: recovered-state front speed")
    print("=" * 72)
    print(f"  beta          = {pred['beta']:g}")
    print(f"  B             = {pred['B']:.4f}")
    print(f"  a_inf         = {pred['a_inf']:.6f}")
    print(f"  c_inf         = {pred['c_inf']:.6f}")
    print(f"  front width   = {pred['width']:.4f}  "
          f"({pred['width_cells']:.1f} cells at h = {H:.6f})")
    print(f"  t_wall        = {pred['t_wall']:.1f}")
    print(f"  F_lev         = {pred['F_lev']:.4f}   "
          f"# 0.5*(1-B), Section 5.6")
    print(f"  fit window    = [{t_lo:.0f}, {t_hi:.0f}]")
    print(f"  Bramson coef  = {pred['bramson']:.4f}")
    curv = (N_DIM - 1) * D_F / R_DELAY_REF
    print(f"  curvature     = {curv:.2e}  "
          f"({100 * curv / pred['c_inf']:.1f}% of c_inf, negligible)")
    print(f"  front data    = {front_dir}")
    print(f"  output        = {out_dir}  "
          f"({', '.join(OUTPUT_FORMATS)})")

    if pred["width_cells"] < 4.0:
        print(f"  [warning] the front spans only {pred['width_cells']:.1f} "
              f"cells; the speed is not resolved at h = {H:.6f}")

    # ---- load ------------------------------------------------------------
    series = []
    ref_rec = None
    col = None

    for A0 in A0_LIST:
        rec = load_front(front_dir, snap_dir, beta, A0, pred["F_lev"])
        if rec is None:
            print(f"  [missing] no front data for beta = {beta:g}, A0 = {A0:g}")
            continue
        if col is None:
            meta = rec.get("meta", {})
            col = primary_column(pred, meta)
            report_grid(meta, pred)
            lv = meta.get("levels", LEVELS)
            if not isinstance(lv, list):
                lv = [lv]
            if col < len(lv):
                print(f"  front level   = {lv[col]:.4f}  (column {col} of "
                      f"{len(lv)})")
        series.append((A0, rec["t"], series_radius(rec, col)))
        if ref_rec is None:
            ref_rec = rec

    if not series:
        raise FileNotFoundError(
            "No front data found.\n"
            f"  looked for: front_beta_{beta:.3f}_A0_*.dat\n"
            f"  in        : {front_dir}\n"
            f"  results root: {RESULTS_DIR}\n"
            f"  snapshot fallback: {snap_dir}/A0_*/field_*.dat\n\n"
            "  If the files are somewhere else, pass the directory directly:\n"
            "    python plot_fig6.py --frontdir /path/to/front_speed\n"
            "  A beta mismatch also looks like this: the file name encodes\n"
            "  beta to three decimals, so --beta must match the run."
        )

    A0_ref, t_ref, r_ref = series[0]

    if t_ref.max() < t_hi:
        print(f"  [warning] data end at t = {t_ref.max():.0f} but the fit "
              f"window needs t = {t_hi:.0f}. Re-run with a longer t_final.")
    if t_ref.size < 20:
        print(f"  [warning] only {t_ref.size} output times; dt_out = 1.0 is "
              f"recommended for a slope fit.")

    # ---- fit -------------------------------------------------------------
    fit, mask = bramson_fit(t_ref, r_ref, pred, t_lo, t_hi)
    rows = []

    if fit is not None:
        dev = 100.0 * (fit["c"] - pred["c_inf"]) / pred["c_inf"]
        print(f"\n  fitted c (Bramson-corrected) = {fit['c']:.6f}  "
              f"({dev:+.2f}% vs c_inf)")
        rows.append({"quantity": "c_fit_bramson", "value": fit["c"],
                     "percent_vs_c_inf": dev})

        c_raw = raw_slope(t_ref, r_ref, t_lo, t_hi)
        if c_raw is not None:
            dev_raw = 100.0 * (c_raw - pred["c_inf"]) / pred["c_inf"]
            print(f"  raw slope (no correction)    = {c_raw:.6f}  "
                  f"({dev_raw:+.2f}% vs c_inf)")
            rows.append({"quantity": "c_fit_raw", "value": c_raw,
                         "percent_vs_c_inf": dev_raw})

    corr = corrected_speed_estimate(t_ref, r_ref, pred, t_lo, t_hi)
    if corr is not None:
        dev_c = 100.0 * (corr["c"] - pred["c_inf"]) / pred["c_inf"]
        print(f"  both corrections removed     = {corr['c']:.6f}  "
              f"({dev_c:+.2f}% vs c_inf,  sd {corr['std']:.6f})")
        rows.append({"quantity": "c_bramson_plus_curvature",
                     "value": corr["c"], "percent_vs_c_inf": dev_c})

    rows.append({"quantity": "c_inf_theory", "value": pred["c_inf"],
                 "percent_vs_c_inf": 0.0})

    sens = level_sensitivity(ref_rec, pred, t_lo, t_hi)
    if sens:
        print("\n  level sensitivity of the fitted speed:")
        for lev, c_lev in sens:
            d = 100.0 * (c_lev - pred["c_inf"]) / pred["c_inf"]
            print(f"    F_lev = {lev:<6g} ->  c = {c_lev:.6f}  ({d:+.2f}%)")
            rows.append({"quantity": f"c_fit_Flev_{lev:g}", "value": c_lev,
                         "percent_vs_c_inf": d})

    # ---- figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
    fig.patch.set_facecolor("white")

    panel_a(axes[0], t_ref, r_ref, pred, fit, mask)
    panel_b(axes[1], t_ref, r_ref, pred, fit, t_lo, t_hi)
    r_ref_level = choose_delay_reference(series)
    delays = panel_c(axes[2], series, r_ref_level)

    if delays:
        print(f"\n  front delay at r = {r_ref_level:g} "
              f"(relative to A0 = {A0_ref:g}):")
        for A0, dt in delays.items():
            print(f"    A0 = {A0:<5g} ->  delta t = {dt:+.1f}")
            rows.append({"quantity": f"delay_A0_{A0:g}", "value": dt,
                         "percent_vs_c_inf": float("nan")})

    fig.tight_layout(rect=[0.01, 0.02, 0.99, 0.98])

    saved = []
    for fmt in OUTPUT_FORMATS:
        path = out_dir / f"fig6.{fmt}"
        kwargs = dict(format=fmt, bbox_inches="tight",
                      facecolor="white", edgecolor="none")
        if fmt == "png":
            kwargs["dpi"] = 600
        fig.savefig(path, **kwargs)
        saved.append(path)

    plt.close(fig)

    csv_path = out_dir / "fig6_front_speed.csv"

    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["quantity", "value", "percent_vs_c_inf"])
        writer.writeheader()
        writer.writerows(rows)

    print("\nSaved:")
    for path in saved:
        print(" -", path)
    print(" -", csv_path)
    print("\nInsert the Bramson-corrected c and its deviation into the two")
    print("\\todo{} marks in Section 5.6 of main.tex, then delete the macro.")
    print("=" * 72)


if __name__ == "__main__":
    main()
