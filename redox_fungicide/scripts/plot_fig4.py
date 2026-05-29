#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 4: Quantitative verification of the asymptotic results.

Panels:
    (a) Fungicide decay against theoretical exponential bounds.
    (b) Redox recovery: |1 - R(.,t)|_infty decays toward zero.
    (c) Fungal-density evolution for several beta values.
    (d) Bifurcation diagram: F_infty versus beta*theta^m.

INPUT
    results/time_series/A_R_F_timeseries/A0_<A0>_beta_<beta>.dat
    results/time_series/theorem_1_5_II/global_summary.txt
        columns:  A0  beta  case_id  beta*theta^m  F_asymptote
                  final_mean_F  final_max_F  |final_mean_F - F_asymptote|

OUTPUT (written to results/figures/)
    fig4.png, fig4.eps

USAGE
    python scripts/plot_fig4.py                 # A0=0.5, reference beta=10
    python scripts/plot_fig4.py 0.5
    python scripts/plot_fig4.py 0.5 10
"""

import sys
import re
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker


# ============================================================
# Repository-relative paths
# ============================================================
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_TS_DIR = REPO_ROOT / "results" / "time_series" / "A_R_F_timeseries"
DEFAULT_GLOBAL_SUM = REPO_ROOT / "results" / "time_series" / "theorem_1_5_II" / "global_summary.txt"
DEFAULT_OUTDIR = REPO_ROOT / "results" / "figures"


# ============================================================
# Model parameters (keep in sync with configs/parameters.h)
# ============================================================
ALPHA = 50.0
ELL = 0.3
EPS = 0.05
THETA = 0.30
M_EXP = 2.0
BETA_CRIT = 1.0 / THETA ** M_EXP


# ============================================================
# Plot style (EPS-safe)
# ============================================================
mpl.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 600,
    "font.size": 12,
    "font.family": "sans-serif",
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "axes.grid": True,
    "axes.axisbelow": False,
    "grid.color": "0.88",
    "grid.alpha": 1.0,
    "grid.linewidth": 0.45,
    "grid.linestyle": "-",
    "lines.linewidth": 1.8,
    "axes.linewidth": 0.9,
})

COLOR_A = "#1f6fbf"
COLOR_R = "#2a8a3e"
COLOR_F = "#c0392b"
COLOR_B = "#1f77b4"
COLOR_CASE_II = "#2c7fb8"
COLOR_CASE_I = "#d7301f"


def blend_with_white(color, amount):
    rgb = np.array(mcolors.to_rgb(color))
    white = np.ones(3)
    return tuple((1.0 - amount) * white + amount * rgb)


FILL_A = blend_with_white(COLOR_A, 0.18)
FILL_R = blend_with_white(COLOR_R, 0.18)
FILL_CASE_II = blend_with_white(COLOR_CASE_II, 0.08)
FILL_CASE_I = blend_with_white(COLOR_CASE_I, 0.08)
MEAN_POINT = blend_with_white(COLOR_B, 0.75)


# ============================================================
# Loaders
# ============================================================
def find_runs_for(ts_dir: Path, A0: float):
    pattern = re.compile(rf"A0_{A0:.3f}_beta_([\d.]+)\.dat")
    runs = []
    for path in ts_dir.glob(f"A0_{A0:.3f}_beta_*.dat"):
        match = pattern.match(path.name)
        if match:
            runs.append((float(match.group(1)), path))
    runs.sort(key=lambda item: item[0])
    return runs


def load_timeseries(path: Path) -> dict:
    data = np.loadtxt(path)
    return {
        "t": data[:, 0],
        "meanA": data[:, 1], "meanR": data[:, 2], "meanF": data[:, 3],
        "minA":  data[:, 4], "maxA":  data[:, 5],
        "minR":  data[:, 6], "maxR":  data[:, 7],
        "minF":  data[:, 8], "maxF":  data[:, 9],
    }


def beta_to_color(beta, beta_min, beta_max):
    if beta < BETA_CRIT:
        frac = (beta - beta_min) / max(BETA_CRIT - beta_min, 1e-12)
        return plt.cm.Blues(0.45 + 0.45 * np.clip(frac, 0.0, 1.0))
    frac = (beta - BETA_CRIT) / max(beta_max - BETA_CRIT, 1e-12)
    return plt.cm.Reds(0.45 + 0.45 * np.clip(frac, 0.0, 1.0))


def set_legend(ax, loc="best", ncol=1):
    leg = ax.legend(
        loc=loc, ncol=ncol,
        frameon=True, facecolor="white", edgecolor="0.75",
        fancybox=False, borderpad=0.45, handlelength=2.2,
    )
    leg.get_frame().set_linewidth(0.6)
    return leg


# ============================================================
# Panels
# ============================================================
def panel_a(ax, ts, A0):
    t = ts["t"]
    upper = A0 * np.exp(-ELL * t)
    lower = A0 * np.exp(-(ALPHA + ELL) * t)

    ax.semilogy(t, ts["maxA"] + 1e-30, color=COLOR_A, lw=2.0,
                label=r"$A_{\max}^{\rm num}(t)$", zorder=4)
    ax.fill_between(t, np.maximum(ts["minA"], 1e-30), ts["maxA"] + 1e-30,
                    color=FILL_A, alpha=1.0, label="min--max", zorder=1)
    ax.semilogy(t, upper + 1e-30, color="k", ls="--", lw=1.1,
                label=r"$A_0e^{-\ell t}$", zorder=5)
    ax.semilogy(t, lower + 1e-30, color="0.45", ls=":", lw=1.1,
                label=r"$A_0e^{-(\alpha+\ell)t}$", zorder=5)

    ax.set_title(r"(a)  Fungicide decay", fontweight="bold")
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"$A(x,t)$")
    ax.set_ylim(1e-12, max(2.0 * A0, 1.2))
    set_legend(ax, loc="lower right")


def panel_b(ax, ts):
    t = ts["t"]
    err_mean = np.abs(1.0 - ts["meanR"])
    err_worst = np.maximum(np.abs(1.0 - ts["minR"]), np.abs(1.0 - ts["maxR"]))
    ref = np.exp(-EPS * t)

    ax.semilogy(t, err_mean + 1e-16, color=COLOR_R, lw=2.0,
                label=r"$|1-\bar R(t)|$", zorder=4)
    ax.semilogy(t, err_worst + 1e-16, color=COLOR_R, ls=":", lw=1.4,
                label=r"$\|1-R(\cdot,t)\|_\infty$", zorder=4)
    ax.semilogy(t, ref, color="k", ls="--", lw=1.1,
                label=r"$e^{-\varepsilon t}$", zorder=3)

    ax.set_title(r"(b)  Redox recovery", fontweight="bold")
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"$|1-R(x,t)|$")
    ax.set_ylim(1e-8, 2.0)
    ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=5))
    ax.yaxis.set_minor_locator(mticker.NullLocator())
    set_legend(ax, loc="lower right")


def panel_c(ax, runs, A0):
    if not runs:
        ax.text(0.5, 0.5, "no time-series data",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    betas = [beta for beta, _ in runs]
    beta_min = min(betas)
    beta_max = max(betas)
    max_y = 0.0

    for beta, path in runs:
        ts = load_timeseries(path)
        t = ts["t"]
        Fmax = ts["maxF"]
        F_inf = max(0.0, 1.0 - beta * THETA ** M_EXP)

        color = beta_to_color(beta, beta_min, beta_max)
        ax.plot(t, Fmax, color=color, lw=1.8,
                label=fr"$\beta={beta:g}$", zorder=3)
        max_y = max(max_y, float(np.max(Fmax)), F_inf)

    ax.axhline(0.0, color="0.25", lw=0.8)
    ax.set_title(r"(c)  Fungal-density evolution", fontweight="bold")
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"$M_\infty(t)$")
    ax.set_ylim(-0.02 * max(max_y, 1e-8), 1.10 * max(max_y, 0.1))
    set_legend(ax, loc=(0.5, 0.55), ncol=2)


def panel_d(ax, A0, global_sum: Path):
    if not global_sum.exists():
        ax.text(0.5, 0.5, f"missing:\n{global_sum}",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    data = np.loadtxt(global_sum)
    if data.ndim == 1:
        data = data[None, :]

    if data.shape[1] >= 7:
        rows_for_A0 = np.isclose(data[:, 0], A0, rtol=0.0, atol=5e-4)
        if np.any(rows_for_A0):
            data = data[rows_for_A0]

    btm = data[:, 3]
    F_theory = data[:, 4]
    F_mean = data[:, 5]
    F_max = data[:, 6]

    order = np.argsort(btm)
    btm = btm[order]
    F_max = F_max[order]

    x_right = max(2.0, float(np.max(btm)) + 0.15)
    x_dense = np.linspace(0.0, x_right, 500)
    F_dense = np.maximum(0.0, 1.0 - x_dense)

    ax.axvspan(0.0, 1.0, color=FILL_CASE_II, alpha=1.0, zorder=0)
    ax.axvspan(1.0, x_right, color=FILL_CASE_I, alpha=1.0, zorder=0)

    ax.plot(x_dense, F_dense, color="k", lw=2.0,
            label=r"$\max(0,1-\beta\theta^m)$", zorder=2)
    ax.plot(btm, F_max, "o", color=COLOR_F, ms=6.5,
            mec="white", mew=0.5,
            label=r"numerical $M_\infty$", zorder=4)

    ax.axvline(1.0, color="0.35", ls="--", lw=1.0, zorder=1)
    ax.text(0.25, 0.96, "Persistence", transform=ax.transAxes,
            ha="center", va="top", fontsize=10,
            color=COLOR_CASE_II, fontweight="bold")
    ax.text(0.74, 0.96, "Extinction", transform=ax.transAxes,
            ha="center", va="top", fontsize=10,
            color=COLOR_CASE_I, fontweight="bold")

    ax.set_title(r"(d)  Bifurcation diagram", fontweight="bold")
    ax.set_xlabel(r"$\mathcal{B}=\beta\theta^m$")
    ax.set_ylabel(r"$F_\infty$")
    ax.set_xlim(0.0, x_right)
    ax.set_ylim(-0.03, 1.05)
    set_legend(ax, loc=(0.6, 0.08))


# ============================================================
# Main
# ============================================================
def main():
    A0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    beta_for_ab = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0

    ts_dir = DEFAULT_TS_DIR
    global_sum = DEFAULT_GLOBAL_SUM
    out_dir = DEFAULT_OUTDIR
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Figure 4: quantitative verification")
    print("=" * 70)
    print(f"A0 = {A0}")
    print(f"reference beta for panels (a),(b) = {beta_for_ab}")
    print(f"time-series dir : {ts_dir}")
    print(f"global summary  : {global_sum}")
    print(f"output dir      : {out_dir}")

    runs = find_runs_for(ts_dir, A0)
    if not runs:
        raise FileNotFoundError(f"No runs found for A0={A0:.3f} in {ts_dir}")
    print(f"found {len(runs)} beta runs: {[beta for beta, _ in runs]}")

    ref_path = None
    ref_beta = None
    for beta, path in runs:
        if abs(beta - beta_for_ab) < 1e-8:
            ref_beta = beta
            ref_path = path
            break
    if ref_path is None:
        ref_beta, ref_path = min(runs, key=lambda item: abs(item[0] - beta_for_ab))
        print(f"beta={beta_for_ab:g} not found; using nearest beta={ref_beta:g}")

    ts_ref = load_timeseries(ref_path)

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.0))
    fig.patch.set_facecolor("white")

    panel_a(axes[0, 0], ts_ref, A0)
    panel_b(axes[0, 1], ts_ref)
    panel_c(axes[1, 0], runs, A0)
    panel_d(axes[1, 1], A0, global_sum)

    fig.tight_layout(rect=[0.02, 0.02, 0.995, 0.995])

    png_path = out_dir / "fig4.png"
    eps_path = out_dir / "fig4.eps"

    fig.savefig(png_path, format="png", dpi=600,
                bbox_inches="tight", facecolor="white", edgecolor="none")
    fig.savefig(eps_path, format="eps",
                bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)

    print(f"saved PNG: {png_path}")
    print(f"saved EPS: {eps_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
