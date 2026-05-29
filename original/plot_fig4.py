"""
Figure 3: Quantitative verification of the asymptotic results.

Panels:
  (a) Fungicide decay and theoretical exponential bounds.
  (b) Redox recovery, shown through the deviation from R_infty = 1.
  (c) Fungal-density evolution for several beta values.
  (d) Bifurcation diagram: F_infty versus beta*theta^m.

Reads from:
    /Users/mkim/Desktop/Research/fungicide_jh/time_series/
        A_R_F_timeseries/A0_<A0>_beta_<beta>.dat
        theorem_1_5_II/global_summary.txt

Saves to:
    /Users/mkim/Desktop/Research/fungicide_jh/figures/fig3.png
    /Users/mkim/Desktop/Research/fungicide_jh/figures/fig3.eps

Usage:
    python fig3.py
    python fig3.py 0.5
    python fig3.py 0.5 10
"""

import sys
import re
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker


# -------------------------------------------------------------
# Paths -- must match main4.cpp
# -------------------------------------------------------------
ROOT = Path("/Users/mkim/Desktop/Research/fungicide_jh")

TS_AFR = ROOT / "time_series" / "A_R_F_timeseries"
TS_THEOREM = ROOT / "time_series" / "theorem_1_5_II"
GLOBAL_SUM = TS_THEOREM / "global_summary.txt"

FIG_ROOT = ROOT / "figures"
FIG_ROOT.mkdir(exist_ok=True)


# -------------------------------------------------------------
# Model parameters from parameters.h
# -------------------------------------------------------------
ALPHA = 50.0
ELL = 0.3
EPS = 0.05
THETA = 0.30
M_EXP = 2.0

BETA_CRIT = 1.0 / THETA**M_EXP


# -------------------------------------------------------------
# Publication style, EPS-safe
# -------------------------------------------------------------
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
    "grid.alpha": 1.0,       # EPS-safe: no transparency
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
    """EPS-safe replacement for alpha blending on a white background.

    EPS/PostScript does not support transparency.  Instead of plotting
    color with alpha=amount, we use the opaque color obtained by blending
    it with white in advance.
    """
    rgb = np.array(mcolors.to_rgb(color))
    white = np.ones(3)
    return tuple((1.0 - amount) * white + amount * rgb)


FILL_A = blend_with_white(COLOR_A, 0.18)
FILL_R = blend_with_white(COLOR_R, 0.18)
FILL_CASE_II = blend_with_white(COLOR_CASE_II, 0.08)
FILL_CASE_I = blend_with_white(COLOR_CASE_I, 0.08)
MEAN_POINT = blend_with_white(COLOR_B, 0.75)


# =============================================================
# Loaders
# =============================================================
def find_runs_for(A0):
    """Find all (beta, dat_path) pairs on disk for this A0."""
    pattern = re.compile(rf"A0_{A0:.3f}_beta_([\d.]+)\.dat")
    runs = []

    for path in TS_AFR.glob(f"A0_{A0:.3f}_beta_*.dat"):
        match = pattern.match(path.name)
        if match:
            runs.append((float(match.group(1)), path))

    runs.sort(key=lambda item: item[0])
    return runs


def load_timeseries(path):
    """Load one A-R-F time-series file."""
    data = np.loadtxt(path)
    return {
        "t": data[:, 0],
        "meanA": data[:, 1],
        "meanR": data[:, 2],
        "meanF": data[:, 3],
        "minA": data[:, 4],
        "maxA": data[:, 5],
        "minR": data[:, 6],
        "maxR": data[:, 7],
        "minF": data[:, 8],
        "maxF": data[:, 9],
    }


def beta_to_color(beta, beta_min, beta_max):
    """Blue scale for beta < beta_c and red scale for beta > beta_c."""
    if beta < BETA_CRIT:
        frac = (beta - beta_min) / max(BETA_CRIT - beta_min, 1e-12)
        return plt.cm.Blues(0.45 + 0.45 * np.clip(frac, 0.0, 1.0))

    frac = (beta - BETA_CRIT) / max(beta_max - BETA_CRIT, 1e-12)
    return plt.cm.Reds(0.45 + 0.45 * np.clip(frac, 0.0, 1.0))


def set_legend(ax, loc="best", ncol=1):
    """Consistent opaque legend for both PNG and EPS."""
    leg = ax.legend(
        loc=loc,
        ncol=ncol,
        frameon=True,
        facecolor="white",
        edgecolor="0.75",
        fancybox=False,
        borderpad=0.45,
        handlelength=2.2,
    )
    leg.get_frame().set_linewidth(0.6)
    return leg


# =============================================================
# Panel (a): fungicide decay
# =============================================================
def panel_a(ax, ts, A0):
    t = ts["t"]

    upper = A0 * np.exp(-ELL * t)
    lower = A0 * np.exp(-(ALPHA + ELL) * t)

    ax.semilogy(
        t,
        ts["maxA"] + 1e-30,
        color=COLOR_A,
        lw=2.0,
        label=r"$A_{\max}^{\rm num}(t)$",
        zorder=4,
    )
    ax.fill_between(
        t,
        np.maximum(ts["minA"], 1e-30),
        ts["maxA"] + 1e-30,
        color=FILL_A,
        alpha=1.0,
        label="min--max",
        zorder=1,
    )
    ax.semilogy(
        t,
        upper + 1e-30,
        color="k",
        ls="--",
        lw=1.1,
        label=r"$A_0e^{-\ell t}$",
        zorder=5,
    )
    ax.semilogy(
        t,
        lower + 1e-30,
        color="0.45",
        ls=":",
        lw=1.1,
        label=r"$A_0e^{-(\alpha+\ell)t}$",
        zorder=5,
    )

    ax.set_title(r"(a)  Fungicide decay", fontweight="bold")
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"$A(x,t)$")
    ax.set_ylim(1e-12, max(2.0 * A0, 1.2))
    set_legend(ax, loc="lower right")


# =============================================================
# Panel (b): redox recovery
# =============================================================
def panel_b(ax, ts):
    t = ts["t"]

    err_mean = np.abs(1.0 - ts["meanR"])
    err_worst = np.maximum(
        np.abs(1.0 - ts["minR"]),
        np.abs(1.0 - ts["maxR"]),
    )
    ref = np.exp(-EPS * t)

    ax.semilogy(
        t,
        err_mean + 1e-16,
        color=COLOR_R,
        lw=2.0,
        label=r"$|1-\bar R(t)|$",
        zorder=4,
    )
    ax.semilogy(
        t,
        err_worst + 1e-16,
        color=COLOR_R,
        ls=":",
        lw=1.4,
        label=r"$\|1-R(\cdot,t)\|_\infty$",
        zorder=4,
    )
    ax.semilogy(
        t,
        ref,
        color="k",
        ls="--",
        lw=1.1,
        label=r"$e^{-\varepsilon t}$",
        zorder=3,
    )

    ax.set_title(r"(b)  Redox recovery", fontweight="bold")
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"$|1-R(x,t)|$")
    ax.set_ylim(1e-8, 2.0)
    
    # y축 라벨 깔끔하게: LogLocator로 주요 틱만 표시
    ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=5))
    ax.yaxis.set_minor_locator(mticker.NullLocator())
    
    set_legend(ax, loc="lower right")


# =============================================================
# Panel (c): F evolution for multiple beta values
# =============================================================
def panel_c(ax, runs, A0):
    if not runs:
        ax.text(
            0.5,
            0.5,
            "no time-series data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
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
        F_inf = max(0.0, 1.0 - beta * THETA**M_EXP)

        color = beta_to_color(beta, beta_min, beta_max)
        ax.plot(
            t,
            Fmax,
            color=color,
            lw=1.8,
            label=fr"$\beta={beta:g}$",
            zorder=3,
        )

        # 점선 제거됨 (175-200 범위의 점선)

        max_y = max(max_y, float(np.max(Fmax)), F_inf)

    ax.axhline(0.0, color="0.25", lw=0.8)
    ax.set_title(r"(c)  Fungal-density evolution", fontweight="bold")
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"$M_\infty(t)$")
    ax.set_ylim(-0.02 * max(max_y, 1e-8), 1.10 * max(max_y, 0.1))
    
    # 레전드: 4행 2열, 우측 밑에서 75% 지점 (0.75, 0.25)
    set_legend(ax, loc=(0.5, 0.55), ncol=2)


# =============================================================
# Panel (d): bifurcation diagram
# =============================================================
def panel_d(ax, A0):
    if not GLOBAL_SUM.exists():
        ax.text(
            0.5,
            0.5,
            f"missing:\n{GLOBAL_SUM}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        return

    data = np.loadtxt(GLOBAL_SUM)
    if data.ndim == 1:
        data = data[None, :]

    # Expected columns:
    # A0, beta, case_id, beta*theta^m, F_asymptote,
    # final_mean_F, final_max_F, error
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
    F_theory = F_theory[order]
    F_mean = F_mean[order]
    F_max = F_max[order]

    x_right = max(2.0, float(np.max(btm)) + 0.15)
    x_dense = np.linspace(0.0, x_right, 500)
    F_dense = np.maximum(0.0, 1.0 - x_dense)

    ax.axvspan(0.0, 1.0, color=FILL_CASE_II, alpha=1.0, zorder=0)
    ax.axvspan(1.0, x_right, color=FILL_CASE_I, alpha=1.0, zorder=0)

    ax.plot(
        x_dense,
        F_dense,
        color="k",
        lw=2.0,
        label=r"$\max(0,1-\beta\theta^m)$",
        zorder=2,
    )
    ax.plot(
        btm,
        F_max,
        "o",
        color=COLOR_F,
        ms=6.5,
        mec="white",
        mew=0.5,
        label=r"numerical $M_\infty$",
        zorder=4,
    )
    # ax.plot(
    #     btm,
    #     F_mean,
    #     "s",
    #     color=MEAN_POINT,
    #     ms=5.8,
    #     mec="white",
    #     mew=0.5,
    #     label=r"numerical $\bar F$",
    #     zorder=3,
    # )

    ax.axvline(1.0, color="0.35", ls="--", lw=1.0, zorder=1)
    ax.text(
        0.25,
        0.96,
        "Persistence",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10,
        color=COLOR_CASE_II,
        fontweight="bold",
    )
    ax.text(
        0.74,
        0.96,
        "Extinction",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10,
        color=COLOR_CASE_I,
        fontweight="bold",
    )

    ax.set_title(r"(d)  Bifurcation diagram", fontweight="bold")
    ax.set_xlabel(r"$\mathcal{B}=\beta\theta^m$")
    ax.set_ylabel(r"$F_\infty$")
    ax.set_xlim(0.0, x_right)
    ax.set_ylim(-0.03, 1.05)
    set_legend(ax, loc=(0.6, 0.08))


# =============================================================
# Main
# =============================================================
def main():
    A0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    beta_for_ab = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0

    print("=" * 70)
    print("Figure 3: quantitative verification")
    print("=" * 70)
    print(f"A0 = {A0}")
    print(f"reference beta for panels (a),(b) = {beta_for_ab}")
    print(f"reading time series from: {TS_AFR}")
    print(f"reading summary from:     {GLOBAL_SUM}")
    print(f"saving figures to:        {FIG_ROOT}")

    runs = find_runs_for(A0)
    if not runs:
        raise FileNotFoundError(f"No runs found for A0={A0:.3f} in {TS_AFR}")

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
    panel_d(axes[1, 1], A0)

    fig.tight_layout(rect=[0.02, 0.02, 0.995, 0.995])

    png_path = FIG_ROOT / "fig3.png"
    eps_path = FIG_ROOT / "fig3.eps"

    fig.savefig(
        png_path,
        format="png",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
    )
    fig.savefig(
        eps_path,
        format="eps",
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
    )

    plt.close(fig)

    print(f"saved PNG: {png_path}")
    print(f"saved EPS: {eps_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()