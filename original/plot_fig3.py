"""
Publication figure: A, R, F dynamics for two beta cases
  Row 1: beta = beta_caseII  (Case II : beta*theta^m < 1, F survives)
  Row 2: beta = beta_caseI   (Case I  : beta*theta^m > 1, F goes extinct)

KEY: F panels share the same y-axis range so the two cases are visually
     comparable -- the extinction case clearly shows F dropping to 0
     against the same scale as the surviving case.

Reads:
    /Users/mkim/Desktop/Research/fungicide_jh/time_series/
        A_R_F_timeseries/A0_<A0>_beta_<beta>.dat

Saves to:
    /Users/mkim/Desktop/Research/fungicide_jh/figures/ARF_two_cases_A0_<A0>.png

Usage:
    python plot_ARF_two_cases.py                  # A0=0.5, beta=8 / 15
    python plot_ARF_two_cases.py 0.5 8 15         # explicit
    python plot_ARF_two_cases.py 0.5 5 20         # wider contrast
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

# -------------------------------------------------------------
# Paths -- must match main4.cpp
# -------------------------------------------------------------
ROOT     = Path("/Users/mkim/Desktop/Research/fungicide_jh")
TS_AFR   = ROOT / "time_series" / "A_R_F_timeseries"
FIG_ROOT = ROOT / "figures"
FIG_ROOT.mkdir(exist_ok=True)

# -------------------------------------------------------------
# Model parameters from parameters.h (kept in sync manually)
# -------------------------------------------------------------
ALPHA = 50.0
ELL   = 0.3
EPS   = 0.05
THETA = 0.30
M_EXP = 2.0
BETA_CRIT = 1.0 / THETA ** M_EXP        # ~ 11.111

# -------------------------------------------------------------
# Plot style
# -------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "font.size": 14,
    "axes.titlesize": 14,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "axes.grid": True,
    "axes.axisbelow": False,
    "grid.color": "0.88",
    "grid.alpha": 1.0,
    "grid.linewidth": 0.45,
    "lines.linewidth": 1.8,
    "axes.linewidth": 0.9,
})

COLOR_A = "#1f6fbf"     # blue
COLOR_R = "#2a8a3e"     # green
COLOR_F = "#c0392b"     # red
ACCENT_CASE_II = "#2c7fb8"
ACCENT_CASE_I  = "#d7301f"


def blend_with_white(color, amount):
    """Return an opaque RGB color equivalent to color drawn with alpha=amount on white.

    EPS/PostScript does not support transparency, so use pre-blended colors
    instead of alpha for publication-safe EPS output.
    """
    rgb = np.array(mcolors.to_rgb(color))
    white = np.ones(3)
    return tuple((1.0 - amount) * white + amount * rgb)

FILL_A = blend_with_white(COLOR_A, 0.18)
FILL_R = blend_with_white(COLOR_R, 0.18)
FILL_F = blend_with_white(COLOR_F, 0.12)
MEAN_F = blend_with_white(COLOR_F, 0.55)


# =============================================================
# Loader
# =============================================================
def load_timeseries(A0, beta):
    path = TS_AFR / f"A0_{A0:.3f}_beta_{beta:.3f}.dat"
    if not path.exists():
        raise FileNotFoundError(f"missing: {path}")
    data = np.loadtxt(path)
    return {
        "t":     data[:, 0],
        "meanA": data[:, 1], "meanR": data[:, 2], "meanF": data[:, 3],
        "minA":  data[:, 4], "maxA":  data[:, 5],
        "minR":  data[:, 6], "maxR":  data[:, 7],
        "minF":  data[:, 8], "maxF":  data[:, 9],
    }


# =============================================================
# Panel plotters
# =============================================================
def plot_A_panel(ax, ts, A0):
    t = ts["t"]
    upper = A0 * np.exp(-ELL * t)
    lower = A0 * np.exp(-(ALPHA + ELL) * t)

    ax.semilogy(t, ts["meanA"] + 1e-30, color=COLOR_A, lw=2.0,
                label=r"$\bar A(t)$")
    ax.fill_between(t, np.maximum(ts["minA"], 1e-30), ts["maxA"] + 1e-30,
                    color=FILL_A, alpha=1.0, label="min~max")
    ax.semilogy(t, upper + 1e-30, color="k", ls="--", lw=1.0,
                label=r"$A_0 e^{-\ell t}$")
    ax.semilogy(t, lower + 1e-30, color="0.45", ls=":", lw=1.0,
                label=r"$A_0 e^{-(\alpha+\ell)t}$")

    ax.set_ylim(1e-12, max(2 * A0, 1.5))
    ax.set_ylabel(r"$A(x,t)$")


def plot_R_panel(ax, ts):
    t = ts["t"]
    ax.plot(t, ts["meanR"], color=COLOR_R, lw=2.0, label=r"$\bar R(t)$")
    ax.fill_between(t, ts["minR"], ts["maxR"], color=FILL_R, alpha=1.0,
                    label="min~max")
    ax.axhline(1.0, color="0.35", ls="--", lw=1.0,
               label=r"$R_\infty=1$")
    ax.axhline(THETA, color="0.45", ls=":", lw=1.0,
               label=r"$\theta$")
    ax.set_ylim(0.0, 1.08)
    ax.set_ylabel(r"$R(x,t)$")


def plot_F_panel(ax, ts, beta, ylim):
    """ylim is shared across both F panels for proper comparison."""
    t = ts["t"]
    F_inf = max(0.0, 1.0 - beta * THETA ** M_EXP)

    ax.plot(t, ts["maxF"], color=COLOR_F, lw=2.2, label=r"$\max_x F(x,t)$")
    ax.plot(t, ts["meanF"], color=MEAN_F, lw=1.4,
            label=r"$\bar F(t)$")
    ax.fill_between(t, ts["minF"], ts["maxF"], color=FILL_F, alpha=1.0)

    ax.axhline(F_inf, color="k", ls="--", lw=1.3,
               label=fr"$F_\infty={F_inf:.2f}$")

    ax.set_ylim(ylim)
    ax.set_ylabel(r"$F(x,t)$")


# =============================================================
# Main
# =============================================================
def main():
    A0          = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    beta_caseII = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
    beta_caseI  = float(sys.argv[3]) if len(sys.argv) > 3 else 15.0

    print(f"A0           = {A0}")
    print(f"Case II beta = {beta_caseII}  (beta*theta^m = {beta_caseII*THETA**M_EXP:.3f})")
    print(f"Case I  beta = {beta_caseI}   (beta*theta^m = {beta_caseI*THETA**M_EXP:.3f})")

    if beta_caseII >= BETA_CRIT:
        print(f"WARNING: beta_caseII={beta_caseII} should be < {BETA_CRIT:.3f}")
    if beta_caseI <= BETA_CRIT:
        print(f"WARNING: beta_caseI={beta_caseI} should be > {BETA_CRIT:.3f}")

    ts_II = load_timeseries(A0, beta_caseII)
    ts_I  = load_timeseries(A0, beta_caseI)

    # ----- Shared y-range for F panels -----
    F_inf_II   = max(0.0, 1.0 - beta_caseII * THETA ** M_EXP)
    F_max_seen = max(ts_II["maxF"].max(), ts_I["maxF"].max(),
                     F_inf_II)
    F_ylim     = (-0.02 * F_max_seen, 1.10 * F_max_seen)

    # ----- Build figure -----
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.2), sharex=True)

    # Row 1: Case II
    plot_A_panel(axes[0, 0], ts_II, A0)
    plot_R_panel(axes[0, 1], ts_II)
    plot_F_panel(axes[0, 2], ts_II, beta_caseII, F_ylim)

    # Row 2: Case I
    plot_A_panel(axes[1, 0], ts_I, A0)
    plot_R_panel(axes[1, 1], ts_I)
    plot_F_panel(axes[1, 2], ts_I, beta_caseI, F_ylim)

    # x-axis labels only on bottom row
    for ax in axes[1, :]:
        ax.set_xlabel(r"time")

    # Per-panel legends (moved to lower left)
    for ax in axes.flat:
        ax.legend(loc="lower right", framealpha=1.0, facecolor="white", edgecolor="0.75", fontsize=11)

    # Column headers
    col_titles = [
        r"(a)  Fungicide $A$",
        r"(b)  Redox state $R$",
        r"(c)  Fungal density $F$",
    ]
    for j, ttl in enumerate(col_titles):
        axes[0, j].set_title(ttl, fontsize=13, fontweight="bold", pad=8)

    # Row labels
    row_labels = [
        (rf"$\beta={beta_caseII:g}$,  "
         rf"$\beta\theta^m={beta_caseII*THETA**M_EXP:.2f}$",
         ACCENT_CASE_II),
        (rf"$\beta={beta_caseI:g}$,  "
         rf"$\beta\theta^m={beta_caseI*THETA**M_EXP:.2f}$",
         ACCENT_CASE_I),
    ]
    for i, (lbl, c) in enumerate(row_labels):
        axes[i, 0].annotate(
            lbl, xy=(-0.30, 0.5), xycoords="axes fraction",
            ha="center", va="center", rotation=90,
            fontsize=12, fontweight="bold", color=c,
        )

    fig.tight_layout(rect=[0.04, 0, 1, 0.99])

    out_eps = FIG_ROOT / "fig2.eps"
    out_png = FIG_ROOT / "fig2.png"

    fig.savefig(out_eps, bbox_inches="tight", format="eps")
    fig.savefig(out_png, bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"\nsaved: {out_eps}")

    print(f"saved: {out_png}")


if __name__ == "__main__":
    main()