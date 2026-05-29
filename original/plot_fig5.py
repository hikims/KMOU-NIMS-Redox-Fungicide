import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

# ============================================================
# Paths
# ============================================================
BASE = "/Users/mkim/Desktop/Research/fungicide_jh/fig5_rerun_A0_0_2p25_5"

SNAP_DIR = os.path.join(BASE, "data")
OUT_DIR = os.path.join(BASE, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

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
# Grid size
# ============================================================
Nx, Ny = 500, 500

# ============================================================
# Rows and columns
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
# Helper functions
# ============================================================
def find_existing_file(folder, candidates):
    """
    Return the first existing file among candidate names.
    """
    for fname in candidates:
        path = os.path.join(folder, fname)
        if os.path.exists(path):
            return path
    return None


def load_F(filepath):
    """
    Load F from field file.
    Expected columns:
    x, y, A, F, R
    """
    data = np.loadtxt(filepath)
    F = data[:, 3].reshape(Nx, Ny)
    return F


def make_physical_initial_patch():
    """
    정확한 t=0 초기조건을 물리 좌표 기준으로 생성.
    원고 설정:
    center = (5,5), radius = 2, F = 0.5 inside, F = 0 outside.
    """
    x = np.linspace(0.0, 10.0, Nx)
    y = np.linspace(0.0, 10.0, Ny)
    X, Y = np.meshgrid(x, y, indexing="ij")

    F0 = np.zeros((Nx, Ny))
    mask = (X - 5.0)**2 + (Y - 5.0)**2 <= 2.0**2
    F0[mask] = 0.5
    return F0


# ============================================================
# Figure
# ============================================================
fig, axes = plt.subplots(
    3,
    4,
    figsize=(12.2, 8.4)
)

fig.patch.set_facecolor("white")

im = None
missing_files = []

F0_exact = make_physical_initial_patch()

# ============================================================
# Plot panels
# ============================================================
for row, (row_label, A0_folder) in enumerate(A0_rows):

    folder = os.path.join(SNAP_DIR, A0_folder)

    for col, (tlabel, candidates) in enumerate(time_files):

        ax = axes[row, col]

        if col == 0:
            F = F0_exact
            source_path = "constructed exact initial patch"

        else:
            path = find_existing_file(folder, candidates)

            if path is None:
                missing_files.append(os.path.join(folder, " OR ".join(candidates)))
                F = None
            else:
                F = load_F(path)
                source_path = path

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
            ax.text(
                0.5,
                0.5,
                "missing",
                ha="center",
                va="center",
                fontsize=15,
                transform=ax.transAxes
            )

        # Column titles only on top row
        if row == 0:
            ax.set_title(
                tlabel,
                fontsize=18,
                pad=10
            )

        # Row labels only on first column
        if col == 0:
            ax.set_ylabel(
                row_label,
                fontsize=18,
                rotation=90,
                labelpad=8
            )

        # Remove axis ticks
        ax.set_xticks([])
        ax.set_yticks([])

        # Thin panel borders
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)

# ============================================================
# Layout
# ============================================================
fig.subplots_adjust(
    left=0.08,
    right=0.86,
    bottom=0.07,
    top=0.92,
    wspace=0.12,
    hspace=0.20
)

# ============================================================
# Shared colorbar
# ============================================================
cax = fig.add_axes([0.89, 0.20, 0.025, 0.60])

if im is not None:
    cbar = fig.colorbar(im, cax=cax)
else:
    norm = Normalize(vmin=0.0, vmax=0.5)
    sm = ScalarMappable(norm=norm, cmap="viridis")
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)

cbar.set_label(
    r"Fungal density $F$",
    fontsize=18,
    labelpad=12
)
cbar.ax.tick_params(labelsize=14)

# ============================================================
# Save
# ============================================================
png_path = os.path.join(OUT_DIR, "fig5.png")
pdf_path = os.path.join(OUT_DIR, "fig5.pdf")

fig.savefig(
    png_path,
    format="png",
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none"
)

fig.savefig(
    pdf_path,
    format="pdf",
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none"
)

plt.close(fig)

# ============================================================
# Output
# ============================================================
print("Saved:")
print(" -", png_path)
print(" -", pdf_path)

if missing_files:
    print("\nMissing files:")
    for f in missing_files:
        print(" -", f)