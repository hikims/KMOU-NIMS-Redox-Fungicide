# Plotting scripts

Python scripts that reproduce every figure in the paper and the supplement
from the simulation output files in `results/`.

| Script | Paper figure | What it shows |
|---|---|---|
| `plot_fig2.py` | Figure 2 | Dose–time map of finite-time practical suppression (β = 8) |
| `plot_fig3.py` | Figure 3 | Two-row A/R/F dynamics: persistence (β = 8) vs. extinction (β = 15) |
| `plot_fig4.py` | Figure 4 | Four-panel verification of the asymptotic results (β sweep) |
| `plot_fig5.py` | Figure 5 | 3 × 4 spatial snapshots of F(x,t) for three doses (β = 10) |
| `plot_figS1.py` | Figure S1 | Estimated empirical threshold vs. multiplicative factor on dA, dF, α |
| `plot_figS2_S3.py` | Figures S2, S3 | Dose–response curves for varying α (S2) and dF (S3) |

All figures are written to `results/figures/`. File names match the paper:
`fig2.{pdf,eps,png}`, `fig3.{eps,png}`, `fig4.{eps,png}`, `fig5.{pdf,png}`,
`figS1.{eps,png}`, `figS2.{eps,png}`, `figS3.{eps,png}`.

---

## Requirements

- Python ≥ 3.8
- `numpy`, `pandas`, `matplotlib`

Install with:

```bash
pip install numpy pandas matplotlib
```

---

## Quick start

From the repository root, after the simulation has produced output files
in `results/`:

```bash
python scripts/plot_fig2.py
python scripts/plot_fig3.py
python scripts/plot_fig4.py
python scripts/plot_fig5.py
python scripts/plot_figS1.py
python scripts/plot_figS2_S3.py
```

Each script auto-locates the repository root from its own file path,
reads from `results/<subfolder>/`, and writes PNG/EPS/PDF files to
`results/figures/`. No environment variables or absolute paths required.

---

## Expected data layout under `results/`

```
results/
├── dose_time_map/
│   └── B1_A0_Tobs_map.txt                       # → fig2
├── time_series/
│   ├── A_R_F_timeseries/
│   │   └── A0_<A0>_beta_<beta>.dat              # → fig3, fig4(a–c)
│   └── theorem_1_5_II/
│       └── global_summary.txt                   # → fig4(d)
├── spatial_snapshots/
│   ├── A0_0.000/
│   │   └── field_t000.dat, field_t010.dat, ...  # → fig5
│   ├── A0_2.250/
│   └── A0_5.000/
├── threshold_parameter/
│   └── summary/
│       ├── threshold_estimates.txt              # → figS1
│       └── threshold_scan_results.txt           # → figS2, figS3
└── figures/                                     # (output)
```

If your output files live elsewhere, every script accepts a custom path
via command-line options where applicable (e.g. `python scripts/plot_fig2.py --data /custom/path/B1_A0_Tobs_map.txt --outdir /custom/out`).

---

## Script-by-script details

### `plot_fig2.py` — Figure 2

**Purpose.** Dose–time heatmaps of the normalized aggregate burden
$\overline{M}(t)$ and the infected-area fraction $\mathcal{I}(t)$, in the
asymptotic persistence regime (β = 8).

**Input.** `results/dose_time_map/B1_A0_Tobs_map.txt`, with columns
`A0  time  total_A  total_F  M_inf  area_percent  min_A  max_A  min_R  max_R  min_F  max_F  log10_total_F  log10_M_inf`.

**Output.**
- `fig2.{pdf,eps,png}` — two-panel figure
- `fig2_Mbar_only.{pdf,eps,png}` — single-panel variant
- `fig2_summary_table.csv` — representative values for A0 = 0, 1, …, 5

**Options.**

```bash
python scripts/plot_fig2.py --data PATH --outdir DIR --contour 2e-2
```

---

### `plot_fig3.py` — Figure 3

**Purpose.** A/R/F time-series in two rows, comparing the persistence case
(β·θ^m < 1, F → 1 − β·θ^m) with the extinction case (β·θ^m > 1, F → 0).
The F panels share a y-axis so the two cases are directly comparable.

**Input.**
`results/time_series/A_R_F_timeseries/A0_<A0>_beta_<beta>.dat`,
with columns
`t  meanA meanR meanF  minA maxA  minR maxR  minF maxF  …`.

**Output.** `fig3.png`, `fig3.eps`.

**Usage.**

```bash
python scripts/plot_fig3.py                # defaults: A0=0.5, beta=8 / 15
python scripts/plot_fig3.py 0.5 8 15       # explicit
python scripts/plot_fig3.py 0.5 5 20       # wider contrast
```

---

### `plot_fig4.py` — Figure 4

**Purpose.** Four-panel quantitative verification of the asymptotic
results:
- (a) fungicide decay against $A_0 e^{-\ell t}$ and $A_0 e^{-(\alpha + \ell)t}$
- (b) redox recovery error $\|1 - R(\cdot, t)\|_\infty$
- (c) fungal-density evolution for several β values
- (d) bifurcation diagram comparing numerical $M_\infty$ with $\max(0, 1 - \beta\theta^m)$

**Input.**
- `results/time_series/A_R_F_timeseries/A0_<A0>_beta_<beta>.dat`
- `results/time_series/theorem_1_5_II/global_summary.txt`,
  with columns
  `A0  beta  case_id  beta*theta^m  F_asymptote  final_mean_F  final_max_F  |final_mean_F − F_asymptote|`

**Output.** `fig4.png`, `fig4.eps`.

**Usage.**

```bash
python scripts/plot_fig4.py                # defaults: A0=0.5, beta_for_(a,b)=10
python scripts/plot_fig4.py 0.5
python scripts/plot_fig4.py 0.5 10
```

---

### `plot_fig5.py` — Figure 5

**Purpose.** A 3 × 4 grid of color maps of $F(x, t)$:
- rows: A0 = 0, 2.25, 5
- columns: t = 0, 10, 50, 100

The t = 0 column is reconstructed analytically from the initial condition
in Section 5.1 of the paper (circular patch, center (5, 5), radius 2,
F = 0.5 inside, F = 0 outside).

**Input.**
`results/spatial_snapshots/A0_<A0>/field_<tag>.dat`, with columns
`x  y  A  F  R`. Accepted name patterns:

| Time | Accepted names |
|---|---|
| t = 0 | reconstructed analytically (no file needed) |
| t = 10 | `field_t010.dat` or `field_0100.dat` |
| t = 50 | `field_t050.dat` or `field_0500.dat` |
| t = 100 | `field_t100.dat` or `field_1000.dat` |

Grid size is hard-coded to 500 × 500. Edit `Nx, Ny` at the top of the
script if you use a different grid.

**Output.** `fig5.png`, `fig5.pdf`.

**Usage.**

```bash
python scripts/plot_fig5.py
```

---

### `plot_figS1.py` — Figure S1

**Purpose.** Estimated empirical fixed-time suppression threshold
$A^*_{0,\mathrm{emp}}(T_\mathrm{obs}; \mathcal{I}_*)$ as a function of
multiplicative factors applied to dA, dF, and α. Filled markers indicate
thresholds found within the tested dose range; crosses indicate parameter
settings for which the criterion was not met up to A0,max = 3.

**Input.**
`results/threshold_parameter/summary/threshold_estimates.txt`, with columns
`group  factor  dA_factor  dF_factor  alpha_factor  A0_lower  A0_upper  threshold_mid`.

**Output.** `figS1.png`, `figS1.eps`.

**Usage.**

```bash
python scripts/plot_figS1.py
```

---

### `plot_figS2_S3.py` — Figures S2, S3

**Purpose.** Dose–response curves of $\overline{M}(T_\mathrm{obs})$ at
$T_\mathrm{obs} = 100$:
- Figure S2: five values of the fungicide consumption rate α
- Figure S3: five values of the fungal diffusivity dF

**Input.**
`results/threshold_parameter/summary/threshold_scan_results.txt`, with
columns
`tag  A0  dA_factor  dF_factor  alpha_factor  dA  dF  alpha  total_F_t20  area_t20  total_F_t50  area_t50  total_F_t100  area_t100`.

**Output.** `figS2.{png,eps}`, `figS3.{png,eps}`.

**Important configuration — `OMEGA_AREA`.**

Near the top of the script:

```python
OMEGA_AREA = 1.0   # see comment below
```

- `OMEGA_AREA = 1.0` if column `total_F_t100` already stores the
  normalized average $\overline{M}$.
- `OMEGA_AREA = 100.0` if it stores the raw integral $M_\mathrm{tot} = \int_\Omega F\,dx$
  (since $|\Omega| = 100$ for $\Omega = [0, 10]^2$).

**Quick sanity check.** At the lowest dose (A0 = 1) the curves should
land in the range $\sim 10^{-2}$ to $10^{-1}$ (close to the persistence
plateau $F_\infty = 0.1$). If they instead come out around $\sim 10^0$ to
$10^1$, change `OMEGA_AREA` to `100.0`.

**Usage.**

```bash
python scripts/plot_figS2_S3.py
```

---

## Notes

- All scripts use the same model parameters as the simulation
  (`THETA = 0.30`, `M_EXP = 2.0`, `ALPHA = 50.0`, `ELL = 0.3`,
  `EPS = 0.05`). If you change these in `configs/parameters.h`, update
  the matching constants at the top of `plot_fig3.py` and `plot_fig4.py`.
- EPS output is rendered without transparency (PostScript does not support
  it), so light fills are pre-blended with white. The PNG output is
  pixel-identical to what you would obtain with alpha blending.
- The scripts do not modify any input files. Re-running them is safe and
  idempotent.
