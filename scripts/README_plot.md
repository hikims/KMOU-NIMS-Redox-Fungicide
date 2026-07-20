# Plotting scripts

Python scripts that reproduce every figure in the paper and the supplement
from the simulation output files in `results/`.

| Script | Paper figure | What it shows |
|---|---|---|
| `plot_fig2.py` | Figure 2 | Dose–time map of finite-time suppression (β = 8) |
| `plot_fig3.py` | Figure 3 | Two-row A/R/F dynamics: persistence (β = 8) vs. extinction (β = 15) |
| `plot_fig4.py` | Figure 4 | Four-panel verification of the recovered-state results (β sweep) |
| `plot_fig5.py` | Figure 5 | 3 × 4 spatial snapshots of F(x,t) for three doses (β = 10) |
| `plot_fig6.py` | Figure 6 | Recovered-state front propagation and the Fisher–KPP speed (β = 10) |
| `plot_figS1.py` | Figure S1 | Infected-area critical dose vs. multiplicative factor on dA, dF, α |
| `plot_figS2_S3.py` | Figures S2, S3 | Dose–response curves for varying α (S2) and dF (S3) |

All figures are written to `results/figures/` as `.pdf`, `.eps` and `.png`.
Fonts are embedded as TrueType (`pdf.fonttype = 42`, `ps.fonttype = 42`), as
required for submission. `plot_fig6.py` additionally writes
`fig6_front_speed.csv` with the fitted speeds and deviations.

---

## Requirements

- Python ≥ 3.8
- `numpy`, `pandas`, `matplotlib`

```bash
pip install numpy pandas matplotlib
```

---

## Quick start

From the repository root, after the drivers have produced output in
`results/`:

```bash
python scripts/plot_fig2.py
python scripts/plot_fig3.py
python scripts/plot_fig4.py
python scripts/plot_fig5.py
python scripts/plot_fig6.py
python scripts/plot_figS1.py
python scripts/plot_figS2_S3.py
```

Each script locates the repository root from its own file path, reads from
`results/<subfolder>/`, and writes to `results/figures/`. No absolute paths
appear anywhere. If the results tree lives outside the repository, point the
scripts at it with one environment variable:

```bash
export REDOX_FUNGICIDE_RESULTS=/scratch/redox_fungicide/results
```

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
├── fig5_rerun_A0_0_2p25_5/
│   └── data/
│       ├── A0_0.000/
│       │   └── field_0000.dat, field_0100.dat,  # → fig5
│       │       field_0500.dat, field_1000.dat
│       ├── A0_2.250/
│       └── A0_5.000/
├── front_speed/
│   ├── front_beta_10.000_A0_0.000.dat           # → fig6
│   ├── front_beta_10.000_A0_2.250.dat
│   ├── front_beta_10.000_A0_5.000.dat
│   └── summary_fig6.txt
├── threshold_parameter/
│   └── summary/
│       └── threshold_scan_results.txt           # → figS1, figS2, figS3
└── figures/                                     # (output)
```

The directory name `theorem_1_5_II` predates the current theorem numbering; it
holds the data for Theorem 3.4 (sharp recovered-state transition).
`plot_fig4.py` also accepts `time_series/recovered_state_transition/`.

Every script accepts an explicit path where applicable, e.g.
`python scripts/plot_figS1.py --data /custom/threshold_scan_results.txt --outdir /custom/out`.

---

## Script-by-script details

### `plot_fig2.py` — Figure 2

Dose–time heatmaps of the normalized aggregate burden $\overline{M}(t)$ and
the infected-area percentage $\mathcal{I}(t)$ in the persistence regime
(β = 8, ℬ = 0.72).

**Input.** `results/dose_time_map/B1_A0_Tobs_map.txt`, columns
`A0  time  total_A  total_F  M_inf  area_percent  min_A  max_A  min_R  max_R  min_F  max_F  log10_total_F  log10_M_inf`.

**Output.** `fig2.{pdf,eps,png}`, `fig2_Mbar_only.{pdf,eps,png}`,
`fig2_summary_table.csv`.

```bash
python scripts/plot_fig2.py --data PATH --outdir DIR --contour 2e-2
```

---

### `plot_fig3.py` — Figure 3

A/R/F time series in two rows: persistence (βθ^m < 1, F → 1 − βθ^m) against
extinction (βθ^m > 1, F → 0). The F panels share a y-axis.

**Input.** `results/time_series/A_R_F_timeseries/A0_<A0>_beta_<beta>.dat`,
columns `t  meanA meanR meanF  minA maxA  minR maxR  minF maxF …`.

**Output.** `fig3.{pdf,eps,png}`.

```bash
python scripts/plot_fig3.py                # A0=0.5, beta=8 / 15
python scripts/plot_fig3.py 0.5 8 15
python scripts/plot_fig3.py 0.5 5 20
```

---

### `plot_fig4.py` — Figure 4

Four-panel verification: (a) fungicide decay against the bounds of
Proposition 3.3, (b) redox recovery error, (c) fungal-density evolution over a
β sweep, (d) bifurcation diagram against max(0, 1 − βθ^m) from Theorem 3.4.

**Input.** the time-series files above, plus
`results/time_series/theorem_1_5_II/global_summary.txt`, columns
`A0  beta  case_id  beta*theta^m  F_asymptote  final_mean_F  final_max_F  |final_mean_F − F_asymptote|`.

**Output.** `fig4.{pdf,eps,png}`.

```bash
python scripts/plot_fig4.py                # A0=0.5, beta_for_(a,b)=10
python scripts/plot_fig4.py 0.5 10
```

---

### `plot_fig5.py` — Figure 5

A 3 × 4 grid of colour maps of F(x,t): rows A0 = 0, 2.25, 5; columns
t = 0, 10, 50, 100. **All four columns, including t = 0, are read from
simulation snapshots**; nothing is reconstructed analytically. Physical
coordinates use cell centres with h = 0.02, and a thin white contour marks the
detection level F = F_thr = 10⁻³ used by the infected-area criterion.

**Input.** `results/fig5_rerun_A0_0_2p25_5/data/A0_<A0>/field_<tag>.dat`,
columns `x  y  A  F  R`. Accepted names per time: `field_0000.dat` /
`field_t000.dat` for t = 0, `field_0100.dat` / `field_t010.dat` for t = 10,
and so on. `results/spatial_snapshots/` is accepted as a fallback layout.

Grid size is fixed at 500 × 500; edit `NX, NY` for a different grid.

**Output.** `fig5.{pdf,eps,png}`.

```bash
python scripts/plot_fig5.py
python scripts/plot_fig5.py --snap-dir results/spatial_snapshots
python scripts/plot_fig5.py --no-contour
```

---

### `plot_fig6.py` — Figure 6

Recovered-state front propagation, verifying the spreading speed of
Theorem 3.5,

    c_inf = 2 sqrt( dF (1 − ℬ) / (1 + θ^m) ),   ℬ = βθ^m.

Panels: (a) equivalent front radius r(t) with the prediction r(0) + c_inf·t and
the Bramson-corrected fit; (b) instantaneous speed against c_inf with the fit
window shaded; (c) front radius for three doses, showing that a finite dose
translates the trajectory in time without changing its asymptotic slope.

Panels (a) and (b) use the A0 = 0 run with R_in ≡ 1, so A ≡ 0 and R ≡ 1 hold
exactly and the F equation is the recovered Fisher–KPP equation from t = 0.
The measured speed therefore contains no redox-recovery transient.

**Input.** `results/front_speed/front_beta_<beta>_A0_<A0>.dat`, a scalar
level-set diagnostic with columns `t  n_lev[0]  n_lev[1]  n_lev[2]`, where
`n_lev[k]` counts grid cells with F ≥ `LEVELS[k]`. Dense field snapshots under
`results/spatial_snapshots/` are accepted as a fallback, but the scalar file is
a few kB against roughly 1.6 GB of fields.

**Output.** `fig6.{pdf,eps,png}` and `fig6_front_speed.csv`.

```bash
python scripts/plot_fig6.py                # beta = 10
python scripts/plot_fig6.py --beta 8       # secondary check
python scripts/plot_fig6.py --scales       # print the Table 3 scales and exit
```

Run `./scripts/front_speed_exec.sh` first; the Figure 5 run is too short
(t_final = 100) and too sparse (four snapshots) to fit a front speed.

---

### `plot_figS1.py` — Figure S1

Estimated infected-area critical dose $\widehat{A}^{\mathcal I}_{0,c}(T_{\rm obs}; F_{\rm thr}, \mathcal I_*)$
as a function of multiplicative factors applied to dA, dF and α. The threshold
is **recomputed from the dose scan** rather than read from a precomputed
estimates file, so the criterion (`--criterion`, `--tobs`) can be varied
without rerunning the simulation. Markers show the midpoint of the first
non-suppressive–suppressive bracket, with vertical error bars for the bracket.
A downward triangle marks a setting already suppressed at the lowest tested
dose; a cross marks one still not suppressed at the largest tested dose.

**Input.** `results/threshold_parameter/summary/threshold_scan_results.txt`,
columns
`tag  A0  dA_factor  dF_factor  alpha_factor  dA  dF  alpha  total_F_t20  area_t20  total_F_t50  area_t50  total_F_t100  area_t100`.

**Output.** `figS1.{pdf,eps,png}`.

```bash
python scripts/plot_figS1.py
python scripts/plot_figS1.py --criterion 1e-8 --tobs 100
```

---

### `plot_figS2_S3.py` — Figures S2, S3

Dose–response curves of $\overline{M}(T_{\rm obs})$ at T_obs = 100:
Figure S2 varies α, Figure S3 varies dF. Both read the same scan file as
Figure S1.

**Configuration — `OMEGA_AREA`.** Set it to `1.0` if `total_F_t100` already
stores the normalized average $\overline{M}$, or to `100.0` if it stores the
raw integral (|Ω| = 100 for Ω = [0,10]²). Sanity check: at the lowest dose the
curves should land around 10⁻² to 10⁻¹, near the persistence plateau
F∞ = 0.1. Values around 10⁰–10¹ mean `OMEGA_AREA` should be `100.0`.

**Output.** `figS2.{pdf,eps,png}`, `figS3.{pdf,eps,png}`.

```bash
python scripts/plot_figS2_S3.py
```

---

## Notes

- The plotting scripts duplicate a few baseline values (`THETA = 0.30`,
  `M_EXP = 2.0`, `ALPHA = 50.0`, `ELL = 0.3`, `EPS = 0.05`, `D_F = 1e-3`) for
  annotation. If these change in `src/parameters.h`, update the constants at
  the top of `plot_fig3.py`, `plot_fig4.py` and `plot_fig6.py`.
- EPS output carries no transparency (PostScript does not support it), so light
  fills are pre-blended with white. The PNG output matches what alpha blending
  would give.
- The scripts do not modify any input file. Re-running them is safe and
  idempotent.
