# Execution scripts

GitHub-friendly wrappers around the C++ drivers in `drivers/` and the shared
solver files in `src/`. They use relative repository paths and write generated
outputs under `results/` by default.

## Location

```text
<repository root>/
├── drivers/
├── src/
├── configs/
├── scripts/        <- these .sh files
├── bin/            <- compiled executables; ignored by git
└── results/        <- generated outputs; ignored by git
```

## Scripts

| Script | Driver target | Paper item |
|---|---|---|
| `convergence_exec.sh` | `main_table_S1_S2_convergence` | Supplementary Tables S1–S2 and invariant diagnostics |
| `dose_time_maps_exec.sh` | `main_fig2_table3_4_dose_time_maps` | Main Figure 2 and Table 2; historical filename kept |
| `front_speed_exec.sh` | `main_fig6_front_speed` | Main Figure 6 and Table 3 |
| `sensitivity_exec.sh` | `main_supp_figS2_S3_sensitivity_curves` | Legacy sensitivity curves; see note below |
| `spatial_dynamics_exec.sh` | `main_fig5_spatial_dynamics` | Main Figure 5 |
| `threshold_exec.sh` | `main_supp_figS1_threshold` | Supplementary Figures S1–S3 |
| `time_exec.sh` | `main_fig3_fig4_verification` | Main Figures 3–4 |

Note on the parameter-dependence scan. `threshold_exec.sh` writes both
`results/threshold_parameter/summary/threshold_estimates.txt` and
`results/threshold_parameter/summary/threshold_scan_results.txt`. The scan file
is the one consumed by `plot_figS1.py`, `plot_figS2_S3.py`, and therefore by
Supplementary Figures S1, S2 and S3. `sensitivity_exec.sh` writes
`results/sensitivity/summary/sensitivity_results.txt`, which no current plotting
script reads; it is retained for the earlier version of the sensitivity study.

## Usage examples

```bash
chmod +x scripts/*.sh
./scripts/time_exec.sh
./scripts/spatial_dynamics_exec.sh
./scripts/front_speed_exec.sh
./scripts/threshold_exec.sh
./scripts/dose_time_maps_exec.sh
./scripts/convergence_exec.sh 0.5
```

Optional environment variables:

```bash
RESULTS_ROOT=/path/to/results MACHINE_LABEL=macstudio ./scripts/dose_time_maps_exec.sh
A0=0.5 BETA_LIST="5 8 10 11 12 15 20" ./scripts/time_exec.sh
TBB_PREFIX=/opt/homebrew/opt/tbb ./scripts/front_speed_exec.sh
```

`front_speed_exec.sh` takes its run parameters as positional arguments rather
than environment variables:

```bash
./scripts/front_speed_exec.sh              # beta = 10, A0 = 0 2.25 5
./scripts/front_speed_exec.sh 8            # secondary check at beta = 8
./scripts/front_speed_exec.sh 10 0 2.25 5  # explicit
```

It delegates the build to the Makefile, so no separate `make` step is needed.
The driver is called once with every dose rather than once per dose: it opens
its summary file in mode `"w"`, so a per-dose loop would leave only the last
dose in `summary_fig6.txt`.

The scripts assume that `Makefile`, `drivers/`, `src/`, and `configs/` are
already present at the repository root.
