# Execution scripts

Place all files in this directory under the repository `scripts/` folder.
The scripts are GitHub-friendly wrappers around the C++ drivers in `drivers/` and the shared solver files in `src/`.
They use relative repository paths and write generated outputs under `results/` by default.

## Location

```text
redox_fungicide_github_package/
├── drivers/
├── src/
├── configs/
├── scripts/        <- put these .sh files here
└── results/        <- generated outputs; ignored by git
```

## Scripts

| Script | Driver target | Paper item |
|---|---|---|
| `convergence_exec.sh` | `main_table_S1_S2_convergence` | Supplementary Tables S1--S2 and invariant diagnostics |
| `dose_time_maps_exec.sh` | `main_fig2_table3_4_dose_time_maps` | Main Table 4 threshold comparison; historical filename kept |
| `sensitivity_exec.sh` | `main_supp_figS2_S3_sensitivity_curves` | Supplementary Figures S2--S3 |
| `spatial_dynamics_exec.sh` | `main_fig5_spatial_dynamics` | Main Figure 5 |
| `threshold_exec.sh` | `main_supp_figS1_threshold` | Supplementary Figure S1 |
| `time_exec.sh` | `main_fig3_fig4_verification` | Main Figures 3--4 |

## Usage examples

```bash
chmod +x scripts/*.sh
./scripts/time_exec.sh
./scripts/spatial_dynamics_exec.sh
./scripts/threshold_exec.sh
./scripts/sensitivity_exec.sh
./scripts/dose_time_maps_exec.sh
./scripts/convergence_exec.sh 0.5
```

Optional environment variables:

```bash
RESULTS_ROOT=/path/to/results MACHINE_LABEL=macstudio ./scripts/dose_time_maps_exec.sh
A0=0.5 BETA_LIST="5 8 10 11 12 15 20" ./scripts/time_exec.sh
```

The scripts assume that `Makefile`, `drivers/`, `src/`, and `configs/` are already present at the repository root.
