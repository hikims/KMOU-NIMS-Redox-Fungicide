# Redox-coupled fungicide–fungal system: reproducibility drivers

This repository contains the C++ driver programs and Python plotting scripts used to generate the numerical data and figures for the manuscript

> *Sharp recovered-state transition and fixed-time critical doses in a redox-coupled fungicide–fungal system*

and its supplementary material.

## Repository layout

```text
drivers/    C++ driver programs, one per manuscript figure or table
src/        Shared solver files: functions.cpp, initial.cpp, headers, parameters.h
configs/    Input case lists for parameter scans
scripts/    Shell drivers (*_exec.sh) and Python plotting scripts (plot_*.py)
docs/       Reproduction manifest and notes
bin/        Compiled executables, created by make, ignored by git
results/    Default output root, created automatically, ignored by git
```

`bin/` and `results/` are not included in the repository; both are created on
first build and first run.

## Reproduction map

Each row is a self-contained pipeline: run the shell driver, then the plotting
script.

| Paper item | Driver | Shell driver | Plotting script |
| --- | --- | --- | --- |
| Main Fig. 2, Table 2 | `main_fig2_table3_4_dose_time_maps.cpp` | `dose_time_maps_exec.sh` | `plot_fig2.py` |
| Main Figs. 3–4 | `main_fig3_fig4_verification.cpp` | `time_exec.sh` | `plot_fig3.py`, `plot_fig4.py` |
| Main Fig. 5 | `main_fig5_spatial_dynamics.cpp` | `spatial_dynamics_exec.sh` | `plot_fig5.py` |
| Main Fig. 6, Table 3 | `main_fig6_front_speed.cpp` | `front_speed_exec.sh` | `plot_fig6.py` |
| Suppl. Figs. S1–S3 | `main_supp_figS1_threshold.cpp` | `threshold_exec.sh` | `plot_figS1.py`, `plot_figS2_S3.py` |
| Suppl. Tables S1–S2 | `main_table_S1_S2_convergence.cpp` | `convergence_exec.sh` | — |

Main Figure 1 is a schematic and is not produced by this repository.

### Data locations read by the plotting scripts

| Plotting script | Input under `results/` |
| --- | --- |
| `plot_fig3.py`, `plot_fig4.py` | `time_series/A_R_F_timeseries/A0_<A0>_beta_<beta>.dat` |
| `plot_fig4.py` (panel d) | `time_series/theorem_1_5_II/global_summary.txt` |
| `plot_fig5.py` | `fig5_rerun_A0_0_2p25_5/data/A0_*/field_*.dat` |
| `plot_fig6.py` | `front_speed/front_beta_<beta>_A0_<A0>.dat` |
| `plot_figS1.py`, `plot_figS2_S3.py` | `threshold_parameter/summary/threshold_scan_results.txt` |

`threshold_exec.sh` writes both `threshold_estimates.txt` and
`threshold_scan_results.txt` into `threshold_parameter/summary/`. The scan file
is the one the plotting scripts read: `plot_figS1.py` recomputes the critical
dose from the scan rather than reading the precomputed estimates, so the
infected-area criterion can be varied without rerunning the simulation. The
legacy driver `main_supp_figS2_S3_sensitivity_curves` writes
`sensitivity/summary/sensitivity_results.txt`, which no current plotting script
reads; `sensitivity_exec.sh` is retained for the earlier sensitivity study.

The directory name `theorem_1_5_II` predates the current theorem numbering; it
holds the data for Theorem 3.4 (sharp recovered-state transition) of the
submitted manuscript. `plot_fig4.py` also accepts the clearer name
`time_series/recovered_state_transition/`.

All plotting scripts write to `results/figures/`, producing `.pdf`, `.eps` and
`.png` for every panel figure. Fonts are embedded as TrueType
(`pdf.fonttype = 42`, `ps.fonttype = 42`) as required for submission.

## Path conventions

Every script resolves its paths relative to the repository, so
`python scripts/<name>.py` behaves the same from any working directory and on
any machine. No absolute or user-specific paths appear anywhere in the
repository.

* C++ drivers take the output root at compile time via `OUTPUT_ROOT`
  (default `./results/`, trailing slash required) and are launched from the
  repository root by the shell drivers.
* Python plotting scripts default to `<repo>/results`. Set `RESULTS_ROOT`
  to read a results tree stored elsewhere. This is the same variable the
  shell drivers accept, so one export covers both halves of the pipeline:

  ```bash
  export RESULTS_ROOT=/scratch/redox_fungicide/results
  python scripts/plot_fig6.py
  ```

## Dependencies

* C++17 compiler
* Intel oneTBB / TBB
* Python 3 with NumPy, pandas and Matplotlib

### Installing TBB

```bash
brew install tbb          # macOS (Homebrew)
sudo apt install libtbb-dev   # Ubuntu/Debian
sudo yum install tbb-devel    # CentOS/RHEL
```

If TBB is installed in a nonstandard location, pass its prefix:

```bash
make all TBB_PREFIX=/opt/homebrew/opt/tbb    # Apple silicon
make all TBB_PREFIX=/usr/local/opt/tbb       # Intel macOS
```

## Build

```bash
make all OUTPUT_ROOT=./results/ MACHINE_LABEL=local
```

or build a single target:

```bash
make main_fig6_front_speed OUTPUT_ROOT=./results/
```

Executables are written to `bin/`. `make clean` removes them. The shell
drivers in `scripts/` build their own target through the same Makefile, so
there is no separate build step before running them.

A failure such as

```text
fatal error: 'tbb/parallel_for.h' file not found
```

means TBB is not installed or not on the include path; install it or pass
`TBB_PREFIX` as above.

## Example runs

```bash
./scripts/dose_time_maps_exec.sh       # Main Fig. 2 and Table 2 data
./scripts/time_exec.sh                 # Main Figs. 3-4 data
./scripts/spatial_dynamics_exec.sh     # Main Fig. 5 snapshots
./scripts/front_speed_exec.sh          # Main Fig. 6 front-speed diagnostic
./scripts/threshold_exec.sh            # Suppl. Fig. S1 threshold scan
./scripts/sensitivity_exec.sh          # Suppl. Figs. S2-S3 sensitivity curves
./scripts/convergence_exec.sh          # Suppl. Tables S1-S2 refinement checks
```

Then produce the figures:

```bash
python scripts/plot_fig2.py
python scripts/plot_fig3.py
python scripts/plot_fig4.py
python scripts/plot_fig5.py
python scripts/plot_fig6.py
python scripts/plot_figS1.py
python scripts/plot_figS2_S3.py
```

### Notes on the Figure 6 run

The front-speed diagnostic records a scalar level-set count rather than field
snapshots, so its entire output is a few kilobytes. It runs to `t_final = 160`
with `dt_out = 1.0`, both longer and denser than the Figure 5 run, because a
front speed cannot be fitted from four snapshots ending at `t = 100`.

Panels (a) and (b) use `A0 = 0` with `R_in = 1`, so that `A = 0` and `R = 1`
hold exactly and the fungal equation is the recovered Fisher–KPP equation from
`t = 0`. This isolates the spreading-speed test from the redox transient. The
positive doses feed only the front-delay panel (c).

Unlike the other drivers, `main_fig6_front_speed` is called once with every
dose rather than once per dose, because it opens its summary file in mode
`"w"`; the shell driver handles this.

## Shared solver files

```text
src/functions.hpp
src/functions.cpp
src/initial.hpp
src/initial.cpp
src/parameters.h
```

Drivers include these as:

```cpp
#include "functions.hpp"
#include "initial.hpp"
#include "parameters.h"
```

Do not edit these unless you intend to change the numerical model, grid
resolution, time step, baseline parameters, or initial conditions. The
plotting scripts duplicate a few baseline values (`alpha`, `ell`, `epsilon`,
`theta`, `m`, `d_F`) for annotation and must be kept in sync with
`src/parameters.h`.

## Output notes

The drivers can create large raw field files. Only code, configuration files,
plotting scripts, and the processed data needed to reproduce the manuscript
figures and tables are kept in the repository. Large raw simulation outputs
are regenerated from the drivers or are available from the corresponding
author on reasonable request.

## Citation and archival

See `CITATION.cff`. The state of the repository at manuscript submission is
tagged; use that tag, or the archived release DOI, rather than the moving tip
of the default branch when reproducing the published results.
