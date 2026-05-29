# Redox-mediated fungicide suppression: reproducibility drivers

This repository contains the C++ driver programs and Python plotting scripts used to generate the numerical data for the manuscript figures, tables, and supplementary figures. The filenames of the uploaded driver programs are preserved.

## Repository layout

```text
drivers/    C++ driver programs for manuscript figures and tables
configs/    Input case lists for parameter scans
src/        Shared solver files: functions.cpp, initial.cpp, headers, parameters.h
scripts/    Convenience shell scripts and Python plotting scripts
docs/       Reproduction manifest and notes
results/    Default output root, created automatically and ignored by git
```

The `results/` directory is not included in the repository. It is created automatically when the drivers are executed.

## Driver-to-paper map

| Paper item                   | Driver filename                             | Main output                                                        |
| ---------------------------- | ------------------------------------------- | ------------------------------------------------------------------ |
| Main Table 4                 | `main_fig2_table3_4_dose_time_maps.cpp`     | `results/TierB/summary/B2_threshold_comparison_<machine>.txt`      |
| Main Figures 3--4            | `main_fig3_fig4_verification.cpp`           | `results/time_series/`                                             |
| Main Figure 5                | `main_fig5_spatial_dynamics.cpp`            | `results/fig5_rerun_A0_0_2p25_5/data/A0_*/field_t*.dat`            |
| Supplementary Figure S1      | `main_supp_figS1_threshold.cpp`             | `results/threshold_parameter/summary/threshold_estimates.txt`      |
| Supplementary Figures S2--S3 | `main_supp_figS2_S3_sensitivity_curves.cpp` | `results/sensitivity/summary/sensitivity_results.txt`              |
| Supplementary Tables S1--S2  | `main_table_S1_S2_convergence.cpp`          | `results/summary/summary_A0_vs_area.txt` and invariant diagnostics |

Note: the historical filename `main_fig2_table3_4_dose_time_maps.cpp` is preserved, but the current driver is primarily the B2 finite-time threshold-comparison driver used for the main threshold table.

## Shared solver files

The following shared C++ solver files are stored in `src/` and are used by the drivers in `drivers/`:

```text
src/functions.hpp
src/functions.cpp
src/initial.hpp
src/initial.cpp
src/parameters.h
```

The driver files include these headers as:

```cpp
#include "functions.hpp"
#include "initial.hpp"
#include "parameters.h"
```

Do not edit these shared files unless you intend to change the numerical model, grid resolution, time step, baseline parameters, or initial conditions.

## Dependencies

The C++ drivers require:

* C++17 compiler
* Intel oneTBB / TBB
* Python 3 with NumPy and Matplotlib for plotting scripts

### Installing TBB

On macOS with Homebrew:

```bash
brew install tbb
```

On Ubuntu/Debian:

```bash
sudo apt install libtbb-dev
```

On CentOS/RHEL:

```bash
sudo yum install tbb-devel
```

If TBB is installed in a nonstandard location, specify the path manually when building:

```bash
make all TBB_PREFIX=/path/to/tbb
```

For example, on Apple Silicon macOS with Homebrew, the TBB prefix is often:

```bash
make all TBB_PREFIX=/opt/homebrew/opt/tbb
```

On Intel macOS with Homebrew, it may be:

```bash
make all TBB_PREFIX=/usr/local/opt/tbb
```

## What was cleaned

* Local absolute paths such as `/Users/mkim/Desktop/Research/fungicide_jh/` were replaced by a portable default output root `./results/`.
* The output root can be overridden at compile time using `OUTPUT_ROOT`, for example `make all OUTPUT_ROOT=/path/to/output/`.
* Input lists were moved to `configs/` with fallback support where appropriate.
* Driver filenames were not changed.
* The file `main_supp_figS1_threshold.cpp` in the uploaded batch was a duplicate of the Figure 5 driver; it has been replaced, under the same filename, with the threshold-scan driver corresponding to Supplementary Figure S1.

## Build

Build all drivers with:

```bash
./build.sh
```

or manually:

```bash
make all OUTPUT_ROOT=./results/ MACHINE_LABEL=local
```

If the build fails with an error such as:

```text
fatal error: 'tbb/parallel_for.h' file not found
```

then TBB is either not installed or its include/library path is not visible to the compiler. Install TBB and rebuild, or specify `TBB_PREFIX` manually:

```bash
make all OUTPUT_ROOT=./results/ MACHINE_LABEL=local TBB_PREFIX=/opt/homebrew/opt/tbb
```

## Example runs

Generate the main threshold-comparison table data:

```bash
./scripts/dose_time_maps_exec.sh
```

Generate Figure 5 snapshot data:

```bash
./scripts/spatial_dynamics_exec.sh
```

Generate Supplementary Figure S1 threshold data:

```bash
./scripts/threshold_exec.sh
```

Generate Supplementary Figures S2--S3 sensitivity data:

```bash
./scripts/sensitivity_exec.sh
```

Generate Main Figures 3--4 verification data:

```bash
./scripts/time_exec.sh
```

Generate Supplementary Tables S1--S2 convergence data:

```bash
./scripts/convergence_exec.sh
```

## Plotting scripts

Python plotting scripts are stored in `scripts/`. They generate the final manuscript figures from the numerical outputs written under `results/`.

```text
plot_fig2.py
plot_fig3.py
plot_fig4.py
plot_fig5.py
plot_figS1.py
plot_figS2_S3.py
```

## Output notes

The drivers can create large raw field files. For GitHub, keep only the code, configuration files, plotting scripts, and processed data needed to reproduce the manuscript figures and tables. Large raw simulation outputs should be archived separately or made available upon reasonable request.

The `results/` directory is ignored by git. It will be created automatically when the drivers are executed.

## Citation

See `CITATION.cff`. Update the DOI after acceptance or after creating a Zenodo release.
