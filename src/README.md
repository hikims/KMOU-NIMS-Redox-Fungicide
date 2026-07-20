# Shared solver files

Shared C++ solver files used by every reproducibility driver in `drivers/`.

```text
functions.hpp    LU factorization helpers, Hill switches, ADI coefficients
functions.cpp
initial.hpp      allocation and initial-condition helpers
initial.cpp
parameters.h     model, grid and time-stepping parameters
```

The drivers include these as:

```cpp
#include "functions.hpp"
#include "initial.hpp"
#include "parameters.h"
```

The Makefile adds `-Isrc`, so no relative paths are needed in the drivers.

Do not edit these files unless you intend to change the numerical model, grid
resolution, time step, baseline parameters, or initial conditions. The Python
plotting scripts duplicate a few of these values for annotation; see the note
at the end of `scripts/README_plot.md`.
